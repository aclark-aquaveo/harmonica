#! python3

from enum import Enum
import os
import subprocess
import random
import string
import shutil
import urllib.request
from zipfile import ZipFile

import pandas as pd

from .tidal_constituents import NOAA_SPEEDS
from .tidal_database import convert_coords, TidalDB


class TidalDBAdcircEnum(Enum):
    """Enum for specifying the type of an ADCIRC database.

    TIDE_NWAT = North West Atlantic database
    TIDE_NEPAC = North East Pacific database
    TIDE_NONE = Enum end - legacy from SMS port.

    """
    TIDE_NWAT = 0  # North West Atlantic Tidal Database
    TIDE_NEPAC = 1  # North East Pacific Tidal Database
    TIDE_NONE = 2  # Used if data not within any ADCIRC database domain


class AdcircDB(TidalDB):
    """The class for extracting tidal data, specifically amplitude and phases, from an ADCIRC database.

    Attributes:
        work_path (str): The path of the working directory.
        exe_with_path (str): The path of the ADCIRC executable.
        db_region (:obj: `TidalDBAdcircEnum`): The type of database.
        cons (:obj:`list` of :obj:`str`): List of the constituents that are valid for the ADCIRC database
        data (:obj:`list` of :obj:`pandas.DataFrame`): List of the constituent component DataFrames with one
            per point location requested from get_components(). Intended return value of get_components().
        grid_no_path (str): Filename of *.grd file to use.
        harm_no_path (str): Filename of *.tdb file to use.
        temp_folder (str): Temporary folder to hold files in while running executables.
        tide_in (str): Temporary 'tides.in' filename with path.
        tide_out (str): Temporary 'tides.out' filename with path.

    """
    def __init__(self, work_path, exe_with_path, db_region):
        """Get the amplitude and phase for the given constituents at the given points.

        Args:
            work_path (str): The path of the working directory.
            exe_with_path (str): The path of the ADCIRC executable.
            db_region (:obj: `TidalDBAdcircEnum`): The db_region of database.

        """
        TidalDB.__init__(self)
        self.work_path = work_path
        self.exe_with_path = exe_with_path
        self.db_region = db_region
        self.cons = ['M2', 'S2', 'N2', 'K1', 'M4', 'O1', 'M6', 'Q1', 'K2']
        self.data = []
        if self.db_region == TidalDBAdcircEnum.TIDE_NWAT:
            self.grid_no_path = 'ec2001.grd'
            self.harm_no_path = 'ec2001.tdb'
        elif self.db_region == TidalDBAdcircEnum.TIDE_NEPAC:
            self.grid_no_path = 'enpac2003.grd'
            self.harm_no_path = 'enpac2003.tdb'
        else:
            self.grid_no_path = ''
            self.harm_no_path = ''
        src_list = list(string.ascii_uppercase + string.digits)
        rand_str = random.choice(src_list)
        self.temp_folder = os.path.join(self.work_path, '_'.join(rand_str))
        # check that the folder does not exist
        while os.path.isdir(self.temp_folder):
            rand_str = random.choice(src_list)
            self.temp_folder = self.temp_folder + rand_str
        self.tide_in = os.path.join(self.temp_folder, 'tides.in')
        self.tide_out = os.path.join(self.temp_folder, 'tides.out')
        self.validate_files()

    def validate_files(self):
        """Check to ensure that the ADCIRC database exists.

        If the database is nonexistent, a new ADCIRC database will be downloaded from Aquaveo.com. Whether the Atlantic
        or Pacific database is download depends on how the object was constructed.

        """
        if not os.path.isdir(self.work_path):
            self.work_path = os.getcwd()

        if not os.path.isfile(self.exe_with_path):  # Make sure the executable exists
            # Download from the Aquaveo website
            if self.db_region == TidalDBAdcircEnum.TIDE_NWAT:  # Download the ADCIRC Atlantic database
                adcirc_db_url = 'http://sms.aquaveo.com/adcircnwattides.zip'
                basename = 'adcircnwattides'
                if not self.exe_with_path:  # check the local directory
                    local_exe = os.path.join(os.getcwd(), "adcircnwattides/adcircnwattides.exe")
                    if os.path.isfile(local_exe):
                        self.exe_with_path = local_exe
                        return
            else:  # Download the ADCIRC Pacific database
                adcirc_db_url = 'http://sms.aquaveo.com/adcircnepactides.zip'
                basename = 'adcircnepactides'
                if not self.exe_with_path:  # check the local directory
                    local_exe = os.path.join(os.getcwd(), "adcircnepactides/adcircnepactides.exe")
                    if os.path.isfile(local_exe):
                        self.exe_with_path = local_exe
                        return
            zip_file = os.path.join(self.work_path, basename + '.zip')
            print("Downloading resource: {}".format(adcirc_db_url))
            with urllib.request.urlopen(adcirc_db_url) as response, open(zip_file, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
            # Unzip the files
            dest_path = os.path.join(self.work_path, basename)
            if not os.path.isdir(dest_path):
                os.mkdir(dest_path)
            print("Unzipping files to: {}".format(dest_path))
            with ZipFile(zip_file, 'r') as unzipper:
                unzipper.extractall(path=dest_path)
            print("Deleting zip file: {}".format(zip_file))
            os.remove(zip_file)  # delete the zip file
            self.exe_with_path = os.path.join(dest_path, basename + ".exe")

    def get_components(self, locs, cons=None, positive_ph=False):
        """Get the amplitude, phase, and speed for the given constituents at the given points.

        Args:
            locs (:obj:`list` of :obj:`tuple` of :obj:`float`): latitude [-90, 90] and longitude [-180 180] or [0 360]
                of the requested points.
            cons (:obj:`list` of :obj:`str`, optional): List of the constituent names to get amplitude and phase for. If
                not supplied, all valid constituents will be extracted.
            positive_ph (bool, optional): Indicate if the returned phase should be all positive [0 360] (True) or
                [-180 180] (False, the default).

        Returns:
            :obj:`list` of :obj:`pandas.DataFrame`: A list of dataframes of constituent information including
                amplitude (meters), phase (degrees) and speed (degrees/hour, UTC/GMT). The list is parallel with locs,
                where each element in the return list is the constituent data for the corresponding element in locs.
                Note that function uses fluent interface pattern.

        """
        # pre-allocate the return value
        constituents = []
        if not cons:
            cons = self.cons  # Get all constituents by default

        # Make sure point locations are valid lat/lon
        locs = convert_coords(locs)
        if not locs:
            return self  # ERROR: Not in latitude/longitude

        self.data = [pd.DataFrame(columns=['amplitude', 'phase', 'speed']) for _ in range(len(locs))]
        for con in cons:
            if self.have_constituent(con):
                con = con.lower()
                constituents.append(con)

        # create the temp directory
        os.mkdir(self.temp_folder)
        os.chdir(self.temp_folder)
        try:
            # write tides.in into the temp directory
            with open(self.tide_in, 'w') as f:
                f.write("{}\n".format(str(len(locs))))
                for pt in locs:
                    # 15.10f
                    f.write("{:15.10f}{:15.10f}\n".format(pt[1], pt[0]))
            # copy the executable and .grd and .tdb file to the temp folder
            temp_exe_with_path = os.path.join(self.temp_folder, os.path.basename(self.exe_with_path))
            temp_grd_with_path = os.path.join(self.temp_folder, self.grid_no_path)
            temp_tdb_with_path = os.path.join(self.temp_folder, self.harm_no_path)

            old_exe_path = os.path.dirname(self.exe_with_path)
            old_grd_with_path = os.path.join(old_exe_path, self.grid_no_path)
            old_tdb_with_path = os.path.join(old_exe_path, self.harm_no_path)

            shutil.copyfile(self.exe_with_path, temp_exe_with_path)
            shutil.copyfile(old_grd_with_path, temp_grd_with_path)
            shutil.copyfile(old_tdb_with_path, temp_tdb_with_path)
            # run the executable
            subprocess.run([temp_exe_with_path])
            # read tides.out from the temp directory
            with open(self.tide_out, 'r') as f:
                all_lines = f.readlines()
                last_name_line = 0
                is_first = True
                used_constituent = False
                con_name = ''
                column_count = 0
                curr_pt = 0
                for count, line in enumerate(all_lines):
                    line = line.strip()
                    if not line:
                        continue
                    elif any(not(c.isdigit() or c.isspace() or c == 'e' or c == 'E'
                                 or c == '.' or c == '-' or c == '+') for c in line):
                        last_name_line = count
                        is_first = True
                        curr_pt = 0
                        continue
                    elif is_first:
                        con_name = all_lines[last_name_line].strip().lower()
                        is_first = False
                        if any(con_name in name for name in constituents):
                            used_constituent = True
                        else:
                            used_constituent = False

                    if used_constituent:
                        if column_count == 0:
                            line_nums = [float(num) for num in line.split()]
                            column_count = len(line_nums)
                        if column_count == 2:
                            amp, pha = map(float, line.split())
                        elif column_count == 4:
                            junk, junk, amp, pha = map(float, line.split())
                        elif column_count == 6:
                            amp, pha, junk, junk, junk, junk = map(float, line.split())
                        elif column_count == 8:
                            junk, junk, amp, pha, junk, junk, junk, junk = map(float, line.split())
                        else:
                            # we have a problem
                            continue
                        self.data[curr_pt].loc[con_name.upper()] = [amp,
                                                                    pha + (360.0 if positive_ph and pha < 0 else 0),
                                                                    NOAA_SPEEDS[con_name.upper()]]
                        curr_pt += 1
        finally:
            # delete the temp directory
            del_files = os.listdir(self.temp_folder)
            for del_file in del_files:
                os.remove(del_file)
            os.chdir(self.work_path)
            os.rmdir(self.temp_folder)

        return self

    def have_constituent(self, a_name):
        """Check if teh given constituent is supported by the ADCIRC tidal database.

        Args:
            a_name (str): The name of the constituent to check.

        Returns:
            True if the constituent is supported by ADCIRC tidal databases, False for unsupported.

        """
        if a_name.upper() in self.cons:
            return True
        else:
            return False
