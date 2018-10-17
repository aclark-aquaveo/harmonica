#! python3

from enum import Enum
import math
import os
import subprocess
import random
import string
import shutil

import pandas as pd
import xarray as xr
import xmsinterp_py

from harmonica import config

from .tidal_database import convert_coords, NOAA_SPEEDS, TidalDB


class TidalDBAdcircEnum(Enum):
    """Enum for specifying the type of an ADCIRC database.

    TIDE_NWAT: North West Atlantic database
    TIDE_NEPAC: North East Pacific database
    TIDE_NONE: Enum end - legacy from SMS port.

    """
    TIDE_NWAT = 0  # North West Atlantic Tidal Database
    TIDE_NEPAC = 1  # North East Pacific Tidal Database
    TIDE_NONE = 2  # Used if data not within any ADCIRC database domain


class AdcircDB(TidalDB):
    """The class for extracting tidal data, specifically amplitude and phases, from an ADCIRC database.

    Attributes:
        exe_with_path (str): The path of the ADCIRC executable.
        db_region (:obj:`harmonica.adcirc_database.TidalDBAdcircEnum`): The type of database.
        grid_no_path (str): Filename of \*.grd file to use.
        harm_no_path (str): Filename of \*.tdb file to use.
        temp_folder (str): Temporary folder to hold files in while running executables.
        tide_in (str): Temporary 'tides.in' filename with path.
        tide_out (str): Temporary 'tides.out' filename with path.

    """
    def __init__(self, resource_dir=None, db_region="adcircnwat"):
        """Constructor for the ADCIRC tidal database extractor.

        Args:
            resource_dir (:obj:`str`, optional): Directory of the ADCIRC resources. If not provided will become a
                subfolder of "data" in the harmonica package location.
            db_region (:obj:`str`, optional): ADCIRC tidal database region. Valid options are 'adcircnwat' and
                'adcircnepac'

        """
        if db_region.lower() == "adcircnwat":
            self.db_region = TidalDBAdcircEnum.TIDE_NWAT
            self.grid_no_path = 'ec2012_v3d_chk.grd'
            self.harm_no_path = 'ec2012_v3d_otis3_fort.53'
        elif db_region.lower() == "adcircnepac":
            self.db_region = TidalDBAdcircEnum.TIDE_NEPAC
            self.grid_no_path = 'enpac2003.grd'
            self.harm_no_path = 'enpac2003.tdb'
        else:
            raise ValueError('unrecognized ADCIRC database region.')
        super().__init__(db_region.lower())
        resource_dir = self.resources.download_model(resource_dir)
        # self.db_with_path = os.path.join(resource_dir, self.resources.model_atts["consts"][0]["M2"])
        self.grid_with_path = os.path.join(resource_dir, self.grid_no_path)
        self.harm_with_path = os.path.join(resource_dir, self.harm_no_path)

        # Build the temp working folder name
        # src_list = list(string.ascii_uppercase + string.digits)
        # rand_str = random.choice(src_list)
        # self.temp_folder = os.path.join(resource_dir, '_'.join(rand_str))
        # # check that the folder does not exist
        # while os.path.isdir(self.temp_folder):
            # rand_str = random.choice(src_list)
            # self.temp_folder = self.temp_folder + rand_str
        # self.tide_in = os.path.join(self.temp_folder, 'tides.in')
        # self.tide_out = os.path.join(self.temp_folder, 'tides.out')

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
            (:obj:`list` of :obj:`pandas.DataFrame`): A list of dataframes of constituent information including
                amplitude (meters), phase (degrees) and speed (degrees/hour, UTC/GMT). The list is parallel with locs,
                where each element in the return list is the constituent data for the corresponding element in locs.
                Note that function uses fluent interface pattern.

        """
        # pre-allocate the return value
        if not cons:
            cons = self.resources.available_constituents()  # Get all constituents by default
        else:
            cons = [con.upper() for con in cons]

        # Make sure point locations are valid lat/lon
        locs = convert_coords(locs)
        if not locs:
            return self  # ERROR: Not in latitude/longitude

        self.data = [pd.DataFrame(columns=['amplitude', 'phase', 'speed']) for _ in range(len(locs))]

        # Step 1: read the file and get geometry:
        con_dsets = self.resources.get_datasets(cons)[0]
        tri_search = xmsinterp_py.geometry.GmTriSearch()
        con_x = con_dsets.x[0].values
        con_y = con_dsets.y[0].values

        mesh_pts = [(float(con_x[idx]), float(con_y[idx]), 0.0) for idx in range(len(con_x))]
        tri_list = con_dsets.element.values.flatten().tolist()
        tri_search.tris_to_search(mesh_pts, tri_list)

        points_and_weights = []
        elem_dset = con_dsets.element[0]
        max_tri = len(elem_dset)
        for pt in locs:
            pt_flip = (pt[1], pt[0])
            tri_idx = tri_search.tri_containing_pt(pt_flip)
            if tri_idx != -1:
                pt_1 = int(tri_list[tri_idx])
                pt_2 = int(tri_list[tri_idx+1])
                pt_3 = int(tri_list[tri_idx+2])
                x1, y1, z1 = mesh_pts[pt_1]
                x2, y2, z2 = mesh_pts[pt_2]
                x3, y3, z3 = mesh_pts[pt_3]
                x = pt_flip[0]
                y = pt_flip[1]
                s1 = abs((x2*y3-x3*y2)-(x*y3-x3*y)+(x*y2-x2*y))
                s2 = abs((x*y3-x3*y)-(x1*y3-x3*y1)+(x1*y-x*y1))
                s3 = abs((x2*y-x*y2)-(x1*y-x*y1)+(x1*y2-x2*y1))
                ta = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                w1 = ((x-x3)*(y2-y3)+(x2-x3)*(y3-y))/ta
                w2 = ((x-x1)*(y3-y1)-(y-y1)*(x3-x1))/ta
                w3 = ((y-y1)*(x2-x1)-(x-x1)*(y2-y1))/ta
                points_and_weights.append(((pt_1, pt_2, pt_3), (w1, w2, w3)))

        for con in cons:
            con_amp_name = con + "_amplitude"
            con_pha_name = con + "_phase"
            con_amp = con_dsets[con_amp_name][0]
            con_pha = con_dsets[con_pha_name][0]
            for i, (pts, weights) in enumerate(points_and_weights):
                amp1 = float(con_amp[pts[0]])
                amp2 = float(con_amp[pts[1]])
                amp3 = float(con_amp[pts[2]])
                pha1 = math.radians(float(con_pha[pts[0]]))
                pha2 = math.radians(float(con_pha[pts[1]]))
                pha3 = math.radians(float(con_pha[pts[2]]))

                self.data[i].loc[con] = [0.0, 0.0, 0.0]
                C1R=amp1*math.cos(pha1)
                C1I=amp1*math.sin(pha1)
                C2R=amp2*math.cos(pha2)
                C2I=amp2*math.sin(pha2)
                C3R=amp3*math.cos(pha3)
                C3I=amp3*math.sin(pha3)
                CTR=C1R*weights[0]+C2R*weights[1]+C3R*weights[2]
                CTI=C1I*weights[0]+C2I*weights[1]+C3I*weights[2]

                new_amp = math.sqrt(CTR*CTR+CTI*CTI)
                self.data[i].loc[con]['amplitude'] = new_amp
                if new_amp == 0.0:
                   new_phase = 0.0
                else:
                   new_phase = math.degrees(math.acos(CTR / new_amp))
                   if CTI < 0.0:
                       new_phase = 360.0 - new_phase
                self.data[i].loc[con]['phase'] = new_phase
                self.data[i].loc[con]['speed'] = NOAA_SPEEDS[con.upper()]

        return self
