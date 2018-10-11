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
        #con_dsets.element
        #con_dsets.x
        #con_dsets.y
        tri_search = xmsinterp_py.geometry.GmTriSearch()
        print(str(len(con_dsets.x[0])))
        con_x = con_dsets.x[0].values
        con_y = con_dsets.y[0].values
        #for idx in range(len(con_dsets.x[0])):
        #    print("idx: {:10} X: {:8f} Y: {:8f}".format(idx, float(con_x[idx]), float(con_y[idx])))
        # for idx, (the_x, the_y) in enumerate(zip(con_x, con_y)):
        #     print("idx: {:10} X: {:8f} Y: {:8f}".format(idx, the_x, the_y))
        # exit()
        mesh_pts = [(float(con_y[idx]), float(con_x[idx]), 0.0) for idx in range(len(con_x))]
        print("Built mesh_pts")
        tri_list = con_dsets.element.values.flatten().tolist()
        print("len: {} {} {}".format(len(mesh_pts), len(con_dsets.element[0]), len(tri_list)))
        # print(str(mesh_pts[:10]))
        tri_search.tris_to_search(mesh_pts, tri_list)
        points_and_weights = []
        elem_dset = con_dsets.element[0]
        max_tri = len(elem_dset)
        for pt in locs:
            tri_idx = tri_search.tri_containing_pt(pt)
            if tri_idx != -1 and tri_idx < max_tri:
                print("Point found! {} {}".format(pt, tri_idx))
                # print("eLen: {}".format(type(elem_dset[tri_idx][0])))
                pt_1 = int(elem_dset[tri_idx][0])
                pt_2 = int(elem_dset[tri_idx][1])
                pt_3 = int(elem_dset[tri_idx][2])
                x1, y1, z1 = mesh_pts[pt_1]
                x2, y2, z2 = mesh_pts[pt_2]
                x3, y3, z3 = mesh_pts[pt_3]
                x = pt[0]
                y = pt[1]
                s1 = abs((x2*y3-x3*y2)-(x*y3-x3*y)+(x*y2-x2*y))
                s2 = abs((x*y3-x3*y)-(x1*y3-x3*y1)+(x1*y-x*y1))
                s3 = abs((x2*y-x*y2)-(x1*y-x*y1)+(x1*y2-x2*y1))
                ta = abs((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))
                w1 = ((x-x3)*(y2-y3)+(x2-x3)*(y3-y))/ta
                w2 = ((x-x1)*(y3-y1)-(y-y1)*(x3-x1))/ta
                w3 = ((y-y1)*(x2-x1)-(x-x1)*(y2-y1))/ta
                points_and_weights.append(((pt_1, pt_2, pt_3), (w1, w2, w3)))
            else:
                print("Not found! {}".format(pt))

        for con in cons:
            con_amp_name = con + "_amplitude"
            con_pha_name = con + "_phase"
            for i, (pts, weights) in enumerate(points_and_weights):
                # print("AMP: {} {}".format(con_amp_name, con_dsets.keys()))
                con_amp = con_dsets[con_amp_name][0]
                con_pha = con_dsets[con_pha_name][0]
                # print("pts0: {} {}".format(pts[0], len(con_amp)))
                amp1 = con_amp[pts[0]]
                amp2 = con_amp[pts[1]]
                amp3 = con_amp[pts[2]]
                pha1 = con_pha[pts[0]]
                pha2 = con_pha[pts[1]]
                pha3 = con_pha[pts[2]]

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
                   new_phase = 180.0 * math.acos(CTR / new_amp) / math.pi
                   if CTI < 0.0:
                       new_phase = 360.0 - new_phase
                self.data[i].loc[con]['phase'] = new_phase
        
        # con_dset_file = self.resources[self.model]["consts"][cons[0]]
        # if config['pre_existing_data_dir']:
            # path = os.path.join(config['pre_existing_data_dir'], self.model, con_dset_file)
        # else:
            # path = os.path.join(config['data_dir'], self.model, con_dset_file)
        # geom = xr.open_dataset(path, group='element')
        # pt_x = xr.open_dataset(path, group='x')
        # pt_y = xr.open_dataset(path, group='y')
        # df = element.to_dataframe()
        # all_grid = pd.read_table(self.grid_with_path, delim_whitespace=True, header=None, skiprows=skiprows,
                               # names=('row_type', 'cmp1', 'cmp2', 'cmp3', 'val'), index_col=1)
        # conns = all_grid[all_grid['row_type'].str.lower() == 'e3t'][['cmp1', 'cmp2', 'cmp3']].values.astype(int) - 1
        # pts = all_grid[all_grid['row_type'].str.lower() == 'nd'][['cmp1', 'cmp2', 'cmp3']].values.astype(float)
        # pts[:, 2] *= -1
        # verts = pd.DataFrame(pts, columns=['x', 'y', 'z'])
        # tris = pd.DataFrame(conns, columns=['v0', 'v1', 'v2'])
        
        # convert to list of tuples:
        # subset = data_set[['data_date', 'data_1', 'data_2']]
        # tuples = [tuple(x) for x in subset.values]
        
        # feed to xmsinterp:
        # element and weights (maybe, if not compute)
        
        # read the harmonic files:
        # use pandas
        # with open(self.harm_with_path) as f:
            # num_cons = int(f.readline())
        # for i in range(num_cons):
            # con_harm = pd.read_table(self.harm_with_path, delim_whitespace=True, header=None, skiprows=lambda x: x < num_cons+2 and x % i,
                                   # names=('cmp1', 'cmp2'))
        
        # find the values at the point locations:
        # ComputeHarmonics
        
        # return as a list
        

        # create the temp directory
        # os.makedirs(self.temp_folder)
        # os.chdir(self.temp_folder)
        # try:
            # # write tides.in into the temp directory
            # with open(self.tide_in, 'w') as f:
                # f.write("{}\n".format(str(len(locs))))
                # for pt in locs:
                    # # 15.10f
                    # f.write("{:15.10f}{:15.10f}\n".format(pt[1], pt[0]))
            # # copy the executable and .grd and .tdb file to the temp folder
            # temp_exe_with_path = os.path.join(self.temp_folder, os.path.basename(self.exe_with_path))
            # temp_grd_with_path = os.path.join(self.temp_folder, self.grid_no_path)
            # temp_tdb_with_path = os.path.join(self.temp_folder, self.harm_no_path)

            # old_exe_path = os.path.dirname(self.exe_with_path)
            # old_grd_with_path = os.path.join(old_exe_path, self.grid_no_path)
            # old_tdb_with_path = os.path.join(old_exe_path, self.harm_no_path)

            # shutil.copyfile(self.exe_with_path, temp_exe_with_path)
            # shutil.copyfile(old_grd_with_path, temp_grd_with_path)
            # shutil.copyfile(old_tdb_with_path, temp_tdb_with_path)
            # # run the executable
            # subprocess.run([temp_exe_with_path])
            # # read tides.out from the temp directory
            # with open(self.tide_out, 'r') as f:
                # all_lines = f.readlines()
                # last_name_line = 0
                # is_first = True
                # used_constituent = False
                # con_name = ''
                # column_count = 0
                # curr_pt = 0
                # for count, line in enumerate(all_lines):
                    # line = line.strip()
                    # if not line:
                        # continue
                    # elif any(not(c.isdigit() or c.isspace() or c == 'e' or c == 'E'
                                 # or c == '.' or c == '-' or c == '+') for c in line):
                        # last_name_line = count
                        # is_first = True
                        # curr_pt = 0
                        # continue
                    # elif is_first:
                        # con_name = all_lines[last_name_line].strip().lower()
                        # is_first = False
                        # if any(con_name in name for name in constituents):
                            # used_constituent = True
                        # else:
                            # used_constituent = False

                    # if used_constituent:
                        # if column_count == 0:
                            # line_nums = [float(num) for num in line.split()]
                            # column_count = len(line_nums)
                        # if column_count == 2:
                            # amp, pha = map(float, line.split())
                        # elif column_count == 4:
                            # junk, junk, amp, pha = map(float, line.split())
                        # elif column_count == 6:
                            # amp, pha, junk, junk, junk, junk = map(float, line.split())
                        # elif column_count == 8:
                            # junk, junk, amp, pha, junk, junk, junk, junk = map(float, line.split())
                        # else:
                            # # we have a problem
                            # continue
                        # self.data[curr_pt].loc[con_name.upper()] = [amp,
                                                                    # pha + (360.0 if positive_ph and pha < 0 else 0),
                                                                    # NOAA_SPEEDS[con_name.upper()]]
                        # curr_pt += 1
        # finally:
            # # delete the temp directory
            # del_files = os.listdir(self.temp_folder)
            # for del_file in del_files:
                # os.remove(del_file)
            # os.chdir(os.path.dirname(self.temp_folder))
            # os.rmdir(self.temp_folder)

        return self
