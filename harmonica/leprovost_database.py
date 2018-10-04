"""LeProvost tidal database extractor

This module contains the tidal database extractor for the LeProvost tidal database.

"""

import math
import os

import pandas as pd

from .tidal_database import convert_coords, NOAA_SPEEDS, TidalDB


class LeProvostDB(TidalDB):
    """Extractor class for the LeProvost tidal database.

    Attributes:
        resource_dir (str): Fully qualified path to a folder containing the LeProvost *.legi files

    """
    def __init__(self, resource_dir=None):
        """Constructor for the LeProvost tidal database extractor.

        Args:
            resource_dir (:obj:`str`, optional): Fully qualified path to a folder containing the LeProvost *.legi files.
                If not provided will be "harmonica/data/leprovost' where harmonica is the location of the pacakage.

        """
        # self.debugger = open("debug.txt", "w")
        super().__init__("leprovost")
        self.resource_dir = self.resources.download_model(resource_dir)

    def get_components(self, locs, cons=None, positive_ph=False):
        """Get the amplitude, phase, and speed of specified constituents at specified point locations.

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
                Empty list on error. Note that function uses fluent interface pattern.

        """
        # If no constituents specified, extract all valid constituents.
        if not cons:
            cons = self.resources.available_constituents()

        # Make sure point locations are valid lat/lon
        locs = convert_coords(locs)
        if not locs:
            return self  # ERROR: Not in latitude/longitude

        self.data = [pd.DataFrame(columns=['amplitude', 'phase', 'speed']) for _ in range(len(locs))]

        deg2rad = 1.0 / 180.0 * math.pi
        rad2deg = 1.0 / deg2rad
        # read the file for each constituent
        for con in cons:
            if self.have_constituent(con):
                con = con.upper()
            else:
                continue
            filename = os.path.join(self.resource_dir, self.resources.model_atts["consts"][0][con])
            with open(filename, 'r') as f:
                # 30 columns in the file
                lon_min, lon_max = map(float, f.readline().split())
                lat_min, lat_max = map(float, f.readline().split())
                d_lon, d_lat = map(float, f.readline().split())
                n_lon, n_lat = map(int, f.readline().split())
                undef_a, undef_p = map(float, f.readline().split())
                if undef_a != undef_p:
                    # there was an error, the undefined values should be the same
                    return []
                else:
                    undef = undef_a
                all_lines = f.readlines()
                g_amp = [[0.0 for _ in range(n_lat)] for _ in range(n_lon)]
                g_pha = [[0.0 for _ in range(n_lat)] for _ in range(n_lon)]
                cur_line = -1
                for lat in range(n_lat):
                    for lon in range(0, n_lon-1, 30):
                        cur_line += 1
                        cur_val = 0
                        line_vals = all_lines[cur_line].split()
                        for k in range(lon, lon+30):
                            g_amp[k][lat] = float(line_vals[cur_val])
                            cur_val += 1
                        cur_line += 1
                        cur_val = 0
                        line_vals = all_lines[cur_line].split()
                        for k in range(lon, lon+30):
                            g_pha[k][lat] = float(line_vals[cur_val])
                            cur_val += 1

            # Extract components for each point for this constituent
            for i, pt in enumerate(locs):
                y_lat = pt[0]
                x_lon = pt[1]
                if x_lon < 0.0:
                    x_lon = x_lon + 360.0
                if x_lon > 180.0:
                    x_lon = x_lon - 360.0
                ixlo = int((x_lon - lon_min) / d_lon) + 1
                xlonlo = lon_min + (ixlo - 1) * d_lon
                ixhi = ixlo + 1
                if ixlo == n_lon:
                    ixhi = 1
                iylo = int((y_lat - lat_min) / d_lon) + 1
                ylatlo = lat_min + (iylo - 1) * d_lat
                iyhi = iylo + 1
                ixlo -= 1
                ixhi -= 1
                iylo -= 1
                iyhi -= 1
                skip = False

                if (ixlo > n_lon or ixhi > n_lon or iyhi > n_lat or iylo > n_lat or ixlo < 0 or ixhi < 0 or iyhi < 0 or
                        iylo < 0):
                    skip = True
                elif (g_amp[ixlo][iyhi] == undef and g_amp[ixhi][iyhi] == undef and g_amp[ixlo][iylo] == undef and
                      g_amp[ixhi][iylo] == undef):
                    skip = True
                elif (g_pha[ixlo][iyhi] == undef and g_pha[ixhi][iyhi] == undef and g_pha[ixlo][iylo] == undef and
                      g_pha[ixhi][iylo] == undef):
                    skip = True

                if skip:
                    self.data[i].loc[con] = [0.0, 0.0, 0.0]
                else:
                    xratio = (x_lon - xlonlo) / d_lon
                    yratio = (y_lat - ylatlo) / d_lat
                    xcos1 = g_amp[ixlo][iyhi] * math.cos(deg2rad * g_pha[ixlo][iyhi])
                    xcos2 = g_amp[ixhi][iyhi] * math.cos(deg2rad * g_pha[ixhi][iyhi])
                    xcos3 = g_amp[ixlo][iylo] * math.cos(deg2rad * g_pha[ixlo][iylo])
                    xcos4 = g_amp[ixhi][iylo] * math.cos(deg2rad * g_pha[ixhi][iylo])
                    xsin1 = g_amp[ixlo][iyhi] * math.sin(deg2rad * g_pha[ixlo][iyhi])
                    xsin2 = g_amp[ixhi][iyhi] * math.sin(deg2rad * g_pha[ixhi][iyhi])
                    xsin3 = g_amp[ixlo][iylo] * math.sin(deg2rad * g_pha[ixlo][iylo])
                    xsin4 = g_amp[ixhi][iylo] * math.sin(deg2rad * g_pha[ixhi][iylo])

                    xcos = 0.0
                    xsin = 0.0
                    denom = 0.0
                    if g_amp[ixlo][iyhi] != undef and g_pha[ixlo][iyhi] != undef:
                        xcos = xcos + xcos1 * (1.0 - xratio) * yratio
                        xsin = xsin + xsin1 * (1.0 - xratio) * yratio
                        denom = denom + (1.0 - xratio) * yratio
                    if g_amp[ixhi][iyhi] != undef and g_pha[ixhi][iyhi] != undef:
                        xcos = xcos + xcos2 * xratio * yratio
                        xsin = xsin + xsin2 * xratio * yratio
                        denom = denom + xratio * yratio
                    if g_amp[ixlo][iylo] != undef and g_pha[ixlo][iylo] != undef:
                        xcos = xcos + xcos3 * (1.0 - xratio) * (1 - yratio)
                        xsin = xsin + xsin3 * (1.0 - xratio) * (1 - yratio)
                        denom = denom + (1.0 - xratio) * (1.0 - yratio)
                    if g_amp[ixhi][iylo] != undef and g_pha[ixhi][iylo] != undef:
                        xcos = xcos + xcos4 * (1.0 - yratio) * xratio
                        xsin = xsin + xsin4 * (1.0 - yratio) * xratio
                        denom = denom + (1.0 - yratio) * xratio

                    xcos = xcos / denom
                    xsin = xsin / denom

                    amp = math.sqrt(xcos * xcos + xsin * xsin)
                    phase = rad2deg * math.acos(xcos / amp)
                    amp /= 100.0
                    if xsin < 0.0:
                        phase = 360.0 - phase
                    phase += (360. if positive_ph and phase < 0 else 0)
                    self.data[i].loc[con] = [amp, phase, NOAA_SPEEDS[con]]

        return self
