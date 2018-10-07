"""LeProvost tidal database extractor

This module contains the tidal database extractor for the LeProvost tidal database.

"""

import math
import numpy

import pandas as pd

from .tidal_database import convert_coords, NOAA_SPEEDS, TidalDB


class LeProvostDB(TidalDB):
    """Extractor class for the LeProvost tidal database.

    Attributes:
        resource_dir (str): Fully qualified path to a folder containing the LeProvost \*.legi files

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
        else:  # Be case-insensitive
            cons = [con.upper() for con in cons]

        # Make sure point locations are valid lat/lon
        locs = convert_coords(locs)
        if not locs:
            return self  # ERROR: Not in latitude/longitude

        self.data = [pd.DataFrame(columns=['amplitude', 'phase', 'speed']) for _ in range(len(locs))]
        dataset_atts = self.resources.model_atts['dataset_atts']

        n_lat = dataset_atts['num_lats']
        n_lon = dataset_atts['num_lons']
        lat_min = -90.0
        lon_min = -180.0
        d_lat = 180.0 / (n_lat - 1)
        d_lon = 360.0 / n_lon

        for d in self.resources.get_datasets(cons):
            nc_names = [x.strip().upper() for x in d.spectrum.values[0]]
            for con in set(cons) & set(nc_names):
                # Extract components for each point for this constituent
                con_idx = nc_names.index(con)
                for i, pt in enumerate(locs):
                    y_lat = pt[0]
                    x_lon = pt[1]
                    if x_lon < 0.0:
                        x_lon = x_lon + 360.0
                    if x_lon > 180.0:
                        x_lon = x_lon - 360.0
                    xlo = int((x_lon - lon_min) / d_lon) + 1
                    xlonlo = lon_min + (xlo - 1) * d_lon
                    xhi = xlo + 1
                    if xlo == n_lon:
                        xhi = 1
                    ylo = int((y_lat - lat_min) / d_lon) + 1
                    ylatlo = lat_min + (ylo - 1) * d_lat
                    yhi = ylo + 1
                    xlo -= 1
                    xhi -= 1
                    ylo -= 1
                    yhi -= 1
                    skip = False

                    # Make sure lat/lon coordinate is in the domain.
                    if (xlo > n_lon or xhi > n_lon or yhi > n_lat or ylo > n_lat or xlo < 0 or xhi < 0 or yhi < 0
                            or ylo < 0):
                        skip = True
                    else:  # Make sure we have at least one neighbor with an active amplitude value.
                        # Read potential contributing amplitudes from the file.
                        xlo_yhi_amp = d.Ha[0][con_idx][xlo][yhi]
                        xlo_ylo_amp = d.Ha[0][con_idx][xlo][ylo]
                        xhi_yhi_amp = d.Ha[0][con_idx][xhi][yhi]
                        xhi_ylo_amp = d.Ha[0][con_idx][xhi][ylo]
                        if (numpy.isnan(xlo_yhi_amp) and numpy.isnan(xhi_yhi_amp) and
                                numpy.isnan(xlo_ylo_amp) and numpy.isnan(xhi_ylo_amp)):
                            skip = True
                        else:  # Make sure we have at least one neighbor with an active phase value.
                            # Read potential contributing phases from the file.
                            xlo_yhi_phase = d.Hg[0][con_idx][xlo][yhi]
                            xlo_ylo_phase = d.Hg[0][con_idx][xlo][ylo]
                            xhi_yhi_phase = d.Hg[0][con_idx][xhi][yhi]
                            xhi_ylo_phase = d.Hg[0][con_idx][xhi][ylo]
                            if (numpy.isnan(xlo_yhi_phase) and numpy.isnan(xhi_yhi_phase) and
                                    numpy.isnan(xlo_ylo_phase) and numpy.isnan(xhi_ylo_phase)):
                                skip = True

                    if skip:
                        self.data[i].loc[con] = [0.0, 0.0, 0.0]
                    else:
                        xratio = (x_lon - xlonlo) / d_lon
                        yratio = (y_lat - ylatlo) / d_lat
                        xcos1 = xlo_yhi_amp * math.cos(math.radians(xlo_yhi_phase))
                        xcos2 = xhi_yhi_amp * math.cos(math.radians(xhi_yhi_phase))
                        xcos3 = xlo_ylo_amp * math.cos(math.radians(xlo_ylo_phase))
                        xcos4 = xhi_ylo_amp * math.cos(math.radians(xhi_ylo_phase))
                        xsin1 = xlo_yhi_amp * math.sin(math.radians(xlo_yhi_phase))
                        xsin2 = xhi_yhi_amp * math.sin(math.radians(xhi_yhi_phase))
                        xsin3 = xlo_ylo_amp * math.sin(math.radians(xlo_ylo_phase))
                        xsin4 = xhi_ylo_amp * math.sin(math.radians(xhi_ylo_phase))

                        xcos = 0.0
                        xsin = 0.0
                        denom = 0.0
                        if not numpy.isnan(xlo_yhi_amp) and not numpy.isnan(xlo_yhi_phase):
                            xcos = xcos + xcos1 * (1.0 - xratio) * yratio
                            xsin = xsin + xsin1 * (1.0 - xratio) * yratio
                            denom = denom + (1.0 - xratio) * yratio
                        if not numpy.isnan(xhi_yhi_amp) and not numpy.isnan(xhi_yhi_phase):
                            xcos = xcos + xcos2 * xratio * yratio
                            xsin = xsin + xsin2 * xratio * yratio
                            denom = denom + xratio * yratio
                        if not numpy.isnan(xlo_ylo_amp) and not numpy.isnan(xlo_ylo_phase):
                            xcos = xcos + xcos3 * (1.0 - xratio) * (1 - yratio)
                            xsin = xsin + xsin3 * (1.0 - xratio) * (1 - yratio)
                            denom = denom + (1.0 - xratio) * (1.0 - yratio)
                        if not numpy.isnan(xhi_ylo_amp) and not numpy.isnan(xhi_ylo_phase):
                            xcos = xcos + xcos4 * (1.0 - yratio) * xratio
                            xsin = xsin + xsin4 * (1.0 - yratio) * xratio
                            denom = denom + (1.0 - yratio) * xratio

                        xcos = xcos / denom
                        xsin = xsin / denom

                        amp = math.sqrt(xcos * xcos + xsin * xsin)
                        phase = math.degrees(math.acos(xcos / amp))
                        amp /= 100.0
                        if xsin < 0.0:
                            phase = 360.0 - phase
                        phase += (360. if positive_ph and phase < 0 else 0)
                        self.data[i].loc[con] = [amp, phase, NOAA_SPEEDS[con]]

        return self
