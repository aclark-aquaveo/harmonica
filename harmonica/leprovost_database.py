"""LeProvost tidal database extractor

This module contains the tidal database extractor for the LeProvost tidal database.

"""

import math
import numpy
import os

import pandas as pd

from .resource import ResourceManager
from .tidal_database import convert_coords, get_complex_components, NOAA_SPEEDS, TidalDB


DEFAULT_LEPROVOST_RESOURCE = 'leprovost'


class LeProvostDB(TidalDB):
    """Extractor class for the LeProvost tidal database.

    """
    def __init__(self, model=DEFAULT_LEPROVOST_RESOURCE):
        """Constructor for the LeProvost tidal database extractor.


        Args:
            model (:obj:`str`, optional): Name of the LeProvost tidal database version. Defaults to the freely
                distributed but outdated 'leprovost' version. See resource.py for additional supported models.

        """
        model = model.lower()  # Be case-insensitive
        if model not in ResourceManager.LEPROVOST_MODELS:  # Check for valid LeProvost model
            raise ValueError("\'{}\' is not a supported LeProvost model. Must be one of: {}.".format(
                model, ", ".join(ResourceManager.LEPROVOST_MODELS).strip()
            ))
        super().__init__(model)

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
        locs = convert_coords(locs, self.model == "fes2014")
        if not locs:
            return self  # ERROR: Not in latitude/longitude

        self.data = [pd.DataFrame(columns=['amplitude', 'phase', 'speed']) for _ in range(len(locs))]
        dataset_atts = self.resources.model_atts.dataset_attributes()

        n_lat = dataset_atts['num_lats']
        n_lon = dataset_atts['num_lons']
        lat_min = -90.0
        lon_min = dataset_atts['min_lon']
        d_lat = 180.0 / (n_lat - 1)
        d_lon = 360.0 / n_lon

        for d in self.resources.get_datasets(cons):
            if self.model == 'leprovost':  # All constituents in one file with constituent name dataset.
                nc_names = [x.strip().upper() for x in d[0].spectrum.data.tolist()]
            else:  # FES2014 has seperate files for each constituent with no constituent name dataset.
                # TODO: Probably need to find a better way to get the constituent name. _file_obj is undocumented, so
                # TODO:     there is no guarantee this functionality will be maintained.
                nc_names = [os.path.splitext(os.path.basename(dset._file_obj._filename))[0].upper() for dset in d]
            for con_idx, con in enumerate(set(cons) & set(nc_names)):
                for i, pt in enumerate(locs):
                    y_lat, x_lon = pt  # lat,lon not x,y
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
                        if self.model == 'leprovost':
                            amp_dset = d[0].amplitude[con_idx]
                            phase_dset = d[0].phase[con_idx]
                        else:
                            amp_dset = d[con_idx].amplitude
                            phase_dset = d[con_idx].phase

                        # Read potential contributing amplitudes from the file.
                        xlo_yhi_amp = amp_dset[yhi][xlo]
                        xlo_ylo_amp = amp_dset[ylo][xlo]
                        xhi_yhi_amp = amp_dset[yhi][xhi]
                        xhi_ylo_amp = amp_dset[ylo][xhi]
                        if (numpy.isnan(xlo_yhi_amp) and numpy.isnan(xhi_yhi_amp) and
                                numpy.isnan(xlo_ylo_amp) and numpy.isnan(xhi_ylo_amp)):
                            skip = True
                        else:  # Make sure we have at least one neighbor with an active phase value.
                            # Read potential contributing phases from the file.
                            xlo_yhi_phase = math.radians(phase_dset[yhi][xlo])
                            xlo_ylo_phase = math.radians(phase_dset[ylo][xlo])
                            xhi_yhi_phase = math.radians(phase_dset[yhi][xhi])
                            xhi_ylo_phase = math.radians(phase_dset[ylo][xhi])
                            if (numpy.isnan(xlo_yhi_phase) and numpy.isnan(xhi_yhi_phase) and
                                    numpy.isnan(xlo_ylo_phase) and numpy.isnan(xhi_ylo_phase)):
                                skip = True

                    if skip:
                        self.data[i].loc[con] = [numpy.nan, numpy.nan, numpy.nan]
                    else:
                        xratio = (x_lon - xlonlo) / d_lon
                        yratio = (y_lat - ylatlo) / d_lat

                        # Get the real and imaginary components from the amplitude and phases in the file. It
                        # would be better if these values were stored in the file like TPXO.
                        complex_comps = get_complex_components(
                            amps=[xlo_yhi_amp, xhi_yhi_amp, xlo_ylo_amp, xhi_ylo_amp],
                            phases=[xlo_yhi_phase, xhi_yhi_phase, xlo_ylo_phase, xhi_ylo_phase],
                        )

                        # Perform bi-linear interpolation from the four cell corners to the target point.
                        xcos = 0.0
                        xsin = 0.0
                        denom = 0.0
                        if not numpy.isnan(xlo_yhi_amp) and not numpy.isnan(xlo_yhi_phase):
                            xcos = xcos + complex_comps[0][0] * (1.0 - xratio) * yratio
                            xsin = xsin + complex_comps[0][1] * (1.0 - xratio) * yratio
                            denom = denom + (1.0 - xratio) * yratio
                        if not numpy.isnan(xhi_yhi_amp) and not numpy.isnan(xhi_yhi_phase):
                            xcos = xcos + complex_comps[1][0] * xratio * yratio
                            xsin = xsin + complex_comps[1][1] * xratio * yratio
                            denom = denom + xratio * yratio
                        if not numpy.isnan(xlo_ylo_amp) and not numpy.isnan(xlo_ylo_phase):
                            xcos = xcos + complex_comps[2][0] * (1.0 - xratio) * (1 - yratio)
                            xsin = xsin + complex_comps[2][1] * (1.0 - xratio) * (1 - yratio)
                            denom = denom + (1.0 - xratio) * (1.0 - yratio)
                        if not numpy.isnan(xhi_ylo_amp) and not numpy.isnan(xhi_ylo_phase):
                            xcos = xcos + complex_comps[3][0] * (1.0 - yratio) * xratio
                            xsin = xsin + complex_comps[3][1] * (1.0 - yratio) * xratio
                            denom = denom + (1.0 - yratio) * xratio
                        xcos = xcos / denom
                        xsin = xsin / denom
                        amp = math.sqrt(xcos * xcos + xsin * xsin)

                        # Compute interpolated phase
                        phase = math.degrees(math.acos(xcos / amp))
                        amp /= 100.0
                        if xsin < 0.0:
                            phase = 360.0 - phase
                        phase += (360. if positive_ph and phase < 0 else 0)
                        speed = NOAA_SPEEDS[con][0] if con in NOAA_SPEEDS else numpy.nan
                        self.data[i].loc[con] = [amp, phase, speed]

        return self
