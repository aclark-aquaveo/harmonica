#! python3

import math

import pandas as pd
import xmsinterp_py

from .tidal_database import convert_coords, NOAA_SPEEDS, TidalDB


class AdcircDB(TidalDB):
    """The class for extracting tidal data, specifically amplitude and phases, from an ADCIRC database.

    """
    def __init__(self):
        """Constructor for the ADCIRC tidal database extractor.

        """
        super().__init__('adcirc')
        self.resources.download_model(None)

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
                c1r = amp1*math.cos(pha1)
                c1i = amp1*math.sin(pha1)
                c2r = amp2*math.cos(pha2)
                c2i = amp2*math.sin(pha2)
                c3r = amp3*math.cos(pha3)
                c3i = amp3*math.sin(pha3)
                ctr = c1r*weights[0]+c2r*weights[1]+c3r*weights[2]
                cti = c1i*weights[0]+c2i*weights[1]+c3i*weights[2]

                new_amp = math.sqrt(ctr*ctr+cti*cti)
                self.data[i].loc[con]['amplitude'] = new_amp
                if new_amp == 0.0:
                    new_phase = 0.0
                else:
                    new_phase = math.degrees(math.acos(ctr / new_amp))
                    if cti < 0.0:
                        new_phase = 360.0 - new_phase
                self.data[i].loc[con]['phase'] = new_phase
                self.data[i].loc[con]['speed'] = NOAA_SPEEDS[con.upper()]

        return self
