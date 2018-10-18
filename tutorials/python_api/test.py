import os

from harmonica.tidal_constituents import Constituents

if __name__ == "__main__":
    # Need to be in (lat, lon), not (x, y)
    all_points = [(39.74, -74.07),
                  (42.32, -70.0),
                  (45.44, -65.0),
                  (43.63, -124.55),
                  (46.18, -124.38)]

    cons = ['M2', 'S2', 'N2', 'K1']
    constituents = Constituents('leprovost')

    # Get astronomical nodal factor data (not dependent on the tidal model)
    nodal_factors = constituents.get_nodal_factor(cons, 15, 30, 8, 2018)

    f = open(os.path.join(os.getcwd(), "tidal_test.out"), "w")
    f.write("Nodal factor:\n")
    f.write(nodal_factors.to_string() + "\n\n")
    f.flush()

    # Get tidal harmonic components for a list of points using the ADCIRC,
    # LeProvost, and TPXO databases.
    f.write("LeProvost components:\n")
    f.flush()
    leprovost_comps = constituents.get_components(all_points, cons)
    for pt in leprovost_comps.data:
        f.write(pt.sort_index().to_string() + "\n\n")

    ad_atlantic_comps = constituents.get_components(all_points, cons, model='adcirc2015')
    f.write("ADCIRC components:\n")
    for pt in ad_atlantic_comps.data:
        f.write(pt.sort_index().to_string() + "\n\n")

    f.write("TPX0 components:\n")
    f.flush()
    tpxo_comps = constituents.get_components(all_points, cons, True, 'tpxo8')
    for pt in tpxo_comps.data:
        f.write(pt.sort_index().to_string() + "\n\n")

    f.write("FES2014 components:\n")
    f.flush()
    fes2014_comps = constituents.get_components(all_points, cons, model='fes2014')
    for pt in fes2014_comps.data:
        f.write(pt.sort_index().to_string() + "\n\n")

    f.close()