import argparse
import os

import harmonica.adcirc_database
import harmonica.leprovost_database
import harmonica.tidal_constituents

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--adcirc_atlantic", default=None,
                        help="path to the ADCIRC Atlantic database")
    parser.add_argument("-p", "--adcirc_pacific", default=None,
                        help="path to the ADCIRC Pacific database")
    parser.add_argument("-l", "--leprovost", default=None,
                        help="path to the LeProvost database folder")
    args = vars(parser.parse_args())

    adcirc_atlantic = args["adcirc_atlantic"]
    adcir_pacific = args["adcirc_pacific"]
    leprovost = args["leprovost"]

    # Need to be in (lat, lon), not (x, y)
    # Invert commented list declarations to test [-180, 180] vs. [0, 360] ranges.
    #atlantic = [(39.74, 285.93),  # Not valid with ADCIRC Pacific database
    #            (42.32, 288.6),
    #            (45.44, 290.62)]
    #pacific = [(43.63, 235.45),  # Not valid with ADCIRC Atlantic database
    #           (46.18, 235.62)]
    atlantic = [(39.74, -74.07),  # Not valid with ADCIRC Pacific database
                (42.32, -71.4),
                (45.44, -69.38)]
    pacific = [(43.63, -124.55),  # Not valid with ADCIRC Atlantic database
               (46.18, -124.38)]
    all_points = []  # All locations valid with LeProvost and TPXO
    all_points.extend(atlantic)
    all_points.extend(pacific)

    good_cons = ['M2', 'S2', 'N2', 'K1']
    # Create an ADCIRC database for the Atlantic locations (default).
    ad_alantic_db = harmonica.adcirc_database.AdcircDB(adcirc_atlantic)
    # Create an ADCIRC database for the Pacific locations.
    ad_pacific_db = harmonica.adcirc_database.AdcircDB(adcir_pacific, "adcircnepac")
    # Create a LeProvost database for all locations.
    leprovost_db = harmonica.leprovost_database.LeProvostDB(leprovost)
    # Create a TPXO database for all locations.
    tpxo_db = harmonica.tidal_constituents.Constituents('tpxo8')

    # Get nodal factor data from the ADCIRC and LeProvost tidal databases
    ad_al_nodal_factor = ad_alantic_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)
    ad_pa_nodal_factor = ad_pacific_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)
    leprovost_nodal_factor = leprovost_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)

    f = open(os.path.join(os.getcwd(), "tidal_test.out"), "w")
    f.write("ADCIRC Atlantic nodal factor:\n")
    f.write(ad_al_nodal_factor.to_string() + "\n\n")
    f.write("ADCIRC Pacific nodal factor:\n")
    f.write(ad_pa_nodal_factor.to_string() + "\n\n")
    f.write("LeProvost nodal factor:\n")
    f.write(leprovost_nodal_factor.to_string() + "\n\n")
    f.flush()

    # Get tidal harmonic components for a list of points using the ADCIRC,
    # LeProvost, and TPXO databases.
    ad_atlantic_comps = ad_alantic_db.get_components(atlantic, good_cons)
    f.write("ADCIRC Atlantic components:\n")
    for pt in ad_atlantic_comps.data:
        f.write(pt.sort_index().to_string() + "\n\n")

    f.write("ADCIRC Pacific components:\n")
    f.flush()
    ad_pacific_comps = ad_pacific_db.get_components(pacific, good_cons)
    for pt in ad_pacific_comps.data:
        f.write(pt.sort_index().to_string() + "\n\n")

    f.write("LeProvost components:\n")
    f.flush()
    leprovost_comps = leprovost_db.get_components(all_points, good_cons)
    for pt in leprovost_comps.data:
        f.write(pt.sort_index().to_string() + "\n\n")

    f.write("TPX0 components:\n")
    f.flush()
    tpxo_comps = tpxo_db.get_components(all_points, good_cons, True)
    for pt in tpxo_comps.data:
        f.write(pt.sort_index().to_string() + "\n\n")

    f.close()