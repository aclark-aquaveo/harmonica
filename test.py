import argparse
import os

import harmonica.adcirc_database
import harmonica.leprovost_database
import harmonica.tidal_constituents
# import TPXODatabase

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--work_dir", default=os.getcwd(),
                        help="directory to save output and temporary files")
    parser.add_argument("-a", "--adcirc_atlantic", default="",
                        help="path to the ADCIRC Atlantic database")
    parser.add_argument("-p", "--adcirc_pacific", default="",
                        help="path to the ADCIRC Pacific database")
    parser.add_argument("-l", "--leprovost", default="",
                        help="path to the LeProvost database folder")
    args = vars(parser.parse_args())
    
    work_dir = args["work_dir"]
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
    all_points = []  # All locations valid with LeProvost
    all_points.extend(atlantic)
    all_points.extend(pacific)

    good_cons = ['M2', 'S2', 'N2', 'K1']
    # bad_con = 'ZZ7'
    ad_alantic_db = harmonica.adcirc_database.AdcircDB(work_dir, adcirc_atlantic,
                                                       harmonica.adcirc_database.TidalDBAdcircEnum.TIDE_NWAT)
    ad_pacific_db = harmonica.adcirc_database.AdcircDB(work_dir, adcir_pacific,
                                                       harmonica.adcirc_database.TidalDBAdcircEnum.TIDE_NEPAC)
    leprovost_db = harmonica.leprovost_database.LeProvostDB(leprovost)
    # tp_db = TPXODatabase.TpxoDB('tpxo8')
    tpx0_db = harmonica.tidal_constituents.Constituents()

    ad_al_nodal_factor = ad_alantic_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)
    ad_pa_nodal_factor = ad_pacific_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)
    leprovost_nodal_factor = leprovost_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)
    # tp_nf = tp_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)

    f = open(os.path.join(work_dir, "tidal_test.out"), "w")
    f.write("ADCIRC Atlantic nodal factor:\n")
    f.write(ad_al_nodal_factor.to_string() + "\n\n")
    f.write("ADCIRC Pacific nodal factor:\n")
    f.write(ad_pa_nodal_factor.to_string() + "\n\n")
    f.write("LeProvost nodal factor:\n")
    f.write(leprovost_nodal_factor.to_string() + "\n\n")
    f.flush()

    all_points1 = ad_alantic_db.get_components(atlantic, good_cons)
    all_points2 = ad_pacific_db.get_components(pacific, good_cons)
    all_points3 = leprovost_db.get_components(all_points, good_cons)
    #all_points4 = tp_db.get_amplitude_and_phase(good_cons, all_points)

    f.write("ADCIRC Atlantic components:\n")
    for pt in all_points1.data:
        f.write(pt.to_string() + "\n\n")
    f.write("ADCIRC Pacific components:\n")
    for pt in all_points2.data:
        f.write(pt.to_string() + "\n\n")
    f.write("LeProvost components:\n")
    for pt in all_points3.data:
        f.write(pt.to_string() + "\n\n")
    f.write("TPX0 components:\n")
    try:
        for pt in all_points:
            components = tpx0_db.get_components(pt, 'tpxo8', good_cons, True)
            f.write(components.data.to_string() + "\n\n")
            f.flush()
    except Exception as e:
        print("Exception thrown during TPX0 extraction: {}".format(e))

    f.close()