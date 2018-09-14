import argparse
import os

import harmonica.adcirc_database
import harmonica.leprovost_database
#import TPXODatabase

def write_nodal(f, nf):
    for n in nf:
        f.write(n.name + ", " +
        str(n.amplitude) + ", " +
        str(n.frequency) + ", " +
        str(n.earth_tide_reduction_factor) + ", " +
        str(n.equilibrium_argument) + ", " +
        str(n.nodal_factor) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--work_dir", default=os.getcwd(),
                        help="directory to save output and temporary files")
    parser.add_argument("-a", "--adcirc_atlantic", default="",
                        help="path to the ADCIRC Atlantic database")
    parser.add_argument("-p", "--adcirc_pacific", default="",
                        help="path to the ADCIRC Pacific database")
    parser.add_argument("-l", "--leprovost",
                        default=os.getcwd(), help="path to the LeProvost database folder")
    args = vars(parser.parse_args())
    
    work_dir = args["work_dir"]
    adcirc_atlantic = args["adcirc_atlantic"]
    adcir_pacific = args["adcirc_pacific"]
    leprovost = args["leprovost"]
    
    atlantic = [(-74.07, 39.74),  # Not valid with ADCIRC Pacific database
                (-71.4, 42.32),
                (-69.38, 45.44)]
    pacific = [(-124.55, 43.63),  # Not valid with ADCIRC Atlantic database
               (-124.38, 46.18)]
    all_points = []  # All locations valid with LeProvost
    all_points.extend(atlantic)
    all_points.extend(pacific)

    good_cons = ['M2', 'S2', 'N2', 'K1']
    # bad_con = 'ZZ7'
    ad_alantic_db = harmonica.adcirc_database.AdcircDB(work_dir, adcirc_atlantic,
                                                       harmonica.adcirc_database.TidalDBAdcircEnum.TIDE_NWAT)
    ad_pacific_db = harmonica.adcirc_database.AdcircDB(work_dir, adcir_pacific,
                                                       harmonica.adcirc_database.TidalDBAdcircEnum.TIDE_NEPAC)
    le_db = harmonica.leprovost_database.LeProvostDB(leprovost)
    # tp_db = TPXODatabase.TpxoDB('tpxo8')

    ad_al_nf = ad_alantic_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)
    ad_pa_nf = ad_pacific_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)
    le_nf = le_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)
    # tp_nf = tp_db.get_nodal_factor(good_cons, 15, 30, 8, 2018)

    f = open(os.path.join(work_dir, "tidal_test.out"), "w")
    f.write("ADCIRC Atlantic nodal factor:\n")
    f.write(ad_al_nf.to_string() + "\n\n")
    f.write("ADCIRC Pacific nodal factor:\n")
    f.write(ad_pa_nf.to_string() + "\n\n")
    f.write("LeProvost nodal factor:\n")
    f.write(le_nf.to_string() + "\n\n")
    f.flush()

    ap1 = ad_alantic_db.get_components(atlantic, good_cons)
    ap2 = ad_pacific_db.get_components(pacific, good_cons)
    ap3 = le_db.get_components(all_points, good_cons)
    #ap4 = tp_db.get_amplitude_and_phase(good_cons, all_points)

    f.write("ADCIRC Atlantic components:\n")
    for pt in ap1:
        f.write(pt.to_string() + "\n\n")
    f.write("ADCIRC Pacific components:\n")
    for pt in ap2:
        f.write(pt.to_string() + "\n\n")
    f.write("LeProvost components:\n")
    for pt in ap3:
        f.write(pt.to_string() + "\n\n")
    # f.write(str(ap4) + "\n\n")

    f.close()