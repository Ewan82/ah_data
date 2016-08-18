import numpy as np
import matplotlib.mlab as mlab


def extract_pcqm(pcqm_csv):
    """
    Reads data from csv pcqm file
    :param csv_file_loc: csv file location
    :return: array of extracted data
    """
    pcqm_dat = mlab.csv2rec(pcqm_csv)
    return pcqm_dat


def update_plot_col(pcqm_dat):
    """ Re-werites plot column for pcqm csv file
    :param pcqm_dat: uneditted pcqm csv file
    :return: updated rec array
    """
    plots = pcqm_dat['plot']
    index = plots[0]
    for x in np.arange(1, len(plots)):
        if plots[x] == -1:
            plots[x] = index
        else:
            index = plots[x]
    return pcqm_dat


def calc_mean_point2plant(sum_p2p, no_pcq):
    """ Calculates mean point to plant distance
    :param sum_p2p: sum of pcqm distances
    :param no_pcq: number of point centered quarters
    :return: mean point to plant distance
    """
    return sum_p2p / no_pcq


def calc_density(meanp2p):
    """ Calculates density of trees per m2
    :param meanp2p: mean point to plant distance
    :return: density of trees per m2
    """
    return 1 / meanp2p**2


def find_density_mass(pcqm_arr, classification='all'):
    """ Calculates tree density and mass per area for pcqm array, given
    a forest classification area.
    :param pcqm_arr: point centered quarter array
    :param classification: section to calculate mass for ('2014 thinned', '2007 thinned' or 'ECN')
    :return: tree density (m-2), tree mass per area (m2)
    """
    if classification == 'all':
        sum_p2p = np.sum(pcqm_arr['distance_m'])
        no_pcq = len(pcqm_arr['plot'])
        dbhs = pcqm_arr['dbh_cm']
        meanp2p = calc_mean_point2plant(sum_p2p, no_pcq)
        density = calc_density(meanp2p)
        calibration = 4.35  # FIND FROM IAN
        mass_area = density*(np.pi*(np.mean(dbhs)/2*100.)**2)*calibration
        ret_val = density, mass_area
    else:
        index = np.where(pcqm_arr['classification'] == classification)[0]
        distances = [pcqm_arr['distance_m'][x] for x in index]
        dbhs = [pcqm_arr['dbh_cm'][x] for x in index]
        sum_p2p = np.sum(distances)
        no_pcq = len(distances)
        meanp2p = calc_mean_point2plant(sum_p2p, no_pcq)
        density = calc_density(meanp2p)
        calibration = 4.35  # FIND FROM IAN
        mass_area = density*(np.pi*(np.mean(dbhs)/2*100.)**2)*calibration
        ret_val = density, mass_area
    return ret_val
