import numpy as np
import matplotlib.mlab as ml


# Transect points
thinned07 = np.concatenate((np.arange(1, 88), np.arange(267, 376)))
thinned14 = np.concatenate((np.arange(88, 267), np.arange(417, 424)))
ecn = np.concatenate((np.arange(376, 417), np.arange(424, 436)))
traps = np.array([10, 144, 281, 227, 325, 373])
ecn_07 = np.concatenate((thinned07, ecn))

# Files
lp80 = 'laik_mean.csv'
hemi = '../hemispherical_photos/processed_hemis/LAI.csv'
trap = 'trap_lai.csv'


def find_lai_lp80(classification, f_name):
    """ Find mean and std of LAI for different regions of the Strait Inclose, Alice Holt
    :param classification: An np array of plot number corresponding to an area of forest
    :param f_name: location of file to extract LAI's from
    :return: array of LAI's, array of plot labels, LAI mean, LAI std
    """
    lai_arr = ml.csv2rec(f_name)
    if classification == 'all':
        lai = lai_arr['lai']
        plot = lai_arr['plot']
    else:
        lai = lai_arr['lai'][np.where(np.in1d(lai_arr['plot'], classification))[0]]
        plot = lai_arr['plot'][np.where(np.in1d(lai_arr['plot'], classification))[0]]
    lai_mean = np.mean(lai)
    lai_std = np.std(lai) / np.sqrt(len(plot))
    return lai, plot, lai_mean, lai_std



