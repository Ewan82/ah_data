import numpy as np
import matplotlib.mlab as mlab


def convert_csv2rec(file_no):
    return mlab.csv2rec('../litter_traps/litterscans/file0'+str(file_no)+'.csv')


def remove_false_data(area_arr, tol=2.0):
    idx = np.where(area_arr < tol)
    return np.delete(area_arr, idx)