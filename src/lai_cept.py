import numpy as np
import netCDF4 as nC
import time as tt
import matplotlib.mlab as mlab
import datetime as dt
import random
import xlrd
import pvlib

# ----------------------------------------------------------------------------------
# NETCDF
# ----------------------------------------------------------------------------------

def create_netcdf_dataset():
    """Creates netcdf dataset for PAR measurements at Alice Holt
    """
    dataset = nC.Dataset('ah_par_test.nc', 'w', format='NETCDF4_CLASSIC')
    time = dataset.createDimension('time', None)
    lat = dataset.createDimension('lat', 435)
    lon = dataset.createDimension('lon', 435)
    times = dataset.createVariable('time', np.float64, ('time',))
    latitudes = dataset.createVariable('latitude', np.float64, ('lat',))
    longitudes = dataset.createVariable('longitude', np.float64, ('lon',))
    above_par = dataset.createVariable('above_par', np.float32, ('time',))
    below_par = dataset.createVariable('below_par', np.float32, ('time', 'lat', 'lon',))
    plot_label = dataset.createVariable('plot_label', np.int32, ('lat', 'lon'))

    dataset.description = 'Record of above and below canopy PAR measurements from the Alice Holt flux site'
    dataset.history = 'Created ' + tt.ctime(tt.time())
    dataset.source = 'Ewan Pinnington, University of Reading'

    latitudes.units = 'degrees north'
    longitudes.units = 'degrees east'
    above_par.units = 'umol m-2 s-1'
    below_par.units = 'umol m-2 s-1'
    plot_label.units = 'None'
    times.units = 'hours since 1900-01-01 00:00:00.0'
    times.calendar = 'gregorian'

    #points_rec = mlab.csv2rec('lonlat_threet.csv')
    #longitudes[:] = points_rec['longitude']
    #latitudes[:] = points_rec['latitude']
    dataset.close()


def add_times2netcdf(dl2e_file, nc_dataset):
    dl2e_arr = extract_par_dl2e(dl2e_file)
    len_n = (convert_dl2e_timestr(dl2e_arr['date'][-1]) - convert_dl2e_timestr(dl2e_arr['date'][0])).seconds
    dataset = nC.Dataset(nc_dataset)
    times = dataset.variables['time']
    dates = []
    for n in range(len_n):
        dates.append(convert_dl2e_timestr(dl2e_arr['date'][0]) + n*dt.timedelta(seconds=1))
    num_dates = nC.date2num(dates, units=times.units, calendar=times.calendar)
    return num_dates

    len_time = len(times)
    for nt in range(len_time, len_n):
        times[nt] = num_dates[nt - len_time]
    dataset.close()


def save_par2netcdf():
    return 'done'


def add_plot_points(nc_dataset):
    dataset = nC.Dataset(nc_dataset)
    plot_label = dataset.variables['plot_label']
    points_rec = mlab.csv2rec('lonlat_threet.csv')
    for x in xrange(435):
        plot_label[x,x] = points_rec['name'][x]
    dataset.close()


def extract_par_lp80(excel_file_loc):
    """
    Reads data from excel lp80 file
    :param excel_file_loc: excel file location
    :return: array of extracted data
    """
    book = xlrd.open_workbook(excel_file_loc)
    sheet = book.sheet_by_index(0)
    lp80_dat = np.array([sheet.col(1)[1:], sheet.col(2)[1:], sheet.col(3)[1:]]).T
    return lp80_dat


def save_lp80_netcdf(lp80_dat, nc_dataset):
    dataset = nC.Dataset(nc_dataset)
    return 'done'

# ----------------------------------------------------------------------------------
# PAR
# ----------------------------------------------------------------------------------


def edit_dl2e_time(time_str):
    """
    Converts time string from DL2e data logger to add year.
    :param time_str: DL2e data logger time string
    :return: updated time string
    """
    return '2015/'+time_str.strip('"')


def edit_dl2e_mv(mv_float):
    """
    Converts DL2e mV reading to PAR umol m-2 s-1
    :param mv_float: mV reading
    :return: PAR value
    """
    return 100.28*float(mv_float) - 0.18


def convert_dl2e_timestr(time_str):
    return dt.datetime.strptime(time_str, "%Y/%d/%m %H:%M:%S")


def all_dl2e_time(time_str):
    time_str_yr = edit_dl2e_time(time_str)
    return dt.datetime.strptime(time_str_yr, "%Y/%d/%m %H:%M:%S   ")


def extract_mv_dl2e(dl2e_file):
    return np.loadtxt(dl2e_file, delimiter=',',
                        skiprows=11, dtype={'names':('date', 'mV'), 'formats':('S19', 'f4')},
                        usecols=(0,3), converters={0: edit_dl2e_time})


def extract_par_dl2e(dl2e_file):
    return np.genfromtxt(dl2e_file, delimiter=',',
                        skip_header=11, names='date, par', dtype='object, f4',
                        usecols=(0,3), converters={0: all_dl2e_time, 3: edit_dl2e_mv})


def nearest_value(array, value):
    return min(array, key=lambda x: abs(x - value))


def idx(array, value):
    return np.where(array==value)[0][0]


def find_nearest_idx(array, value):
    return idx(array, nearest_value(array, value))


def extract_lai_fpar(above_par_dat, below_par_dat):
    above_par_ra = mlab.csv2rec(above_par_dat)
    below_par_ra = mlab.csv2rec(below_par_dat)
    points_ra = mlab.csv2rec('lonlat_threet.csv')
    plot = below_par_ra['plot']
    date = below_par_ra['date']
    below_par = below_par_ra['par']
    lats = np.array(points_ra['latitude'].tolist()*2)
    lons = np.array(points_ra['longitude'].tolist()*2)
    above_par = []
    fapar = []
    for time in enumerate(date):
        par_idx = find_nearest_idx(above_par_ra['date'], time[1])
        above_par.append(np.mean((above_par_ra['par'][par_idx-1], above_par_ra['par'][par_idx],
                                 above_par_ra['par'][par_idx+1])))
        if above_par_ra['par'][par_idx] < below_par[time[0]]:
            fapar.append(0)
        else:
            fapar.append((above_par_ra['par'][par_idx] - below_par[time[0]]) /
                     above_par_ra['par'][par_idx])
    above_par = np.array(above_par)
    fapar = np.array(fapar)
    newra = np.column_stack((date, plot, lats, lons, above_par, below_par, fapar))
    new_ra = np.core.records.fromarrays(newra.transpose(),
                                        dtype=[('date', 'object'),
                                               ('plot', 'i'), ('lat', 'f'),
                                               ('lon', 'f'), ('above_par', 'f'),
                                               ('below_par', 'f'), ('fapar', 'f')])
    return new_ra


def create_lai_arr(par_csv):
    par_ra = mlab.csv2rec(par_csv)
    above_par = par_ra['above_par']
    below_par = par_ra['below_par']
    lats = par_ra['lat']
    lons = par_ra['lon']
    times = par_ra['date']
    plots = par_ra['plot']
    lai = []
    zenith = []
    for time in enumerate(times):
        index = time[0]
        zenith.append(find_solar_zenith(lons[index], lats[index], time[1]))
        try:
            lai.append(calc_lai(lons[index], lats[index], time[1], above_par[index], below_par[index]))
        except ValueError:
            lai.append(float('NaN'))
            print index

    newra = np.column_stack((times, plots, lats, lons, above_par, below_par, zenith, lai))
    new_ra = np.core.records.fromarrays(newra.transpose(),
                                    dtype=[('date', 'object'),
                                           ('plot', 'i'), ('lat', 'f'),
                                           ('lon', 'f'), ('above_par', 'f'),
                                           ('below_par', 'f'), ('solar_zenith', 'f'),
                                           ('lai', 'f')])
    return new_ra


def save_ra2csv(rec_arr, fname):
    mlab.rec2csv(rec_arr, fname)

# ----------------------------------------------------------------------------------
# PCQM plots
# ----------------------------------------------------------------------------------


def draw_rand_plots(no_plots):
    plots = random.sample(np.arange(1,436,2), no_plots)
    plots.sort()
    points_ra = mlab.csv2rec('lonlat_threet.csv')
    lats = []
    lons = []
    mens = []

    for point in enumerate(plots):
        plot_idx = find_nearest_idx(points_ra['name'], point[1])
        lats.append(points_ra['latitude'][plot_idx])
        lons.append(points_ra['longitude'][plot_idx])
        mens.append(points_ra['m_plot'][plot_idx])

    plots = np.array(plots)
    lats = np.array(lats)
    lons = np.array(lons)
    mens = np.array(mens)

    newra = np.column_stack((plots, lats, lons, mens))
    new_ra = np.core.records.fromarrays(newra.transpose(),
                                        dtype=[('plot', 'i'),
                                               ('lat', 'f'), ('lon', 'f'),
                                               ('m_plot', 'i')])
    return new_ra

def draw_plots_inc_mens():
    all_plots = np.arange(1,436).tolist()
    men_hemi_lst = ah_mens_plots+ah_other_plots
    men_hemi_lst.sort()

    for x in men_hemi_lst:
        if x in all_plots:
            all_plots.remove(x)
    for x in men_hemi_lst:
        if x+1 in all_plots:
                all_plots.remove(x+1)
    for x in men_hemi_lst:
        if x-1 in all_plots:
            all_plots.remove(x-1)

    rm_lst = []
    for x in xrange(len(all_plots)-1):
        if all_plots[x+1]-all_plots[x] == 1:
            rm_lst.append(all_plots[x])
    for x in rm_lst:
        if x in all_plots:
            all_plots.remove(x)

    return all_plots

def twenty_m_points():
    temp_plot_lst = draw_plots_inc_mens()
    temp_plots = random.sample(np.array(temp_plot_lst), 25)
    plots = temp_plots + ah_mens_plots + ah_other_plots
    plots.sort()
    points_ra = mlab.csv2rec('lonlat_threet.csv')
    lats = []
    lons = []
    mens = []

    for point in enumerate(plots):
        plot_idx = find_nearest_idx(points_ra['name'], point[1])
        lats.append(points_ra['latitude'][plot_idx])
        lons.append(points_ra['longitude'][plot_idx])
        mens.append(points_ra['m_plot'][plot_idx])

    plots = np.array(plots)
    lats = np.array(lats)
    lons = np.array(lons)
    mens = np.array(mens)

    newra = np.column_stack((plots, lats, lons, mens))
    new_ra = np.core.records.fromarrays(newra.transpose(),
                                        dtype=[('plot', 'i'),
                                               ('lat', 'f'), ('lon', 'f'),
                                               ('m_plot', 'i')])
    return new_ra


ah_mens_plots = [1, 23, 43, 63, 88, 108, 116, 124, 162, 185, 186, 207, 236, 281,
                 311, 324, 325, 339, 363, 384, 404, 420, 435]

ah_other_plots = [5, 10, 15, 28, 33, 38, 48, 53, 59, 68, 73, 79, 84, 93, 98, 103,
                  113, 121, 129, 134, 139, 144, 149, 154, 159, 167, 172, 177, 191,
                  196, 201, 212, 217, 222, 227, 232, 241, 246, 251, 256, 261, 266,
                  271, 276, 286, 291, 296, 301, 306, 316, 330, 336, 345, 349, 354,
                  359, 368, 373, 378, 389, 394, 400, 409, 414, 425, 430]


# ----------------------------------------------------------------------------------
# LP-80 LAI calculations
# ----------------------------------------------------------------------------------


def find_solar_zenith(lon, lat, time):
    """ Calculates solar zenith angle
    :param lon: longitude value
    :param lat: latitude value
    :param time: time at which angle to be calculated
    :return: solar zenith angle
    """
    return pvlib.solarposition.get_solarposition(time, lat, lon)['zenith'].item()


def calc_r(zen_ang, above_par):
    """ Calculates r value for LP80 decagon ceptometer
    :param above_par: above canopy PAR reading
    :return: r value for given zenith angle
    """
    if (zen_ang*np.pi / 180.) > 1.5:
        raise ValueError('Nighttime measurement invalid')
    r = above_par / (2550.*np.cos(zen_ang*np.pi/180.))
    if r > 0.82:
        r = 0.82
    elif r < 0.2:
        r = 0.2
    else:
        r = r
    return r


def calc_fb(r):
    """ Calculates beam fraction for LP80 ceptometer
    :param r: r values from calc_r fn.
    :return: beam fraction value
    """
    fb = 1.395 + r*(-14.43 + r*(48.57 + r*(-59.024 + r*24.835)))
    return fb


def calc_extinc_coef_decagon(chi, zen_ang):
    """ Calculates extinction coefficient for canopy (K), using decagon formula
    :param chi: leaf angle distribution
    :param zen_ang: solar zenith angle
    :return: K
    """
    k = np.sqrt(chi**2 + np.tan(zen_ang*np.pi/180.)**2) / (chi + 1.744*(chi + 1.182)**-0.733)
    return k


def calc_extinc_coef_campbell(chi, zen_ang):
    """ Calculates extinction coefficient for canopy (K), using campbell formula
    :param chi: leaf angle distribution
    :param zen_ang: solar zenith angle
    :return: K
    """
    k = np.sqrt(chi**2 + 1/np.tan(zen_ang*np.pi/180.)**2) / \
        (1.47 + 0.45*chi + 0.1223*chi**2 - 0.013*chi**3 + 0.000509*chi**4)
    return k


def calc_leaf_absorptivity_term(a=0.9):
    """ Calculates leaf absorptivity term in LP80 LAI calculation
    :param a: leaf absorptivity (assumed to be 0.9 by decagon)
    :return: A, LP80 leaf absorptivity term
    """
    A = 0.283 + 0.785*a - 0.159*a**2
    return A


def calc_tau(above_par, below_par):
    """ Calculates tau, the ratio of PAR measured below the canopy to PAR above the canopy
    :param above_par: above canopy PAR reading
    :param below_par: below canopy PAR reading
    :return: tau
    """
    tau = below_par / above_par
    if tau >= 1.:
        raise ValueError('below par greater than above par, cannot be!!')
    return tau


def calc_lai(lon, lat, time, above_par, below_par):
    """ Calculates LAI given appropriate data
    :param lon: longitude of sampled point
    :param lat: latitude of sampled point
    :param time: time of measurement as a datetime object
    :param above_par: corresponding above canopy PAR measurement
    :param below_par: below canopy PAR measurement from LP80
    :return: LAI (m^2/m^2)
    """
    zenith = find_solar_zenith(lon, lat, time)
    r = calc_r(zenith, above_par)
    fb = calc_fb(r)
    # fb = 0.
    k = calc_extinc_coef_campbell(1., zenith)
    # k = 1./2.*np.cos(zenith*np.pi/180.)
    A = calc_leaf_absorptivity_term()
    tau = calc_tau(above_par, below_par)
    lai = ((1. - 1./(2.*k))*fb - 1.)*np.log(tau) / A*(1 - 0.47*fb)
    return lai


def calc_lai_test(lon, lat, time, above_par, below_par):
    zenith = find_solar_zenith(lon, lat, time)
    r = calc_r(zenith, above_par)
    #fb = calc_fb(r)
    fb = 0.4
    #k = calc_extinc_coef_campbell(1., zenith)
    k = 1./2.*np.cos(zenith)
    A = calc_leaf_absorptivity_term()
    tau = calc_tau(above_par, below_par)
    lai = ((1. - 1./(2.*k))*fb - 1.)*np.log(tau) / A*(1 - 0.47*fb)
    return r, fb, k, A, tau, np.log(tau), lai
