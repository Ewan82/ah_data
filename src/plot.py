import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import seaborn as sns
import lai_da as ld


def plot_lai_comparison(transect):
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    palette = sns.color_palette("colorblind", 11)

    south_t = np.arange(1, 186)
    mid_t = np.arange(186, 325)
    north_t = np.arange(325, 436)
    lp80_lai = ml.csv2rec('laik_mean.csv')
    hemi_lai = ml.csv2rec('../hemispherical_photos/processed_hemis/LAI.csv')
    trap_lai = ml.csv2rec('../littertrap_loc.csv')
    if transect == 'south':
        ax.plot(south_t, lp80_lai['lai'][lp80_lai['plot'] < 186], color=palette[0], label='Ceptometer')
        ax.plot(hemi_lai['plot'][hemi_lai['plot'] < 186], hemi_lai['lai'][hemi_lai['plot'] < 186], 'o',
                color=palette[1], label='Hemi. photos')
        ax.plot(trap_lai['plot'][trap_lai['plot'] < 186], trap_lai['lai'][trap_lai['plot'] < 186], 'o',
                color=palette[2], label='Litter traps')
    elif transect == 'mid':
        ax.plot(mid_t, lp80_lai['lai'][(185 < lp80_lai['plot']) & (lp80_lai['plot'] < 325)], color=palette[0]
                , label='Ceptometer')
        ax.plot(hemi_lai['plot'][[(185 < hemi_lai['plot']) & (hemi_lai['plot'] < 325)]],
                hemi_lai['lai'][[(185 < hemi_lai['plot']) & (hemi_lai['plot'] < 325)]], 'o', color=palette[1],
                label='Hemi. photos')
        ax.plot(trap_lai['plot'][[(185 < trap_lai['plot']) & (trap_lai['plot'] < 325)]],
                trap_lai['lai'][[(185 < trap_lai['plot']) & (trap_lai['plot'] < 325)]], 'o',color=palette[2],
                label='Litter traps')
        plt.xlim(mid_t[-1], mid_t[0])
    if transect == 'north':
        ax.plot(north_t, lp80_lai['lai'][lp80_lai['plot'] > 324], color=palette[0], label='Ceptometer')
        ax.plot(hemi_lai['plot'][hemi_lai['plot'] > 324], hemi_lai['lai'][hemi_lai['plot'] > 324], 'o',
                color=palette[1], label='Hemi. photos')
        ax.plot(trap_lai['plot'][trap_lai['plot'] > 324], trap_lai['lai'][trap_lai['plot'] > 324], 'o',
                color=palette[2], label='Litter traps')
    ax.set_xlabel('Plot label')
    ax.set_ylabel('LAI (m2/ m2)')
    plt.legend()
    return ax, fig


def plot_lai_comparison2():
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 8))
    palette = sns.color_palette("colorblind", 11)

    south_t = np.arange(1, 186)
    mid_t = np.arange(186, 325)
    north_t = np.arange(325, 436)
    lp80_lai = ml.csv2rec('laik_mean.csv')
    hemi_lai = ml.csv2rec('../hemispherical_photos/processed_hemis/LAI.csv')
    trap_lai = ml.csv2rec('../littertrap_loc.csv')
    ax.plot(south_t, lp80_lai['lai'][lp80_lai['plot'] < 186], color=palette[0], label='Ceptometer')
    ax.plot(hemi_lai['plot'][hemi_lai['plot'] < 186], hemi_lai['lai'][hemi_lai['plot'] < 186], 'o',
            color=palette[1], label='Hemi. photos')
    ax.plot(trap_lai['plot'][trap_lai['plot'] < 186], trap_lai['lai'][trap_lai['plot'] < 186], 'o',
            color=palette[2], label='Litter traps')
    ax.axvline(186, color='k', linestyle='--')
    ax.text(5, 0.25, 'South transect')

    ax.plot(mid_t, lp80_lai['lai'][(185 < lp80_lai['plot']) & (lp80_lai['plot'] < 325)], color=palette[0])
    ax.plot(hemi_lai['plot'][[(185 < hemi_lai['plot']) & (hemi_lai['plot'] < 325)]],
            hemi_lai['lai'][[(185 < hemi_lai['plot']) & (hemi_lai['plot'] < 325)]], 'o', color=palette[1])
    ax.plot(trap_lai['plot'][[(185 < trap_lai['plot']) & (trap_lai['plot'] < 325)]],
            trap_lai['lai'][[(185 < trap_lai['plot']) & (trap_lai['plot'] < 325)]], 'o',color=palette[2])
    ax.axvline(325, color='k', linestyle='--')
    ax.text(190, 0.25, 'Middle transect')

    ax.plot(north_t, lp80_lai['lai'][lp80_lai['plot'] > 324], color=palette[0])
    ax.plot(hemi_lai['plot'][hemi_lai['plot'] > 324], hemi_lai['lai'][hemi_lai['plot'] > 324], 'o',
            color=palette[1])
    ax.plot(trap_lai['plot'][trap_lai['plot'] > 324], trap_lai['lai'][trap_lai['plot'] > 324], 'o',
            color=palette[2])
    ax.text(329, 0.25, 'North transect')
    ax.set_xlabel('Plot label')
    ax.set_ylabel('LAI (m2/ m2)')
    plt.legend(loc='upper left')
    plt.xlim(south_t[0], north_t[-1])
    return ax, fig


def plot_lai_comparison_thin(classification):
    sns.set_context(rc={'lines.linewidth': 0.8, 'lines.markersize': 6})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 8))
    palette = sns.color_palette("colorblind", 11)

    lai_c, plot_c, lai_m_c, lai_s_c = ld.find_lai_lp80(classification, ld.lp80)
    lai_h, plot_h, lai_m_h, lai_s_h = ld.find_lai_lp80(classification, ld.hemi)
    lai_t, plot_t, lai_m_t, lai_s_t = ld.find_lai_lp80(classification, ld.trap)

    lai_h_fill = np.ones(len(lai_c))*float('NaN')
    lai_h_fill[plot_c.searchsorted(plot_h)] = lai_h
    lai_t_fill = np.ones(len(lai_c))*float('NaN')
    lai_t_fill[plot_c.searchsorted(plot_t)] = lai_t

    xlist = np.arange(len(lai_c))
    ax.plot(xlist, lai_c, color=palette[0], label='Ceptometer')
    ax.plot(xlist, lai_h_fill, 'o', color=palette[1], label='Hemi. photos')
    ax.plot(xlist, lai_t_fill, 'o', color=palette[2], label='Litter traps')
    ax.axhline(lai_m_c, color=palette[0], linestyle='--')
    ax.axhline(lai_m_h, color=palette[1], linestyle='--')
    ax.axhline(lai_m_t, color=palette[2], linestyle='--')

    ax.set_xlabel('Number of LAI observations')
    ax.set_ylabel('LAI (m2/ m2)')
    plt.legend(loc='upper left')
    plt.ylim(0, 8)
    plt.xlim(0, xlist[-1]+1)
    return ax, fig