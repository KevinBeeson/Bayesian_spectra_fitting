import os
import matplotlib
if os.environ.get('DISPLAY') is None:
    # enables figure saving on clusters with ssh conection
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.io import fits
from glob import glob
from scipy.interpolate import splrep, splev
from scipy.signal import savgol_filter, argrelextrema

GALAH_BLUE  = [4718, 4903]
GALAH_GREEN = [5649, 5873]
GALAH_RED   = [6481, 6739]
GALAH_IR    = [7590, 7890]


def move_to_dir(path):
    if not(os.path.isdir(path)):
        os.mkdir(path)
    os.chdir(path)


def plot_tsne_results(axis1, axis2, colour_table, table_rows, suffix='', prefix='', ps=0.5):
    for col in table_rows:
        plot_name = prefix+'plot_' + col + suffix + '.png'
        if not(os.path.isfile(plot_name)):
            plt.close()
            plt.scatter(axis1, axis2, marker='o', lw=0, cmap='jet', c=colour_table[col], s=ps)  # ,vmin=-2, vmax=2)
            plt.title(col)
            plt.colorbar()
            plt.tight_layout()
            plt.savefig(plot_name, dpi=500)
            plt.close()


def get_spectra_dr52(object, bands=[1,2,3,4], root='', read_sigma=False, remove_nan=False,
                     extension=4, individual=False):
    if individual:
        subfolder = 'all'
    else:
        subfolder = 'com'
    fits_path = root + object[0:6] + '/standard/'+subfolder+'/' + object
    fits_path_2 = root + object[0:6] + '/standard/'+subfolder+'/' + object
    # determine path of spectra
    if not(os.path.isfile(fits_path + '1.fits')):
        fits_path = fits_path_2
        if not(os.path.isfile(fits_path + '1.fits')):
            print ('Spectra file not found')
            if read_sigma:
                return np.array([]), np.array([]), np.array([])
            else:
                return np.array([]), np.array([])
    # read selected spectrum bands
    spect_all = list([])
    if read_sigma:
        sigma_all = list([])
    wvls_all = list([])

    for i_b in range(len(bands)):
        band = bands[i_b]
        fits_path_band = fits_path + str(band) + '.fits'
        fits_data = fits.open(fits_path_band, memmap=False)
        if len(fits_data) < (extension+1):
            return list([]), list([])
        data_len = len(fits_data[extension].data)
        spect_all.append(fits_data[extension].data)
        if read_sigma:
            noise_data_len = len(fits_data[1].data)
            sigma_all.append(fits_data[1].data[noise_data_len-data_len:noise_data_len])
        header = fits_data[extension].header
        # calculate wavelengths of observed spectra
        wvls_all.append(header.get('CRVAL1') + np.float64(range(0, data_len)) * header.get('CDELT1'))
        # print header.get('CDELT1')
        fits_data.close()
    if read_sigma:
        return spect_all, wvls_all, sigma_all
    else:
        return spect_all, wvls_all


def plot_spectrum(spec, wvl, out_path='spectra.png', title=None, emis_lines=None, emis_lines_names=False, grad=False):
    data_shape = spec.shape
    n_plots = data_shape[0]
    if not(os.path.isfile(out_path)):
        print ('Plotting spectrum data.')
        plt.close()
        fig = plt.figure(1, figsize=(12, 3*n_plots))
        for id_plot in np.arange(n_plots)+1:
            fig_sub = plt.subplot(n_plots, 1, id_plot)  # rows, cols, id
            if id_plot == 1 and title is not None:
                fig_sub.set_title(title)
            if emis_lines is not None:
                for line in emis_lines:
                    plt.axvline(x=line, color='black', linewidth=0.6)
            plt.plot(wvl[id_plot-1], spec[id_plot-1], linewidth=.8, color='red')
            if grad:
                plt.ylabel('Gradient')
                plt.ylim([-0.5, 0.5])
            else:
                plt.ylabel('Flux')
                plt.ylim([0.5, 1.1])
            # determine spectral range
            wvl_mean = np.nanmean(wvl[id_plot-1])  #use only valid data for the band detemination
            if wvl_mean < 5000:
                plt.xlim([4710, 4755])
                # plt.xlim([4710, 4900])
            if wvl_mean > 5000 and wvl_mean < 6000:
                # plt.xlim([5650, 5880])
                plt.xlim([5815, 5835])
            if wvl_mean > 6000 and wvl_mean < 7000:
                # plt.xlim([6480, 6740])
                plt.xlim([6570, 6590])
            if wvl_mean > 7000:
                plt.xlim([7590, 7890])
        plt.xlabel('Wavelength')
        plt.tight_layout()
        plt.savefig(out_path, dpi=400)
        plt.close()


def spectra_resample(spectra, wvl_orig, wvl_target, k=3):
    """
    :param spectra:
    :param wvl_orig:
    :param wvl_target:
    :param k:
    :return:
    """
    idx_finite = np.isfinite(spectra)
    min_wvl_s = np.nanmin(wvl_orig[idx_finite])
    max_wvl_s = np.nanmax(wvl_orig[idx_finite])
    bspline = splrep(wvl_orig[idx_finite], spectra[idx_finite], k=k)
    idx_target = np.logical_and(wvl_target >= min_wvl_s,
                                wvl_target <= max_wvl_s)
    new_flux = splev(wvl_target[idx_target], bspline)
    nex_flux_out = np.ndarray(len(wvl_target))
    nex_flux_out.fill(np.nan)
    nex_flux_out[idx_target] = new_flux
    return nex_flux_out


def _evaluate_norm_fit(orig, fit, idx, sigma_low, sigma_high):
    # diffence to the original data
    diff = orig - fit
    std_diff = np.nanstd(diff[idx])
    # select data that will be fitted
    idx_outlier = np.logical_or(diff < (-1. * std_diff * sigma_low),
                                diff > (std_diff * sigma_high))
    return np.logical_and(idx, ~idx_outlier)


def spectra_normalize(wvl, spectra_orig, steps=5, sigma_low=2., sigma_high=2.5, window=15, order=5, n_min_perc=5.,
                      return_fit=False, return_idx=False, sg_filter=False, func='cheb', fit_on_idx=None, fit_mask=None):
    # perform sigma clipping before the next fitting cycle
    idx_fit = np.logical_and(np.isfinite(wvl), np.isfinite(spectra_orig))
    spectra = np.array(spectra_orig)

    if fit_mask is not None:
        idx_fit = np.logical_and(idx_fit, fit_mask)

    if fit_on_idx is not None:
        idx_fit = np.logical_and(idx_fit, fit_on_idx)
        steps = 1  # no clipping performed, one iteration, forced fitting on selected pixels
    else:
        # filter noisy original spectra, so it is easier to determine continuum levels
        if sg_filter:
            spectra = savgol_filter(spectra_orig, window_length=15, polyorder=5)
        init_fit = np.nanmedian(spectra)
        idx_fit = _evaluate_norm_fit(spectra, init_fit, idx_fit, sigma_low*2.5, sigma_high*2.5)
    data_len = np.sum(idx_fit)
    n_fit_points_prev = np.sum(idx_fit)
    for i_f in range(steps):  # number of sigma clipping steps
        # print i_f
        if func == 'cheb':
            chb_coef = np.polynomial.chebyshev.chebfit(wvl[idx_fit], spectra[idx_fit], order)
            cont_fit = np.polynomial.chebyshev.chebval(wvl, chb_coef)
        if func == 'legen':
            leg_coef = np.polynomial.legendre.legfit(wvl[idx_fit], spectra[idx_fit], order)
            cont_fit = np.polynomial.legendre.legval(wvl, leg_coef)
        if func == 'poly':
            poly_coef = np.polyfit(wvl[idx_fit], spectra[idx_fit], order)
            cont_fit = np.poly1d(poly_coef)(wvl)
        if func == 'spline':
            # if i_f == 1:
            #     chb_coef = np.polynomial.chebyshev.chebfit(wvl[idx_fit], spectra[idx_fit], 5)
            #     cont_fit = np.polynomial.chebyshev.chebval(wvl, chb_coef)
            #     idx_fit = _evaluate_norm_fit(spectra, cont_fit, idx_fit, sigma_low, sigma_high)
            spline_coef = splrep(wvl[idx_fit], spectra[idx_fit], k=order, s=window)
            cont_fit = splev(wvl, spline_coef)
            print (i_f, 'points:', n_fit_points_prev, 'knots:', len(spline_coef[0]))
        idx_fit = _evaluate_norm_fit(spectra, cont_fit, idx_fit, sigma_low, sigma_high)
        n_fit_points = np.sum(idx_fit)
        if 100.*n_fit_points/data_len < n_min_perc:
            break
        if n_fit_points == n_fit_points_prev:
            break
        else:
            n_fit_points_prev = n_fit_points
    if return_fit:
        if return_idx:
            return cont_fit, idx_fit
        else:
            return cont_fit
    else:
        return spectra_orig / cont_fit


def spectra_logspace(flx, wvl):
    """
    :param flx:
    :param wvl:
    :return:
    """
    wvl_new = np.logspace(np.log10(wvl[0]), np.log10(wvl[-1]), num=len(wvl))
    return np.interp(wvl_new, wvl, flx), wvl_new


def spectra_linspace(flx, wvl):
    """
    :param flx:
    :param wvl:
    :return:
    """
    wvl_new = np.linspace(wvl[0], wvl[-1], num=len(wvl))
    return np.interp(wvl_new, wvl, flx), wvl_new


def spectra_rvshift(flx, wvl, rv, linerize=True):
    """
    :param flx:
    :param wvl:
    :return:
    """
    wvl_new = wvl * (1. + rv / 299792.458)
    if linerize:
        flx_new = np.interp(np.linspace(wvl_new[0], wvl_new[-1], num=len(wvl_new)), wvl_new, flx)
        wvl_new = np.linspace(wvl_new[0], wvl_new[-1], num=len(wvl_new))
        return flx_new, wvl_new
    else:
        return flx, wvl_new