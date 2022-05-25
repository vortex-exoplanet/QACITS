from qacits.util.bin_images import bin_images
from qacits.util.psf_flux import get_psf_flux, get_all_di
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


def calibrate_qacits(psf_ON, psf_OFF, img_sampling, calib_tt, cx=None, cy=None, 
        radii={'inner':(0,1.7),'outer':(1.7,2.3),'full':(0, 2.7)},
        nbin=0, ratio=0, plot_fig=True, verbose=False):

    """
    Args:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
        img_sampling (float):
            image sampling in pix per lambda/D
        calib_tt : 1D array
            when in model calibration mode, the true tip-tilt amplitude must be 
            provided in order to perform the fit of the models.
    """

    # get flux from off-axis PSF frame (photutils aperture photometry)
    psf_flux = get_psf_flux(psf_OFF, img_sampling, cx=cx, cy=cy, verbose=verbose)

    # bin + normalize on-axis PSF cube
    psf_ON = bin_images(psf_ON, nbin)
    psf_ON /= psf_flux

    # compute the differential intensities in the 3 regions
    all_di_mod, _ = get_all_di(psf_ON, radii, img_sampling, ratio=ratio, cx=cx, cy=cy)

    #***** model calibration mode ********************************************
    tt_fit_lim = {'inner':(0., 0.1), 'outer':(0., 0.5), 'full':(0.2, 0.5)}
    if plot_fig is True:
        colors = {'inner':[0.,0.3,.7],'outer':[.7,0.,0.3],'full':[0.,.7,0.5]}
        plt.figure(num=1, figsize=(12,9))
        plt.clf()
        fig, ax = plt.subplots(nrows=3,ncols=2,num=1)
        fig.subplots_adjust(hspace=0)
    coeffs = {}
    for i, region in enumerate(['inner', 'outer', 'full']):
        ind_x = np.where((calib_tt > tt_fit_lim[region][0]) & 
                         (calib_tt < tt_fit_lim[region][1]))[0]
        x = calib_tt[ind_x]
        yy  = all_di_mod[region]
        y = yy[ind_x]
        if region == 'full':
            y  = np.abs(y)**(1/3) # full estimator 
        coeff = linregress(x,y).slope
        fit_coeff = calib_tt*coeff
        if region == 'full':
            coeff = coeff**3
            fit_coeff = fit_coeff**3
        error = (yy - fit_coeff)/yy*100
        coeffs[region] = coeff

        if plot_fig is True:
            ax[i,0].set_xlabel(r'True tip-tilt [$\lambda/D$]')
            ax[i,0].set_ylabel('Normalized Diff. Intensity')
            ax[i,0].plot(calib_tt, yy, 'o', color=colors[region], alpha=.9, markersize=2,
                        label=region+r' - r = {0:.1f} to {1:.1f} $\lambda/D$'
                        .format(radii[region][0], radii[region][1]))
            ax[i,0].plot(calib_tt, fit_coeff, color=colors[region], alpha=.6, linestyle='--', 
                        label=r'Fit coeff. [{0:.2f}-{1:.2f}] $\lambda/D$ = {2:.3f}'
                        .format(tt_fit_lim[region][0],tt_fit_lim[region][1],coeff))
            ax[i,0].grid(color='.8',linestyle='--')
            ax[i,0].set_xlim(0.,)
            ax[i,0].legend()
            #error = ((yy - fit_coeff) / psf_flux) * 100
            ax[i,1].set_xlabel(r'True tip-tilt [$\lambda/D$]')
            ax[i,1].set_ylabel('Model Error [%]')
            ax[i,1].plot(calib_tt, error, 
                        'o', markersize=2, color=colors[region], alpha=.6)
            ax[i,1].set_ylim(-20., 20.)
            ax[i,1].set_xlim(0.,)
            ax[i,1].grid(color='.8',linestyle='--')

    if verbose is True:
        print('\nModel calibration results:'+
                '\nInner slope = {0:.3f}\nOuter slope = {1:.3f}\nFull coeff  = {2:.3f}'
                .format(*coeffs.values()))
    
    return coeffs