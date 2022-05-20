import numpy as np


def calibrate_qacits(calib_tt=None):


    calib_tt=None, model_calibration=False, , 




    #***** model calibration mode ********************************************
    tt_fit_lim = {'inner':(0., 0.1),'outer':(0., 0.5),'full':(0.2, 0.5)}
    coeffs = {}
    for i, region in enumerate(tt_fit_lim):
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
        error = ((yy - fit_coeff) / yy) * 100
        #error = ((yy - fit_coeff) / psf_flux) * 100
        ax[i,1].set_xlabel(r'True tip-tilt [$\lambda/D$]')
        ax[i,1].set_ylabel('Model Error [%]')
        ax[i,1].plot(calib_tt, error, 
                    'o', markersize=2, color=colors[region], alpha=.6)
        ax[i,1].set_ylim(-20., 20.)
        ax[i,1].set_xlim(0.,)
        ax[i,1].grid(color='.8',linestyle='--')
        coeffs[region] = coeff

    qacits_params['inner_slope'] = coeffs['inner']
    qacits_params['outer_slope'] = coeffs['outer']
    qacits_params['full_coeff'] = coeffs['full']
    print('\nModel calibration results:'+
            '\nInner slope = {0:.3f}\nOuter slope = {1:.3f}\nFull coeff  = {2:.3f}'
            .format(*coeffs.values()))