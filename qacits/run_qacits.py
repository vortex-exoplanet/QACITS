from qacits.util.bin_images import bin_images
from qacits.util.psf_flux import get_psf_flux, get_all_di
import numpy as np


def run_qacits(psf_ON, psf_OFF, img_sampling, cx=None, cy=None, force='outer',
        coeffs={'inner':1, 'outer':1, 'full':1},
        radii={'inner':(0,1.7), 'outer':(1.7,2.3), 'full':(0, 2.7)},
        nbin=0, ratio=0, phase_tolerance=60, modul_tolerance=0.33,
        small_tt_regime=0.3, large_tt_regime=0.2, verbose=False, **qacits_params):

    """
    Pointing error estimation using the QACITS model with calibrated linear 
    coefficients (see calibrate_qacits.py), based on the QACITS method for a 
    Vortex coronagraph of charge 2.

    Args:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
        img_sampling (float):
            image sampling in pix per lambda/D
        cx (float):
            x position of the sub-image center [pix], defaults to the image center
        cy (float):
            y position of the sub-image center [pix], defaults to the image center
        force (str):
            force the QACITS estimator to use a specific estimator, defaults to 'outer'
        coeffs (dict of float):
            linear coefficients in the QACITS model

    Return:
        full_estimate_output (float ndarray):
            full estimate output; first two columns are tip-tilt estimate
    """

    # get flux from off-axis PSF frame (photutils aperture photometry)
    psf_flux = get_psf_flux(psf_OFF, img_sampling/2, cx=cx, cy=cy, verbose=verbose)

    # bin + normalize on-axis PSF cube
    psf_ON = bin_images(psf_ON, nbin)
    psf_ON /= psf_flux

    # compute the differential intensities in the 3 regions
    all_di_mod, all_di_arg = get_all_di(psf_ON, radii, img_sampling, ratio=ratio, cx=cx, cy=cy)
    
    # Pointing error estimation mode
    # ------------------------------
    
    ncube = len(psf_ON)
    inner_slope = coeffs['inner']
    outer_slope = coeffs['outer']
    full_coeff = coeffs['full']
    phase_tolerance = phase_tolerance
    modul_tolerance = modul_tolerance
    small_tt_regime = small_tt_regime
    large_tt_regime = large_tt_regime

    #-- inner region: linear
    inner_est      = np.zeros((ncube, 2))
    inner_est[:,0] = all_di_mod['inner'] / coeffs['inner']
    inner_est[:,1] = all_di_arg['inner'] + np.pi
    #-- outer region: linear
    outer_est      = np.zeros((ncube, 2))
    outer_est[:,0] = all_di_mod['outer'] / coeffs['outer']
    outer_est[:,1] = all_di_arg['outer']
    #-- full region: cubic
    full_est      = np.zeros((ncube, 2))
    full_est[:,0] = np.abs(all_di_mod['full']/coeffs['full'])**(1/3)
    full_est[:,1] = all_di_arg['full']    

    #-- Estimator selection: 
    final_est = np.zeros((ncube, 2))
    test_output = np.zeros((ncube, 3))
    if force == 'inner':
        final_est = inner_est
    elif force == 'outer':
        final_est = outer_est
    elif force == 'full':
        final_est = full_est
    else :
        for i in range(ncube):
            # modulus to be trusted for choosing tt regime
            outer_modulus = outer_est[i,0]
            
            # build complex phasors
            inner_phasor = inner_est[i,0] * np.exp(1j * inner_est[i,1]) 
            outer_phasor = outer_est[i,0] * np.exp(1j * outer_est[i,1])
            full_phasor  = full_est[i,0]  * np.exp(1j * full_est[i,1])
            
            # test estimate agreement
            #-- phase agreement: IN/OUT
            test_phasor = np.exp(1j * inner_est[i,1]) * np.exp(-1j * outer_est[i,1])
            inout_test_phase = np.arctan2(np.imag(test_phasor), np.real(test_phasor))
            in_out_phase_agreement = (np.abs(inout_test_phase) < phase_tolerance/180*np.pi)
            #-- phase agreement: OUT/FULL
            test_phasor = np.exp(1j * full_est[i,1]) * np.exp(-1j * outer_est[i,1])
            fullout_test_phase = np.arctan2(np.imag(test_phasor), np.real(test_phasor))
            full_out_phase_agreement = (np.abs(fullout_test_phase) < phase_tolerance/180*np.pi)
            full_out_modul_agreement = (np.abs(full_est[i,0]-outer_est[i,0]) < full_est[i,0]*modul_tolerance)
            if verbose is True :
                print(i, inner_est[i,0], outer_est[i,0], full_est[i,0], in_out_phase_agreement, 
                        np.abs(test_phasor)*180./np.pi, in_out_phase_agreement)
                print('{0:03d} -- '.format(i) +
                       '\n \t IN   {0:.3f} l/D {1:.1f} deg'.format(inner_est[i,0],inner_est[i,1]*180/np.pi)+
                       '\n \t OUT  {0:.3f} l/D {1:.1f} deg'.format(outer_est[i,0],outer_est[i,1]*180/np.pi)+
                       '\n \t FULL {0:.3f} l/D {1:.1f} deg'.format(full_est[i,0],full_est[i,1]*180/np.pi))
            
            ### SMALL TIPTILT REGIME
            if outer_modulus < small_tt_regime :
                if verbose is True :
                    print('\n \t > IN-OUT phase agreement is {}'.format(in_out_phase_agreement))
                if in_out_phase_agreement == True:
                    meanphasor = (inner_phasor + outer_phasor)/2.
                    if verbose is True :
                        print('\t => small1: in+out')
                else :
                    meanphasor = outer_phasor 
                    if verbose is True :
                        print('\t => small2: out')

            ### LARGE TIPTILT REGIME
            else :
                if verbose is True :
                    print('\n \t > FULL-OUT phase agreement is {}'.format(full_out_phase_agreement)+
                    '\n \t > FULL-OUT modulus agreement is {}'.format(full_out_modul_agreement))
                if full_out_phase_agreement == True:
                    if full_out_modul_agreement == True:
                        meanphasor = ( outer_phasor + full_phasor )/2.
                        if verbose is True :
                            print('\t => large1: full+out')
                    else :
                        if full_est[i,0] < 1.:
                            meanphasor = full_phasor
                        else:
                            meanphasor = np.exp(1j * full_est[i,1]) # set the maximal estimate to 1 lbd/D
                        if verbose is True :
                            print('\t => large2: full')
                else :
                    if full_est[i,0] < 1.:
                        meanphasor = full_phasor
                    else:
                        meanphasor = np.exp(1j * full_est[i,1]) # set the maximal estimate to 1 lbd/D
                    if verbose is True :
                        print('\t => large3: full')

            if verbose is True:
                final_phase = np.arctan2(np.imag(meanphasor), np.real(meanphasor))*180./np.pi
                print('\t    final estimator mod = {0:.3f} l/D phase = {1:.1f} deg'.format(np.abs(meanphasor), final_phase))
            final_est[i,0] = np.abs(meanphasor)
            final_est[i,1] = np.arctan2(np.imag(meanphasor),np.real(meanphasor))
            test_output[i,:] = np.array([inout_test_phase,fullout_test_phase,np.abs(full_est[i,0]-outer_est[i,0])])
    
    #-- Final estimator in X,Y: 
    final_est_xy = np.zeros_like(final_est)
    final_est_xy[:,0] = final_est[:,0] * np.cos(final_est[:,1])
    final_est_xy[:,1] = final_est[:,0] * np.sin(final_est[:,1])

    full_estimate_output = np.ndarray((ncube, 11))
    full_estimate_output[:,0:2] = final_est_xy
    full_estimate_output[:,2:4] = inner_est
    full_estimate_output[:,4:6] = outer_est
    full_estimate_output[:,6:8] = full_est
    full_estimate_output[:,8:] = test_output

    return full_estimate_output