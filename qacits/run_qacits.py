from qacits.util.bin_images import bin_images
from qacits.util.circle_mask import circle_mask
from qacits.util.get_delta_i import get_delta_i
from photutils import CircularAperture, aperture_photometry


import numpy as np


def run_qacits(psf_ON, psf_OFF, fwhm, inner_slope=1, outer_slope=1, full_coeff=1,
        radii={'inner':(0,1.7),'outer':(1.7,2.3),'full':(0, 2.7)}, force=None,
        nbin=0, ratio=0, phase_tolerance=60, modul_tolerance=0.33,
        small_tt_regime=0.3, large_tt_regime=0.2, verbose=False, **params):

    """
    Args:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
        fwhm (float):
            full width half max, in pix per lambda/D
    """

  

    ### Compute the differential intensities in the 3 regions #################

    # get flux from off-axis PSF frame (photutils aperture photometry)
    nimg = psf_OFF.shape[-1]
    aper = CircularAperture((nimg//2, nimg//2), r=fwhm/2)
    psf_flux = aperture_photometry(psf_OFF, aper)['aperture_sum'].data
    if verbose is True:
        print('photutils aperture photometry: psf_flux=%s'%np.round(psf_flux,5))

    # bin + normalize on-axis PSF cube
    psf_ON = bin_images(psf_ON, nbin)
    psf_ON /= psf_flux
    ncube, nx, ny = psf_ON.shape

    # create the corresponding masks
    # compute the differential intensities

    masks = {}
    for region in radii:
        masks[region]= (circle_mask(psf_ON[0], radii[region][1] * fwhm, 
                              cx=nx//2, cy=ny//2) * 
                   (1 - circle_mask(psf_ON[0], radii[region][0] * fwhm, 
                              cx=nx//2, cy=ny//2)))
    
    all_dixy = {}
    for region in masks:
        all_dixy1 = get_delta_i(psf_ON, radii[region][1] * fwhm,
                                cx=nx//2, cy=ny//2)
        if radii[region][0] !=0 :
            all_dixy0 = get_delta_i(psf_ON, radii[region][0] * fwhm,
                                    cx=nx//2, cy=ny//2)
        else :
            all_dixy0 = 0.
        all_dixy[region] = all_dixy1 - all_dixy0

    
    #-- Debias the diff. int. computed on full area from the linear component
    #   (estimated from the outer area)
    all_dixy['full'] += (ratio * all_dixy['outer'])

    # Transform to mod-arg (modulus-argument)
    all_di_mod = {}
    all_di_arg = {}
    for region in ['inner','outer','full']:
        all_di_mod[region] = np.sqrt(all_dixy[region][:,0]**2 + all_dixy[region][:,1]**2)
        all_di_arg[region] = np.arctan2(all_dixy[region][:,1], all_dixy[region][:,0])

    #-- inner region: linear
    inner_est      = np.zeros((ncube, 2))
    inner_est[:,0] = all_di_mod['inner'] / inner_slope
    inner_est[:,1] = all_di_arg['inner'] + np.pi
    #-- outer region: linear
    outer_est      = np.zeros((ncube, 2))
    outer_est[:,0] = all_di_mod['outer'] / outer_slope
    outer_est[:,1] = all_di_arg['outer']
    #-- full region: cubic
    full_est      = np.zeros((ncube, 2))
    full_est[:,0] = np.abs(all_di_mod['full']/full_coeff)**(1./3.)
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
    
    ### Final estimator in X,Y ################################################
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