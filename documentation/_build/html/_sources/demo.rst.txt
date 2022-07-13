.. _demo:

Basic usage
######################

This is intended as a very brief overview of the steps necessary to get ``QACITS`` running.
This quick step guide also exists in a Jupyter notebook version: `demo.ipynb <https://github.com/vortex-exoplanet/QACITS/blob/main/notebooks/demo.ipynb>`_.


Import packages
===========================

.. code-block:: python

    from qacits import calibrate_qacits, run_qacits
    import os
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

Load demo data
===========================

QACITS is supposed to run in closed loop during the full science data acquisition 
sequence. To predict pointing sensing performance, the user can also run QACITS 
on simulated data sets. For the purpose of this demonstration, we use a sample of 
short simulated observing sequences (8 seconds each).

.. code-block:: python

    os.chdir(os.path.normpath(os.path.expandvars('$HOME/QACITS')))
    psf_OFF = fits.getdata('data/offaxis_PSF_L_CVC.fits')
    psf_ON_calib = fits.getdata('data/onaxis_PSF_L_CVC_calib.fits.gz')
    psf_ON_jitter = fits.getdata('data/onaxis_PSF_L_CVC_jitter.fits.gz')
    tt_lamD_calib = fits.getdata('data/point_8s_100ms_L_calib.fits')
    tt_lamD_jitter = fits.getdata('data/point_8s_100ms_L_jitter.fits')
    img_sampling = 3.98 # pix per lambda/D
    lamD = 21.78 # mas per lambda/D
    dit = 0.1 # seconds


QACITS default parameters
===========================

.. code-block:: python

    qacits_params = dict(
        force = None,
        coeffs = {'inner':1, 'outer':1, 'full':1},
        radii = {'inner':(0,1.7), 'outer':(1.7,2.3), 'full':(0, 2.7)},
        nbin = 0,
        ratio = 0,
        phase_tolerance = 60,
        modul_tolerance = 0.33,
        small_tt_regime = 0.3,
        large_tt_regime = 0.2, 
    )


Calibrate QACITS with known tip-tilt
=======================================

The linear coefficients in the QACITS model linking the differential intensity 
measured in the coronagraphic image and the pointing offsets that we want to 
determine need to be optimized because they depend on the shape of the telescope 
pupil and of the Lyot stop. During the calibration process, we measure the 
differential intensities in two regions: the inner region ranging from 0 to 
1.7 λ/D, and the outer region ranging from 1.7 to 2.3 λ/D. In addition to the 
“inner” and “outer” estimators, we also assess the “full” estimator, based on 
the cubic law relation between differential intensities and pointing error, over 
the full region (0 - 2.7 λ/D). For the remainder of our analysis, we only 
consider the outer estimator because it provides small model errors (<10%) over 
a large range of offsets (up to 0.5 λ/D).

.. code-block:: python

    qacits_params['coeffs'] = calibrate_qacits(psf_ON_calib, psf_OFF, img_sampling, tt_lamD_calib, plot_fig=True, verbose=True)
    qacits_params['force'] = 'outer'


Run QACITS: pointing error estimation
=======================================

Now that the QACITS algorithm has been calibrated, we can continue our 
performance analysis using the linear coefficient from the outer estimator, and 
perform simulations of realistic pointing error estimation.

.. code-block:: python

    tiptilt_estimate = run_qacits(psf_ON_jitter, psf_OFF, img_sampling, **qacits_params)

The first two columns of the output array correspond to the final estimated x 
and y tilt. They can be converted to milliarcseconds (mas) and compared to the 
read tip-tilt, in order to infer the RMS error value. 

.. code-block:: python

    # convert to mas
    tt_est = -tiptilt_estimate[:,0:2]*lamD
    tt_true = tt_lamD_jitter*lamD
    # rms error
    dist = np.sqrt((tt_est[:,0]-tt_true[:,0])**2 + (tt_est[:,1]-tt_true[:,1])**2)
    rms = np.sqrt(np.mean(dist**2))
