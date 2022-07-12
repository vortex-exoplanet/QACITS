.. _QACITS:

QACITS package
######################

Running QACITS consists of two phases. First, we determine the linear 
coefficients in the QACITS model during the calibration phase, using known 
tip-tilt offsets. Then, we use QACITS with the calculated coefficients to 
perform real-time pointing error estimation on realistic data sets.

calibrate_qacits
-----------------
Calibration function for computing the linear coefficients in the QACITS model
using known tip-tilt offsets, based on the QACITS method for a Vortex 
coronagraph of charge 2.

    Args:
        psf_ON (float ndarray):
            cube of on-axis PSFs
        psf_OFF (float ndarray):
            off-axis PSF frame
        img_sampling (float):
            image sampling in pix per lambda/D
        tt_lamD (2D float ndarray):
            true x and y tip-tilt values in lambda/D used to fit the model
    Return:
        coeffs (dict of float):
            linear coefficients in the QACITS model

run_qacits
-----------------
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
        coeffs (dict of float):
            linear coefficients in the QACITS model
    Return:
        full_estimate_output (float ndarray):
            full estimate output; first two columns are tip-tilt estimate