.. _QACITS:

Software design
######################

Main modules
========================

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

Utilities
========================

QACITS provides a few utilities for image processing and photometry, stored as 
standalone modules in the ``util`` folder.


bin_images
--------------
Returns a cube of images averaged by bins of nbin images.
    
Args:
    cube (float ndarray):
        single image or image cube of ncube frames
    nbin (int):
        number of images in the returned cube
        
        - if 0: no frame averaging, the returned cube is a copy of the input
        - if 1: returns the average of the whole cube
        - else: the returned cube is made of nbin images, each being the average of ncube/nbin images of the input cube

Returns:
    cube_binned (float ndarray):
        binned image cube


get_psf_flux
--------------
Computes the aperture photometry of the PSF core for a given radius
generally corresponding to half the FWHM. Based on photometry routines
from the photutils package (if available). 

Args:
    cube (float ndarray):
        single image or image cube of ncube frames
    radius (float):
        radius [pix] of the PSF core

Returns:
    psf_flux (float):
        flux in the PSF core

get_di_xy
-------------
Computes the differential intensities along the x and y axes. Based on 
photometry routines from the photutils package (if available). 
    
Args:
    cube (float ndarray):
        single image or image cube of ncube frames
    radius (float):
        radius [pix] of the region of interest (full, inner, or outer area)
    cx (float, optional):
        x position of the sub-image center [pix], defaults to the image center
    cy (float, optional):
        y position of the sub-image center [pix], defaults to the image center

Returns:
    di_xy (float ndarray):
        2D element containing the differential intensities measured along the
        x and y axes

get_all_di
------------
Computes the differential intensities for all regions (full, inner, or outer 
area) along the x and y axes, then converts the outputs to modulus and argument.
Based on photometry routines from the photutils package (if available).
    
Args:
    cube (float ndarray):
        single image or image cube of ncube frames
    radii (dict):
        dictionary containing the radii in lambda/D for each region of interest 
    img_sampling (float):
        image sampling in pix per lambda/D

Returns:
    all_di_mod (dict):
        dictionary containing the modulus of the differential intensities
        for each region of interest 
    all_di_arg (dict):
        dictionary containing the argument of the differential intensities
        for each region of interest 

Archive: qacits_vlt_package_v4_ehuby
------------------------------------------