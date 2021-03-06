.. _utils:

Utilities
####################

QACITS provides a few utilities for image processing and photometry, stored as 
standalone modules in the ``util`` folder.

Image processing
========================

bin_images
--------------

Returns a cube of images averaged by bins of nbin images.
    
    Args:
        cube (float ndarray):
            single image or image cube of ncube frames
        nbin (int):
            number of images in the returned cube.
            
            - if 0: no frame averaging, the returned cube is a copy of the input.
            - if 1: returns the average of the whole cube.
            - else: the returned cube is made of nbin images, each being the average of ncube/nbin images of the input cube.
    
    Returns:
        cube_binned (float ndarray):
            binned image cube.


Photometry
========================

get_psf_flux
--------------

get_di_xy
-------------
Computes the differential intensities along the x and y axes, based on the
photometry routines in the module photutils (if available). 

    Args:
        cube : array_like
            Input 2D array.
        radius : float
            radius of the region of interest (full, inner, or outer area).
        cx : float, optional
            x position of the sub-image center [pix], can be fractional.
            If not specified, the center is defined at the center of the image.
        cy : float, optional
            y position of the sub-image center [pix], can be fractional.
            If not specified, the center is defined at the center of the image.

    Returns:
        di_xy : array_like
            2D element containing the differential intensities measured along the
            x and y axes.

get_all_di
------------
Compute the differential intensities for all regions.


Archive: qacits_vlt_package_v4_ehuby
=============================================

