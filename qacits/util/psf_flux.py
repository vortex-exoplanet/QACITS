import numpy as np
try:
    # import photutils
    from photutils import aperture
    _exact_default_ = True
except:
    _exact_default_ = False

def get_psf_flux(img, radius, cx=None, cy=None, exact=_exact_default_, verbose=False):
    """ 
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
    """

    ny, nx = img.shape
    if cx == None :
        cx = (nx - 1)/2
    if cy == None :
        cy = (ny - 1)/2

    if exact is False:
        x, y = np.meshgrid(np.arange(nx)-cx, np.arange(ny)-cy)
        r = np.abs(x + 1j*y)
        psf_flux = np.sum(img*(r <= radius))
    else:
        aper = aperture.CircularAperture((cx, cy), radius)
        psf_flux = aperture.aperture_photometry(img, aper, method='exact')['aperture_sum'][0]

    if verbose is True:
        print('psf_flux = %s (photutils is %s)'%(np.round(psf_flux, 5), exact))

    return psf_flux


def get_di_xy(cube, radius, cx=None, cy=None, exact=_exact_default_):
    """ 
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
    """

    cube = np.array(cube, ndmin=3)
    ncube, ny, nx = cube.shape
    if cx == None :
        cx = (nx - 1)/2
    if cy == None :
        cy = (ny - 1)/2

    di_xy = []
    if radius == 0:
        di_xy = np.zeros((ncube, 2))

    elif exact is False:
        x, y = np.meshgrid(np.arange(nx)-cx, np.arange(ny)-cy)
        r = np.abs(x + 1j*y)
        for img in cube:
            img_mask = img*(r <= radius)
            Sxy = np.sum(img_mask)
            Sx = np.cumsum(np.sum(img_mask, axis=0))
            Sy = np.cumsum(np.sum(img_mask, axis=1))
            Ix = Sxy - 2*np.interp(cx, np.arange(nx)+.5, Sx)
            Iy = Sxy - 2*np.interp(cy, np.arange(ny)+.5, Sy)
            di_xy.append([Ix, Iy])
    else:
        aper = aperture.CircularAperture((cx, cy), radius)
        y, x = np.indices((ny,nx))
        cx1 = np.floor(cx)-1
        cy1 = np.floor(cy)-1
        for img in cube:
            Sxy = aperture.aperture_photometry(img, aper, method='exact')['aperture_sum'][0]
            x_flux1 = aperture.aperture_photometry(img*(x<=cx1), aper, method='exact')['aperture_sum'][0]
            x_flux2 = aperture.aperture_photometry(img*(x<=(cx1+1)), aper, method='exact')['aperture_sum'][0]
            x_flux3 = aperture.aperture_photometry(img, aper, method='exact')['aperture_sum'][0]
            Sx = [x_flux1, x_flux2, x_flux3]
            y_flux1 = aperture.aperture_photometry(img*(y<=cy1), aper, method='exact')['aperture_sum'][0]
            y_flux2 = aperture.aperture_photometry(img*(y<=(cy1+1)), aper, method='exact')['aperture_sum'][0]
            y_flux3 = aperture.aperture_photometry(img, aper, method='exact')['aperture_sum'][0]
            Sy = [y_flux1, y_flux2, y_flux3]
            Ix = Sxy - 2*np.interp(cx, [cx1+.5, cx1+1.5, nx], Sx)
            Iy = Sxy - 2*np.interp(cy, [cy1+.5, cy1+1.5, ny], Sy)
            di_xy.append([Ix, Iy])

    return np.float32(di_xy)


def get_all_di(cube, radii, img_sampling, ratio=0, cx=None, cy=None, exact=_exact_default_):
    """ 
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
    """

    all_dixy = {}
    for region in ['inner', 'outer', 'full']:
        r0, r1 = radii[region]
        all_dixy0 = get_di_xy(cube, r0*img_sampling, cx=cx, cy=cy)
        all_dixy1 = get_di_xy(cube, r1*img_sampling, cx=cx, cy=cy)
        all_dixy[region] = all_dixy1 - all_dixy0

    #-- Debias the diff. int. computed on full area from the linear component
    #   (estimated from the outer area)
    all_dixy['full'] += (ratio * all_dixy['outer'])

    # Transform to mod-arg (modulus-argument)
    all_di_mod = {}
    all_di_arg = {}
    for region in ['inner', 'outer', 'full']:
        all_di_mod[region] = np.sqrt(all_dixy[region][:,0]**2 + all_dixy[region][:,1]**2)
        all_di_arg[region] = np.arctan2(all_dixy[region][:,1], all_dixy[region][:,0])

    return all_di_mod, all_di_arg