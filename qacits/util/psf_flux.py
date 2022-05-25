import numpy as np
try:
    import photutils
    _exact_default_ = True
except:
    _exact_default_ = False


def get_psf_flux(img, radius, cx=None, cy=None, exact=_exact_default_, verbose=False):

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
        aper = photutils.CircularAperture((cx, cy), radius)
        psf_flux = photutils.aperture_photometry(img, aper, method='exact')['aperture_sum'][0]

    if verbose is True:
        print('psf_flux = %s (photutils = %s)'%(np.round(psf_flux, 5), exact))

    return psf_flux


def get_di_xy(cube, radius, cx=None, cy=None, exact=_exact_default_):
    """ 
    Computes the differential intensities along the x and y axes, based on the
    photmetry routines in the module photutils (if available). 
    
    Parameters
    ----------
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

    Returns
    -------
    di_xy : array_like
        2D element containing the differential intensities measured along the
        x and y axes.
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
        aper = photutils.CircularAperture((cx, cy), radius)
        y, x = np.indices((ny,nx))
        cx1 = np.floor(cx)-1
        cy1 = np.floor(cy)-1
        for img in cube:
            Sxy = photutils.aperture_photometry(img, aper, method='exact')['aperture_sum'][0]
            x_flux1 = photutils.aperture_photometry(img*(x<=cx1), aper, method='exact')['aperture_sum'][0]
            x_flux2 = photutils.aperture_photometry(img*(x<=(cx1+1)), aper, method='exact')['aperture_sum'][0]
            x_flux3 = photutils.aperture_photometry(img, aper, method='exact')['aperture_sum'][0]
            Sx = [x_flux1, x_flux2, x_flux3]
            y_flux1 = photutils.aperture_photometry(img*(y<=cy1), aper, method='exact')['aperture_sum'][0]
            y_flux2 = photutils.aperture_photometry(img*(y<=(cy1+1)), aper, method='exact')['aperture_sum'][0]
            y_flux3 = photutils.aperture_photometry(img, aper, method='exact')['aperture_sum'][0]
            Sy = [y_flux1, y_flux2, y_flux3]
            Ix = Sxy - 2*np.interp(cx, [cx1+.5, cx1+1.5, nx], Sx)
            Iy = Sxy - 2*np.interp(cy, [cy1+.5, cy1+1.5, ny], Sy)
            di_xy.append([Ix, Iy])

    return np.float32(di_xy)


def get_all_di(cube, radii, img_sampling, ratio=0, cx=None, cy=None, exact=_exact_default_):
    """
    Compute the differential intensities for all regions.
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