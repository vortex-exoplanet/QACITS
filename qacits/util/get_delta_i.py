import numpy as np
try:
    import photutils
    _exact_default_ = True
except:
    _exact_default_ = False


def get_delta_i(cube, cx=None, cy=None, radius=None, exact=_exact_default_):
    """ 
    Computes the differential intensities along the x and y axes, based on the
    photmetry routines in the module photutils (if available). 
    
    Parameters
    ----------
    cube : array_like
        Input 2D array.
    cx : float, optional
        x position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.
    cy : float, optional
        y position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.
    radius : float, optional
        radius of the region of interest (full, inner, or outer area).
        Must be specified to use photutils.

    Returns
    -------
    delta_i : array_like
        2D element containing the differential intensities measured along the
        x and y axes.
    """

    cube = np.array(cube, ndmin=3)
    ny, nx = cube[0].shape
    if cx == None :
        cx = nx // 2
    if cy == None :
        cy = ny // 2

    delta_i = []

    if exact is False or radius is None:
        for img in cube:
            Sx = np.cumsum(np.sum(img, axis=0))
            Sy = np.cumsum(np.sum(img, axis=1))
            Ix = Sx[-1] - 2. * np.interp(cx, np.arange(nx)+.5, Sx)
            Iy = Sy[-1] - 2. * np.interp(cy, np.arange(ny)+.5, Sy)
            delta_i.append(np.array([Ix, Iy]))

    else:
        y, x = np.indices((ny,nx))
        cx1 = np.floor(cx)-1
        cy1 = np.floor(cy)-1
        aper = photutils.CircularAperture((cx, cy), radius)
        for img in cube:
            obj_flux = photutils.aperture_photometry(img*(x<=cx1), aper, method='exact')
            x_flux1 = obj_flux['aperture_sum'][0]
            obj_flux = photutils.aperture_photometry(img*(x<=(cx1+1)), aper, method='exact')
            x_flux2 = obj_flux['aperture_sum'][0]
            obj_flux = photutils.aperture_photometry(img, aper, method='exact')
            x_flux3 = obj_flux['aperture_sum'][0]
            Sx = [x_flux1,x_flux2,x_flux3]
            obj_flux = photutils.aperture_photometry(img*(y<=cy1), aper, method='exact')
            y_flux1 = obj_flux['aperture_sum'][0]
            obj_flux = photutils.aperture_photometry(img*(y<=(cy1+1)), aper, method='exact')
            y_flux2 = obj_flux['aperture_sum'][0]
            obj_flux = photutils.aperture_photometry(img, aper, method='exact')
            y_flux3 = obj_flux['aperture_sum'][0]
            Sy = [y_flux1,y_flux2,y_flux3]
            Ix = Sx[-1] - 2. * np.interp(cx, [cx1+.5,cx1+1.5,nx], Sx)
            Iy = Sy[-1] - 2. * np.interp(cy, [cy1+.5,cy1+1.5,ny], Sy)
            delta_i.append(np.array([Ix, Iy]))        

    return np.array(delta_i)