import numpy as np

def circle_mask(image, radius_pix, cx=None, cy=None):
    """ 
    Creates a circular mask with values of 1 or 0, of the same dimensions as the
    input image, and defined by a radius and center coordinates. 
    Values are 1 inside the circle, 0 outside.

    Parameters
    ----------
    image : array_like
        Input 2D array.
    radius_pix : float
        Radius of the mask in pixels.
    cx : float, optional
        x position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.
    cy : float, optional
        y position of the sub-image center [pix], can be fractional.
        If not specified, the center is defined at the center of the image.

    Returns
    -------
    circle_mask : array_like
        2D array containing 1 and 0 values only.
    """
    ny, nx = image.shape

    if cx == None :
        cx = (nx-1.) / 2.
    if cy == None :
        cy = (ny-1) / 2.

    cx2 = np.round(cx*2.)/2.
    cy2 = np.round(cy*2.)/2.

    gridy, gridx = np.indices(image.shape)
    gridx = gridx - cx2
    gridy = gridy - cy2
    gridr = np.sqrt(gridx**2. + gridy**2.)

    circle_mask = gridr < radius_pix

    return circle_mask