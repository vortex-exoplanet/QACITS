import numpy as np

def bin_images(cube, nbin):
    """
    Returns a cube of images averaged by bins of nbin images.
    
    Parameters
    ----------
    cube (float ndarray):
        single image or image cube of ncube frames
    nbin (int):
        number of images in the returned cube.
        - if 0: no frame averaging, the returned cube is a copy of the input.
        - if 1: returns the average of the whole cube.
        - else: the returned cube is made of nbin images, 
            each being the average of ncube/nbin images of the input cube.
    
    Returns
    -------
    cube_binned (float ndarray):
        binned image cube.
    """

    # image cube must be 3D numpy array
    cube = np.array(cube, ndmin=3)
    ncube, ny, nx = cube.shape
    assert nbin <= ncube, 'nbin must be <= ncube'

    # case 0:
    if nbin == 0:
        cube_binned = cube.copy()
    # case 1: all images averaged
    elif nbin == 1:
        cube_binned = np.array(np.mean(cube, axis=0), ndmin=3)
    # else: bin the images
    elif ncube > 1:
        cube_binned = np.zeros((nbin, ny, nx))
        bin_width = ncube // nbin
        i0 = ncube - bin_width * nbin
        for i in range(nbin):
            cube_binned[i] = np.mean(cube[i0+bin_width*i:i0+bin_width*(i+1),:,:], axis=0)
    else :
        cube_binned = cube.copy()

    return cube_binned