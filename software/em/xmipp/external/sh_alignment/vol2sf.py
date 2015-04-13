#!/usr/bin/env python

'''
Created on Sep 29, 2010

@author: chen
'''

import numpy as np

def vol2sf(vol, r, b, center=None):
    """Transfer a volume into a serial of spherical functions. The volume will be decomposed into a serials of concentric shells.
    It will sample the volume on equiangular 2B*2B grid. (Trilinear interpolation) For more detail, see "http://www.cs.dartmouth.edu/~geelong/sphere/".
    
    Parameters
    ----------
    vol: Target volume
         numpy.ndarray

    r: The radius of shell you want to get from the volume (in voxel)
       Integer

    b: Bandwidth, which determines the sampling on the sphere
       Integer

    center: The center of the circle (By default is the center of the volume)
            list [x, y, z]
    
    Returns
    -------
    f(the_0, phi_0) ... f(the_0, phi_2B-1), f(the_2B-1, phi_0) ... f(the_2B-1, phi_2B-1)
    List
    """
    if r > vol.shape[0]/2 or r > vol.shape[1]/2 or r > vol.shape[2]/2:
        raise RuntimeError("Given radius is larger than the volume!")
    elif r <= 0:
        raise RuntimeError("Radius should be larger than the 0!")
    else:
        pass
    
    if center:
        m_x, m_y, m_z = center
    else:
        # the middle point of the volume
        m_x = vol.shape[0]/2; m_y = vol.shape[1]/2; m_z = vol.shape[2]/2

    the = np.pi*(2*np.arange(2*b)+1) / (4*b)
    phi = np.pi*np.arange(2*b) / b

    # C order mesh
    phi, the = np.meshgrid(phi, the)
    the = the.flatten()
    phi = phi.flatten()

    # compute the coordinates
    x = r*np.cos(phi)*np.sin(the) + m_x
    y = r*np.sin(phi)*np.sin(the) + m_y
    z = r*np.cos(the) + m_z

    # interpolate
    from scipy.ndimage import map_coordinates
    res = map_coordinates(vol, [x, y, z], order=2)
    
    return res


def fvol2sf(vol, r, b):
    """Transfer a volume in Fourier space into a serial of spherical functions.
    Makes use of the Hermitian symmetry to speed up.

    Parameters
    ----------
    vol: The volume of Fourier coefficients. It should be full and zero frequency in the center!
         numpy.ndarray

    r: Radius in k-space
       Integer

    b: Bandwidth, which determines the sampling on the sphere
       Integer

    Returns
    -------
    f(the_0, phi_0) ... f(the_0, phi_2B-1), f(the_2B-1, phi_0) ... f(the_2B-1, phi_2B-1)
    List
    """
    if r > vol.shape[0]/2 or r > vol.shape[1]/2 or r > vol.shape[2]/2:
        raise RuntimeError("Given radius is larger than the volume!")
    elif r <= 0:
        raise RuntimeError("Radius should be larger than the 0!")
    else:
        pass
    
    # zero frequency position
    m_x = vol.shape[0]/2; m_y = vol.shape[1]/2; m_z = vol.shape[2]/2

    # only half
    the = np.pi*(2*np.arange(b)+1) / (4*b)
    phi = np.pi*np.arange(2*b) / b

    # C order mesh
    phi, the = np.meshgrid(phi, the)
    the = the.flatten()
    phi = phi.flatten()

    # compute the coordinates
    x = r*np.cos(phi)*np.sin(the) + m_x
    y = r*np.sin(phi)*np.sin(the) + m_y
    z = r*np.cos(the) + m_z

    # use spline intepolation on the real/imaginary parts
    # can be done better, but now it suffices
    vol_r = np.real(vol)
    vol_i = np.imag(vol)

    from scipy.ndimage import map_coordinates
    res_r_a = map_coordinates(vol_r, [x, y, z], order=2)
    res_i_a = map_coordinates(vol_i, [x, y, z], order=2)

    # fill in the other half
    ind = np.arange(2*b**2)
    cind = (2*b-1-(ind+2*b**2)/(2*b))*2*b + np.mod(np.mod(ind, 2*b)+b, 2*b)

    res_r = np.zeros((4*b**2,))
    res_r[:2*b**2] = res_r_a
    res_r[2*b**2:] = res_r_a[cind]

    res_i = np.zeros((4*b**2,))
    res_i[:2*b**2] = res_i_a
    res_i[2*b**2:] = -res_i_a[cind]
    
    return [res_r, res_i]


def fourier_sf_shift(sf, r, shape, dx, dy, dz):
    """Do a Fourier shift of a spherical function to save some interpolation time.

    Parameters
    ----------
    sf: A (complex) spherical function in Fourier space
        numpy.ndarray

    r: Radius in k-space
       Integer

    shape: Original volume shape
           numpy.shape

    dx: Shift in X dim in real space
        Integer / float

    dy: Shift in Y dim in real space
        Integer / float

    dz: Shift in Z dim in real space
        Integer / float

    b: Bandwidth, which determines the sampling on the sphere
       Integer

    Returns
    -------
    A shifted (complex) spherical function in Fourier space.
    """
    b = int(len(sf)**0.5/2)
    assert len(sf) == 4*b**2

    the = np.pi*(2*np.arange(2*b)+1) / (4*b)
    phi = np.pi*np.arange(2*b) / b

    # C order mesh
    phi, the = np.meshgrid(phi, the)
    the = the.flatten()
    phi = phi.flatten()

    # compute the coordinates
    xx = r*np.cos(phi)*np.sin(the)
    yy = r*np.sin(phi)*np.sin(the)
    zz = r*np.cos(the)

    sx = shape[0]
    sy = shape[1]
    sz = shape[2]

    shift = np.exp(-2j*np.pi/sx*xx*dx) * np.exp(-2j*np.pi/sy*yy*dy) * np.exp(-2j*np.pi/sz*zz*dz)

    return sf*shift


def sf2file(sf, filename, conv=True):
    """Save the spherical function as a file.
    The file format is is same as "http://www.cs.dartmouth.edu/~geelong/sphere/".

    Parameters
    ----------
    sf: Spherical function
        List

    filename: File name to save as
              String
    """
    f = open(filename, 'w')
    
    for v in sf:
        f.write("%f\n" % v)
        if conv is False:
            f.write("0\n") # imaginary part is always 0
    
    f.close()


if __name__ == '__main__':
    import sys, getopt
    usage = './scriptname -v volume -r radius -b bandwidth -c flag_for_convolution -o output_filename'
    
    if len(sys.argv) == 1:
        print usage
        sys.exit()        
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hcv:r:b:o:", ["help"])
    except getopt.GetoptError:
        print 'Command not right. Exit!'
        sys.exit()
    
    conv = False
    for o,a in opts:
        if o in ("-h", "--help"):
            print usage
            sys.exit()
        if o in ("-v"):
            vol_name = a
        if o in ("-r"):
            radius = int(a)
        if o in ("-b"):
            bw = int(a)
        if o in ("-c"):
            conv = True
        if o in ("-o"):
            filename = a
    
    from tompy.io import read
    vol = read(vol_name)
    sf = vol2sf(vol, radius, bw)
    sf2file(sf, filename, conv)


