#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy

def normalize(v):
    """Normalize the data.

    @param v: input volume.

    @return: Normalized volume.
    """
    m = np.mean(v)
    v = v-m
    s = np.std(v)
    v = v/s
    return v

def bandpass(v, low=0, high=-1, sigma=0):
    """Do a bandpass filter on a given volume.

    @param v: input volume.
    @param low: low frequency in k-space.
    @param high: high frequency in k-space.
    @param sigma: smoothness factor.

    @return: bandpass filtered volume.
    """
    assert low >= 0
    from tools import create_sphere
    if high == -1:
        high = np.min(v.shape)/2
    assert low < high

    if low == 0:
        mask = create_sphere(v.shape, high, sigma)
    else:
        # BUG! TODO
        # the sigma
        mask = create_sphere(v.shape, high, sigma) - create_sphere(v.shape, low, sigma)

    from transform import fourier_filter
    res = fourier_filter(v, mask, True)

    return res


def median3d(data, size=3):
    """Median filter.

    @param data: data to be filtered.
    @param size: size of the median filter.

    @return: filtered image.
    """
    from scipy.ndimage.filters import median_filter
    d = median_filter(data, size)
    return d


def gaussian3d(data, sigma=0):
    """Gaussian filter.

    @param data: data to be filtered.
    @param sigma: sigma of Gaussian.

    @return: filtered data.
    """
    from scipy.ndimage.filters import gaussian_filter
    d = gaussian_filter(data, sigma)
    return d


class Wedge(object):
    """Defines the missing wedge filter in Fourier space."""
    def __init__(self):
        super(Wedge, self).__init__()

    def apply(self, data):
        raise NotImplementedError('Abstract method! Please overwrite it.')

    def toSphericalFunc(self, bw, radius):
        raise NotImplementedError('Abstract method! Please overwrite it.')


class GeneralWedge(Wedge):
    """General wedge."""
    def __init__(self, wedge_vol, half=True, isodd=False):
        """Initialize a general wedge with given Fourier volume.

        @param half: if true, the given volume is in half and zero frequency at the corner.
                     Otherwise, the given volume is full and zero frequency in the center.
        @param isodd: when half is true, this field tells the z dimension of the full Fourier volume is odd-/even-sized.
        """
        self.set_wedge_volume(wedge_vol, half)

    def set_wedge_volume(self, wedge_vol, half=True, isodd=False):
        if half:
            self._volume = wedge_vol

            # human understandable version with 0-freq in the center
            from transform import fourier_reduced2full, fftshift
            self._whole_volume = fftshift(fourier_reduced2full(self._volume, isodd))
        else:
            self._whole_volume = wedge_vol

            from transform import fourier_full2reduced, ifftshift
            self._volume = fourier_full2reduced(ifftshift(self._whole_volume))

    def apply(self, data, rotation=None):
        if rotation is not None: # rotate the wedge first
            assert len(rotation) == 3
            from transform import rotate3d, fourier_full2reduced, ifftshift
            filter_vol = rotate3d(self._whole_volume, rotation[0], rotation[1], rotation[2])
            filter_vol = fourier_full2reduced(ifftshift(filter_vol))
        else:
            filter_vol = self._volume

        from transform import fourier_filter
        res = fourier_filter(data, filter_vol, False)

        return res

    def toSphericalFunc(self, bw, radius):
        assert(bw<=128)

        # start sampling
        from vol2sf import vol2sf
        sf = vol2sf(self._whole_volume, radius, bw)
        
        return sf


class SingleTiltWedge(Wedge):
    """Missing wedge of single tilt geometry. Assume Y axis is the rotation axis."""
    def __init__(self, start_ang, end_ang):
        super(SingleTiltWedge, self).__init__()
        self.start_ang = start_ang
        self.end_ang = end_ang
        self._volume_shape = None # store the wedge volume shape (whole!)
        self._volume = None # store the wedge volume in k-space (only half!)

        self._bw = None # store the bandwidth of the spherical function
        self._sf = None # store the spherical function
    
    def _create_wedge_volume(self, size):
        from transform import fftshift, fourier_full2reduced
        # if no missing wedge
        if self.start_ang == -90 and self.end_ang == 90:
            filter_vol = np.ones(size)
            self._volume = fourier_full2reduced(filter_vol)
            self._volume_shape = size
            return

        filter_vol = np.ones(size)
        x, z = scipy.mgrid[0.:size[0], 0.:size[2]]
        x -= size[0]/2
        ind = np.where(x) # find the non-zeros
        z -= size[2]/2

        angles = np.zeros(z.shape)
        angles[ind] = np.arctan(z[ind]/x[ind]) * 180 / np.pi

        angles = np.reshape(angles, (size[0],1,size[2]))
        angles =  np.repeat(angles, size[1], axis=1)

        filter_vol[angles > -self.start_ang] = 0
        filter_vol[angles < -self.end_ang] = 0

        filter_vol[size[0]/2, :, :] = 0
        filter_vol[size[0]/2, :, size[2]/2] = 1

        # create a sphere and multiple it with the wedge
        from tools import create_sphere
        mask = create_sphere(size)
        filter_vol *= mask

        # shift and cut in half
        filter_vol = fftshift(filter_vol)
        self._volume = fourier_full2reduced(filter_vol)

        self._volume_shape = size

    def apply(self, data, rotation=None):
        """
        @param rotation: apply rotation to the wedge first
        """
        # if no missing wedge
        if self.start_ang == -90 and self.end_ang == 90:
            return data

        if self._volume is not None and np.array_equal(self._volume_shape, data.shape):
            pass
        else:
            self._create_wedge_volume(data.shape)

        if rotation is not None: # rotate the wedge first
            assert len(rotation) == 3
            from transform import rotate3d, fourier_reduced2full, fourier_full2reduced, fftshift, ifftshift
            isodd = self._volume_shape[2] % 2
            filter_vol = fftshift(fourier_reduced2full(self._volume, isodd))
            filter_vol = rotate3d(filter_vol, rotation[0], rotation[1], rotation[2])
            filter_vol = fourier_full2reduced(ifftshift(filter_vol))
        else:
            filter_vol = self._volume

        from transform import fourier_filter
        res = fourier_filter(data, filter_vol, False)

        return res

    def toSphericalFunc(self, bw, radius=None, threshold=0.5):
        """Convert the wedge from k-space to a spherical function.
        
        @param bw: bandwidth of the spherical function.
        @param radius: radius in k-space. For general Wedge, not used for SingleTiltWedge.
        @param threshold: threshold, above which the value there would be set to 1.

        @return: a spherical function in numpy.array
        """
        assert(bw<=128)

        # if no missing wedge
        if self.start_ang == -90 and self.end_ang == 90:
            self._sf = np.ones((4*bw**2,))
            return self._sf

        r = 45 # this radius and the volume size should be sufficient for sampling b <= 128
        if self._volume is None:
            self._create_wedge_volume((100,100,100))
        
        if self._bw == bw and self._sf is not None:
            return self._sf
        else:
            self._bw = bw
        
        from transform import fourier_reduced2full, fftshift
        isodd = self._volume_shape[2] % 2
        filter_vol = fftshift(fourier_reduced2full(self._volume, isodd))

        # start sampling
        from math import pi, sin, cos
        res = []
        
        for j in xrange(2*bw):
            for k in xrange(2*bw):
                the = pi*(2*j+1)/(4*bw) # (0,pi)
                phi = pi*k/bw # [0,2*pi)
                
                # this part actually needs interpolation
                x = int(cos(phi)*sin(the)*r+50)
                y = int(sin(phi)*sin(the)*r+50)
                z = int(cos(the)*r+50)
                
                # if the value is bigger than the threshold, we include it
                if filter_vol[x,y,z] > threshold:
                    res.append(1.0)
                else:
                    res.append(0.0)
        
        # store it so that we don't have to recompute it next time
        self._sf = np.array(res)
        
        return self._sf


