#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy

def NCC(v1, v2):
    """Compute the Normalized Cross Correlation between the two volumes.

    @param v1: volume 1.
    @param v2: volume 2.

    @return: NCC.
    """
    from filter import normalize
    vv1 = normalize(v1)
    vv2 = normalize(v2)
    score = np.sum(vv1*vv2)/vv1.size

    return score


################################
#     FLCF Relevant Stuff      #
################################

def mean_under_mask(vol, mask):
    # number of non-zeros inside the mask
    p = np.sum(mask)
    m = vol * mask
    
    return np.sum(m)/p


def std_under_mask(vol, mask, m):
    # number of non-zeros inside the mask
    p = np.sum(mask)
    
    squareV = vol**2
    
    res = mean_under_mask(squareV, mask) - m**2
    
    return np.abs(res)**0.5


def mean_vol_under_mask(volume, mask):
    """Calculate the mean volume under the mask (Both should have the same size).

    @param volume: input volume.
    @param mask: mask.
    
    @return: the calculated mean volume under mask.
    """
    p = np.sum(mask)

    # do the (circular) convolution
    from transform import rfft, irfft, fftshift
    size = volume.shape
    # somehow this should be conjugated
    res = fftshift(irfft(rfft(volume) * np.conjugate(rfft(mask)), size)) / p

    return res


def std_vol_under_mask(volume, mask, meanV):
    """Calculate the STD volume under the mask.

    @param volume: input volume.
    @param mask: mask.
    @param meanV: mean volume under mask.
    
    @return: the calculated STD volume under mask.
    """
    # calculate the square of the volume
    copyV = volume**2
    copyMean = meanV**2

    result = mean_vol_under_mask(copyV, mask) - copyMean
    
    result = np.abs(result)**0.5

    return result


def FLCF(volume, template, mask=None, stdV=None):
    '''Fast local correlation function
    
    @param volume: target volume
    @param template: template to be searched. It can have smaller size then target volume.
    @param mask: template mask. If not given, a default sphere mask will be used.
    @param stdV: standard deviation of the target volume under mask, which do not need to be calculated again when the mask is identical.
    
    @return: the local correlation function
    '''
    if volume.shape[0] < template.shape[0] or volume.shape[1] < template.shape[1] or volume.shape[2] < template.shape[2]:
        raise Exception('Template size is bigger than the target volume!')

    # generate the mask 
    if mask is None:
        from tools import create_sphere
        mask = create_sphere(template.shape)
    else:
        if template.shape[0]!=mask.shape[0] and template.shape[1]!=mask.shape[1] and template.shape[2]!=mask.shape[2]:
            raise Exception('Template and mask sizes are not the same!')

    # normalize the template under mask
    meanT = mean_under_mask(template, mask)
    stdT = std_under_mask(template, mask, meanT)
    
    temp = (template - meanT)/stdT
    temp = temp * mask

    # construct both the template and the mask which has the same size as target volume
    from tools import paste_in_center
    tempV = temp
    if volume.shape[0] != temp.shape[0] or volume.shape[1] != temp.shape[1] or volume.shape[2] != temp.shape[2]:
        tempV = np.zeros(volume.shape)
        tempV = paste_in_center(temp, tempV)
    
    maskV = mask
    if volume.shape[0] != mask.shape[0] or volume.shape[1] != mask.shape[1] or volume.shape[2] != mask.shape[2]:
        maskV = np.zeros(volume.shape)
        maskV = paste_in_center(mask, maskV)
    
    # calculate the mean and std of volume
    meanV = mean_vol_under_mask(volume, maskV)
    stdV = std_vol_under_mask(volume, maskV, meanV)
    
    from transform import rfft, irfft, fftshift
    size = volume.shape
    fT = rfft(tempV)
    fT = np.conjugate(fT)
    result = fftshift(irfft(fT*rfft(volume), size))/stdV
    
    return result/np.sum(mask)


def find_peak_position(score, peak_prior=None):
    if peak_prior is not None:
        sc = score * peak_prior
    else:
        sc = score

    idx = np.where(sc == sc.max())

    # return the first occurance of the maximum
    return idx[0][0], idx[1][0], idx[2][0]


