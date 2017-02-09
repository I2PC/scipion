# **************************************************************************
# *
# * Authors:    Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from django.http import HttpResponse

from pyworkflow.gui import getPILImage
import xmipp

from pyworkflow.em.constants import *


#===============================================================================
# UTILS METHODS
#===============================================================================

def proccessModeFilter(mode, value):
    # Order : low - high - decay
    
    # FIX to use in case value is not an array
    if (type(value) is int) or (type(value) is float):
        assign = 0
    else: 
        assign = 1
        

    if mode == FILTER_LOW_PASS:
        print "filter low pass"
#        value[0] = 0.
        value[0] = 0.01
    
    elif mode == FILTER_HIGH_PASS:
        print "filter high pass"
        value[1] = 0.01
    
    elif mode == FILTER_BAND_PASS:
        print "filter band pass"
    
    elif mode == FILTER_LOW_PASS_NO_DECAY:
        print "filter low pass no decay"
        value = [0.01, value, 0.01]
        
    elif mode == FILTER_NO_DECAY:
        print "filter band pass no decay"
        value = [value[0], value[1], 0.01]
        
    return value

def setValueOnRange(valueList, min, max):
    # Check if a value is on a range and if not set it
    for i, val in enumerate(valueList):
        if val < min:
            val = min
        if val > max:
            val = max
        valueList[i] = val

    return valueList

def validateMaskRadius(value, xdim, radius):
    """ Validation done for some radius length """

    # MASK RADIUS
    if radius == 1:
        if value > xdim/2 or value == -1:
            value = xdim/2
    
    # MASK RADII
    elif radius == 2:
        if value[0] > xdim/2 or value[0] == -1:
            value[0] = xdim/2
        if value[1] > xdim/2 or value[1] == -1:
            value[1] = xdim/2
    
    return value
    

def validateSet(setOf):
    """Validation for a set of micrographs."""
    if setOf is None:
        res = "errorInput"
    else:
        res = 1
    return res


def validateParticles(particles):
    """Validation for a set of particles."""
    if particles is None:
        res = "errorInput"
    elif particles.getSize() == 0:
        res = "errorEmpty"
    else:
        res = 1
        
    return res


def get_image_psd(request):
    """
    Function to get the computing psd image
    """
    imagePath = request.GET.get('image', None)
    downsample = request.GET.get('downsample', None)
    dim = request.GET.get('dim', None)
    
    # create a xmipp image empty
    imgXmipp = xmipp.Image()

    #    ** DEBUG **
#    print str(imagePath)
#    print float(downsample) 
#    print int(dim)

    # compute the PSD image
    xmipp.fastEstimateEnhancedPSD(imgXmipp, str(imagePath), 
                                  float(downsample), int(dim), 2)
    
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
       
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response


def get_image_bandpass(request):
    """
    Function to get the computing image with a fourier filter applied
    """
    imagePath = request.GET.get('image', None)
    lowFreq = request.GET.get('lowFreq', 0.)
    highFreq = request.GET.get('highFreq', 1.0)
    decay = request.GET.get('decayFreq', 0.)
    dim = request.GET.get('dim', None)
    
#    ** DEBUG **
#    print str(imagePath)
#    print float(lowFreq)
#    print float(highFreq)
#    print float(decay)
#    print int(dim)
    
    # create a xmipp image empty
    imgXmipp = xmipp.Image()
    
    # compute the Fourier Filter in the image
    xmipp.bandPassFilter(imgXmipp, str(imagePath), float(lowFreq), float(highFreq), float(decay), int(dim))
    
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
        
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    
    return response


def get_image_gaussian(request):
    """
    Function to get the computing image with a gaussian filter applied
    """
    imagePath = request.GET.get('image', None)
    freqSigma = request.GET.get('sigmaFreq', None)
    dim = request.GET.get('dim', None)
    
    # create a xmipp image empty
    imgXmipp = xmipp.Image()
    
    # compute the Gaussian Filter in the image
    xmipp.gaussianFilter(imgXmipp, str(imagePath), float(freqSigma), int(dim))
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
        
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response

def get_image_mask(request):
    """
    Function to get the computing image with a mask filter applied
    """
    imagePath = request.GET.get('image', None)
    dim = request.GET.get('dim', None)
    
    # create a xmipp image empty
    imgXmipp = xmipp.Image()
    
    # from PIL import Image
    img = getPILImage(imgXmipp, dim)
        
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response
