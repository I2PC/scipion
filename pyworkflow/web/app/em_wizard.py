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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import basename
import xmipp
from views_util import * 
from views_protocol import updateProtocolParams
from pyworkflow.manager import Manager
from pyworkflow.project import Project
from django.shortcuts import render_to_response
from pyworkflow.gui import getImage, getPILImage
from django.http import HttpResponse
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.constants import *

from pyworkflow.em import SetOfImages, SetOfMicrographs, Volume, SetOfParticles, SetOfVolumes, ProtCTFMicrographs
from pyworkflow.em.wizard import EmWizard

from pyworkflow.em.packages.xmipp3.convert import xmippToLocation, locationToXmipp
from pyworkflow.em.packages.spider.convert import locationToSpider
from pyworkflow.em.packages.spider.wizard import filter_spider
from pyworkflow.em.packages.xmipp3.constants import *


# Imports for web wizards
from pyworkflow.web.app.wizards.xmipp_wizard import *
from pyworkflow.web.app.wizards.spider_wizard import *
from pyworkflow.web.app.wizards.relion_wizard import *


#===============================================================================
#    Wizard base function (to call the others)
#===============================================================================

def wizard(request):
    # Get the Wizard Name
    requestDict = getattr(request, "POST")
#    functionName = requestDict.get("wizName")

    # Get and instance the wizard class
    className = requestDict.get("wizClassName")
    wizClass = globals().get(className, None)()

    # Get the protocol object
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)

    print "======================= in wizard: " + str(wizClass)
    
    # Obtain the parameters for the wizard
    return wizClass._run(protocol, request)


#===============================================================================
#    Wizard common resources (to build the context)
#===============================================================================

def wiz_base(request, context):
    context_base = {'general_style': getResourceCss('general'),
               'wizard_style': getResourceCss('wizard'),
               'jquery_ui_style': getResourceCss('jquery_ui'),
               'font_awesome': getResourceCss('font_awesome'),
               'jquery': getResourceJs('jquery'),
               'jquery_ui': getResourceJs('jquery_ui'),
               'jquery_ui_touch': getResourceJs('jquery_ui_touch'),
               'wizard_utils': getResourceJs('wizard_utils'),
               'raphael':getResourceJs('raphael'),
               'projectName': request.session['projectName'],
               'loading' : getResourceIcon('loading'),
               }
    
    context.update(context_base)
    return context


#===============================================================================
#    Methods for utils
#===============================================================================

def proccessModeFilter(mode, value):
    # Order : low - high - decay
    
    if mode == 0:
        # low pass
        value[0] = 0.
    elif mode == 1:
        # high pass
        value[1] = 1.0
    elif mode== 2:
        #band pass
        pass
    
#     if mode == 0:
#        # low pass
#        lowFreq = 0.
#        highFreq = value[1]
#        decay = value[2]
#    elif mode == 1:
#        # high pass
#        lowFreq = value[0]
#        highFreq = 1.0
#        decay = value[2]
#    elif mode== 2:
#        #band pass
#        highFreq = value[0]
#        lowFreq = value[1]
#        decay = value[2]
#    
    return value

def validateMaskRadius(value, xdim, radius):
    if radius == 1:
        if value > xdim :
            value = xdim
        elif value == -1 :
            value = xdim/2
    elif radius == 2:
        if value[0] > xdim :
            value[0] = xdim
        elif value[1] > xdim :
            value[1] = xdim
        elif value[0] == -1 :
            value[0] = xdim/2
        elif value[1] == -1 :
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
#        else:
#            res = parts
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

    # compute the PSD image
    xmipp.fastEstimateEnhancedPSD(imgXmipp, str(imagePath), float(downsample), int(dim), 2)
    
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

