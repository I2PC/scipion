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

from xmipp_wizard import *
#from spider_wizard import *
#from relion_wizard import *


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
    
    # Obtain the parameters for the wizard
    functionName, params = wizClass._run(protocol)
    function = globals().get(functionName, None)
    
    print "======================= in wizard: " + functionName
    
    if function is None:
        pass  # redirect to error: wizard not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        return function(protocol, params, request)

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
               'projectName': request.session['projectName']
               }
    
    context.update(context_base)
    return context

#===============================================================================
#    Wizard classes (to execute differents wizards)
#===============================================================================

def wiz_downsampling(protocol, params, request):
    
    # Get params
    micrographs = params['input'].get()
    label = params['label']
    value = params['value']
    
    res = validateSet(micrographs)
    
    if res is not 1:
        return HttpResponse(res)
    else:           
        mics = [mic.clone() for mic in micrographs]
        for m in mics:
            m.basename = basename(m.getFileName())
             
        context = {'objects': mics,
                   label : value
                   }

        context = wiz_base(request, context)
        return render_to_response('wizards/wiz_downsampling.html', context)


def wiz_ctf(protocol, params, request):
    
    # Get params
    micrographs = params['input'].get()
    label = params['label']
    value = params['value']

    res = validateSet(micrographs)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        mics = [mic.clone() for mic in micrographs]
        for m in mics:
            m.basename = basename(m.getFileName())    
        
        context = {'objects': mics,
                   label[0]: value[0],
                   label[1] : value[1],
                   }
        
        context = wiz_base(request, context)
        return render_to_response('wizards/wiz_ctf.html', context)

def wiz_particle_mask_radius(protocol, params, request):
    
    # Get params
    particles = params['input'].get()
    label = params['label']
    value = params['value']
    
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            xdim = getImageXdim(request, parts[0].text)
    
            mask_radius = value
             
            if mask_radius > xdim :
                mask_radius = xdim
            elif mask_radius == -1 :
                mask_radius = xdim/2
                
            context = {'objects': parts,
                       label: mask_radius,
                       'xdim':xdim,
                       }

            context = wiz_base(request, context)
            return render_to_response('wizards/wiz_particle_mask_radius.html', context)

def wiz_particles_mask_radii(protocol, params, request):
    
    # Get params
    particles = params['input'].get()
    label = params['label']
    value = params['value']
    
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            xdim = getImageXdim(request, parts[0].text)
                
            context = {'objects': parts,
                       label[0]: value[0],
                       label[1]: value[1],
                       'xdim':xdim,
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_particles_mask_radii.html', context)

def wiz_volume_mask_radius(protocol, params, request):
    
    # Get params
    volumes = params['input'].get()
    label = params['label']
    value = params['value']
    
    res = validateSet(volumes)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        vols = []
        if isinstance(volumes, Volume):
            vols.append(volumes)
        else: 
            vols = [vol.clone() for vol in volumes]
        
        for v in vols:
            v.basename = basename(v.getFileName())
                    
        xdim = getImageXdim(request, vols[0].getFileName())
    
        mask_radius = value
         
        if mask_radius > xdim :
            mask_radius = xdim
        elif mask_radius == -1 :
            mask_radius = xdim/2
            
        context = {'objects': vols,
                   label: mask_radius,
                   'xdim': xdim,
                   }
        
        context = wiz_base(request, context)
        
        return render_to_response('wizards/wiz_volume_mask_radius.html', context)

def wiz_volumes_mask_radii(protocol, params, request):

    # Get params
    volumes = params['input'].get()
    label = params['label']
    value = params['value']
    
    res = validateSet(volumes)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        vols = []
        if isinstance(volumes, Volume):
            vols.append(volumes)
        else: 
            vols = [vol.clone() for vol in volumes]
        
        for v in vols:
            v.basename = basename(v.getFileName())
                    
        xdim = getImageXdim(request, vols[0].getFileName())
            
        context = {'objects': vols,
                   label[0]: value[0],
                   label[1]: value[1],
                   'xdim': xdim,
                   }
        
        context = wiz_base(request, context)
        
        return render_to_response('wizards/wiz_volumes_mask_radii.html', context)

def wiz_filter_particle(protocol, params, request):
    
    # Get params
    particles = params['input'].get()
    label = params['label']
    value = params['value']
    mode = params['mode']
    
    if mode == 0:
        # low pass
        lowFreq = 0.
        highFreq = value[1]
        decay = value[2]
    elif mode == 1:
        # high pass
        lowFreq = value[0]
        highFreq = 1.0
        decay = value[2]
    elif mode== 2:
        #band pass
        highFreq = value[0]
        lowFreq = value[1]
        decay = value[2]
    
    res = validateParticles(particles)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            context = {'objects': parts,
                       'mode': mode,
                       label[0]: lowFreq,
                       label[1]: highFreq,
                       label[2]: decay,
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_filter_particles.html', context)
        
def wiz_filter_volume(protocol, params, request):
    
    # Get params
    volumes = params['input'].get()
    label = params['label']
    value = params['value']
    mode = params['mode']
    
    if mode == 0:
        # low pass
        lowFreq = 0.
        highFreq = value[1]
        decay = value[2]
    elif mode == 1:
        # high pass
        lowFreq = value[0]
        highFreq = 1.0
        decay = value[2]
    elif mode== 2:
        #band pass
        highFreq = value[0]
        lowFreq = value[1]
        decay = value[2]
    
    res = validateSet(volumes)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        vols = []
        if isinstance(volumes, Volume):
            vols.append(volumes)
        else: 
            vols = [vol.clone() for vol in volumes]
        
        for v in vols:
            v.basename = basename(v.getFileName())
                
        if len(vols) == 0:
            return HttpResponse("errorIterate")
        else:
            context = {'objects': vols,
                       'mode': mode,
                       label[0]: lowFreq,
                       label[1]: highFreq,
                       label[2]: decay,
                       'unit': UNIT_PIXEL
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_filter_volumes.html', context)


def wiz_particle_gaussian(protocol, params, request):
    
    # Get params
    particles = params['input'].get()
    label = params['label']
    value = params['value']
    
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles.clone(),100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            context = {'objects': parts,
                       label: value,
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_gaussian_particle.html', context)
        
        
def wiz_volume_gaussian(protocol, params, request):
    
    # Get params
    volumes = params['input'].get()
    label = params['label']
    value = params['value']
    
    res = validateSet(volumes)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        vols = []
        if isinstance(volumes, Volume):
            vols.append(volumes)
        else: 
            vols = [vol.clone() for vol in volumes]
        
        for v in vols:
            v.basename = basename(v.getFileName())
        
        if len(vols) == 0:
            return HttpResponse("errorIterate")
        else:
            context = {'objects': vols,
                       label: value,
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_gaussian_vol.html', context)
    
    
#===============================================================================
#    Methods for utils
#===============================================================================

def getParticleSubset(particles, num):
    """
    Method to prepare the particles to use a wizard
    """
    particleList = []
    for i, part in enumerate(particles):
        
        # Cloning the particle
        particle = part.clone() 
        
        if i == num: # Only load up to NUM particles
            break
        particle.text = particle.getFileName()
        particle.basename = basename(particle.text)
        
        index = particle.getIndex()
        if index:
            particle.text = "%03d@%s" % (index, particle.text)
            particle.basename = "%03d@%s" % (index, particle.basename)
        
        particleList.append(particle)

    return particleList

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
        
    response = HttpResponse(mimetype="imag($(this))e/png")
    img.save(response, "PNG")
    return response


def wiz_filter_spider(protocol, request):
    particles = protocol.inputParticles.get()
    
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            context = {'objects': parts,
                       'filterType': protocol.filterType.get(),
                       'filterMode': protocol.filterMode.get(),
                       'usePadding': protocol.usePadding.get(),
                       'protocol': protocol,
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_filter_spider.html', context)
       
       
#========================================================================
# SPIDER
#========================================================================

def get_image_filter_spider(request):
    """
    Function to get the computing image with a spider filter applied
    """
    pars={}
    
    imagePath = request.GET.get('image', None)
    dim = request.GET.get('dim', None)
    filterType = int(request.GET.get('filterType', None))
    pars["filterType"] = filterType
    pars["filterMode"] = int(request.GET.get('filterMode', None))
    pars["usePadding"] = request.GET.get('usePadding', None)    
    pars["op"]="FQ"
    
    # Copy image to filter to Tmp project folder
    outputName = os.path.join("Tmp", "filtered_particle")
    outputPath = outputName + ".spi"
    
    outputLoc = (1, outputPath)
    ih = ImageHandler()
    ih.convert(xmippToLocation(imagePath), outputLoc)
    outputLocSpiStr = locationToSpider(1, outputName)
    
    # check values
    if filterType < 2:
        pars["filterRadius"] = request.GET.get('radius', None)
    else: 
        pars["highFreq"] = float(request.GET.get('highFreq', None))
        pars["lowFreq"] = float(request.GET.get('lowFreq', None))
        if filterType == 2:
            pars["temperature"] = request.GET.get('temperature', None)    

    filter_spider(outputLocSpiStr, outputLocSpiStr, **pars)
    
    # Get output image and update filtered image
    img = xmipp.Image()
    locXmippStr = locationToXmipp(1, outputPath)
    img.read(locXmippStr)
    
    # from PIL import Image
    img = getPILImage(img, dim)
    
    response = HttpResponse(mimetype="image/png")
    img.save(response, "PNG")
    return response

def wiz_spider_particle_mask_radius(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= "radius"
    protParams['value']= protocol.radius.get()
#    return em_wizard.wiz_particle_mask_radius(protocol, protParams, request)
    return wiz_particle_mask_radius(protocol, protParams, request)


def wiz_spider_particle_mask_radii(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= ["innerRadius", "outerRadius"]
    protParams['value']= [protocol.innerRadius.get(), protocol.outerRadius.get()]
#    return em_wizard.wiz_particles_mask_radii(protocol, protParams, request)
    return wiz_particles_mask_radii(protocol, protParams, request)


#===============================================================================
# RELION
#===============================================================================

def wiz_relion_bandpass(protocol, request):
    particles = protocol.inputParticles.get()
    
    highFreq = protocol.iniLowPassFilter.get()
    
    res = validateParticles(particles)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles.clone(),100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            context = {'objects': parts,
                       'lowFreq': 0,
                       'highFreq': highFreq,
                       'decayFreq': 0,
                       'unit': UNIT_PIXEL
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_relion_bandpass.html', context)
    
