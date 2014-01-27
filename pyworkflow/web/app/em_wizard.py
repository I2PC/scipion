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

from pyworkflow.em.packages.xmipp3.convert import xmippToLocation, locationToXmipp
from pyworkflow.em.packages.spider.convert import locationToSpider
from pyworkflow.em.packages.spider.wizard import filter_spider
from pyworkflow.em.packages.xmipp3.constants import *

def wizard(request):
    # Get the Wizard Name
    requestDict = getattr(request, "POST")
    functionName = requestDict.get("wizName")
    function = globals().get(functionName, None)
    
    print "======================= in wizard: " + functionName
    
    # Get the protocol object
    project, protocol = loadProtocolProject(request)
    updateProtocolParams(request, protocol, project)
    
    if function is None:
        pass  # redirect to error: wizard not found
    elif not callable(function):
        pass  # redirect to error: name is not a function
    else:
        return function(protocol, request)
    
def wiz_base(request, context):
    context_base = {'general_style': getResourceCss('general'),
               'wizard_style': getResourceCss('wizard'),
               'jquery_ui_style': getResourceCss('jquery_ui'),
               'font_awesome': getResourceCss('font_awesome'),
               'jquery': getResourceJs('jquery'),
               'jquery_ui': getResourceJs('jquery_ui'),
               'wizard_utils': getResourceJs('wizard_utils'),
               'raphael':getResourceJs('raphael'),
               'projectName': request.session['projectName']
               }
    
    context.update(context_base)
    
    return context


def wiz_downsampling(protocol, request):
    micrographs = protocol.inputMicrographs.get()
    res = validateSet(micrographs)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        for m in micrographs:
            m.basename = basename(m.getFileName())
        
#        mics = [mic for mic in micrographs]
#        for m in mics:
#            m.basename = basename(m.getFileName())
            
        context = {'objects': micrographs,
                   'downFactor': protocol.downFactor.get()
                   }
        
        context = wiz_base(request, context)

        return render_to_response('wizards/wiz_downsampling.html', context)

def wiz_ctf(protocol, request):
    micrographs = protocol.inputMicrographs.get()
    
    res = validateSet(micrographs)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        mics = [mic for mic in micrographs]
        for m in mics:
            m.basename = basename(m.getFileName())    
        
        context = {'objects': mics,
                   'high_res' : protocol.highRes.get(),
                   'low_res': protocol.lowRes.get(),
                   }
        
        context = wiz_base(request, context)
        
        return render_to_response('wizards/wiz_ctf.html', context)

def wiz_particle_mask(protocol, request):
    particles = protocol.inputParticles.get()
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            xdim = getImageXdim(request, parts[0].text)
    
            mask_radius = protocol.maskRadius.get()
             
            if mask_radius > xdim :
                mask_radius = xdim
            elif mask_radius == -1 :
                mask_radius = xdim/2
                
            context = {'objects': parts,
                       'maskRadius': mask_radius,
                       'xdim':xdim,
                       }

            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_particle_mask.html', context)

def wiz_particle_mask_radii(protocol, request):
    particles = protocol.inputParticles.get()
    res = validateParticles(particles) 
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
        if len(parts) == 0:
            return HttpResponse("errorIterate")
        else:
            xdim = getImageXdim(request, parts[0].text)
            inner_radius = protocol.innerRadius.get()
            outer_radius = protocol.outerRadius.get()
                
            context = {'objects': parts,
                       'innerRadius': inner_radius,
                       'outerRadius': outer_radius,
                       'xdim':xdim,
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_particle_mask_radii.html', context)

def wiz_volume_mask(protocol, request):
    volumes = protocol.input3DReferences.get()
    res = validateSet(volumes)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        vols = [vol for vol in volumes]
        
        for v in vols:
            v.basename = basename(v.getFileName())
                    
        xdim = getImageXdim(request, vols[0].getFileName())
    
        mask_radius = protocol.maskRadius.get()
         
        if mask_radius > xdim :
            mask_radius = xdim
        elif mask_radius == -1 :
            mask_radius = xdim/2
            
        context = {'objects': vols,
                   'maskRadius': mask_radius,
                   'xdim': xdim,
                   }
        
        context = wiz_base(request, context)
        
        return render_to_response('wizards/wiz_volume_mask.html', context)

def wiz_volume_mask_radii(protocol, request):
    volumes = protocol.input3DReferences.get()
    res = validateSet(volumes)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        vols = [vol for vol in volumes]
        
        for v in vols:
            v.basename = basename(v.getFileName())
                    
        xdim = getImageXdim(request, vols[0].getFileName())
    
        inner_radius = protocol.innerRadius.get()
        outer_radius = protocol.outerRadius.get()
            
        context = {'objects': vols,
                   'innerRadius': inner_radius,
                   'outerRadius': outer_radius,
                   'xdim': xdim,
                   }
        
        context = wiz_base(request, context)
        
        return render_to_response('wizards/wiz_volume_mask_radii.html', context)

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

def wiz_bandpass(protocol, request):
    particles = protocol.inputParticles.get()
    mode = protocol.fourierMode.get()
    
    if mode == 0:
        # low pass
        lowFreq = 0.
        highFreq = protocol.highFreq.get()
        decay = protocol.freqDecay.get()
    elif mode == 1:
        # high pass
        lowFreq = protocol.lowFreq.get()
        highFreq = 1.0
        decay = protocol.freqDecay.get()
    elif mode== 2:
        #band pass
        highFreq = protocol.highFreq.get()
        lowFreq = protocol.lowFreq.get()
        decay = protocol.freqDecay.get()
    
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
                       'lowFreq': lowFreq,
                       'highFreq': highFreq,
                       'decayFreq': decay,
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_bandpass.html', context)

def wiz_gaussian(protocol, request):
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
                       'freqSigma': protocol.freqSigma.get(),
                       }
            
            context = wiz_base(request, context)
            
            return render_to_response('wizards/wiz_gaussian.html', context)
        
def wiz_relion_bandpass(protocol, request):
    particles = protocol.inputParticles.get()
    
    highFreq = protocol.iniLowPassFilter.get()
    
    res = validateParticles(particles)
    
    if res is not 1:
        return HttpResponse(res)
    else:
        parts = getParticleSubset(particles,100)
        
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
        

def getParticleSubset(particles, num):
    """
    Method to prepare the particles to use a wizard
    """
    particleList = []
    for i, particle in enumerate(particles):
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
    """
    Validation for a set of micrographs
    """
    if setOf is None:
        res = "errorInput"
    else:
        res = 1
    return res

def validateParticles(particles):
    """
    Validation for a set of particles
    """
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
    
