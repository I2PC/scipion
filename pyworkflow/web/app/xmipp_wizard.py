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

import em_wizard
from pyworkflow.em.packages.xmipp3 import *


def wiz_xmipp_downsampling(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputMicrographs
    protParams['label']= "downFactor"
    protParams['value']= protocol.downFactor.get()
    return em_wizard.wiz_downsampling(protocol, protParams, request)
    
    
def wiz_xmipp_ctf(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputMicrographs
    protParams['label']= ["lowRes", "highRes"]
    protParams['value']= [protocol.lowRes.get(), protocol.highRes.get()]
    return em_wizard.wiz_ctf(protocol, protParams, request)


def wiz_xmipp_particle_mask_radius(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= "radius"
    protParams['value']= protocol.radius.get()
    return em_wizard.wiz_particle_mask_radius(protocol, protParams, request)


def wiz_xmipp_volume_mask_radius(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputVolumes
    protParams['label']= "radius"
    protParams['value']= protocol.radius.get()
    return em_wizard.wiz_volume_mask_radius(protocol, protParams, request)


def wiz_xmipp_particle_mask_radii(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= ["innerRadius", "outerRadius"]
    protParams['value']= [protocol.innerRadius.get(), protocol.outerRadius.get()]
    return em_wizard.wiz_particles_mask_radii(protocol, protParams, request)


def wiz_xmipp_filter_particle(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= ["lowFreq", 
                          "highFreq", 
                          "freqDecay"]
    protParams['value']= [protocol.lowFreq.get(), 
                          protocol.highFreq.get(), 
                          protocol.freqDecay.get()]
    protParams['mode']= protocol.fourierMode.get()
    return em_wizard.wiz_filter_particle(protocol, protParams, request)


def wiz_xmipp_gaussian_particle(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= "freqSigma"
    protParams['value']= protocol.freqSigma.get()
    return em_wizard.wiz_particle_gaussian(protocol, protParams, request)


def wiz_xmipp_filter_volume(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputParticles
    protParams['label']= ["lowFreq", 
                          "highFreq", 
                          "freqDecay"]
    protParams['value']= [protocol.lowFreq.get(), 
                          protocol.highFreq.get(), 
                          protocol.freqDecay.get()]
    protParams['mode']= protocol.fourierMode.get()
    return em_wizard.wiz_filter_volume(protocol, protParams, request)


def wiz_xmipp_gaussian_volume(protocol, request):
    protParams = {}
    protParams['input']= protocol.inputVolumes
    protParams['label']= "freqSigma"
    protParams['value']= protocol.freqSigma.get()
    return em_wizard.wiz_volume_gaussian(protocol, protParams, request)



    
    
    
