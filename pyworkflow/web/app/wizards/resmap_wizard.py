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

import os
from os.path import basename
import xmipp

from pyworkflow.em.wizard import EmWizard
from django.shortcuts import render_to_response
from django.http import HttpResponse

from pyworkflow.em.packages.resmap.wizard import ResmapPrewhitenWizard, INSTRUCTIONS
from pyworkflow.web.app.em_wizard import *
from tools import *
from pyworkflow.web.app.views_base import base_wiz
from pyworkflow.web.app.views_util import savePlot, loadProject, getResourceCss, getResourceJs
from pyworkflow.gui.plotter import Plotter
from pyworkflow.em.packages.resmap import ProtResMap
from pyworkflow.em.packages.resmap.wizard import myCreatePrewhiteningFigure


class ResmapPrewhitenWizardWeb(ResmapPrewhitenWizard):
    _environments = [WEB_DJANGO]
    
    def _run(self, protocol, request):
        
        if not protocol.inputVolume.get().hasValue():
            return HttpResponse("errorInput")
        elif protocol.useSplitVolume == True:
            if not protocol.splitVolume.get().hasValue():
                return HttpResponse("errorInput")
        else:
            plotUrl, min_ang = getPlotResMap(request, protocol)

            ang = protocol.prewhitenAng.get()
            if ang < min_ang:
                ang = min_ang
            
            context = {
                       # Params
                       'inputId': protocol.inputVolume.get().getObjId(),
                       'ang': round(ang, 2),
                       'ramp': protocol.prewhitenRamp.get(),
                       # Extra Params
                       'pValue': protocol.pVal.get(),
                       'minRes': protocol.minRes.get(), 
                       'maxRes': protocol.maxRes.get(),
                       'stepRes': protocol.stepRes.get(),
                       # Others
                       'plot': plotUrl,
                       'min_ang': round(min_ang, 2),
                       'messi_js': getResourceJs('messi'),
                       'messi_css': getResourceCss('messi'),
                       }
            
            if protocol.useSplitVolume == True and protocol.splitVolume.get().hasValue():
                context.update({'useSplit' : 1,
                                'splitId' : protocol.splitVolume.getObjId(),
                                })
           
            if protocol.applyMask == True and protocol.maskVolume.get().hasValue():
                context.update({'useMask' : 1,
                                'maskId' : protocol.maskVolume.getObjId(),
                                }) 
            
            context = base_wiz(request, context)
            return render_to_response('wizards/wiz_resmap.html', context)

#===============================================================================
# RESMAP REQUEST UPDATE 
#===============================================================================

def get_resmap_plot(request):
    # LOAD Project
    projectName = request.session['projectName']
    project = loadProject(request)
    
    # Create protocol
    newProtocol = project.newProtocol(ProtResMap)
    
    # GET and SET Parameters
    newAng = request.GET.get('ang')
    newProtocol.prewhitenAng.set(float(newAng))
  
    newRamp = request.GET.get('ramp')
    newProtocol.prewhitenRamp.set(float(newRamp))
    
    inputId = request.GET.get('inputId')
    inputVolume = project.mapper.selectById(inputId)
    newProtocol.inputVolume.set(inputVolume)
    
    useSplit = request.GET.get('useSplit', None)
    if useSplit is not None:
        newProtocol.useSplitVolume = True
        splitId = request.GET.get('splitId', None)
        splitVolume = project.mapper.selectById(splitId)
        newProtocol.splitVolume.set(splitVolume.get())
        
    useMask = request.GET.get('useMask', None)
    if useSplit is not None:
        newProtocol.applyMask = True
        maskId = request.GET.get('maskId',  None)
        inputMask = project.mapper.selectById(maskId)
        newProtocol.maskVolume.set(inputMask.get())
        
    # Extra Params
    pValue = request.GET.get('pValue', None)
    if pValue is not None:
        newProtocol.pVal.set(float(pValue))
        
    minRes = request.GET.get('minRes', None)
    if minRes is not None:
        newProtocol.minRes.set(float(minRes))
        
    maxRes = request.GET.get('maxRes', None)
    if maxRes is not None:
        newProtocol.maxRes.set(float(maxRes))
        
    stepRes = request.GET.get('stepRes', None)
    if stepRes is not None:
        newProtocol.stepRes.set(float(stepRes))

    plotUrl, _ = getPlotResMap(request, newProtocol)
    
    return HttpResponse(plotUrl, mimetype='application/javascript')   

#===============================================================================
# RESMAP UTILS 
#===============================================================================

def getPlotResMap(request, protocol):
     #1st step: Convert input volumes
    results = _beforePreWhitening(protocol, os.getcwd())
    min_ang = 2.1* results['vxSize']
    #2nd step: Generate plot
    plotter = _runPreWhiteningWeb(protocol, results)
    #3rd step: Save to png image
    
#     plotUrl = "/" + savePlot(request, plotter)
    plotUrl = savePlot(request, plotter) 
    return plotUrl, min_ang

def _beforePreWhitening(protocol, dir):
    from pyworkflow.em.convert import ImageHandler
    # Convert input volumes
    ih = ImageHandler()
    inputVolume = protocol.inputVolume.get()
    path = join(dir, 'volume1.map')
    print path
    ih.convert(inputVolume, path)
    if protocol.useSplitVolume:
        ih.convert(protocol.splitVolume.get(), join(dir, 'volume2.map'))
    
    return protocol.runResmap(dir, wizardMode=True)

     
def _runPreWhiteningWeb(protocol, results):
    
    newAng = protocol.prewhitenAng.get()
    newRamp = protocol.prewhitenRamp.get()
    
    figsize = (18, 9)
    gridsize = [0, 0]
    plotter = Plotter(*gridsize, figsize=figsize, windowTitle="ResMap")
    
    myCreatePrewhiteningFigure(results, plotter.getFigure(),
                             newAng, newRamp
                            )
    
    return plotter
