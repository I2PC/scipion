# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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


from pyworkflow.em.viewer import CommandView, Viewer, DESKTOP_TKINTER
from pyworkflow.protocol.params import (LabelParam, StringParam, EnumParam,
                                        IntParam, LEVEL_ADVANCED)
from pyworkflow.viewer import ProtocolViewer
from pyworkflow.em.constants import (COLOR_CHOICES, COLOR_OTHER, COLOR_JET, 
                                     COLOR_TERRAIN, COLOR_GIST_EARTH, 
                                     COLOR_GIST_NCAR, COLOR_GNU_PLOT,
                                      COLOR_GNU_PLOT2, AX_X, AX_Y, AX_Z)
from pyworkflow.gui.plotter import Plotter
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.em import ImageHandler, ChimeraView
from pyworkflow.em.data import Volume
#import numpy as np
#import matplotlib.pyplot as plt
from matplotlib import cm
from convert import getEnviron
from pyworkflow.em.viewer import DataView, LocalResolutionViewer
from protocol_blocres import BsoftProtBlocres, FN_RESOLMAP, FN_HALF1


#------------------------ Some views and  viewers ------------------------
        

class BsoftVolumeView(CommandView):
    def __init__(self, inputFile, **kwargs):

        # Small trick to handle .vol Spider volumes
        if inputFile.endswith('.vol'):
            inputFile += ":spi"

        CommandView.__init__(self, 'bshow "%s" &' % inputFile,
                             env=getEnviron(), **kwargs)

             
class BsoftViewer(Viewer):
    _environments = [DESKTOP_TKINTER]
    _targets = [Volume]
    
    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    def visualize(self, obj, **kwargs):
        cls = type(obj)

        fn = obj.getFileName()
        BsoftVolumeView(fn).show()

class BsoftPlotter(EmPlotter):
    """ Class to create several plots with Bsoft utilities"""
    def __init__(self, x=1, y=1, mainTitle="", **kwargs):
        EmPlotter.__init__(self, x, y, mainTitle, **kwargs)
 
binaryCondition = ('(colorMap == %d) ' % (COLOR_OTHER))

class BsoftViewerBlocres(LocalResolutionViewer):
    """
    Visualization tools for blocres results.

    blocres is a bsoft package for computing the local resolution of 3D
    density maps studied in structural biology, primarily by cryo-electron
    microscopy (cryo-EM).    
    """
    _label = 'viewer blocres'
    _targets = [BsoftProtBlocres]
    _environments = [DESKTOP_TKINTER]
    
    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, *args, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        
        form.addParam('doShowVolumeSlices', LabelParam,
                      label="Show resolution slices")
        
        form.addParam('doShowOriginalVolumeSlices', LabelParam,
                      label="Show original volume slices")

        form.addParam('doShowResHistogram', LabelParam,
                      label="Show resolution histogram")
        
        group = form.addGroup('Colored resolution Slices and Volumes')
        group.addParam('colorMap', EnumParam, choices=COLOR_CHOICES,
                      default=COLOR_JET,
                      label='Color map',
                      help='Select the color map to apply to the resolution map. '
                      'http://matplotlib.org/1.3.0/examples/color/colormaps_reference.html.')
        
        group.addParam('otherColorMap', StringParam, default='jet',
                      condition = binaryCondition,
                      label='Customized Color map',
                      help='Name of a color map to apply to the resolution map.'
                      ' Valid names can be found at '
                      'http://matplotlib.org/1.3.0/examples/color/colormaps_reference.html')
        group.addParam('sliceAxis', EnumParam, default=AX_Z,
                       choices=['x', 'y', 'z'],
                       display=EnumParam.DISPLAY_HLIST,
                       label='Slice axis')
        group.addParam('doShowVolumeColorSlices', LabelParam,
                       label="Show colored slices")
        group.addParam('doShowOneColorslice', LabelParam, 
                       expertLevel=LEVEL_ADVANCED, 
                      label='Show selected slice')
        group.addParam('sliceNumber', IntParam, default=-1,
                       expertLevel=LEVEL_ADVANCED, 
                       label='Show slice number')
        group.addParam('doShowChimera', LabelParam,
                      label="Show Resolution map in Chimera")


    def _getVisualizeDict(self):
        self.protocol._createFilenameTemplates()
        return {'doShowOriginalVolumeSlices': self._showOriginalVolumeSlices,
                'doShowVolumeSlices': self._showVolumeSlices,
                'doShowVolumeColorSlices': self._showVolumeColorSlices,
                'doShowOneColorslice': self._showOneColorslice,
                'doShowResHistogram': self._plotHistogram,
                'doShowChimera': self._showChimera,
                }
       
    def _showVolumeSlices(self, param=None):
        out_cm = DataView(self.protocol._getFileName(FN_RESOLMAP))
        
        return [out_cm]


    def _showOriginalVolumeSlices(self, param=None):
        out_cm = DataView(self.protocol.inputVolume.get().getFileName())

        return [out_cm]


    def _showVolumeColorSlices(self, param=None):
        imageFile = self.protocol._getFileName(FN_RESOLMAP)
        imgData, min_Res, max_Res = self.getImgData(imageFile)

        xplotter = BsoftPlotter(x=2, y=2, mainTitle="Local Resolution Slices "
                                                     "along %s-axis."
                                                     %self._getAxis())
        #The slices to be shown are close to the center. Volume size is divided in 
        # 9 segments, the fouth central ones are selected i.e. 3,4,5,6
        for i in xrange(3,7): 
            sliceNumber = self.getSlice(i, imgData)
            a = xplotter.createSubPlot("Slice %s" % (sliceNumber), '', '')
            matrix = self.getSliceImage(imgData, sliceNumber, self._getAxis())
            plot = xplotter.plotMatrix(a, matrix, min_Res, max_Res,
                                       cmap=self.getColorMap(),
                                       interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]


    def _showOneColorslice(self, param=None):
        imageFile = self.protocol._getFileName(FN_RESOLMAP)
        imgData, min_Res, max_Res = self.getImgData(imageFile)

        xplotter = BsoftPlotter(x=1, y=1, mainTitle="Local Resolution Slices "
                                                     "along %s-axis."
                                                     %self._getAxis())
        sliceNumber = self.sliceNumber.get()
        if sliceNumber < 0:
            x ,_ ,_ ,_ = ImageHandler().getDimensions(imageFile)
            sliceNumber = x/2
        else:
            sliceNumber -= 1
        a = xplotter.createSubPlot("Slice %s" % (sliceNumber+1), '', '')
        matrix = self.getSliceImage(imgData, sliceNumber, self._getAxis())
        plot = xplotter.plotMatrix(a, matrix, min_Res, max_Res,
                                       cmap=self.getColorMap(),
                                       interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]


    def _plotHistogram(self, param=None):
        imageFile = self.protocol._getFileName(FN_RESOLMAP)
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        imgList = imgData.flatten()
        imgListNoZero = filter(lambda x: x > 0, imgList)
        nbins = 30
        plotter = EmPlotter(x=1,y=1,mainTitle="  ")
        plotter.createSubPlot("Resolution histogram",
                              "Resolution (A)", "# of Counts")
        fig = plotter.plotHist(imgListNoZero, nbins)
        return [plotter]
    
    
    def _showChimera(self, param=None):
        fnResVol = self.protocol._getFileName(FN_RESOLMAP)
        fnOrigMap = self.protocol._getFileName(FN_HALF1)
        cmdFile = self.protocol._getExtraPath('chimera_resolution_map.cmd')
        sampRate = self.protocol.resolution_Volume.getSamplingRate()
        self.createChimeraScript(cmdFile, fnResVol, fnOrigMap, sampRate)
        view = ChimeraView(cmdFile)
        return [view]
        

    def _getAxis(self):
        return self.getEnumText('sliceAxis')

    
    def getColorMap(self):
        if (COLOR_CHOICES[self.colorMap.get()] is 'other'): 
            cmap = cm.get_cmap(self.otherColorMap.get())
        else:
            cmap = cm.get_cmap(COLOR_CHOICES[self.colorMap.get()])
        if cmap is None:
            cmap = cm.jet
        return cmap





