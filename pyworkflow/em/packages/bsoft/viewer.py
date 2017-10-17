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
from pyworkflow.protocol.params import LabelParam, StringParam, EnumParam
from pyworkflow.viewer import ProtocolViewer
from pyworkflow.em.constants import *
from pyworkflow.gui.plotter import Plotter
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.em import ImageHandler
from pyworkflow.em.data import Volume
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from convert import getEnviron
from pyworkflow.em.viewer import ChimeraView, DataView, LocalResolutionViewer
from protocol_blocres import BsoftProtBlocres, FN_RESOLMAP


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

class BsoftViewerBlocres(LocalResolutionViewer, BsoftProtBlocres):
    """
    Visualization tools for blocres results.

    blocres is a bsoft package for computing the local resolution of 3D
    density maps studied in structural biology, primarily by cryo-electron
    microscopy (cryo-EM).    
    """
    _label = 'viewer blocres'
    _targets = [BsoftProtBlocres]
    _environments = [DESKTOP_TKINTER]
    
    @staticmethod
    def getColorMapChoices():
        return plt.colormaps()
   
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

    def _getVisualizeDict(self):
        self.protocol._createFilenameTemplates()
        return {'doShowOriginalVolumeSlices': self._showOriginalVolumeSlices,
                'doShowVolumeSlices': self._showVolumeSlices,
                'doShowVolumeColorSlices': self._showVolumeColorSlices,
                'doShowResHistogram': self._plotHistogram,
                }
       
    def _showVolumeSlices(self, param=None):
        cm = DataView(self.protocol.outputVolume.getFileName())
        
        return [cm]
    
    def _showOriginalVolumeSlices(self, param=None):
        cm = DataView(self.protocol.inputVolume.get().getFileName())
        cm2 = DataView(self.protocol.inputVolume2.get().getFileName())

        return [cm]
    
    def _showVolumeColorSlices(self, param=None):
        import os
        imageFile = self.protocol._getFileName(FN_RESOLMAP)
        imgData, min_Res, max_Res = self.getImgData(imageFile)

        xplotter = BsoftPlotter(x=2, y=2, mainTitle="Local Resolution Slices "
                                                     "along %s-axis."
                                                     %self._getAxis())
        #The slices to be shown are close to the center. Volume size is divided in 
        # 9 segments, the fouth central ones are selected i.e. 3,4,5,6
        for i in xrange(3,7): 
            slice = self.getSlice(i, imgData)
            a = xplotter.createSubPlot("Slice %s" % (slice), '', '')
            matrix = self.getSliceImage(imgData, i, self._getAxis())
            plot = xplotter.plotMatrix(a, matrix, min_Res, max_Res,
                                       cmap=self.getColorMap(),
                                       interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]


    def _plotHistogram(self, param=None):
        imageFile = self.protocol._getFileName(FN_RESOLMAP)
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        Res = imgData[imgData!=self.protocol.fill.get()]
        minres = np.amin(Res)
        maxres = np.amax(Res)
        steps = 30
        stepSize = (maxres - minres)/steps
        x_axis = []
        y_axis = []
    
        for idx in range(0,steps):
            x_axis_aux = minres + idx*stepSize
            x_axis.append(x_axis_aux)
            auxRes = Res[Res<(minres+(idx+1)*stepSize)]
            auxRes2 = auxRes[auxRes>=(minres+(idx*stepSize))]
            length = len(auxRes2)
            y_axis_aux = int(length)
            y_axis.append(y_axis_aux)
        fig = plt.figure()        
        plt.bar(x_axis, y_axis, width = stepSize)    
        plt.title("Resolutions Histogram")
        plt.xlabel("Resolution (A)")
        plt.ylabel("Counts")    

        return [Plotter(figure = fig)]

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





