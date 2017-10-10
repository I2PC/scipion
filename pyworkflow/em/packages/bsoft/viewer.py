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
from pyworkflow.em import ImageHandler
from pyworkflow.em.data import Volume
import numpy as np
from protocol_blocres import BsoftProtBlocres, OUTPUT_RESOLUTION_FILE
import matplotlib.pyplot as plt
from matplotlib import cm
from convert import getEnviron
from pyworkflow.em.viewer import ChimeraView, DataView, LocalResolutionViewer


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
                      help='Name of a color map to apply to the resolution map. Valid names can be found at '
                            'http://matplotlib.org/1.3.0/examples/color/colormaps_reference.html')
        group.addParam('sliceAxis', EnumParam, default=AX_Z,
                       choices=['x', 'y', 'z'],
                       display=EnumParam.DISPLAY_HLIST,
                       label='Slice axis')

        group.addParam('doShowVolumeColorSlices', LabelParam,
              label="Show colored slices")


        
    def _getVisualizeDict(self):
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
        imageFile = self.protocol._getExtraPath(OUTPUT_RESOLUTION_FILE)
        imgData2, min_Res, max_Res = self.getImgData(imageFile)

        fig, im = self._plotVolumeSlices('blocres slices', imgData2,
                                         min_Res, max_Res, self.getColorMap(), dataAxis=self._getAxis())
        cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(im, cax=cax)
        cbar.ax.invert_yaxis()

        return [Plotter(figure = fig)]


    def _plotHistogram(self, param=None):
        imageFile = self.protocol._getExtraPath(OUTPUT_RESOLUTION_FILE)
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

    def _plotVolumeSlices(self, title, volumeData, vminData, vmaxData, cmap, **kwargs):
        """ Helper function to create plots of volumes slices. 
        Params:
            title: string that will be used as title for the figure.
            volumeData: numpy array representing a volume from where to take the slices.
            cmap: color map to represent the slices.
        """
        # Get some customization parameters, by providing with default values
        titleFontSize = kwargs.get('titleFontSize', 14)
        titleColor = kwargs.get('titleColor','#104E8B')
        sliceFontSize = kwargs.get('sliceFontSize', 10)
        sliceColor = kwargs.get('sliceColor', '#104E8B')
        size = kwargs.get('n', volumeData.shape[0])
        origSize = kwargs.get('orig_n', size)
        dataAxis = kwargs.get('dataAxis', 'z')
    
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        f.suptitle(title, fontsize=titleFontSize, color=titleColor, fontweight='bold')
    
        def showSlice(ax, index):
            sliceTitle = 'Slice %s' % int(index*size/9)
            ax.set_title(sliceTitle, fontsize=sliceFontSize, color=sliceColor)
            return ax.imshow(self.getSliceImage(volumeData, index, dataAxis), vmin=vminData, vmax=vmaxData,
                             cmap=self.getColorMap(), interpolation="nearest")
        
        im = showSlice(ax1, 3)
        showSlice(ax2, 4)
        showSlice(ax3, 5)
        showSlice(ax4, 6)
        
        return f, im 
    
   
    def getColorMap(self):
        if (COLOR_CHOICES[self.colorMap.get()] is 'other'): 
            cmap = cm.get_cmap(self.otherColorMap.get())
        else:
            cmap = cm.get_cmap(COLOR_CHOICES[self.colorMap.get()])
        if cmap is None:
            cmap = cm.jet
        return cmap





