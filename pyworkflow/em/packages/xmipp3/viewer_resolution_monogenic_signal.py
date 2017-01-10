# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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

from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.em.viewer import ChimeraView, ObjectView, DataView
from protocol_resolution_monogenic_signal import XmippProtMonoRes, OUTPUT_RESOLUTION_FILE
from pyworkflow.em.metadata import MetaData, MDL_X, MDL_COUNT
from pyworkflow.em import ImageHandler
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


class XmippMonoResViewer(ProtocolViewer):
    """
    Visualization tools for MonoRes results.
    
    MonoRes is a Xmipp packagefor computing the local resolution of 3D
    density maps studied in structural biology, primarily by cryo-electron
    microscopy (cryo-EM).
    """
    _label = 'viewer MonoRes'
    _targets = [XmippProtMonoRes]      
    _environments = [DESKTOP_TKINTER]
   
    
    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, *args, **kwargs)


    def _defineParams(self, form):
        form.addSection(label='Visualization')
        
        form.addParam('doShowVolumeSlices', LabelParam,
                      label="Show volume slices")
        
        form.addParam('doShowVolumeColorSlices', LabelParam,
              label="Show colored slices")
        
        form.addParam('doShowChimera', LabelParam,
                      label="Show Resolution map in Chimera")

        form.addParam('doShowResHistogram', LabelParam,
                      label="Show resolution histogram")

        
    def _getVisualizeDict(self):
        return {'doShowVolumeSlices': self._showVolumeSlices,
                'doShowVolumeColorSlices': self._showVolumeColorSlices,
                'doShowResHistogram': self._plotHistogram,
                'doShowChimera': self._showChimera,
                }
       
    def _showVolumeSlices(self, param=None):
        cm = DataView(self.protocol.outputVolume.getFileName())
        
        return [cm]
    
    def _showVolumeColorSlices(self, param=None):
        imageFile = self.protocol._getExtraPath(OUTPUT_RESOLUTION_FILE)
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        max_Res = np.amax(imgData)

        #  This is to generate figures for the paper
        # min_Res = np.amin(imgData)
        # imgData2 = imgData
        imgData2 = np.ma.masked_where(imgData < 0.01, imgData, copy=True)
        
        min_Res = np.amin(imgData2)
        fig, im = self._plotVolumeSlices('MonoRes slices', imgData2, 
                                                  min_Res, max_Res, self.getColormap())
        cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cbar = fig.colorbar(im, cax=cax)
        cbar.ax.invert_yaxis()

        return [Plotter(figure=fig)]

    def _plotHistogram(self, param=None):
        md = MetaData()
#         path_ = self.protocol._getPath('extra/hist.xmd')
        md.read(self.protocol._getPath('extra/hist.xmd'))
        x_axis = []
        y_axis = []

        i = 0
        for idx in md:
            x_axis_ = md.getValue(MDL_X, idx)
            if i==0:
                x0 = x_axis_
            elif i==1:
                x1 = x_axis_
            y_axis_ = md.getValue(MDL_COUNT, idx)

            i+=1
            if (y_axis_== 0):
                continue
            x_axis.append(x_axis_)
            y_axis.append(y_axis_)
        delta = x1-x0
        for ii in range(len(x_axis)):
            x_axis[ii] = x_axis[ii]-0.5*delta
        plt.bar(x_axis, y_axis, width = delta)
        plt.title("Resolutions Histogram")
        plt.xlabel("Resolution (A)")
        plt.ylabel("Counts")
        
        return [plt.show()]

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
    
        def getSlice(slice):
            if dataAxis == 'y':
                return volumeData[:,slice,:]
            elif dataAxis == 'x':
                return volumeData[:,:,slice]
            else:
                return volumeData[slice,:,:]
    
        def showSlice(ax, index):
            sliceTitle = 'Slice %s' % int(index*size/9)
            slice = int(index*origSize/9)
            ax.set_title(sliceTitle, fontsize=sliceFontSize, color=sliceColor)
            return ax.imshow(getSlice(slice), vmin=vminData, vmax=vmaxData,
                             cmap=cmap, interpolation="nearest")
        
        im = showSlice(ax1, 3)
        showSlice(ax2, 4)
        showSlice(ax3, 5)
        showSlice(ax4, 6)
        
        return f, im 

    def _showChimera(self, param=None):
        cmdFile = self.protocol._getPath('Chimera_resolution.cmd')
        view = ChimeraView(cmdFile)
        return [view]

    @staticmethod
    def getColormap():

        c = mcolors.ColorConverter().to_rgb
        rvb = XmippMonoResViewer.make_colormap(
             [c('blue'),
              c('cyan'), 0.25, c('cyan'),
              c('green'), 0.50, c('green'),
              c('yellow'), 0.75, c('yellow'),
              c('red')])
        return rvb

    @staticmethod
    def make_colormap(seq):
        """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
        seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
        cdict = {'red': [], 'green': [], 'blue': []}
        for i, item in enumerate(seq):
            if isinstance(item, float):
                r1, g1, b1 = seq[i - 1]
                r2, g2, b2 = seq[i + 1]
                cdict['red'].append([item, r1, r2])
                cdict['green'].append([item, g1, g2])
                cdict['blue'].append([item, b1, b2])
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)
