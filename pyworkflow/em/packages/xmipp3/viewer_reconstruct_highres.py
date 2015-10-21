# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
# *           Slavica Jonic                (jonic@impmc.upmc.fr)
# * Ported to Scipion:
# *           J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Nov 2014
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

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, Viewer
from pyworkflow.em.viewer import ObjectView, DataView
import pyworkflow.em.showj as showj
from protocol_reconstruct_highres import XmippProtReconstructHighRes
from plotter import XmippPlotter
from glob import glob
from os.path import exists, join
from pyworkflow.em.packages.xmipp3.convert import getImageLocation

class XmippReconstructHighResViewer(Viewer):
    """ Visualize the output of protocol reconstruct highres """
    _label = 'viewer reconstruct highres'
    _targets = [XmippProtReconstructHighRes]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    def _visualize(self, prot, **kwargs):
        views = []
        
        # FSC
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
        xplotter = XmippPlotter(windowTitle="FSC")
        a = xplotter.createSubPlot("FSC", "Frequency (1/A)", "Zscore")
        fnIterDirs = sorted(glob(prot._getExtraPath("Iter*")))
        legends = []
        lastFullIter = -1
        for it in range(0,len(fnIterDirs)):
            fnDir = prot._getExtraPath("Iter%03d"%it)
            fnFSC = join(fnDir,"fsc.xmd")
            if exists(fnFSC):
                lastFullIter=it
                show = True
                legends.append('iter %d' % it)
                self._plotFSC(a, fnFSC)
                xplotter.showLegend(legends)
        a.grid(True)
        views.append(xplotter)
        
        # Angle file
        fnLastAngles=join(prot._getExtraPath("Iter%03d"%lastFullIter),"angles.xmd")
        if hasattr(prot, "outputParticles"):
            obj = prot.outputParticles
            fn = obj.getFileName()
            labels = 'id enabled _filename _xmipp_zScore _xmipp_cumulativeSSNR '
            labels += '_ctfModel._defocusU _ctfModel._defocusV _xmipp_shiftX _xmipp_shiftY _xmipp_tilt _xmipp_continuousX _xmipp_continuousY _xmipp_scale _xmipp_maxCC _xmipp_weight'
            labels += " _xmipp_cost _xmipp_weightContinuous2 _xmipp_angleDiff _xmipp_weightJumper _xmipp_weightSSNR"
            views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={showj.ORDER: labels, 
                                                      showj.VISIBLE: labels, 
                                                      showj.MODE: showj.MODE_MD,
                                                      showj.RENDER:'_filename'}))
        else:
            if lastFullIter>0:
                views.append(DataView(fnLastAngles, viewParams={showj.MODE: showj.MODE_MD}))

        # Volume
        if hasattr(prot, "outputVolume"):
            obj = prot.outputVolume
            fn = getImageLocation(obj)
            views.append(ObjectView(self._project, obj.strId(), fn, viewParams={showj.RENDER: 'image', showj.SAMPLINGRATE: obj.getSamplingRate()}))
        else:
            if lastFullIter>0:  
                fnLastVolume=join(prot._getExtraPath("Iter%03d"%lastFullIter),"volumeAvg.mrc")
                views.append(ObjectView(self._project, None, fnLastVolume, viewParams={showj.RENDER: 'image'}))
        
        # Jumper weights                                    
        if lastFullIter>0:
            if prot.weightJumper and lastFullIter>1:
                import xmipp
                xplotter = XmippPlotter(windowTitle="Jumper weight")
                a = xplotter.createSubPlot("Jumper weight", "Weight", "Count")
                xplotter.plotMdFile(fnLastAngles,xmipp.MDL_WEIGHT_JUMPER,xmipp.MDL_WEIGHT_JUMPER,nbins=100)
                views.append(xplotter)

        return views

    def _plotFSC(self, a, fnFSC):
        import xmipp
        md = xmipp.MetaData(fnFSC)
        resolution_inv = [md.getValue(xmipp.MDL_RESOLUTION_FREQ, id) for id in md]
        frc = [md.getValue(xmipp.MDL_RESOLUTION_FRC, id) for id in md]
        self.maxFrc = max(frc)
        self.minInv = min(resolution_inv)
        self.maxInv = max(resolution_inv)
        a.plot(resolution_inv, frc)
        a.xaxis.set_major_formatter(self._plotFormatter)
        a.set_ylim([-0.1, 1.1])

    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999
        if value:
            inv = 1/value
        return "1/%0.2f" % inv
