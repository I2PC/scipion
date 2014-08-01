# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

from os import remove
from os.path import exists
from protocol_ctf_discrepancy import XmippProtCTFDiscrepancy
from pyworkflow.em import data
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.em.viewer import ObjectView, MODE, MODE_MD, ORDER, VISIBLE
from pyworkflow.protocol.params import FloatParam, IntParam, HiddenBooleanParam
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO, \
    ProtocolViewer
import collections
import numpy as np

#class XmippCTFDiscrepancyViewer(Viewer):
class XmippCTFDiscrepancyViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer CTF Discrepancy'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtCTFDiscrepancy]
    _memory = False
    resolutionThresholdOLD = -1
    #temporary metadata file with ctf that has some resolution greathan than X
    tmpMetadataFile='viewersTmp.sqlite'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        #group = form.addGroup('Overall results')
        form.addParam('resolutionThreshold', FloatParam, default=0., 
                      label='resolution threshold (A)',
                      help='Select only CTF consistent at this resolution (in A).')
        form.addParam('visualizeTable', HiddenBooleanParam, default=False,
                      label="Visualize Table.")
        form.addParam('visualizeMatrix', HiddenBooleanParam, default=False,
                      label="Visualize Matrix.")
        form.addParam('visualizeHistogram', IntParam, default=10,
                      label="Visualize Histogram (Bin size)")



    def _getVisualizeDict(self):
        return {
                 'visualizeTable': self._visualizeTable
                ,'visualizeMatrix': self._visualizeMatrix
                ,'visualizeHistogram': self._visualizeHistogram
                }        

    def _calculateAuxiliaryFile(self):
        """create new ctf_set with ctf that satisfies the
        constraint and persist it
        """
        #metadata file with protocol output
        self.sourceFile = self.protocol.outputCTF.getFileName()
        #get temporary fileName for metadata file
        self.targetFile = self.protocol._getTmpPath(self.tmpMetadataFile)
        resolutionThreshold = self.resolutionThreshold.get()
        print "TODO: this should be closer to the mapper. Here it does not make any sense. ROB"
        #TODO check if this is necessary
        if exists(self.targetFile):
            remove(self.targetFile)
        #metadata with selected CTFs
        _setOfCTFsConst  = data.SetOfCTF(filename=self.targetFile)
        #object read metadata file
        ctfs  = data.SetOfCTF(filename = self.sourceFile)
        #condition to be satisfized for CTFs
        for ctf in ctfs:
            #print "ctf", ctf.printAll()
            if ctf.resolution.get()>resolutionThreshold:
                _setOfCTFsConst.append(ctf)
        #new file with selected CTFs
        _setOfCTFsConst.write()
        #check if empty
        if _setOfCTFsConst.getSize()<1:
            print "WARNING: Empty set of CTFs."
        #close the mapper, if not the object cannot be reused (Although it should be able)
        _setOfCTFsConst.close()

    def isComputed(self):
        _resolutionThreshold = self.resolutionThreshold.get()
        if self.resolutionThresholdOLD != _resolutionThreshold:
            self._calculateAuxiliaryFile()
            self.resolutionThresholdOLD = _resolutionThreshold

    def _visualizeTable(self, e=None):
        """ From here call all visualizations
        """
        self.isComputed()
        views = []

        #display metadata with selected variables
        labels = 'id enabled _micObj._filename method1 method2 resolution _defocusU _defocusV _defocusAngle' 
        views.append(ObjectView(self._project.getName(), 
                                self.protocol.strId(), 
                                self.targetFile,
                                viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels}))
        return views 

    def _visualizeHistogram(self, e=None):
        views = []
        self.isComputed()
        numberOfBins = self.visualizeHistogram.get()
        numberOfBins = min(numberOfBins, self.setOfCTFsConst.getSize())
        plotter = EmPlotter()
        plotter.createSubPlot("Resolution Discrepancies histogram", 
                      "Resolution", "# of Micrographs")
        resolution = [ctf.resolution.get() for ctf in self.setOfCTFsConst]
        plotter.plotHist(resolution, nbins=numberOfBins)
        return views.append(plotter)

    def _visualizeMatrix(self,e=None):
        views = []
        self.isComputed()
        inputCTFs=self.protocol.inputCTFs
        _matrix = np.zeros(shape=(len(inputCTFs), len(inputCTFs)))
        _matrix[0][0]=1
        _matrix[1][0]=2
        _matrix[0][1]=3
        _matrix[1][1]=4
        
        ticksLablesMajor=[]

        self.methodNames = collections.OrderedDict()
        for i, ctf in enumerate(inputCTFs):
            protocol = self.protocol.getMapper().getParent(ctf.get())
            name = "(%d) %s " % (i+1, protocol.getClassLabel())
            ticksLablesMajor.append(name)
        plotter = EmPlotter(fontsize=14)
        resolution=2
        plotter.createSubPlot("# micrographs with resolution\n lower than %d"%(resolution), 
                      "Resolution", "# of Micrographs")
#        plotter.plotMatrix(_matrix,cmap='seismic_r'
        plotter.plotMatrix(_matrix,cmap='Greens'
                        , xticksLablesMajor=ticksLablesMajor
                        , yticksLablesMajor=ticksLablesMajor)
        return views.append(plotter)

