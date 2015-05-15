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


from protocol_ctf_discrepancy import XmippProtCTFDiscrepancy
from pyworkflow.em import data
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.em.viewer import ObjectView
from pyworkflow.em.showj import MODE, MODE_MD, ORDER, VISIBLE
from pyworkflow.protocol.params import FloatParam, IntParam, LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.utils.path import cleanPath
import numpy as np



class XmippCTFDiscrepancyViewer(ProtocolViewer):
    """ This protocol computes the maximum resolution up to which two
     CTF estimations would be ``equivalent'', defining ``equivalent'' as having
      a wave aberration function shift smaller than 90 degrees
    """
    _label = 'viewer CTF Discrepancy'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtCTFDiscrepancy]
    _memory = False
    resolutionThresholdOLD = -1
    #temporary metadata file with ctf that has some resolution greathan than X
    tmpMetadataFile='viewersTmp.sqlite'

    def _defineParams(self, form):
        
        self.averagesFile, self.pairsFile = self.protocol._getAnalyzeFiles()
        
        form.addSection(label='Visualization')
        #group = form.addGroup('Overall results')
        form.addParam('resolutionThreshold', FloatParam, default=999999.,
                      label='Resolution threshold (A)',
                      help='Select only CTF consistent at this resolution (in A).')
        form.addParam('visualizePairs', LabelParam,
                      label="Visualize comparison table.",
                      help="List with resolution at which the CTF estimated by a pair of methods"
                           " is no longer equivalent."  )
        form.addParam('visualizeAverage', LabelParam,
                      label="Visualize average table.",
                      help="Show a table with the averaged CTFs."  )     
        form.addParam('visualizeMatrix', LabelParam,
                      label="Visualize comparison matrix.",
                      help="Number of micrographs that have a CTF estimation"
                           " -given by two methods- that are equivalent at resolution=threshold")
        form.addParam('visualizeHistogram', IntParam, default=10,
                      label="Visualize Histogram (Bin size)",
                      help="Histogram of the resolution at which two methods are equivalent")

    def _getVisualizeDict(self):
        return {
                 'visualizePairs': self._visualizePairs,
                 'visualizeAverage': self._visualizeAverages,
                 'visualizeMatrix': self._visualizeMatrix,
                 'visualizeHistogram': self._visualizeHistogram
                }        

    def _calculateAuxiliaryFile(self):
        """create new ctf_set with ctf that satisfies the
        constraint and persist it
        """
        try:
            self.setOfCTFsConst
        except AttributeError:
            pass
        else:
            #TODO close the mapper, if not the object cannot be reused (Although it should be able)
            self.setOfCTFsConst.close()

        #metadata file with protocol output
        #get temporary fileName for metadata file
        self.targetFile = self.protocol._getTmpPath(self.tmpMetadataFile)
        resolutionThreshold = self.resolutionThreshold.get()
        print "TODO: this should be closer to the mapper. Here it does not make any sense. ROB"
        #TODO check if this is necessary
        cleanPath(self.targetFile)
        
        #metadata with selected CTFs
        self.setOfCTFsConst  = data.SetOfCTF(filename=self.targetFile)
        #object read metadata file
        ctfs  = data.SetOfCTF(filename=self.pairsFile)
        #condition to be satisfized for CTFs
        for ctf in ctfs:
            if ctf.resolution < resolutionThreshold:
                self.setOfCTFsConst.append(ctf)
        #new file with selected CTFs
        self.setOfCTFsConst.write()
        #check if empty
        if self.setOfCTFsConst.getSize() < 1:
            print "WARNING: Empty set of CTFs."
        #TODO aggregation function to compute the minimum resolution by micrograph
        #TODO aggregation function to compute the average resolution by micrograph
        #self._setOfCTFsConst.close()

    def isComputed(self):
        _resolutionThreshold = self.resolutionThreshold.get()
        if self.resolutionThresholdOLD != _resolutionThreshold:
            self._calculateAuxiliaryFile()
            self.resolutionThresholdOLD = _resolutionThreshold

    def _visualizePairs(self, e=None):
        """ From here call all visualizations
        """
        self.isComputed()
        views = []

        #display metadata with selected variables
        labels = 'id enabled _micObj._filename method1 method2 resolution _defocusU _defocusV _defocusAngle' 
        views.append(ObjectView(self._project,
                                self.protocol.strId(), 
                                self.targetFile,
                                viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels}))
        return views 

    def _visualizeAverages(self, e=None):
        """ From here call all visualizations
        """
        views = []

        #display metadata with selected variables
        labels = '_micObj._filename averageDefocusU averageDefocusV averageDefocusAngle averageResolution' 
        views.append(ObjectView(self._project,
                                self.protocol.strId(), 
                                self.averagesFile,
                                viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels}))
        return views 
    
    def _visualizeHistogram(self, e=None):
        views = []
        self.isComputed()
        numberOfBins = self.visualizeHistogram.get()
        numberOfBins = min(numberOfBins, self.setOfCTFsConst.getSize())
        plotter = EmPlotter()
        plotter.createSubPlot("Resolution Discrepancies histogram", 
                      "Resolution (A)", "# of Comparisons")
        resolution = [ctf.resolution.get() for ctf in self.setOfCTFsConst]
        plotter.plotHist(resolution, nbins=numberOfBins)
        return views.append(plotter)

    def extract(self, string, start='(', stop=')'):
        return string[string.index(start)+1:string.index(stop)]

    def _visualizeMatrix(self,e=None):
        views = []
        self.isComputed()
        inputCTFs=self.protocol.inputCTFs
        shape=len(inputCTFs)#number of methods
        _matrix = np.zeros(shape=(shape, shape))
        ticksLablesMajorX=[None] * shape
        ticksLablesMajorY=[None] * shape
        for ctf in self.setOfCTFsConst:
            m1 = int(self.extract(ctf.getAttributeValue('method1')))-1
            m2 = int(self.extract(ctf.getAttributeValue('method2')))-1
            _matrix[m1][m2] += 1
            _matrix[m2][m1] = _matrix[m1][m2]
            #rather inefficient but since  outputCTF is not ordered...
            if ticksLablesMajorX[m1] is None:
                ticksLablesMajorX[m1] = "(%d)"%(m1+1)
                ticksLablesMajorY[m1] = ctf.getAttributeValue('method1')
            if ticksLablesMajorX[m2] is None:
                ticksLablesMajorX[m2] = "(%d)"%(m2+1)
                ticksLablesMajorY[m2] = ctf.getAttributeValue('method2')

        plotter = EmPlotter(fontsize=14)
        resolution=self.resolutionThreshold.get()
        plotter.createSubPlot("# micrographs with resolution (A)\n lower than %d"%(resolution),
                      "Method #", "Method")
#        plotter.plotMatrix(_matrix,cmap='seismic_r'
        plotter.plotMatrix(_matrix,cmap='Greens'
                        , xticksLablesMajor=ticksLablesMajorX,rotationX=0
                        , yticksLablesMajor=ticksLablesMajorY)
        views.append(plotter)
        return views
        #return views.append(plotter)

