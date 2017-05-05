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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from protocol_refmac import CCP4ProtRunRefmac
from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer, Viewer
from pyworkflow.gui.text import _open_cmd
from pyworkflow.em.data import EMSet, EMObject
from pyworkflow.object import Float, String
from pyworkflow.em.viewer import ObjectView, TableView
from tkMessageBox import showerror

from pyworkflow.em.viewer import ImageView, ChimeraView
import os


def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except:
        print("Error:", msg)


class SingletonParseFile(object):
    """class that parse de log file. It is defined as a singleton so
       it may be called from all routines that need to parse data but the actual
       parsing happends only once"""
    _instance = None

    def __init__(self, fileName, tkParent=None, lastIteration=0):
        if not self.__initialized:
            self.__initialized = True
            object.__init__(self)
            self.headerDict = {}#parsed header goes here
            self.dataDict = {}#parsed data goes here
            self.tkParent = tkParent
            self.fileName = fileName
            self._parsefile(lastIteration)

    def __new__(cls, *args, **kwargs):
        """Singleton magic happend here """
        if not cls._instance:
            cls._instance = super(SingletonParseFile, cls).__new__(cls, *args, **kwargs)
            cls._instance.__initialized = False

        return cls._instance

    def _parseFinalResults(self, filePointer):
        headerList = []
        dataList = []
        stop = False
        while 1:
            line = filePointer.readline()
            if line.strip() == '$TEXT:Result: $$ Final results $$':  # detect final results
                break
            if not line:
                stop = True
                break
        if stop:
            self.msg = 'Can not find "Final result" information in log file: %s'% self.fileName
        else:
            # finalResultsDict={'header':}
            #parse header
            headerList.append(" ")
            line = filePointer.readline()
            words = line.strip().split()
            headerList.extend([words[0],words[1]])
            #parse data: another 4 lines
            for i in range(4):
                row = []
                line = filePointer.readline()
                words = line.strip().split()
                #the first column has 2 words
                row.extend([words[0]+" "+ words[1], words[2], words[3]])
                dataList.append(tuple(row))
            #TODO: remove debug lines
        print("FinalResults", headerList)
        print("FinalResults", dataList)
        return headerList, dataList

    def retrievefinalResults(self):
        return self.headerDict['finalResults'], self.dataDict['finalResults']

    def _parseLastIteration(self, filePointer, iteration):
        print("iteration",iteration)
        headerList = ["variable","value"]
        dataList = []
        stop = False
        print("line","$GRAPHS:Cycle   %d. M(Fom) v. resln :N:1,3,5,7,8,9,10:"%iteration)
        while 1:
            line = filePointer.readline()
            if line.strip() == '$GRAPHS:Cycle   %d. M(Fom) v. resln :N:1,3,5,7,8,9,10:'%iteration:  # detect final results
                break
            if not line:
                stop = True
                break
        print("stop1",stop)
        if stop:
            self.msg = 'Can not find "Last Iteration" information in log file: %s' % self.fileName
        else:
            # find three lines with $$
            counter=4
            while counter!=0:
                line = filePointer.readline()
                if line.find("$$")!=-1:
                    counter -= 1
                if not line:
                    stop = True
                    break
            print("stop2", stop)

            for i in range(14):
                row = []
                line = filePointer.readline()
                print "LINE", line
                words = line.strip().split("=")
                # the first column has 2 words
                row.extend([words[0].strip(), words[1].strip()])
                dataList.append(tuple(row))
            # TODO: remove debug lines
            print("LastIteration", headerList)
            print("LastIteration", dataList)
        return headerList, dataList

    def retrievelastIteration(self):
        return self.headerDict['lastIteration'], self.dataDict['lastIteration']

    def _parsefile(self, lastIteration=0):
        """ call the different functions that parse the data in the right order"""
        with open(self.fileName,"r") as filePointer:
            headerList,dataList = self._parseLastIteration(filePointer, lastIteration)
            self.headerDict['lastIteration'] = headerList
            self.dataDict['lastIteration']   = dataList
            headerList,dataList = self._parseFinalResults(filePointer)
            self.headerDict['finalResults'] = headerList
            self.dataDict['finalResults']   = dataList

class CCP4ProtRunRefmacViewer(ProtocolViewer):
    """ Viewer for CCP4 program refmac
    """
    _label = 'Refmac Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [CCP4ProtRunRefmac]

    # ROB: do we need this memory for something?
    # _memory = False
    # temporary metadata file with ctf that has some resolution greathan than X
    # tmpMetadataFile = 'viewersTmp.sqlite'

    def _defineParams(self, form):
        form.addSection(label='Visualization of Refmac results')
        # group = form.addGroup('Overall results')
        form.addParam('displayMask', LabelParam,
                      label="PDB based Mask",
                      help="Display Masked map")
        form.addParam('showFinalResults', LabelParam,
                      label="Final Results Table",
                      help="Table of Final Results from refine.log file")
        form.addParam('showLogFile', LabelParam,
                      label="Show log file",
                      help="open refmac log file in a text editor")
        form.addParam('showLastIteration', LabelParam,
                      label="Results Table (last iteration)",
                      help="Table stored in log file summarizing the last iteration")
        form.addParam('displayRFactorPlot', LabelParam,
                      label="R-factor vs. iteration",
                      help="Plot R-factor as a function of the iteration")
        form.addParam('displayFOMPlot', LabelParam,
                      label="FOM vs. iteration",
                      help="Plot Figure Of Merit as a function of the iteration")
        form.addParam('displayLLPlot', LabelParam,
                      label="-LL vs. iteration",
                      help="Plot Log likelihood as a function of the iteration")
        form.addParam('displayLLfreePlot', LabelParam,
                      label="-LLfree vs. iteration",
                      help="Plot Log likelihood as a function of the iteration")
        form.addParam('displayGeometryPlot', LabelParam,
                      label="Geometry vs. iteration",
                      help="""Plot Geometry as a function of the iteration:
Geometry includes rmsBOND (root mean square bond lengths)
zBOND (zscore of the deviation of bond lengths)
rmsANGL (root mean square bond angles)
zANGL (zscore of the deviation of bond angles)
and rmsCHIRAL (root mean square of chiral index""")

    def _getVisualizeDict(self):
        return {
            'showFinalResults': self._visualizeFinalResults,
            'showLastIteration': self._visualizeLastIteration,
            'displayMask': self._visualizeMask,
            'displayRFactorPlot': self._visualizeRFactorPlot,
            'displayFOMPlot': self._visualizeFOMPlot,
            'displayLLPlot': self._visualizeLLPlot,
            'displayLLfreePlot': self._visualizeLLfreePlot,
            'displayGeometryPlot': self._visualizeGeometryPlot,
            'showLogFile': self._visualizeLogFile
        }

    def _visualizeMask(self):
        pass

    def _visualizeFinalResults(self, e=None):


        """
        views = []
        labels = '_1 _2'
        emSet = EMSet(filename="/tmp/kk.sqlite")
        emObject = EMObject()
        emObject._1 = String('first parameter')
        emObject._2 = Float(12.)
        emSet.append(emObject)
        emObject = EMObject()
        emObject._1 = String('second parameter')
        emObject._2 = Float(22.)
        emSet.append(emObject)
        emSet.write()
        views.append(ObjectView(self._project,
                                self.protocol.strId(),
                                "/tmp/kk.sqlite",
                                viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels}))
        return views
"""
        #Selection of lines from 'refine.log' file that include Refmac final results.
        #These lines will be saved in outputLines list

        singletonParseFile = SingletonParseFile(self.protocol._getExtraPath(self.protocol.refineLogFileName),
                                                self.getTkRoot(), self.protocol.nRefCycle.get() + 1)
        headerList, dataList = singletonParseFile.retrievefinalResults()
        if not dataList:
            errorWindow(self.getTkRoot(),singletonParseFile.msg)
            return

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg="Values for a good fitted 3D map. R factor ~ 0.3, Rms BondLength ~ 0.02.",
                  title= "Refmac: Final Results Summary",
                  height=len(dataList), width=200,padding=40)

    def _visualizeLogFile(self, e=None):
        """Show refmac log file."""
        refineLogFileName = self.protocol._getExtraPath(self.protocol.refineLogFileName)
        _open_cmd(refineLogFileName, self.getTkRoot())

    def _visualizeLastIteration(self, e=None):

        # Selection of lines from 'refine.log' file that include Refmac final results.
        # These lines will be saved in outputLines list
        singletonParseFile = SingletonParseFile(self.protocol._getExtraPath(self.protocol.refineLogFileName),
                                                self.protocol.nRefCycle.get() + 1)
        dataList=singletonParseFile.dataDict['lastIteration']
        TableView(headerList=singletonParseFile.headerDict['lastIteration'],
                  dataList=dataList,
                  mesg=" ",
                  title= "Refmac: Last Iteration summary",
                  height=len(dataList), width=200,padding=40)

        """
        f = self.protocol._getExtraPath(self.protocol.refineLogFileName)
        outputLines = []
        with open(f) as input_data:
            for line in input_data:
                if line.strip() == '$TABLE: Cycle   '+self.protocol.nRefCycle.get()+1+'. FSC and  Fom(<cos(DelPhi)>-acentric, centric, overall v resln:':
                    break
            for line in input_data:
                if line.strip() == 'Things for loggraph, R factor and others vs cycle':
                    break
                outputLines.append(line)
        outputLines = outputLines[25:38]
        # Creation of two lists (headerList and dataList) with the first line and the remaining lines of outputLines, respectively
        headerList = []
        dataList = []
        for i in range(0, len(outputLines)):
            if i == 0:
                headerList.extend([' ', outputLines[i].strip().split()[0], outputLines[i].strip().split()[1]])
            else:
                dataList.extend([outputLines[i].strip().split()[0] + ' ' + outputLines[i].strip().split()[1],
                                 outputLines[i].strip().split()[2],
                                 outputLines[i].strip().split()[3]])
        # Conversion of dataList in a list of 3-element tuples
        it = iter(dataList)
        dataList = zip(it, it, it)
        # Arrangement of the final table
        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg="This list includes a summary of Refmac execution final results",
                  title="Refmac: Final Results Summary",
                  height=len(dataList))
"""
    def _visualizeRFactorPlot(self, e=None):
        """
        xplotter = Plotter(windowTitle="SAXS Curves")
        a = xplotter.createSubPlot('SAXS curves', 'Angstrongs^-1', 'log(SAXS)', yformat=False)
        a.plot(x[:, 0], numpy.log(x[:, 1]))
        a.plot(x[:, 0], numpy.log(x[:, 2]))
        if obj.experimentalSAXS.empty():
            xplotter.showLegend(['SAXS in solution', 'SAXS in vacuo'])
        else:
            xplotter.showLegend(['Experimental SAXS', 'SAXS from volume'])
        xplotter.show()
        """
        pass

    def _visualizeFOMPlot(self, e=None):
        pass

    def _visualizeLLPlot(self, e=None):
        pass

    def _visualizeLLfreePlot(self, e=None):
        pass

    def _visualizeGeometryPlot(self, e=None):
        pass