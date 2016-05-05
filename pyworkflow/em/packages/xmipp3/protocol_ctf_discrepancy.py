# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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

import collections
from itertools import izip

from pyworkflow.utils.path import cleanPath, removeBaseExt
from pyworkflow.object import Set, Float, String, Object
import pyworkflow.protocol.params as params

from pyworkflow.em.metadata import Row, MetaData
import pyworkflow.em as em
import convert 
import xmipp


class XmippProtCTFDiscrepancy(em.ProtCTFMicrographs):
    """
    Protocol to estimate the agreement between different estimation of the CTF
    for the same set of micrographs. The algorithm assumes that two CTF are consistent
    if the phase (wave aberration function) of the two CTFs are closer than 90 degrees.
    The reported resolution is the resolution at which the two CTF phases differ in 90 degrees.
    """
    _label = 'ctf discrepancy'
    
    def __init__(self, **args):
        em.ProtCTFMicrographs.__init__(self, **args)
        self._freqResol = {}

    def _defineParams(self, form):
        form.addSection(label='Input')
        # Read N ctfs estimations
        form.addParam('inputCTFs', params.MultiPointerParam, pointerClass='SetOfCTF',
                      label="input CTFs",
                      help='Select the first set of CTFs to compare')        

        form.addParallelSection(threads=4, mpi=0)       
        
#--------------------------- INSERT steps functions --------------------------------------------  
                                
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self.setOfCTF = self.inputCTFs[0].get()
        #self.setOfCTF2 = self.inputCTFs[1].get()
        self.methodNames = self._loadMethods()
        
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # For each ctf pair insert the steps to process it
        # check same size, same micrographs
        for method1, method2, ctf in self._iterMethods():
            stepId = self._insertFunctionStep('_computeCTFDiscrepancyStep' 
                                              , ctf.getObjId()
                                              , method1
                                              , method2
                                              ,prerequisites=[]) # Make estimation steps independent between them
            deps.append(stepId)
        # Insert step to create output objects       
        self._insertFunctionStep('createAnalyzeFilesStep', prerequisites=deps)
    
    def _computeCTFDiscrepancyStep(self, ctfId, method1, method2):
        #TODO must be same micrographs
        #convert to md
        mdList = [MetaData(), MetaData()]
        ctfList = [self.inputCTFs[method1].get()[ctfId], self.inputCTFs[method2].get()[ctfId]]
        ctfRow = Row()
        
        for md, ctf in izip(mdList, ctfList):
            objId = md.addObject()
            convert.ctfModelToRow(ctf, ctfRow)
            convert.micrographToRow(ctf.getMicrograph(), ctfRow, alignType=convert.ALIGN_NONE)
            ctfRow.writeToMd(md, objId)
        self._freqResol[(method1, method2, ctfId)] = xmipp.errorMaxFreqCTFs2D(*mdList)

    def createAnalyzeFilesStep(self):
        """ This method will create two sqlite files that 
        will be used in by the viewer in the 'Analyze result' action.
        ctfSet: a table with average ctf values between all methods.
        ctfSetPair: a table where each row is a pair between two methods
        """
        ctfSetFn, ctfSetPairFn = self._getAnalyzeFiles()
        
        cleanPath(ctfSetFn, ctfSetPairFn)
        
        ctfSet     = Set(filename=ctfSetFn)
        ctfSetPair = em.SetOfCTF(filename=ctfSetPairFn)
        #import pdb
        #pdb.set_trace()
        minimumResolution   = {}
        maximumResolution   = {}
        averageResolution   = {}
        averageDefocusU     = {}
        averageDefocusV     = {}
        averageDefocusAngle = {}
        
        for ctf in self.setOfCTF:
            micFileName = self._getMicName(ctf)
            minimumResolution[micFileName]   = 0.
            maximumResolution[micFileName]   = 0.
            averageResolution[micFileName]   = 0.
            averageDefocusU[micFileName]     = 0.
            averageDefocusV[micFileName]     = 0.
            averageDefocusAngle[micFileName] = 0.
        
        for method1, method2, ctfId in self._freqResol:
            ctf = em.CTFModel()
            ctf1 = self.inputCTFs[method1].get()[ctfId]
            ctf2 = self.inputCTFs[method2].get()[ctfId]
            ctf.setDefocusU(    (ctf1.getDefocusU() + ctf2.getDefocusU())/2. )
            ctf.setDefocusV(    (ctf1.getDefocusV() + ctf2.getDefocusV())/2. )
            ctf.setDefocusAngle((ctf1.getDefocusAngle() + ctf2.getDefocusAngle())/2. )
            ctf.setMicrograph(   ctf1.getMicrograph())

            #Clean objId since we can have ctf from the same micrograph
            # and by default it is set to micrograph id
            ctf.cleanObjId()
            resolution = self._freqResol[(method1, method2, ctfId)]
            ctf.resolution = Float(resolution)
            ctf.method1 = String(self.methodNames[method1])
            ctf.method2 = String(self.methodNames[method2])
            # save the values of defocus for each micrograph in a list
            ctfSetPair.append(ctf)
            
            micFileName = self._getMicName(ctf1)
            
            if (not micFileName in minimumResolution or
                not micFileName in maximumResolution):
                pass
            
            
            if resolution < minimumResolution[micFileName]:
                minimumResolution[micFileName] = resolution
            
            if resolution > maximumResolution[micFileName]:
                maximumResolution[micFileName] = resolution
            
            averageResolution[micFileName]    += resolution
            averageDefocusU[micFileName]      += ctf.getDefocusU()
            averageDefocusV[micFileName]      += ctf.getDefocusV()
            averageDefocusAngle[micFileName]  += ctf.getDefocusAngle()
            
        size = float(len(self.setOfCTF))

        for ctf in self.setOfCTF:
            ctfAvg = Object()
            micFileName = self._getMicName(ctf)
            ctfAvg._micObj = ctf.getMicrograph()
            ctfAvg.averageDefocusU     = Float(averageDefocusU[micFileName] / size)
            ctfAvg.averageDefocusV     = Float(averageDefocusV[micFileName]     / size)
            ctfAvg.averageDefocusAngle = Float(averageDefocusAngle[micFileName] / size)
            ctfAvg.averageResolution   = Float(averageResolution[micFileName]   / size)
            ctfSet.append(ctfAvg)
            
        ctfSetPair.write()
        ctfSet.write()

    def _citations(self):
        return ['Marabini2014a']
    
    def _summary(self):
        message = []
        for i, ctf in enumerate(self.inputCTFs):
            protocol = self.getMapper().getParent(ctf.get())
            message.append("Method %d %s" % (i+1, protocol.getClassLabel()))
        #TODO size de la cosa calculada
        ####message.append("Comparered <%d> micrograph" % (size,'micrographs'))
        return message    
    
    def _methods(self):
        pass#nothing here
    
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        #same micrographs in both CTF??
        errors = [ ] 
        # Add some errors if input is not valid
        return errors
    
    def _iterMethods(self):
        for method1 in self.methodNames:
            for method2 in self.methodNames:
                if method1 < method2:
                    for ctf in self.setOfCTF:
                        yield method1, method2, ctf
                        
    def _loadMethods(self):
        """ Load the methods names for the protocols
        that produced the selected input CTFs.
        """
        methodNames = collections.OrderedDict()
        
        for i, ctf in enumerate(self.inputCTFs):
            protocol = self.getMapper().getParent(ctf.get())
            methodNames[i] = "(%d) %s " % (i+1, protocol.getClassLabel())

        return methodNames
    
    def _getMicName(self, ctf):
        """ Get the micName to be used as Key from
        a given ctf object. 
        """
        return removeBaseExt(ctf.getMicrograph().getFileName())
    
    def _getAnalyzeFiles(self):
        """ Return the name of the analyze result files. """
        return (self._getExtraPath('ctf_set.sqlite'),
                self._getExtraPath('ctf_set_pairs.sqlite'))
        
        