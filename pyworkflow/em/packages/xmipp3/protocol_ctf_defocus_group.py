# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
This file implements CTF defocus groups using xmipp 3.1
"""
from pyworkflow.em import *  
from pyworkflow.utils import *  
from convert import createXmippInputImages, writeSetOfDefocusGroups
import xmipp
from math import pi

# TODO: change the base class to a more apropiated one
class XmippProtCTFDefocusGroup(ProtProcessParticles):
    """
    Create a list of defocus values that delimite 
    defocus groups.
    """
    _label = 'defocus group'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        """ Define the parameters that will be input for the Protocol.
        This definition is also used to generate automatically the GUI.
        """
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, 
                      pointerClass='SetOfParticles', pointerCondition='hasCTF',
                      label="Input particles with CTF", 
                      help='Select the input particles. \n '
                           'they should have information about the CTF (hasCTF=True)')
        form.addParam('ctfGroupMaxDiff', FloatParam, default=1,
                      label='Error for grouping', validators=[GE(1.,'Error must be greater than 1')],
                      help='Maximum error when grouping, the higher the more groups'
                           'This is a 1D program, only defocus U is used\n '
                           'the frequency at which the phase difference between the CTF\n'
                           'belonging to 2 particles is equal to Pi/2 is computed \n '
                           'If this difference is less than 1/(2*factor*sampling_rate)\n' 
                           'then images are placed in different groups')          
        
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        """ In this function the steps that are going to be executed should
        be defined. Two of the most used functions are: _insertFunctionStep or _insertRunJobStep
        """
        #TODO: when aggregation functions are defined in Scipion set
        # this step can be avoid and the protocol can remove Xmipp dependencies
        imgsFn          = createXmippInputImages(self, self.inputParticles.get())
        ctfGroupMaxDiff = self.ctfGroupMaxDiff.get()
        
        #verifyFiles = []
        self._insertFunctionStep('createOutputStep', imgsFn, ctfGroupMaxDiff)
    
    #--------------------------- STEPS functions --------------------------------------------       
    def createOutputStep(self, imgsFn, ctfGroupMaxDiff):
        """ Create defocus groups and generate the output set """
        fnScipion = self._getPath('defocus_groups.sqlite')
        fnXmipp   = self._getPath('defocus_groups.xmd')
        setOfDefocus = SetOfDefocusGroup(filename=fnScipion)
        df = DefocusGroup()
        mdImages    = xmipp.MetaData(imgsFn)
        mdGroups    = xmipp.MetaData()
        mdGroups.aggregateMdGroupBy(mdImages, xmipp.AGGR_COUNT, 
                                               [xmipp.MDL_CTF_DEFOCUSU,
                                                xmipp.MDL_CTF_DEFOCUS_ANGLE,
                                                xmipp.MDL_CTF_SAMPLING_RATE,
                                                xmipp.MDL_CTF_VOLTAGE,
                                                xmipp.MDL_CTF_CS,
                                                xmipp.MDL_CTF_Q0], 
                                                xmipp.MDL_CTF_DEFOCUSU, 
                                                xmipp.MDL_COUNT)
        mdGroups.sort(xmipp.MDL_CTF_DEFOCUSU)
        
        mdCTFAux = xmipp.MetaData()
        idGroup  = mdGroups.firstObject()
        idCTFAux = mdCTFAux.addObject()
        mdCTFAux.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, 0., idCTFAux);
        mdCTFAux.setValue(xmipp.MDL_CTF_SAMPLING_RATE, mdGroups.getValue(xmipp.MDL_CTF_SAMPLING_RATE,idGroup), idCTFAux)
        mdCTFAux.setValue(xmipp.MDL_CTF_VOLTAGE,       mdGroups.getValue(xmipp.MDL_CTF_VOLTAGE,idGroup),       idCTFAux)
        mdCTFAux.setValue(xmipp.MDL_CTF_CS,            mdGroups.getValue(xmipp.MDL_CTF_CS,idGroup),            idCTFAux)
        mdCTFAux.setValue(xmipp.MDL_CTF_Q0,            mdGroups.getValue(xmipp.MDL_CTF_Q0,idGroup),            idCTFAux)
        resolutionError= mdGroups.getValue(xmipp.MDL_CTF_SAMPLING_RATE,idGroup)

        counter = 0
        minDef  = mdGroups.getValue(xmipp.MDL_CTF_DEFOCUSU,idGroup)
        maxDef  = minDef
        avgDef  = minDef

        mdCTFAux.setValue(xmipp.MDL_CTF_DEFOCUSU, minDef, idCTFAux)
        for idGroup in mdGroups:
            defocusU = mdGroups.getValue(xmipp.MDL_CTF_DEFOCUSU,idGroup)
            mdCTFAux.setValue(xmipp.MDL_CTF_DEFOCUSV, defocusU, idCTFAux);
            counter += 1
            resolution = xmipp.errorMaxFreqCTFs(mdCTFAux,pi/2.)
            if  resolution > resolutionError/ctfGroupMaxDiff:
                avgDef /=  counter
                df.cleanObjId()
                df.setDefocusMin(minDef)
                df.setDefocusMax(maxDef)
                df.setDefocusAvg(avgDef)
                setOfDefocus.append(df)
                counter = 0
                minDef = defocusU
                avgDef = defocusU
                mdCTFAux.setValue(xmipp.MDL_CTF_DEFOCUSU, defocusU, idCTFAux);
            else:
                avgDef  += defocusU
                
            maxDef = defocusU
        ###
        setOfDefocus.printAll()
        ###
        writeSetOfDefocusGroups(setOfDefocus, fnXmipp)
        self._defineOutputs(outputDefocusGroups=setOfDefocus)
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        errors = [ ] 
        # Add some errors if input is not valid
        return errors
    
    def _citations(self):
        cites = []            
        return cites
    
    def _summary(self):
        summary = []

        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
    