# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin  (jmdelarosa@cnb.csic.es)
# *              Yuxiang Chen 
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************

import os

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as em  

from convert import readSetOfVolumes



class ProtPyTomImport(em.ProtImport):
    """ Subtomogram averaging using pytom autofocus """
    _label = 'frm 3d align'

    #--------------------------- DEFINE param functions --------------------------------------------
    
    def _defineParams(self, form):
        
        form.addSection('Input')
        form.addParam('inputXml', params.FileParam,
              label='Input particles XML file',
              help="Select the particles xml file thats have \n"
                   "the result of some PyTom program.")
        
        form.addParam('samplingRate', params.FloatParam,
           label=pwutils.properties.Message.LABEL_SAMP_RATE)
        
        form.addParam('copyFiles', params.BooleanParam, default=False, 
                      expertLevel=params.LEVEL_ADVANCED, 
                      label="Copy files?",
                      help="By default the files are not copied into the\n"
                           "project to avoid data duplication and to save\n"
                           "disk space. Instead of copying, symbolic links are\n"
                           "created pointing to original files. This approach\n"
                           "has the drawback that if the project is moved to\n"
                           "another computer, the links need to be restored.\n")

        form.addParallelSection(threads=0, mpi=0)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep') 
        
    #--------------------------- STEPS functions --------------------------------------------        

    def createOutputStep(self):
        volSet = self._createSetOfVolumes()
        volSet.setSamplingRate(self.samplingRate.get())
        readSetOfVolumes(self.inputXml.get(), volSet)
        
        self._defineOutputs(outputVolumes=volSet)

    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        """ Check that some preconditions are met before launching 
        the auto-focus classification run. 
        """
        errors = []
        return errors
        
    def _summary(self):
        summary = []
        return summary
    
    def _citations(self):
        return ['Chen2014']
    
    def _methods(self):
        return []
    
    #--------------------------- UTILS functions --------------------------------------------   
    
    def _getScript(self, *paths):
        return os.path.join(os.environ['PYTOM_HOME'], *paths)
    
