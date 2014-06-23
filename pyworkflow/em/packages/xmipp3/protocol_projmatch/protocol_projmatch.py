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
This sub-package implement projection matching using xmipp 3.1
"""

from pyworkflow.em import ProtRefine3D, ProtClassify3D

from projmatch_form import _defineProjectionMatchingParams

       
        
class XmippProtProjMatch(ProtRefine3D, ProtClassify3D):
    """ 3D reconstruction and classification using multireference projection matching"""

    _label = 'projection matching'

    #--------------------------- DEFINE param functions --------------------------------------------   
        
    def _defineParams(self, form):
        """ Since the form definition is very very large,
        we have do it in a separated function.
        """
        _defineProjectionMatchingParams(self, form)
         
         
    #--------------------------- INSERT steps functions --------------------------------------------  
    
    def _insertAllSteps(self):
        self._insertInitSteps()
        self._insertIterationSteps()
        self._insertEndSteps()
        
    def _insertInitSteps(self):
        """ Insert some initialization steps required
        before executing the steps per iteration. 
        """
        pass
                                                    
    def _insertIterationSteps(self):
        """ Insert several steps needed per iteration. """
        pass
                                                    
    def _insertEndSteps(self):
        """ Insert steps after the iteration steps 
        have been performed.
        """ 
        self._insertFunctionStep('createOutputStep')
                 
                                       
    #--------------------------- STEPS functions --------------------------------------------       
    
    def createOutputStep(self):
        print "output generated..........."


    #--------------------------- INFO functions -------------------------------------------- 
    
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return self._summary()  # summary is quite explicit and serve as methods
    
    #--------------------------- UTILS functions --------------------------------------------
    
