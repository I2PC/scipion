# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              
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

import sys
import os
from pyworkflow.em.protocol.protocol_import import ProtImportImages, ProtImportFiles
import pyworkflow.protocol.params as params
from os.path import basename
from h5py import File
from pyworkflow.utils.path import removeExt, replaceExt

#from pyworkflow.em.convert import ImageHandler
#from pyworkflow.em import Protocol

#from pyworkflow.em.data import Volume
#from pyworkflow.em.protocol import ProtReconstruct3D
#from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
#from pyworkflow.utils import getFloatListFromValues
#from pyworkflow.utils.path import cleanPattern, cleanPath, copyFile
#import xmipp
#from pyworkflow.object import Float, String
#from math import sqrt
#from plotter import XmippPlotter


class ProtImportTiltSeries(ProtImportImages):
    """    
    This prtocol is to import tilt seies and related info included......
    """
    _label = 'import tilt series'  
    
    _outputClassName = 'TiltSeries'    
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):               
        #filesCondition = self._getFilesCondition()
        
        form.addSection(label='Import')

        form.addParam('filesPath', params.PathParam, 
                      #condition=filesCondition,
                      label="Files directory",
                      help="Directory with the files you want to import.")        
        form.addParam('copyFiles', params.BooleanParam, default=False, 
                      expertLevel=params.LEVEL_ADVANCED, 
                      label="Copy files?",
                      help="By default the files are not copied into the\n"
                           "project to avoid data duplication and to save\n"
                           "disk space. Instead of copying, symbolic links are\n"
                           "created pointing to original files. This approach\n"
                           "has the drawback that if the project is moved to\n"
                           "another computer, the links need to be restored.\n")
        
    
        
        #self._defineAcquisitionParams(form)
        
    #def _defineAcquisitionParams(self, form): ########################### in ghesmat dar akhar hamrah ba focal series kamel shavad
    #    """ Define acquisition parameters, it can be overriden
    #    by subclasses to change what parameters to include.
    #    """
    #    group = form.addGroup('Acquisition info')
    #    group.addParam('haveDataBeenPhaseFlipped', params.BooleanParam,
    #                   default=False,
    #                  label='Have data been phase-flipped?',
    #                  help='Set this to Yes if the images have been ctf-phase '
    #                       'corrected.')
    #    group.addParam('acquisitionWizard', params.LabelParam, important=True,
    #                   condition='importFrom != %d' % self.IMPORT_FROM_FILES,
    #                   label='Use the wizard button to import acquisition.',
    #                   help='Depending on the import Format, the wizard\n'
    #                        'will try to import the acquisition values.\n'
    #                        'If not found, required ones should be provided.')
    #    group.addParam('voltage', params.FloatParam, default=200,
    #               label=Message.LABEL_VOLTAGE, 
    #               help=Message.TEXT_VOLTAGE)
    #    group.addParam('sphericalAberration', params.FloatParam, default=2,
    #               label=Message.LABEL_SPH_ABERRATION, 
    #               help=Message.TEXT_SPH_ABERRATION)
    #    group.addParam('amplitudeContrast', params.FloatParam, default=0.1,
    #                  label=Message.LABEL_AMPLITUDE,
    #                  help=Message.TEXT_AMPLITUDE)
    #    group.addParam('magnification', params.IntParam, default=50000,
    #               label=Message.LABEL_MAGNI_RATE, 
    #               help=Message.TEXT_MAGNI_RATE)
    #    return group      
                      
        form.addParallelSection(threads=0, mpi=0)
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):                
        
        #self._insertFunctionStep(funcName, self.getPattern(),
        #                         self.voltage.get(),
        #                         self.sphericalAberration.get(),
        #                        self.amplitudeContrast.get(),
        #                         self.magnification.get())
        
        self._insertFunctionStep('importHdf5Step')         
    #--------------------------- STEPS functions --------------------------------------------
    
    def importHdf5Step(self):
        """ reading hdf5 file to use in xmipp_xray_import as an input """
        
        
        fnIn = self.filesPath.get()
        copyOrLink = self.getCopyOrLink()
        dst = self._getExtraPath(basename(fnIn))
        copyOrLink(fnIn, dst)
                
        fhHdf5 = File(dst, 'r')
        if not "TomoNormalized" in fhHdf5:
            print "Input data need to be normalized ..."
            self.runJob("xmipp_xray_import", 
                    "--mistral %s --oroot %s" % (dst, removeExt(dst)))
            imgSet = removeExt(dst) + '.mrc'
            print "imgSet=", imgSet
            
             
        else:
            print "Input data were already normalized ..."
            
                
        
        
        
        
    
    
    #--------------------------- INFO functions --------------------------------        
    def _validate(self):
        errors = []
        
        if not self.filesPath.get():
            errors.append("The path can not be empty!!!")
        
        return errors
        
    def _summary(self):######ba tabe e khorooji check shavad,,,mesle validateoverfitting##################################################### bayad eslah shavad
        summary = []
        #outputSet = self._getOutputSet()
        if not self.outPutReady:
            summary.append("Output " + self._outputClassName + " are not ready yet.") 
            if self.copyFiles:
                summary.append("*Warning*: You select to copy files into your project.\n"
                               "This will make another copy of your data and may take \n"
                               "more time to import. ")
        else:
            summary.append("TiltSeries are imported from %s" % (self.getFileName()))
            
        return summary
    
    def _methods(self):########################################################################  bayad eslah shavad
        methods = []
        outputSet = self._getOutputSet()
        if outputSet is not None:
            methods.append("Kino......")
            
        return methods    
    #--------------------------- UTILS functions -------------------------------
    def _getOutputName(self):
        # We assume that the import output is always a 'SetOfSomething'
        return self._outputClassName.replace('SetOf', 'output')
        
    def _getOutputSet(self):
        return getattr(self, self._getOutputName(), None)
    
        
    #--------------------------- BASE methods to be overriden ------------------        
    #def _getFilesCondition(self):
    #    """ Return an string representing the condition
    #    when to display the files path and pattern to grab
    #    files.
    #    """
    #    return '(importFrom == %d)' % self.IMPORT_FROM_FILES
