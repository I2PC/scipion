# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This sub-package contains the XmippCtfMicrographs protocol
"""


from pyworkflow.em import *  
from pyworkflow.utils import *
from pyworkflow.em.packages.xmipp3.data import *
from pyworkflow.utils.path import makePath, replaceBaseExt, join, basename


class XmippProtCTFMicrographs(ProtCTFMicrographs):
    """Protocol to perform CTF estimation on a set of micrographs in the project"""
    
    __prefix = join('%(micrographDir)s','xmipp_ctf')
    _templateDict = {
        # This templates are relative to a micrographDir
        'micrographs': 'micrographs.xmd',
        'prefix': __prefix,
        'ctfparam': __prefix +  '.ctfparam',
        'psd': __prefix + '.psd',
        'enhanced_psd': __prefix + '_enhanced_psd.xmp',
        'ctfmodel_quadrant': __prefix + '_ctfmodel_quadrant.xmp',
        'ctfmodel_halfplane': __prefix + '_ctfmodel_halfplane.xmp'
#        'ctffind_ctfparam': join('%(micrographDir)s', 'ctffind.ctfparam'),
#        'ctffind_spectrum': join('%(micrographDir)s', 'ctffind_spectrum.mrc')
        }

    _definition = DefCTFMicrographs()
        
    def __init__(self, **args):
        
        Protocol.__init__(self, **args)
        
    def _defineSteps(self):
        ''' insert the steps to perform ctf estimation on a set of micrographs
        '''
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get() 
        
        self.params = {'kV': self.inputMics.microscope.voltage.get(),
                       'Cs': self.inputMics.microscope.sphericalAberration.get(),
                       'sampling_rate': self.inputMics.samplingRate.get(),
                       'ctfmodelSize': 256,
                       'Q0': self.ampContrast.get(),
                       'min_freq': self.lowRes.get(),
                       'max_freq': self.highRes.get(),
                       'pieceDim': self.windowSize.get(),
                       'defocus_range': (self.maxDefocus.get()-self.minDefocus.get())*10000/2,
                       'defocusU': (self.maxDefocus.get()+self.minDefocus.get())*10000/2
                       }
        
        # For each micrograph insert the steps to process it
        for fn, micrographDir in self.__iterMicrographs():
            
            # CTF estimation with Xmipp
            self._insertFunctionStep('estimateCTF', fn, micrographDir)
                    
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')

    def estimateCTF(self, mic, micDir):
        ''' Run the estimate CTF program '''
        
        # Create micrograph dir under extra directory
        makePath(micDir)
            
        # CTF estimation with Xmipp
        args = "--micrograph " + mic + \
               " --oroot " + self._getFilename('prefix', micrographDir=micDir) + \
               " --fastDefocus"
                 
        for par, val in self.params.iteritems():
            args+= " --" + par + " " + str(val)
                
        runJob(None, 'xmipp_ctf_estimate_from_micrograph', args)
            
    
    def __iterMicrographs(self):
        """Iterate over micrographs and yield
        micrograph name and a directory to process """
        for mic in self.inputMics:
            fn = mic.getFileName()
            micrographDir = self._getExtraPath(removeExt(basename(fn))) 
            yield (fn, micrographDir)       

    def createOutput(self):
        # Create micrographs metadata with CTF information
        mdOut = self._getPath(self._getFilename('micrographs'))
        
        md = MetaData()
        for fn, micrographDir in self.__iterMicrographs():
            objId = md.addObject()
            md.setValue(MDL_MICROGRAPH, fn, objId)
            ctfparam = self._getFilename('ctfparam', micrographDir=micrographDir)
            labels = [MDL_PSD, MDL_PSD_ENHANCED, MDL_CTF_MODEL, MDL_IMAGE1, MDL_IMAGE2]
            if exists(ctfparam): # Get filenames
                keys = ['psd', 'enhanced_psd', 'ctfparam', 'ctfmodel_quadrant', 'ctfmodel_halfplane']
                values = [self._getFilename(key, micrographDir=micrographDir) for key in keys]
            else: # No files
                values = ['NA' for i in range(len(labels))]
    
            # Set values in metadata
            for label, value in zip(labels, values):
                md.setValue(label, value, objId)
        
        if not md.isEmpty():
            md.sort(MDL_MICROGRAPH)
            md.write(mdOut)
            
        dirOut,fnOut=os.path.split(mdOut)
        runJob(None,"xmipp_ctf_sort_psds","-i %s -o %s/aux_%s"%(mdOut,dirOut,fnOut))
        runJob(None,"mv","-f %s/aux_%s %s"%(dirOut,fnOut,mdOut))

        # Create the SetOfMicrographs object on the database
        self.outputMicrographs = XmippSetOfMicrographs(filename=mdOut)     
        self.outputMicrographs.copyInfo(self.inputMics)
        # This property should only be set by CTF estimation protocols
        self.outputMicrographs._ctf.set(True)
        
        self._defineOutputs(outputMicrographs=self.outputMicrographs)
