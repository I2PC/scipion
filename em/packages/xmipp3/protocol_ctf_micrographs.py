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
from pyworkflow.em.packages.xmipp3.data import *
from pyworkflow.utils.path import makePath, basename, join, exists
from pyworkflow.utils import runJob


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

    def createOutput(self):
        # Create micrographs metadata with CTF information
        mdOut = FileName("Micrographs@" + 
                         self._getPath(self._getFilename('micrographs')))        
        micSet = XmippSetOfMicrographs(str(mdOut))
        
        # Label to set values in metadata
        labels = [MDL_PSD, MDL_PSD_ENHANCED, MDL_CTF_MODEL, MDL_IMAGE1, MDL_IMAGE2]
        # Filenames key related to ctf estimation
        keys = ['psd', 'enhanced_psd', 'ctfparam', 'ctfmodel_quadrant', 'ctfmodel_halfplane']
        
        for fn, micrographDir in self._iterMicrographs():
            mic = XmippMicrograph(fn)
            ctfparam = self._getFilename('ctfparam', micrographDir=micrographDir)
            if exists(ctfparam): # Get filenames
                values = [self._getFilename(key, micrographDir=micrographDir) for key in keys]
            else: # No files
                values = ['NA'] * len(labels)
    
            # Set values to the micrograph
            mic.setValue(*zip(labels, values))
            micSet.append(mic)
            
        micSet.sort()
        micSet.write()
            
        auxMdOut = FileName("Micrographs@" + 
                            self._getTmpPath(self._getFilename('micrographs')))
        runJob(None,"xmipp_ctf_sort_psds","-i %s -o %s" % (mdOut, auxMdOut))
        runJob(None,"mv","-f %s %s" % (auxMdOut.removeBlockName(),
                                       mdOut.removeBlockName()))

        # Create the SetOfMicrographs object on the database
        micSet.copyInfo(self.inputMics)
        # This property should only be set by CTF estimation protocols
        micSet._ctf.set(True)       
        self._defineOutputs(outputMicrographs=micSet)
