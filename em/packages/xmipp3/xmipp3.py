# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
This sub-package will contains Xmipp3.0 specific protocols
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_ORIGINAL, MDL_MICROGRAPH_TILTED, MDL_MICROGRAPH_TILTED_ORIGINAL


class XmippProtImportMicrographs(ProtImportMicrographs):
    pass


class XmippProtParticlePicking(ProtParticlePicking):
    pass


class XmippProtDownsampleMicrographs(ProtDownsampleMicrographs):
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        self.inputMicrographs = SetOfMicrographs(args.get('inputMicrographs'))
        self.downfactor = Integer(args.get('downfactor', 2))    
        self.IOTable = {}
        
    def defineSteps(self):
        '''for each micrograph call the downsampling function
        '''
        for mic in self.inputMicrographs:
            fn = mic.getFileName()
            fnOut = self.getPath(os.path.basename(f))
            self.insertFunctionStep('downsampleMicrograph', fn, fnOut,
                                    self.downfactor.get())
            self.IOTable[fn] = fnOut
        
        # create output objects       
        self.insertFunctionStep('createOutput')

    def downsampleMicrograph(self, mic, micOut, downfactor):
        """Downsample a micrograph with a downsampling factor.
        Accepts a micrograph path, the output path and a downsampling factor."""
        
        #TODO: Set log
        
        runJob(log,"xmipp_transform_downsample", "-i %(mic)s -o %(micOut)s --step %(downFactor)f --method fourier" % locals())
        
    def createOutput(self):
        self.downsampledmics = SetOfMicrographsXmipp()     
        self.downsampledmics.microscope.voltage.set(self.inputMicrographs.microscope.getVoltage())
        self.downsampledmics.microscope.aberration.set(self.inputMicrographs.microscope.getSphericalAberration())
        self.downsampledmics.sampling.set(self.inputMicrographs.sampling.get())
        
        #Create the xmipp metadata micrographs.xmd  
        mdOut = self.getPath("micrographs.xmd")
         
        md = xmipp.MetaData()      
        for i, v in IOTable:
            objId = famMD.addObject()
            MD.setValue(xmipp.MDL_MICROGRAPH,v,objId)
            MD.setValue(xmipp.MDL_MICROGRAPH_ORIGINAL,fnMicrograph,objId)
            #TODO: Handle Tilted micrographs
#            if tiltPairs:
#                MD.setValue(xmipp.MDL_MICROGRAPH_TILTED,IOTable[fnMicrographTilted],objId)
#                MD.setValue(xmipp.MDL_MICROGRAPH_TILTED_ORIGINAL,fnMicrographTilted,objId)
        MD.write("micrographs"+"@"+mdOut)
        
        self.downsampledmics.setFileName(mdOut)

        self.defineOutputs(micrograph=self.downsampledmics)
        pass

class XmippProtAlign(ProtAlign):
    pass


class XmippProtClassify(ProtClassify):
    pass


class XmippProtAlignClassify(ProtAlignClassify):
    pass