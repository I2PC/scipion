# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
This sub-package contains the XmippProtNormalizeStrain protocol
"""
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.protocol.params import MultiPointerParam
from pyworkflow.utils.path import createLink
import os
import xmipp


class XmippProtNormalizeStrain(ProtAnalysis3D):
    """
    Normalize the local strain and rotations amongst several runs
    """
    _label = 'normalize strain'
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputRuns', MultiPointerParam, pointerClass='XmippProtVolumeStrain',
                      label="input strain calculations",
                      help='Select the runs of strain calculations to be normalized')        

#--------------------------- INSERT steps functions --------------------------------------------  
                                
    def _insertAllSteps(self):
        self._insertFunctionStep('normalize',"localrot")
        self._insertFunctionStep('normalize',"strain")
            
    def normalize(self,what):
        # Get overall minimum and maximum
        V=xmipp.Image()
        minAll=1e38
        maxAll=-1e38
        for prot in self.inputRuns:
            protId=prot.get().getObjId()
            protDir=prot.get()._getPath('')
            fnVol=os.path.join(protDir,"extra","result_%s.vol"%what)
            V.read(fnVol)
            _, _, minVal, maxVal = V.computeStats()
            minAll = min(minAll,minVal)
            maxAll = max(maxAll,maxVal)
        
        # Write the Chimera file
        for prot in self.inputRuns:
            protId=prot.get().getObjId()
            protDir=prot.get()._getPath('')
            fnRoot=os.path.relpath(os.path.join(protDir,"extra","result"),self._getPath(''))
            scriptFile = self._getPath('%d_result_%s_chimera.cmd'%(protId,what))
            fhCmd = open(scriptFile, 'w')
            fhCmd.write("open %s\n" % (fnRoot+"_final.vol"))
            fhCmd.write("open %s\n" % (fnRoot+"_%s.vol"%what))
            fhCmd.write("vol #1 hide\n")
            fhCmd.write("scolor #0 volume #1 cmap rainbow cmapRange %f,%f reverseColors True\n"%(minAll,maxAll))
            fhCmd.close()
