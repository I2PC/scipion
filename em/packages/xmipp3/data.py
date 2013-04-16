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
This modules contains basic hierarchy
for specific Xmipp3 EM data objects
"""

from pyworkflow.em import *   
from xmipp import MetaData, MDL_MICROGRAPH
    
class XmippSetOfMicrographs(SetOfMicrographs):
    
    def __iter__(self):
        """Iterate over the set of micrographs in the MetaData"""
        md = MetaData(self.getFileName())
        for objId in md:    
            m = Micrograph()
            m.setFileName(md.getValue(MDL_MICROGRAPH, objId))       
            yield m
        
    def convertFrom(self, setOfMics):
        """Method to convert a set of micrographs to be used in Xmipp"""
        md = MetaData()

        for mic in setOfMics:
            objId = md.addObject()
            md.setValue(MDL_MICROGRAPH, mic.getFileName(), objId)

        md.write(self.getFileName())
            
        # finalname = self._getTmpPath(replaceBaseExt(fn, "tmp.mrc"))
        
        
    
