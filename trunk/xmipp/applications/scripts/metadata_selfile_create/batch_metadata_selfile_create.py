#!/usr/bin/env xmipp_python
"""/***************************************************************************
 *
 * Authors:     Roberto Marabini
 *
 * Universidad Autonoma de Madrid
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
"""
#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import XmippScript

class ScriptCreateMetadata(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Create a metadata from a file pattern.')
        ## params
        self.addParamsLine(' -p <pattern>          : Pattern to match')
        self.addParamsLine('   alias --pattern;')
        self.addParamsLine(' -o <metadata>         : Output metadata')
        self.addParamsLine('[ -s]                  : Check if the images are stacks')
        self.addParamsLine('   alias --isstack;')    
        ## examples
        self.addExampleLine('   xmipp_metadata_selfile_create -p "Images/*xmp" -o all_images.sel -isstack')
            
    def run(self):
        pattern = self.getParam('--pattern')
        isStack = self.checkParam('-s')
        import glob
        files = glob.glob(pattern)
        files.sort()
    
        from xmipp import MetaData, FileName, SingleImgSize, MDL_IMAGE, MDL_ENABLED
        mD = MetaData()
        inFile = FileName()
        
        nSize = 1
        for file in files:
            fileAux=file
            if isStack:
                if file.endswith(".mrc"):
                    fileAux=file+":mrcs"
                x, x, x, nSize = SingleImgSize(fileAux)
            if nSize != 1:
                counter = 1
                for jj in range(nSize):
                    inFile.compose(counter, fileAux)
                    objId = mD.addObject()
                    mD.setValue(MDL_IMAGE, inFile, objId)
                    mD.setValue(MDL_ENABLED, 1, objId)
                    counter += 1
            else:
                objId = mD.addObject()
                mD.setValue(MDL_IMAGE, fileAux, objId)
                mD.setValue(MDL_ENABLED, 1, objId)
                
        outFile = self.getParam('-o')
        mD.write(outFile)

if __name__ == '__main__':
    ScriptCreateMetadata().tryRun()
