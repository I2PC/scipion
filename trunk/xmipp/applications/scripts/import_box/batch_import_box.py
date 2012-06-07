#!/usr/bin/env xmipp_python
"""/***************************************************************************
 *
 * Authors:     J.M. de la Rosa Trevin    (jmdelarosa@cnb.csic.es)
 *             
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

from os.path import basename, exists, join
from glob import glob
from xmipp import MetaData, MDL_MICROGRAPH, Image, HEADER
from protlib_xmipp import XmippScript
from protlib_filesystem import replaceFilenameExt
from protlib_import import convertBox


class ScriptImportBox(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Import an EMAN picking .box file into Xmipp coordinates')
        self.addUsageLine('This program take the left-upper coner coordinates of the EMAN boxer')
        self.addUsageLine('And convert to particle center coordinate in Xmipp metadata file')
        self.addUsageLine('The input is the EMAN coordinate file and the output the Xmipp metadata file')
        ## params
        self.addParamsLine(' -i <micMetadata> <boxFolder>        : Input should be the micrographs metadata and the EMAN .box folder')
        self.addParamsLine(' -o <outputFolder>                   : Output folder where to generate output .pos files')
        self.addParamsLine('[--family <family=DefaultFamily>]   : Family name in Xmipp to which import particles')
        self.addParamsLine('[--particle_size <bs>] : box size of particles, by default use the one on EMAN')
        ## examples
        
        self.addExampleLine('Import .box files from folder Input/BoxFiles/ to particle picking run_001', False)
        self.addExampleLine('xmipp_import_box -i ParticlePicking/Manual/run_001/micrographs.xmd Input/BoxFiles/ -o ParticlePicking/Manual/run_001/')
      
    def run(self):
        micFn = self.getParam('-i', 0)
        boxFolder = self.getParam('-i', 1)
        family =  self.getParam('--family')
        ysize = particleSize = None
#        
#        if self.checkParam('--invert'):
#            ysize = self.getIntParam('--invert')
        
        if self.checkParam('--particle_size'):
            particleSize = self.getIntParam('--particle_size')
            
        outputFolder = self.getParam('-o')        
        micMd = MetaData(micFn)
        img = Image()
        
        for objId in micMd:
            micFile = micMd.getValue(MDL_MICROGRAPH, objId)
            boxFile = join(boxFolder, replaceFilenameExt(basename(micFile), '.box'))
            if exists(boxFile):
                img.read(micFile, HEADER)
                ysize = img.getDimensions()[1]
                posFile = join(outputFolder, replaceFilenameExt(basename(micFile), '.pos'))
                convertBox(boxFile, posFile, ysize, family, particleSize)
                
if __name__ == '__main__':
    ScriptImportBox().tryRun()
