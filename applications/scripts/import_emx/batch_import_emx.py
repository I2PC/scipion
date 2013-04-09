#!/usr/bin/env xmipp_python
"""/***************************************************************************
 *
 * Authors:     Roberto Marabini          (roberto@cnb.csic.es)
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

from os.path import basename, splitext
from protlib_xmipp import XmippScript
from xmipp import MetaData, MDL_CTF_MODEL, MD_APPEND, MD_OVERWRITE, FileName
from emxLib.emxLib import ctfMicXmippToEmx
from emx.emxmapper import *
from emx.emx import *
from emxLib.emxLib import ctfMicEMXToXmipp, coorEMXToXmipp

class ScriptImportEMX(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine("Convert  from EMX metadata files");
        self.addParamsLine(' -i <text_file_emx>              : Input metadata file ');
        self.addParamsLine(' [--mode <mode=micCTF>]          : information to extract')
        self.addParamsLine("         where <mode>");
        self.addParamsLine("             micCTF              : import micrograph ctf");
        self.addParamsLine("             Coordinates         : import particle coordinates (so far only works for a single micrograph)");
        self.addParamsLine("     alias -m;");
        self.addExampleLine("Import information from EMX file to Xmipp", False);
        self.addExampleLine("xmipp_import_emx -i particlePicking.emx -m micCTF ");
      
    def run(self):
#        baseFilename = 'EMXread.emx'
#        scriptDir = os.path.dirname(__file__)
#        goldStandard = os.path.join(scriptDir,baseFilename)
#        emxData   = EmxData()
#        mapper    = EmxMapper(emxData, goldStandard, globals())
#        mapper.convertToEmxData(emxData)

        emxFileName  = self.getParam('-i')
        mode         = self.getParam('-m')
        #object to store emx data
        emxData   = EmxData()
        #object to read/write emxData
        mapper = XmlMapper(emxData)
        #read file
        mapper.readEMXFile(emxFileName)
#        xmlMapperW.writeEMXFile(fileName)
        #create xmd files with mic CTF information and auxiliary files
        if mode == 'micCTF':
            ctfMicEMXToXmipp(emxData,MICROGRAPH)
        elif mode == 'Coordinates':
            coorEMXToXmipp(emxData,PARTICLE,emxFileName)

if __name__ == '__main__':
    ScriptImportEMX().tryRun()
