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

from os.path import basename, splitext,exists
from os import remove
from protlib_xmipp import XmippScript
from xmipp import MetaData, MDL_CTF_MODEL, MD_APPEND, MD_OVERWRITE, FileName
from emx.emxmapper import *
from emx.emx import *
from protlib_emx import *


class ScriptImportEMX(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine("Convert  TO emx metadata files");
        self.addParamsLine(' -i <metadataXMIPP>             : Input xmipp metadata file ');
        self.addParamsLine(' [-o <metadataEMX>]               : Output emx  metadata file ');
#        self.addParamsLine(" [--binaryFile <fileName_mrc>]   : One input binary file ");
#        self.addParamsLine("     alias -b;");
        self.addParamsLine(' [--mode <mode=micCTF>]         : information to extract')
        self.addParamsLine("         where <mode>");
        self.addParamsLine("             alignment          : export particle shift and rotations");
        self.addParamsLine("             coordinates        : export particle coordinates (so far only works for a single image)");
        self.addParamsLine("             micCTF             : export micrograph ctf");
        self.addParamsLine("             micCTFChallenge    : export micrograph ctf to challenge format"); 
        self.addParamsLine("             particles          : export particles");       
        self.addParamsLine("     alias -m;");
#        self.addParamsLine(' [--amplitudeContrast <Q=0.1>]   : amplitudeContrast, mandatory when mode=micCTF');
#        self.addParamsLine("     alias -a;");
        self.addExampleLine("Export information from Metadata XmippFile file to EMX", False);
        self.addExampleLine("xmipp_export_emx -i microgaph.xmd -a 0.1");
        self.addExampleLine("input is the file micrograph.xmd created by screen micrograph protocol");
      
    def run(self):
        xmdFileName = FileName(self.getParam('-i'))
        emxFileName = self.getParam('-o')
        if (emxFileName=='metadataXMIPP.emx'):
            emxFileName = xmdFileName.withoutExtension()+'.emx'
        mode         = self.getParam('-m')

        if exists(emxFileName):
           remove(emxFileName)
        emxData = EmxData()
        mapper  = XmlMapper(emxData)
        #create emx files with mic CTF information
        if mode   == 'micCTF':
            ctfMicXmippToEmx(emxData,xmdFileName)
        elif mode == 'micCTFChallenge':
            ctfMicXmippToEmxChallenge(emxData,xmdFileName)
        elif mode == 'coordinates':
            coorrXmippToEmx(emxData,xmdFileName)
        elif mode == 'alignment':
            alignXmippToEmx(emxData,xmdFileName)
        mapper.writeEMXFile(emxFileName)
                    
        
if __name__ == '__main__':
    ScriptImportEMX().tryRun()
