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
from xmipp import MetaData, MDL_CTF_MODEL, MD_APPEND, MD_OVERWRITE
from protlib_import import convertCtfparam
from emx_data_model import MICROGRAPH, PARTICLE, EmxData
from emx_reader import EmxXmlReader
from lib_emx import ctfMicXmippFromEmx

class ScriptImportEMX(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine("Convert  from EMX metadata files");
        self.addParamsLine(' -i <text_file_emx>              : Input metadata file ');
        self.addParamsLine(" [--binaryFile <fileName_mrc>]   : One input binary file ");
        self.addParamsLine("     alias -b;");
#        self.addParamsLine(' [--path <val=emxImport>]       : directory for output files')
#        self.addParamsLine("     alias -p;");
        self.addParamsLine(' [--mode <mode=micCTF>]          : information to extract')
        self.addParamsLine("         where <mode>");
        self.addParamsLine("             micCTF              : extract micrograph ctf");
        self.addParamsLine("     alias -m;");
        #self.addKeyWords("import emx");
        self.addExampleLine("Import information from EMX file to Xmipp", False);
        self.addExampleLine("xmipp_import_emx -i particlePicking.emx -b mic.ctf ");
      
    def run(self):
        emxFileName  = self.getParam('-i')
        binFileName  = self.getParam('-b')
        mode         = self.getParam('-m')
        
        #emx class to store emx data
        emxData = EmxData()
        #emx class for reading emx data
        reader       = EmxXmlReader()
        reader.read(emxFileName,emxData)
        
        #create xmd files with mic CTF information and auxiliary files
        if mode == 'micCTF':
            ctfMicXmippFromEmx(emxData,emxFileName,oroot)

        #detect binary data type micrograph/particle
        #type =emxData.findObjectType(binFileName)
        #declare a 
#        if type == MICROGRAPH:
#            pass
#        elif type == PARTICLE:
#            pass
#        else:
#            raise Exception("Unknown object type associated with file=%s" % dictFK )
            
        
if __name__ == '__main__':
    ScriptImportEMX().tryRun()
