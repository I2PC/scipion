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

from os.path import basename, splitext
from protlib_xmipp import XmippScript
from xmipp import MetaData, MDL_CTF_MODEL, MD_APPEND, MD_OVERWRITE, FileName,XmippError
from emxLib.emxLib import ctfMicXmippToEmx
from emx.emxmapper import *
from emx.emx import *
from emxLib.emxLib import ctfMicEMXToXmipp, coorEMXToXmipp, alignEMXToXmipp

class ScriptImportEMX(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine("Convert  from EMX metadata files")
        self.addParamsLine(' -i <text_file_emx>              : Input metadata file ')
        self.addParamsLine(' [--mode <mode=micCTF>]          : information to extract')
        self.addParamsLine("         where <mode>")
        self.addParamsLine("             alignment           : export particle shift and rotations")
        self.addParamsLine("             coordinates         : import particle coordinates (so far only works for a single micrograph)")
        self.addParamsLine("             micCTF              : import micrograph ctf")
        self.addParamsLine("     alias -m;")
        self.addParamsLine(' [--notValidate]                 : do not validate xml against schema');
        self.addParamsLine('                                 : validation may require net connection');
        self.addParamsLine('                                 : and program xmllint installed');
        self.addParamsLine('     alias -n;');
        self.addParamsLine(' [--onlyValidate]                : do not import');
        self.addParamsLine('     alias -y;');
        self.addParamsLine(' [--schema <schema=default>]     : validate again this schema ');
        self.addParamsLine("     alias -s;")
        self.addExampleLine(' default= http://sourceforge.net/p/emexchange/code/ci/master/tree/trunk/resourcesEmx/schemas/emx.xsd?format=raw>');
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
        doNotValidate= self.checkParam('-n')
        onlyValidate = self.checkParam('-y')
        schema       = self.getParam('-s')
        #validate emx file
        if not doNotValidate:
            try:
                if schema == 'default':
                    (code, _out,_err)=validateSchema(emxFileName)
                else:
                    (code, _out,_err)=validateSchema(emxFileName,schema)
                if code !=0:
                   raise Exception (_err) 
            except Exception, e:
                print >> sys.stderr, "XMIPP_ERROR -1:", str(e)
                exit(-1)
            if onlyValidate:
                print emxFileName, " was successfully validated"
                exit(0)
        #object to store emx data
        emxData   = EmxData()
        #object to read/write emxData
        mapper = XmlMapper(emxData)
	try:
		#read file
		mapper.readEMXFile(emxFileName)
	#        xmlMapperW.writeEMXFile(fileName)
		#create xmd files with mic CTF information and auxiliary files
		if mode == 'micCTF':
		    ctfMicEMXToXmipp(emxData)
		elif mode == 'coordinates':
		    coorEMXToXmipp(emxData,PARTICLE,emxFileName)
		elif mode == 'alignment':
		    xmdFileName = emxFileName.replace(".emx",".xmd")
		    alignEMXToXmipp(emxData,PARTICLE,xmdFileName)
	except Exception, e:
                print >> sys.stderr, "XMIPP_ERROR -1:", str(e)
                exit(-1)

if __name__ == '__main__':
    ScriptImportEMX().tryRun()
    exit(0)
