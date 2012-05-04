#!/usr/bin/env xmipp_python
"""/***************************************************************************
 *
 * Authors:     Roberto Marabini          (roberto@cnb.csic.es)
 *              J.M. de la Rosa Trevin    (jmdelarosa@cnb.csic.es)
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
from xmipp import MetaData, MDL_CTFMODEL, MD_APPEND, MD_OVERWRITE
from protlib_import import convertCtfparam


class ScriptImportCtfdat(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Import a ctfparam file from Xmipp2.4 version')
        self.addUsageLine('The old style ctfparam file is converted to a new')
        self.addUsageLine('metadata file and each ctfparam is also converted')
        self.addUsageLine('Two files will be generated: ')
        self.addUsageLine('$oroot.ctfdat with new ctfdata')
        self.addUsageLine('$oroot.ctfparam with all ctfparam in different blocks')
        self.addUsageLine('New Xmipp3.0 CTF definition can be found:')
        self.addUsageLine('http://i2pc.cnb.csic.es/emx/CargarDictionaryFormat.htm?type=Convention#ctf')
        ## params
        self.addParamsLine(' -i <old_ctfdat>          : Old Xmipp2.4 CTF list file (.ctfdat)')
        self.addParamsLine(' --oroot <oroot>          : Output files prefix')
        ## examples
        self.addExampleLine('   xmipp_import_ctfparam -i old.ctfparam -o new.ctfparam')
      
    def run(self):
        oldCtfdat = self.getParam('-i')
        oroot = self.getParam('--oroot')
        md = MetaData()
        md.readPlain(oldCtfdat, "image CTFModel")        
        ctfparamsDict = {}
        mode = MD_OVERWRITE
        
        for objId in md:
            oldCtfparam = md.getValue(MDL_CTFMODEL, objId)
            if oldCtfparam in ctfparamsDict:
                newCtfparam = ctfparamsDict[oldCtfparam]
            else:
                block = splitext(basename(oldCtfparam))[0]
                newCtfparam = "%(block)s@%(oroot)s.ctfparam" % locals()
                md2 = convertCtfparam(oldCtfparam)
                md2.write(newCtfparam, mode)
                ctfparamsDict[oldCtfparam] = newCtfparam
                mode = MD_APPEND
            md.setValue(MDL_CTFMODEL, newCtfparam, objId)
            
        md.write("%s.ctfdat" % oroot)
        
if __name__ == '__main__':
    ScriptImportCtfdat().tryRun()
