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

import os
from protlib_xmipp import XmippScript
from protlib_import import convertCtfparam


class ScriptImportCtfparam(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Import a ctfparam file from Xmipp2.4 version')
        self.addUsageLine('The old style ctfparam file is converted to a new')
        self.addUsageLine('metadata file and each')
        self.addUsageLine('New Xmipp3.0 CTF definition can be found:')
        self.addUsageLine('http://i2pc.cnb.csic.es/emx/CargarDictionaryFormat.htm?type=Convention#ctf')
        ## params
        self.addParamsLine(' -i <old_ctfparam>          : Old Xmipp2.4 CTF definition file (.ctfparam)')
        self.addParamsLine(' -o <new_ctfparam>          : New CTF definition file')
        ## examples
        self.addExampleLine('   xmipp_import_ctfparam -i old.ctfparam -o new.ctfparam')
      
    def run(self):
        oldCtf = self.getParam('-i')
        newCtf = self.getParam('-o')        
        md = convertCtfparam(oldCtf)
        md.write(newCtf)
        
if __name__ == '__main__':
    ScriptImportCtfparam().tryRun()
