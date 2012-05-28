#!/usr/bin/env xmipp_python
"""/***************************************************************************
 *
 * Authors:     Carlos Oscar Sorzano
 *              J. M. de la Rosa Trevin
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

import os
from protlib_xmipp import XmippScript
from protlib_utils import reportError 
from protlib_filesystem import getXmippPath

class ScriptCompile(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine('Compile a C++ program using Xmipp libraries')
        ## params
        self.addParamsLine(' -i <cpp_file>          : C++ file to compile')
        self.addParamsLine('   alias --input;')
        self.addParamsLine(' [--debug]              : Compile with debugging flags')
        ## examples
        self.addExampleLine('Compile myprogram.cpp', False)
        self.addExampleLine('xmipp_compile myprogram.cpp')

    def run(self):        
        fn = self.getParam('-i')
        from os.path import splitext
        [fnBase,ext]=splitext(fn)
        if ext!=".cpp":
            reportError(fn+" is not a .cpp file")
        command='g++ -z -LD_MSG:off=15,134,85 ';
        if self.checkParam("--debug"):
            command +="-g -pg";
        xmippDir=getXmippPath();
        if xmippDir=="":
           reportError("$XMIPP_HOME is not set in the environment")
        command+=" -o "+fnBase+" "+fn+" -O -D_LINUX "+\
                 "-L"+xmippDir+"/lib "+\
                 "-I"+xmippDir+"/libraries "+\
                 "-I"+xmippDir+" "+\
                 "-lXmippClassif -lXmippData -lXmippExternal -lXmippInterface -lXmippRecons -lfftw3 -lfftw3_threads -lsqlite3 -ltiff -ljpeg"
        os.system(command)

if __name__ == '__main__':
    ScriptCompile().tryRun()
