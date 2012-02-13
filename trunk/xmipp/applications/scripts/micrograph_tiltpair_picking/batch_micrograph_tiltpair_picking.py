#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import ScriptAppIJ
from xmipp import Program, FileName
import xmipp

class ScriptTiltPairPicking(ScriptAppIJ):
    def __init__(self):
        ScriptAppIJ.__init__(self, 'xmipp.particlepicker.tiltpair.Main')
        
    def defineOtherParams(self):
        self.addParamsLine(" -o <directory>                                       : Output directory for load/save session data without updating model.");
        self.addParamsLine("    alias --output;                                   ");
        
    def readOtherParams(self):
        input = self.getParam('-i')
        output = self.getParam('-o')

        self.args = "%(input)s %(output)s" % locals()
        
        
if __name__ == '__main__':
    ScriptTiltPairPicking().tryRun()
