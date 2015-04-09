#!/usr/bin/env python

import os
from protlib_xmipp import ScriptAppIJ
from xmipp import Program, FileName
import xmipp

class ScriptTiltPairPicking(ScriptAppIJ):
    def __init__(self):
        ScriptAppIJ.__init__(self, 'xmipp.viewer.particlepicker.tiltpair.TiltPairPickerRunner')
        
    def defineOtherParams(self):
        self.addParamsLine(" -o <directory>                                       : Output directory for load/save session data without updating model.");
        self.addParamsLine("    alias --output;                                   ");
        self.addParamsLine(" --mode <pick_mode=manual>                            : Mode in wich Particle Picker will be used");
        self.addParamsLine("    where <pick_mode>");
        self.addParamsLine("       manual                                         : Enables manual mode. User will pick particles manually.");
        self.addParamsLine("       readonly                                       : Enables readonly mode. User can see the picked particles, but cannot modify them.");
        
    def readOtherParams(self):
        input = self.getParam('-i')
        output = self.getParam('-o')
        mode = self.getParam('--mode')

        self.args = "--input %(input)s --output %(output)s --mode %(mode)s" % locals()
        
        
if __name__ == '__main__':
    ScriptTiltPairPicking().tryRun()
