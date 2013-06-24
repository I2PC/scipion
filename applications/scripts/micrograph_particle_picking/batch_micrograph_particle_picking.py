#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import ScriptAppIJ
from xmipp import Program, FileName
import xmipp
#512 is too low memory, should be 1024 the default.
class ScriptTrainingPicking(ScriptAppIJ):
    def __init__(self):
        ScriptAppIJ.__init__(self, 'xmipp.viewer.particlepicker.training.Main')
        
    def defineOtherParams(self):
        self.addParamsLine(" -o <directory>                                       : Output directory for load/save session data without updating model.");
        self.addParamsLine("    alias --output;                                   ");
        self.addParamsLine(" --mode <pick_mode=manual>                            : Mode in wich Particle Picker will be used");
        self.addParamsLine("    where <pick_mode>");
        self.addParamsLine("       manual  <thr=1> <fast=True> <incore=False>     : Enables manual mode. User will pick particles manually and can switch to supervised mode if desired.");
        self.addParamsLine("                                                      : Particles from manual mode will be used to train software and switch to supervised mode");
        self.addParamsLine("                                                      : Autopicker will use number of threads, fast and incore modes provided");
        self.addParamsLine("       review <file>                                  : Enables review mode. User reviews/corrects particles set provided on specified file");
        self.addParamsLine("                                                      : without updating model. ");
        self.addParamsLine("       readonly                                       : Enables readonly mode. User can see the picked particles, but cannot modify them.");
        
    def readOtherParams(self):
        input = self.getParam('-i')
        output = self.getParam('-o')
        mode = self.getParam('--mode')

        if (not os.path.exists(output)):
            os.makedirs(output)
        manual = (mode == 'manual')
        if manual:
            
            numberOfThreads = self.getIntParam('--mode', 1)
            fastMode = self.getParam('--mode', 2)
            incore = self.getParam('--mode', 3)
            
        review = (mode == 'review')
        if review:
            file = self.getParam('--mode', 1)

        self.args = "%(input)s %(output)s %(mode)s" %locals()
        if manual:
            self.args += " %(numberOfThreads)d %(fastMode)s %(incore)s" %locals()
        if review:
            self.args += " %(file)s" %locals()
        
if __name__ == '__main__':
    ScriptTrainingPicking().tryRun()

       

