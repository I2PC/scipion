#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import XmippScript, estimateMemory
from protlib_filesystem import getXmippPath

class ScriptParticlePicking(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine("Particles picking utility.")
        ## params
        self.addParamsLine(" -i <metadata>                                        : Input metadata containing micrographs information");
        self.addParamsLine("    alias --input;")
        self.addParamsLine(" -o <directory>                                       : Output directory for load/save session data");
        self.addParamsLine("    alias --output;")
        self.addParamsLine(' == Picking mode == ')
        self.addParamsLine(" --mode <pick_mode=manual>                            : Mode in wich Particle Picker will be used");
        self.addParamsLine("    where <pick_mode>");
        self.addParamsLine("       manual                                         : Enables manual mode. User will pick particles manually.");
        self.addParamsLine("       supervised <thr=1> <fast=True> <incore=False>  : Enables supervised mode. User will use autopicking. Then review/correct particles selected."); 
        self.addParamsLine("                                                      : Particles from manual mode can be used to train software and switch to supervised mode");
        self.addParamsLine("                                                      : Autopicker will use number of threads and fast and incore modes provided");
        self.addParamsLine("       review <file>                                  : Enables review mode. User reviews/corrects particles set provided on file");
        self.addParamsLine("                                                      : without updating model. ");
        self.addParamsLine(' == Java options == ')
        self.addParamsLine(' [-m <mem="">]                                        : Memory amount for JVM');
        self.addParamsLine('    alias --memory;');
        
    def run(self):
        input = self.getParam('-i')
        output = self.getParam('-o')
        plugins_dir = getXmippPath("external/imagej/plugins/")
        ij_jar = getXmippPath("external/imagej/ij.jar")
        memory = self.getParam('-m')
        mode = self.getParam('--mode')
        if len(memory) == 0:
            memory=str(2*estimateMemory(input))+"m"
            print "No memory size provided. Estimated: " + memory
        supervised = (mode == 'supervised')
        if supervised:
            numberOfThreads = self.getIntParam('--mode', 1)
            fastMode=self.getParam('--mode', 2)
            incore=self.getParam('--mode', 3)
        review = (mode == 'review')
        if review:
            file = self.getParam('--mode', 1)

        jar = "Xmipp_PP.jar"
        cmd = "java -Xmx%(memory)s -Dplugins.dir=%(plugins_dir)s -cp %(plugins_dir)s*:%(ij_jar)s: particlepicker.training.Main %(input)s %(output)s %(mode)s" % locals()
        if supervised:
            cmd+=" %(numberOfThreads)d %(fastMode)s %(incore)s"%locals()
        if review:
            cmd+=" %(file)s"%locals()
        print(cmd)
        os.system(cmd)
    
if __name__ == '__main__':
    ScriptParticlePicking().tryRun()
