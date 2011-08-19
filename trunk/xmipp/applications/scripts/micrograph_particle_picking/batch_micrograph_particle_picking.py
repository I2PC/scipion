#!/usr/bin/env python

import os
from protlib_xmipp import XmippScript
from protlib_filesystem import getXmippPath

class ScriptParticlePicking(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine("Particles picking utility.")
        ## params
        self.addParamsLine(" -i <metadata>          : Input metadata containing micrographs information")
        self.addParamsLine("   alias --input;")
        self.addParamsLine(" -o <directory>            : Output directory for load/save picking results")
        self.addParamsLine("   alias --output;")
        self.addParamsLine('  [--memory <mem="1024m">]              : Memory ammount for JVM');
        self.addParamsLine('         alias -m;');
        self.addParamsLine(' == Automatic picking == ')
        self.addParamsLine(' [--auto <thr=1> <fast=True> <incore=False>]  : Activates auto mode with a given number of threads')    
        self.addParamsLine('                                              : fast and incore modes')    
    
    def run(self):
        input = self.getParam('-i')
        output = self.getParam('-o')
        plugins_dir = getXmippPath("external/imagej/plugins/*")
        memory = self.getParam('--memory')
        if len(memory) == 0:
            memory = "1024m"
            print "No memory size provided. Using default: " + memory
        auto = self.checkParam('--auto')
        if auto:
            numberOfThreads = self.getIntParam('--auto',0)
            fastMode=self.getParam('--auto',1)
            incore=self.getParam('--auto',2)
        jar = "Xmipp_PP.jar"
        cmd = "java -Xmx%(memory)s -cp %(plugins_dir)s: Xmipp.Main %(input)s %(output)s" % locals()
        if auto:
            cmd+=" %(numberOfThreads)d %(fastMode)s %(incore)s"%locals()
        os.system(cmd)
    
if __name__ == '__main__':
    ScriptParticlePicking().tryRun()
