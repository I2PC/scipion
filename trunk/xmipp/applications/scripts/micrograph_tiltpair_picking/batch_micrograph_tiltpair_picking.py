#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import XmippScript, estimateMemory
from protlib_filesystem import getXmippPath

class ScriptParticlePicking(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)
        
    def defineParams(self):
        self.addUsageLine("Tilt pair picking utility.")
        ## params
        self.addParamsLine(" -i <micrograph_pairs_file>                           : Input metadata containing micrographs information");
        self.addParamsLine("    alias --input;")
        self.addParamsLine(" -o <directory>                                       : Output directory for load/save session data without updating model.");
        self.addParamsLine("    alias --output;                                   ");
        self.addParamsLine(' == Java options == ')
        self.addParamsLine(' [-m <mem="">]                                        : Memory amount for JVM');
        self.addParamsLine('    alias --memory;');
        
            
    
    def run(self):
        input = self.getParam('-i')
        output = self.getParam('-o')
        plugins_dir = getXmippPath("external/imagej/plugins/*")
        ij_jar = getXmippPath("external/imagej/ij.jar")
        memory = self.getParam('-m')
        if len(memory) == 0:
            memory=str(3*estimateMemory(input))+"m"
            print "No memory size provided. Using default: " + memory

        jar = "Xmipp_PP.jar"
        cmd = "java -Xmx%(memory)s -Dplugins.dir=%(plugins_dir)s -cp %(plugins_dir)s:%(ij_jar)s: particlepicker.tiltpair.Main %(input)s %(output)s" % locals()
        os.system(cmd)
    
if __name__ == '__main__':
    ScriptParticlePicking().tryRun()
