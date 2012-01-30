#!/usr/bin/env xmipp_python

import os
#from protlib_xmipp import XmippScript, estimateMemory
#from protlib_filesystem import getXmippPath
#
#class ScriptParticlePicking(XmippScript):
#    def __init__(self):
#        XmippScript.__init__(self)
#        
#    def defineParams(self):
#        self.addUsageLine("Tilt pair picking utility.")
#        ## params
#        self.addParamsLine(" -i <micrograph_pairs_file>                           : Input metadata containing micrographs information");
#        self.addParamsLine("    alias --input;")
#        self.addParamsLine(" -o <directory>                                       : Output directory for load/save session data without updating model.");
#        self.addParamsLine("    alias --output;                                   ");
#        self.addParamsLine(' == Java options == ')
#        self.addParamsLine(' [-m <mem="">]                                        : Memory amount for JVM');
#        self.addParamsLine('    alias --memory;');
#        
#            
#    
#    def run(self):
#        input = self.getParam('-i')
#        output = self.getParam('-o')
#        xmipp_path = getXmippPath();
#        plugins_dir = getXmippPath("external/imagej/plugins/")
#        ij_jar = getXmippPath("external/imagej/ij.jar")
#        memory = self.getParam('-m')
#        if len(memory) == 0:
#            memory=str(3*estimateMemory(input))+"m"
#            print "No memory size provided. Using default: " + memory
#        cmd = "java -Xmx%(memory)s -Dplugins.dir=%(plugins_dir)s -cp %(plugins_dir)s*:%(ij_jar)s:%(xmipp_path)s/java/lib/*: particlepicker.tiltpair.Main %(input)s %(output)s" % locals()
#        print cmd
#        os.system(cmd)
#    
#if __name__ == '__main__':
#    ScriptParticlePicking().tryRun()
from protlib_xmipp import ScriptAppIJ
from xmipp import Program, FileName
import xmipp

class ScriptParticlePicking(ScriptAppIJ):
    def __init__(self):
        ScriptAppIJ.__init__(self, 'xmipp.particlepicker.tiltpair.Main')
        
    def defineOtherParams(self):
#        self.addParamsLine(" -i <micrograph_pairs_file>                           : Input metadata containing micrographs information");
#        self.addParamsLine("    alias --input;")
        self.addParamsLine(" -o <directory>                                       : Output directory for load/save session data without updating model.");
        self.addParamsLine("    alias --output;                                   ");
        
    def readOtherParams(self):
        input = self.getParam('-i')
        output = self.getParam('-o')
        self.args = "%(input)s %(output)s" % locals()
        
        
if __name__ == '__main__':
    ScriptParticlePicking().tryRun()