#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import XmippScript
from protlib_utils import runJavaJar

class ScriptStitchingJ(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)

    def defineParams(self):
        self.addUsageLine("Images stitching utility.")
        ## params
        self.addParamsLine('  --memory <memory>			: Maximum memory. ');
        self.addParamsLine('         alias -m;');
        self.addParamsLine('  --input <filename1> <filename2>	: Images to stitch. ');
        self.addParamsLine('         alias -i;');
        self.addParamsLine('  --output <filename>		: Filename to store results. ');
        self.addParamsLine('         alias -o;');
        self.addParamsLine('  --parameters <filename>		: Parameters for SIFT features extraction and images alignment (java properties format)');
        self.addParamsLine('         alias -p;');

    def readParams(self):
	self.args = ""
	self.mem = ""

        if self.checkParam('--memory'):
		self.mem = self.getParam('--memory')
	else:
	        self.mem = "512m"
        	print "No memory size provided. Using default: " + self.mem

        if self.checkParam('--input'):
        	self.args += " -i %s" % self.getParam('--input')
        if self.checkParam('--output'):
        	self.args += " -o %s" % self.getParam('--output') 
        if self.checkParam('--parameters'):
        	self.args += " -p %s" % self.getParam('--parameters') 

    def run(self):
	runJavaJar(self.mem, "external/Stitching/Stitching.jar", self.args, False)

if __name__ == '__main__':
        ScriptStitchingJ().tryRun()
