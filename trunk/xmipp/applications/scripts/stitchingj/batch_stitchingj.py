#!/usr/bin/env xmipp_python

from protlib_xmipp import XmippScript
from protlib_utils import runJavaJar

class ScriptStitchingJ(XmippScript):
    def __init__(self):
        XmippScript.__init__(self)

    def defineParams(self):
        self.addUsageLine("Images stitching utility.")
        ## params
        self.addParamsLine('  [--memory <memory="512m">]	: Maximum memory.');
        self.addParamsLine('         alias -m;');
        self.addParamsLine('  --input <...>			: Images to stitch.');
        self.addParamsLine('         alias -i;');
        self.addParamsLine('  --output <filename>		: Filename to store result image.');
        self.addParamsLine('         alias -o;');
        self.addParamsLine('  --parameters <filename>		: Parameters for SIFT features extraction and images alignment (java properties format)');
        self.addParamsLine('         alias -p;');
        self.addParamsLine('  [--stack <filename>]		: Filename to store result as stack.');
        self.addParamsLine('         alias -s;');
        self.addParamsLine('  [-x <initial_x_coordinate>]		: Initial x position for second image.');
        self.addParamsLine('  [-y <initial_y_coordinate>]		: Initial y position for second image.');

    def readParams(self):
        self.args = ""
        self.mem = ""

        if self.checkParam('--memory'):
		    self.mem = self.getParam('--memory')
        else:
            self.mem = "512m"
        	print "No memory size provided. Using default: " + self.mem

        if self.checkParam('--input'):
        	#self.args += " -i %s" % self.getParam('--input')
		self.args = "-i %s" % ' '.join(self.getListParam('--input'))
        if self.checkParam('--output'):
        	self.args += " -o %s" % self.getParam('--output')
        if self.checkParam('--parameters'):
        	self.args += " -p %s" % self.getParam('--parameters')
        if self.checkParam('--stack'):
        	self.args += " -stack %s" % self.getParam('--stack')
        if self.checkParam('-x'):
        	self.args += " -x %s" % self.getParam('-x')
        if self.checkParam('-y'):
        	self.args += " -y %s" % self.getParam('-y')

    def run(self):
	   runJavaJar(self.mem, "external/Stitching/Stitching.jar", self.args, False)

if __name__ == '__main__':
        ScriptStitchingJ().tryRun()
