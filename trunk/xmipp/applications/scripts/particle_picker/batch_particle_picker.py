#!/usr/bin/env python

import os
from protlib_xmipp import XmippScript
from protlib_filesystem import getXmippPath

class ScriptParticlePicker(XmippScript):
	def __init__(self):
		XmippScript.__init__(self)
		
	def defineParams(self):
		self.addUsageLine("Particles picking utility.")
		## params
		self.addParamsLine(" -i <metadata>          : Input metadata containing micrographs information")
		self.addParamsLine("   alias --input;")
		self.addParamsLine(" -o <directory>	        : Output directory for load/save picking results")
		self.addParamsLine("   alias --output;")
		self.addParamsLine('  [--memory <mem="1024m">]			  : Memory ammount for JVM');
		self.addParamsLine('		 alias -m;');	
	
	def run(self):
		input = self.getParam('-i')
		output = self.getParam('-o')
		plugins_dir = getXmippPath("external/imagej/plugins/*")
		memory = self.getParam('--memory')
		if len(memory) == 0:
			memory = "1024m"
			print "No memory size provided. Using default: " + memory
		jar = "Xmipp_PP.jar"
		cmd = " java -Xmx%(memory)s -cp %(plugins_dir)s: model.MyRunnable %(input)s %(output)s " % locals()
		os.system(cmd)
	
if __name__ == '__main__':
	ScriptParticlePicker().tryRun()

