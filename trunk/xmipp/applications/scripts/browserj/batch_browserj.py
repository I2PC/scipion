#!/usr/bin/env xmipp_python

from protlib_utils import runImageJPlugin
from protlib_xmipp import XmippScript

class ScriptBrowserJ(XmippScript):
	def defineParams(self):
		self.addParamsLine('  [--memory <mem="512m">]        : Memory ammount for JVM');
		self.addParamsLine('         alias -m;');
		self.addParamsLine('  --dir <directory="">           : List of params ');
		self.addParamsLine('         alias -d;');
			
	def readParams(self):
		self.memory = self.getParam('--memory')
		if self.memory == "512m":
			print "No memory size provided. Using default: " + self.memory
		workdir = self.getParam('--dir')
		self.args = ""
		if len(workdir) > 0:
			self.args += "-dir " + workdir
		
	def run(self):
		runImageJPlugin(self.memory, "XmippBrowser.txt", self.args)
			
if __name__ == '__main__':
	ScriptBrowserJ().tryRun()		
