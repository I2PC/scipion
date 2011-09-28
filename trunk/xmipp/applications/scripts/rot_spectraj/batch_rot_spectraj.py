#!/usr/bin/env xmipp_python

from protlib_xmipp import ScriptPluginIJ

class ScriptRotSpectraJ(ScriptPluginIJ):
	def __init__(self):
		ScriptPluginIJ.__init__(self, "XmippRotSpectraViewer.txt")

	def defineParams(self):
	        self.addParamsLine('  [--memory <mem="512m">]              : Memory ammount for JVM');
        	self.addParamsLine('         alias -m;');
		self.addParamsLine('  [--vectors <vectorsfile>]           : Vectors file ');
		self.addParamsLine('         alias -f;');
		self.addParamsLine('  [--classes <classesfile>]                            : Classes file');
		self.addParamsLine('         alias -c;');
		self.addParamsLine('  [--data <datafile>]                            : Vectros data file');
		self.addParamsLine('         alias -d;');

	def readParams(self):
		self.memory = self.getParam('--memory')
		if self.memory == "512m":
			print "No memory size provided. Using default: " + self.memory

		self.args = ""

		if self.checkParam('--vectors'):
			self.args += " --vectors %s" % self.getParam('--vectors')
		if self.checkParam('--classes'):
			self.args += " --classes %s" % self.getParam('--classes') 
		if self.checkParam('--data'):
			self.args += " --data %s" % self.getParam('--data') 

if __name__ == '__main__':
	ScriptRotSpectraJ().tryRun()

