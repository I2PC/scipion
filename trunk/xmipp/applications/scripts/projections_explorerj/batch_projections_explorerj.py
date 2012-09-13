#!/usr/bin/env xmipp_python

from protlib_xmipp import ScriptPluginIJ

class ScriptProjExplorer(ScriptPluginIJ):
	def __init__(self):
		ScriptPluginIJ.__init__(self, "XmippProjectionsExplorer.txt")
		
	def defineOtherParams(self):
		self.addParamsLine('  [--angles <anglesfile>]           : Associated euler angles files');
		self.addParamsLine('         alias -a;');		
			
	def readOtherParams(self):
		if self.checkParam('--angles'):
			self.args += " --angles %s" % self.getParam('--angles')
		
if __name__ == '__main__':
	ScriptProjExplorer().tryRun()
