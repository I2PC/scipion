#!/usr/bin/env python

from pyworkflow.em.packages.xmipp3 import ScriptPluginIJ

class ScriptTomo(ScriptPluginIJ):
	def __init__(self):
		ScriptPluginIJ.__init__(self, "XmippTomo.txt")
		
	def defineOtherParams(self):
		pass	
		# self.addParamsLine('  [--angles <anglesfile="">]           : Associated euler angles files');
		# self.addParamsLine('         alias -a;');		
			
	def readOtherParams(self):
		pass
		# if self.checkParam('--angles'):
	# 		self.args += " --angles %s" % self.getParam('--angles')
		
if __name__ == '__main__':
	ScriptTomo().tryRun()
