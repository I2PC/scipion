#!/usr/bin/env xmipp_python

from protlib_xmipp import ScriptPluginIJ

class ScriptShowJ(ScriptPluginIJ):
	def __init__(self):
		ScriptPluginIJ.__init__(self, "XmippBrowser.txt")
		
	def defineOtherParams(self):
		self.addParamsLine('  [--mode <mode_value=image>]           : List of params ');
		self.addParamsLine('     where <mode_value> image table')
		self.addParamsLine('         alias -o;');
		self.addParamsLine('  [--poll ]                            : Keep checking for changes on input files');
		self.addParamsLine('         alias -p;');
		self.addParamsLine('  [--rows <rows>]                            : number of rows in table');
		self.addParamsLine('         alias -r;');
		self.addParamsLine('  [--columns <columns>]                            : number of columns in table');
		self.addParamsLine('         alias -c;');
			
	def readOtherParams(self):
		if self.checkParam('--mode'):
			self.args += " --mode %s" % self.getParam('--mode')
		if self.checkParam('--poll'):
			self.args += " --poll"
		if self.checkParam('--rows'):
			self.args += " --rows %s" % self.getParam('--rows') 
		if self.checkParam('--columns'):
			self.args += " --columns %s" % self.getParam('--columns') 
		
if __name__ == '__main__':
	ScriptShowJ().tryRun()

