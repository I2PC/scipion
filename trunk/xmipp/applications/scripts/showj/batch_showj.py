#!/usr/bin/env xmipp_python

import os
from protlib_xmipp import ScriptAppIJ
from xmipp import Program, FileName
import xmipp

class ScriptShowJ(ScriptAppIJ):
	def __init__(self):
		ScriptAppIJ.__init__(self, 'xmipp.viewer.Viewer')
		
	def defineOtherParams(self):
		self.addParamsLine('  [--mode <mode_value=image>]           : List of params ')
		self.addParamsLine('     where <mode_value> image gallery metadata rotspectra')
		self.addParamsLine('         alias -d;')
		self.addParamsLine('  [--poll]                            : Keeps checking for changes on input files  (for image mode only!)')
		self.addParamsLine('         alias -p;')
		self.addParamsLine('  [--render]	: Activates images rendering (for metadata mode only)')
		self.addParamsLine('         alias -e;')
		self.addParamsLine('  [--rows <rows>]                            : number of rows in table')
		self.addParamsLine('         alias -r;')
		self.addParamsLine('  [--columns <columns>]                            : number of columns in table')
		self.addParamsLine('         alias -c;')
		self.addParamsLine('  [--zoom <zoom>]                            : zoom for images')
		self.addParamsLine('         alias -z;')
		self.addParamsLine('  [--view <axis="z">]                        : Viewer position (for volumes only)')
		self.addParamsLine('     where <axis> z y x z_pos y_pos x_pos')
		self.addParamsLine('  [--debug] : debug')
		
	
	def readOtherParams(self):	
		if self.checkParam('--mode'):
			self.args += " --mode %s" % self.getParam('--mode') 	
		if self.checkParam('--poll'):
			self.args += " --poll"
		if self.checkParam('--render'):
			self.args += " --render"
		if self.checkParam('--rows'):
			self.args += " --rows %s" % self.getParam('--rows') 
		if self.checkParam('--columns'):
			self.args += " --columns %s" % self.getParam('--columns') 
		if self.checkParam('--zoom'):
			self.args += " --zoom %s" % self.getParam('--zoom') 
		if self.checkParam('--debug'):
			self.args += " --debug"
		if self.checkParam('--view'):
			self.args += " --view %s" % self.getParam('--view')
		
if __name__ == '__main__':
	ScriptShowJ().tryRun()

