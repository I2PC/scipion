#!/usr/bin/env xmipp_python

from protlib_xmipp import ScriptPluginIJ

class ScriptMicrographViewerJ(ScriptPluginIJ):
	def __init__(self):
		ScriptPluginIJ.__init__(self, "XmippMicrographViewer.txt")

if __name__ == '__main__':
	ScriptMicrographViewerJ().tryRun()