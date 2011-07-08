#!/usr/bin/env python

import optparse
from protlib_utils import runImageJPlugin

def command_line_options():
	""" add command line options here"""
	_usage="""Usage: xmipp_showj [--mem <memory_ammount>] [-i file1 [-i file2 [..]]] [--poll]"""
	parser = optparse.OptionParser(_usage)
	parser.add_option("-m", "--memory",  dest="memory", default="", help="Memory ammount for JVM")        
	parser.add_option("-i", "--input", dest="inputFile", help="input volume file to show")
	parser.add_option("-a", "--angles", dest="anglesFile", default="", help="associated euler angles file")

	(options, args) = parser.parse_args()

	return (options.memory, options.inputFile, options.anglesFile)

memory, inputFile, anglesFile = command_line_options();

args = ""
if inputFile:
	args += " -vol " + inputFile
	if anglesFile:
		args += " -angles " + anglesFile

runImageJPlugin(memory, "xmippExplorer.txt", args)
