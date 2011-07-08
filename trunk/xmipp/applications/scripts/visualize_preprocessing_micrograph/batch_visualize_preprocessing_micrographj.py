#!/usr/bin/env python

import optparse
from protlib_utils import runImageJPlugin

def command_line_options():
	""" add command line options here"""
	_usage="""Usage: visualize_preprocessing_micrographj [--mem <memory_ammount>] [-i file1 [-i file2 [..]]>"""
	parser = optparse.OptionParser(_usage)
	parser.add_option("-m", "--memory",  dest="memory", default="", help="Memory ammount for JVM")        
	parser.add_option("-i", "--input", action="append", dest="inputFiles", help="input files to show")

	(options, args) = parser.parse_args()

	return (options.memory,options.inputFiles)

memory, files = command_line_options();

args = ""
if files:
	args = "-i"
	for i in range(len(files)):
		args += " " + files[i]

runImageJPlugin(memory, "xmippVisualizeMicrograph.txt", args)
