#!/usr/bin/env python

import optparse
from protlib_utils import runImageJPlugin

def command_line_options():
	""" add command line options here"""
	_usage="""Usage: xmipp_showj [--mem <memory_ammount>] [-i file1 [-i file2 [..]]] [--poll]"""
	parser = optparse.OptionParser(_usage)
	parser.add_option("-m", "--memory",  dest="memory", default="", help="Memory ammount for JVM")        
	parser.add_option("-i", "--input", action="append", dest="inputFiles", help="input files to show")
	parser.add_option("-o", "--mode", dest="mode", default="image", help="Mode to open files: image, table")
	parser.add_option("-p", "--poll", action="store_true", dest="poll", default=False, help="Keep checking for changes on input files")

	(options, args) = parser.parse_args()

	return (options.memory,options.inputFiles, options.mode, options.poll)

memory, files, mode, poll = command_line_options();

args = ""
if files:
	args += "-i"
	for i in range(len(files)):
		args += " " + files[i]

if mode:
	args += " --mode " + mode

if poll:
	args += " --poll"

runImageJPlugin(memory, "xmippBrowser.txt", args)