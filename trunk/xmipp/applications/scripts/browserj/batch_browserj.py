#!/usr/bin/env python

import optparse
from protlib_utils import runImageJPlugin

def command_line_options():
	""" add command line options here"""
	_usage="""Usage: xmipp_browserj [--mem <memory_ammount>] [--dir dir]"""
	parser = optparse.OptionParser(_usage)
	parser.add_option("-m", "--memory",  dest="memory", default="512m", help="Memory ammount for JVM")        
	parser.add_option("-d", "--dir", default="", dest="workdir", help="work directory")

	(options, args) = parser.parse_args()

	return (options.memory,options.workdir)

memory, workdir = command_line_options();

if memory == "512m":
	print "No memory size provided. Using default: " + memory

args = ""
if len(workdir) > 0:
	args += "-dir " + workdir

runImageJPlugin(memory, "xmippBrowser.txt", args)
