import os, glob, sys, optparse

def command_line_options():
	""" add command line options here"""
	_usage="""Usage: xmipp_showj [--mem <memory_ammount>] [-i file1 [-i file2 [..]]] [--poll]"""
	parser = optparse.OptionParser(_usage)
	parser.add_option("-m", "--memory",  dest="memory", default="512m", help="Memory ammount for JVM")        
	parser.add_option("-i", "--input", dest="inputFile", help="input volume file to show")
	parser.add_option("-a", "--angles", dest="anglesFile", default="", help="associated euler angles file")

	(options, args) = parser.parse_args()

	return (options.memory, options.inputFile, options.anglesFile)

memory, inputFile, anglesFile = command_line_options();

if memory == "512m":
	print "No memory size provided. Using default: " + memory

if inputFile:
	angles_str = ""
	if anglesFile:
		angles_str = "-angles " + anglesFile

	imagej_home = "../../../external/imagej"
	plugins_dir = imagej_home + "/plugins/"
	macros_dir = imagej_home + "/macros/"
	imagej_jar = imagej_home + "/ij.jar"
	macro = macros_dir + "xmippExplorer.txt"
	cmd = "echo 'java -Xmx%s -Dplugins.dir=%s -jar %s -macro %s \"%s %s\"'" % (memory, plugins_dir, macro, imagej_jar, inputFile, angles_str)
	#$JVM/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar -macro $IMAGEJ_HOME/macros/xmippExplorer.txt "-vol $VOLFILE $ANGLESFILE"
	os.system(cmd)
