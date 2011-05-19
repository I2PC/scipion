import os, glob, sys, optparse

def command_line_options():
	""" add command line options here"""
	_usage="""Usage: xmipp_showj [--mem <memory_ammount>] [-i file1 [-i file2 [..]]] [--poll]"""
	parser = optparse.OptionParser(_usage)
	parser.add_option("-m", "--memory",  dest="memory", default="512m", help="Memory ammount for JVM")        
	parser.add_option("-d", "--dir", default="", dest="workdir", help="work directory")

	(options, args) = parser.parse_args()

	return (options.memory,options.workdir)

memory, workdir = command_line_options();

if memory == "512m":
	print "No memory size provided. Using default: " + memory

workdir_str = ""
if workdir:
	workdir_str = "-dir " + workdir

imagej_home = "../../../external/imagej"
plugins_dir = imagej_home + "/plugins/"
macros_dir = imagej_home + "/macros/"
imagej_jar = imagej_home + "/ij.jar"
macro = macros_dir + "xmippBrowser.txt"
cmd = "echo 'java -Xmx%s -Dplugins.dir=%s -jar %s -macro %s \"%s\"'" % (memory, plugins_dir, macro, imagej_jar, workdir_str)
#$JVM/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar -macro $IMAGEJ_HOME/macros/xmippBrowser.txt "-dir $WORKDIR"
os.system(cmd)
