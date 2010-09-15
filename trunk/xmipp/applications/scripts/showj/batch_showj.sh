#!/usr/bin/env sh
SCRIPTPATH=`readlink -f $0`
CWD=`dirname $SCRIPTPATH`
XMIPP_BASE="$CWD/.."
JAVA_HOME="$XMIPP_BASE/external/java"
JVM="$JAVA_HOME/jvm"

if [ ! -L "$JVM" ]
then
	echo "JVM not installed... fixing it..."
	$JAVA_HOME/install_jvm.sh
fi

if [ ! -L "$JVM" ]
then
	echo "JVM is missing, so program can't be run. Check your XMIPP installation."
else
	# GETOPTS
	if [ $# -gt 0 ]
	then
	    MEM=$1
	else
	    MEM=512m
	    echo "No memory size provided. Using default: $MEM"
	    echo "Usage: xmipp_showj <Memory size>. Example: xmipp_showj 1024m"
	fi

	export LD_LIBRARY_PATH=$XMIPP_BASE/lib
	IMAGEJ_HOME=$XMIPP_BASE/external/imagej
	$JVM/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar
fi
