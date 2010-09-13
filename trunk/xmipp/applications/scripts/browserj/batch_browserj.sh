#!/usr/bin/env sh
XMIPP_BASE="$HOME/xmipp"
CWD="$XMIPP_BASE/external/java"
JAVA_HOME="$CWD/jvm"

if [ ! -L "$JAVA_HOME" ]
then
	echo "JVM not installed... fixing it..."
	$XMIPP_BASE/external/java/install_jvm.sh
fi

if [ ! -L "$JAVA_HOME" ]
then
	echo "JVM is missing, so program can't be run. Check your XMIPP installation."
else
	# GETOPTS
	ERROR=0
	if [ $# -gt 0 ]
	then
	    MEM=$1
	else
	    MEM=512m
	    ERROR=1
	fi

	if [ $# -gt 1 ]
	then
	    WORKDIR="$2"
	else
	    WORKDIR=`pwd`
	    ERROR=`expr $ERROR + 2`
	fi

	if [ "$ERROR" != "0" ]
	then
    		if [ "$ERROR" != "2" ]
		then
			echo "No memory size provided. Using default: $MEM"
		fi
		if [ "$ERROR" != "1" ]
		then
			echo "Not enough arguments."
		fi
		echo "Usage: xmipp_browserj <Memory size> <work_directory>. Example: xmipp_browserj 1024m $HOME"
	fi

	export LD_LIBRARY_PATH=$XMIPP_BASE/lib
	IMAGEJ_HOME=$XMIPP_BASE/external/imagej
	$JAVA_HOME/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar -macro $IMAGEJ_HOME/macros/xmippBrowser.txt "$WORKDIR"
fi
