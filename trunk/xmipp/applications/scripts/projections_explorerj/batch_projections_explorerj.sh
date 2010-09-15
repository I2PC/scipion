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
	ERROR=0
	if [ $# -gt 3 ]
	then
	    if [ $# -gt 4 ]
	    then
		    MEM=$1
		    MACRO_ARGS="$2 $3 $4 $5"
	    else
		    ERROR=1
		    MEM=512m
		    MACRO_ARGS="$1 $2 $3 $4"
	    fi
	else
	    ERROR=2
	fi

	if [ "$ERROR" != "0" ]
	then
		echo "Usage: xmipp_projection_explorer <Memory size> -vol <volume_file> -angles <angles_file>. Example: xmipp_projection_explorer 1024m -vol file.vol -angles angles.txt"
		echo "No memory size provided. Using default: $MEM"
	fi

	if [ "$ERROR" = "2" ]
	then
		echo "Not enough arguments."
	else
		export LD_LIBRARY_PATH=$XMIPP_BASE/lib
		IMAGEJ_HOME=$XMIPP_BASE/external/imagej
		$JVM/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar -macro $IMAGEJ_HOME/macros/xmippExplorer.txt "$MACRO_ARGS"
	fi
fi
