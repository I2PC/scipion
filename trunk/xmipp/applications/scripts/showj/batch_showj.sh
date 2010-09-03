#!/usr/bin/env sh
XMIPP_PROTOCOLS=`which xmipp_protocols`
XMIPP_BASE=`dirname $XMIPP_PROTOCOLS`/..
JAVA_HOME=$XMIPP_BASE/external/java/jvm

if [ ! -e "$JAVA_HOME" ]
then
	echo "JVM not installed... fixing it..."
	$XMIPP_BASE/external/java/install_jvm.sh
fi

if [ ! -e "$JAVA_HOME" ]
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
	$JAVA_HOME/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar
fi
