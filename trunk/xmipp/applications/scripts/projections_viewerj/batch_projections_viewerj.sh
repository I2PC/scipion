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
	if [ $# -le 1 ]
	then
	    echo "Not enough arguments."
	    echo "Usage: xmipp_projections_viewerj <Memory size> -vol <volume_file>. Example: xmipp_projections_viewerj 1024m -vol file.vol"
	else
		export LD_LIBRARY_PATH=$XMIPP_BASE/lib
		IMAGEJ_HOME=$XMIPP_BASE/external/imagej
		$JAVA_HOME/bin/java -Xmx$1 -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar -macro $IMAGEJ_HOME/macros/xmippProjectionsViewer.txt "$2"
	fi
fi
