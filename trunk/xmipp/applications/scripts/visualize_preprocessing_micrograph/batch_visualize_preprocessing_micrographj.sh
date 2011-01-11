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
	while [ $# -gt 0 ]
	do
		case "$1" in
			--mem)MEM=$2;READFILE=0;shift;;
			-i)READFILE=1;;	# Activates FILE input mode
			---)shift; break;;
			-*)READFILE=0;echo >&2 \
				"Unknown parameter: $1"
				exit 1;;
			*)test "$READFILE" = "1" && FILE="$FILE $1";;
		esac
		shift
	done

	if test -z "$MEM" || test -z "$FILE"
	then
		SHOW_HELP=1
	fi

	if [ -z "$MEM" ]
	then
		MEM=512m	# Default memory value.
		echo "No memory size provided. Using default: $MEM"
	fi

	if [ -n "$FILE" ]
	then
		FILE="-i$FILE"
	else
		echo "No file provided."
	fi

	if [ "$SHOW_HELP" = "1" ]
	then
		echo "Usage: visualize_preprocessing_micrographj [--mem <memory_ammount>] <-i file1 [file2 [..]]>"
	fi

	if [ ! -z "$FILE" ]
	then
		export LD_LIBRARY_PATH=$XMIPP_BASE/lib
		IMAGEJ_HOME=$XMIPP_BASE/external/imagej
		$JVM/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar -macro $IMAGEJ_HOME/macros/xmippVisualizeMicrograph.txt "$FILE"
	fi
fi
