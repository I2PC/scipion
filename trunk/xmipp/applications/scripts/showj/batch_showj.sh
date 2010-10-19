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
			--mem)MEM=$2;READIMG=0;READSEL=0;READVOL=0;shift;;
			--img)READIMG=1;READSEL=0;READVOL=0;;	# Activates IMG input mode
			--sel)READIMG=0;READSEL=1;READVOL=0;;	# Activates SEL input mode
			--vol)READIMG=0;READSEL=0;READVOL=1;;	# Activates VOL input mode
			--poll)POLL=-poll;READIMG=0;READSEL=0;READVOL=0;;	# Sets polling
			---)shift; break;;
			-*)READIMG=0;READSEL=0;READVOL=0;echo >&2 \
				"Unknown parameter: $1"
				exit 1;;
			*)test "$READIMG" = "1" && IMG="$IMG $1";test "$READSEL" = "1" && SEL="$SEL $1";test "$READVOL" = "1" && VOL="$VOL $1";;
		esac
		shift
	done

	if test -z "$MEM" || test -z "$IMG" || test -z "$SEL" || test -z "$VOL"
	then
		SHOW_HELP=1
	fi

	if [ -z "$MEM" ]
	then
		MEM=512m	# Default memory value.
		echo "No memory size provided. Using default: $MEM"
	fi

	if [ -n "$IMG" ]
	then
		IMG="-img$IMG"
	fi

	if [ -n "$SEL" ]
	then
		SEL="-sel$SEL"
	fi

	if [ -n "$VOL" ]
	then
		VOL="-vol$VOL"
	fi

	if [ "$SHOW_HELP" = "1" ]
	then
		echo "Usage: xmipp_showj [--mem <memory_ammount>] [--img|vol|sel <file1 [file2 [..]]>] [--poll]"
	fi

	export LD_LIBRARY_PATH=$XMIPP_BASE/lib
	IMAGEJ_HOME=$XMIPP_BASE/external/imagej
	$JVM/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar -macro $IMAGEJ_HOME/macros/xmippBrowser.txt "$IMG $SEL $VOL $POLL"
fi
