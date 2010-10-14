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
			--mem)MEM=$2;READIMGS=0; READSEL=0;shift;;
			--img)READIMGS=1;;	# Activates IMGS input mode
			--sel)READSEL=1;;
			--poll)POLL=-poll;READIMGS=0;READSEL=0;;	# Sets polling
			---)shift; break;;
			-*)READIMGS=0; READSEL=0;echo >&2 \
				"Unknown parameter: $1"
				exit 1;;
			*)test "$READIMGS" = "1" && IMGS="$IMGS $1";test "$READSEL" = "1" && SEL="$SEL $1";;
		esac
		shift
	done
echo "MEM: $MEM"
echo "IMGS: $IMGS"
echo "SEL: $SEL"
echo "POLL: $POLL"

	if test -z "$MEM" || test -z "$IMGS"
	then
		SHOW_HELP=1
	fi

	if [ -z "$MEM" ]
	then
		MEM=512m	# Default memory value.
		echo "No memory size provided. Using default: $MEM"
	fi

	if [ -n "$IMGS" ]
	then
		IMGS="-img$IMGS"
	fi

	if [ -n "$SEL" ]
	then
		SEL="-sel$SEL"
	fi

	echo "SEL: $SEL"

	if [ "$SHOW_HELP" = "1" ]
	then
		echo "Usage: xmipp_showj [--mem <memory_ammount>] [--img <file1 [file2 [..]]>] [--poll]"
	fi

	export LD_LIBRARY_PATH=$XMIPP_BASE/lib
	IMAGEJ_HOME=$XMIPP_BASE/external/imagej
	$JVM/bin/java -Xmx$MEM -Dplugins.dir=$IMAGEJ_HOME/plugins/ -jar $IMAGEJ_HOME/ij.jar -macro $IMAGEJ_HOME/macros/xmippBrowser.txt "$IMGS $SEL $POLL"
fi
