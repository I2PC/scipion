#!/usr/bin/env sh
LINUX="java-6-sun-1.6.0.20_linux/"
LINUX64="jdk1.6.0_14"
MACOSX="java-6-sun-1.6.0.20_macosx/"

SCRIPT=`readlink -f $0`
SCRIPTPATH=`dirname $SCRIPT`
XMIPP_BASE="$SCRIPTPATH/../.."
JAVA_HOME="$XMIPP_BASE/external/java"
JVM="$SCRIPTPATH/jvm"

echo $JVM

# Which is our OS?
case "$(uname -s)" in
Darwin)
	platform=$MACOSX;;
Linux)
	case "$(uname -m)" in
		x86_64) platform=$LINUX64;;
		*) platform=$LINUX;;
	esac;;
esac

# Creates a symbolic link to the java virtual machine
rm $JVM
ln "$JAVA_HOME/$platform" -s $JVM

