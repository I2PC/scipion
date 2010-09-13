#!/usr/bin/env sh
LINUX="java-6-sun-1.6.0.20_linux/"
LINUX64="jdk1.6.0_14"
MACOSX="java-6-sun-1.6.0.20_macosx/"

XMIPP_BASE="$HOME/xmipp"
CWD="$XMIPP_BASE/external/java"
JVM="$CWD/jvm"

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

# Sets JAVA_HOME...
JAVA_HOME="$CWD/$platform"

# Creates a symbolic link to the java virtual machine
test -e $JVM && rm $JVM
ln $JAVA_HOME -s $JVM

