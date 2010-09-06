#!/usr/bin/env sh
LINUX="java-6-sun-1.6.0.20_linux/"
LINUX64="java-6-sun-1.6.0.20_linux64/"
MACOSX="java-6-sun-1.6.0.20_macosx/"

XMIPP_PROTOCOLS=`which xmipp_protocols`
XMIPP_BASE=`dirname $XMIPP_PROTOCOLS`/..
JVM="$XMIPP_BASE"/external/java/"jvm"

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
test -z $JAVA_HOME && JAVA_HOME="$XMIPP_BASE"/external/java/"$platform"

# Creates a symbolic link to the java virtual machine
test -e $JVM && rm $JVM
ln $JAVA_HOME -s $JVM

