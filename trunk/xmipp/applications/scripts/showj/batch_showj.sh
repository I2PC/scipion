#!/usr/bin/env sh
XMIPP_PROTOCOLS=`which xmipp_protocols`
XMIPP_BASE=`dirname $XMIPP_PROTOCOLS`/..
# GETOPTS
MEM=512m
java -Xmx$MEM -jar $XMIPP_BASE/external/imagej/ij.jar

