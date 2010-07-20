#!/usr/bin/env sh
XMIPP_PROTOCOLS=`which xmipp_protocols`
XMIPP_BASE=`dirname $XMIPP_PROTOCOLS`/..
# GETOPTS
if [ $# -gt 0 ]
then
    MEM=$1
else
    MEM=512m
    echo "No memory size provided. Using default: $MEM"
fi

java -Xmx$MEM -Dplugins.dir=$XMIPP_BASE/external/imagej/plugins -jar $XMIPP_BASE/external/imagej/ij.jar

