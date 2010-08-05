#!/usr/bin/env sh
XMIPP_PROTOCOLS=`which xmipp_protocols`
XMIPP_BASE=`dirname $XMIPP_PROTOCOLS`/..
# GETOPTS
ERROR=0
if [ $# -gt 0 ]
then
    MEM=$1
else
    MEM=512m
    ERROR=1
fi

if [ $# -gt 4 ]
then
    MACRO_ARGS="$2 $3 $4 $5"
else
    ERROR=`expr $ERROR + 2`
fi

if [ "$ERROR" != "0" ]
then
    if [ "$ERROR" != "2" ]
    then
        echo "No memory size provided. Using default: $MEM"
    fi
    if [ "$ERROR" != "1" ]
    then
        echo "Not enough arguments."
    fi
    echo "Usage: xmipp_projection_explorer <Memory size> -vol <volume_file> -angles <angles_file>"
fi

java -Xmx$MEM -Dplugins.dir=$XMIPP_BASE/external/imagej/plugins/ -Dmacros.dir=$XMIPP_BASE/external/imagej/macros/ -jar $XMIPP_BASE/external/imagej/ij.jar -macro xmippExplorer.txt "$MACRO_ARGS"
