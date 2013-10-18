#!/bin/bash
# Update Xmipp version number in the files that require it (see FILE variable below)

if [ $# -ne 1 ]
then
    echo Syntax: $0 new_xmipp_version
    exit 1
fi

NEW_VERSION=$1

FILE=Doxyfile
if [ ! -f $FILE ]
then
	echo Please run this script from XMIPP_HOME
	exit 2
fi

sed -i "s/PROJECT_NUMBER\(.*\)=\(.*\)[0-9]\.[0-9]/PROJECT_NUMBER\1=\2${NEW_VERSION}/g" $FILE

FILE=protocols/compile_gui.py
sed  -i "s/XMIPP_VERSION = \"[0-9]\.[0-9]\"/XMIPP_VERSION = \"${NEW_VERSION}\"/g" $FILE

FILE=protocols/about_gui.py
sed  -i "s/XMIPP_VERSION = \"[0-9]\.[0-9]\"/XMIPP_VERSION = \"${NEW_VERSION}\"/g" $FILE

FILE=VERSION.txt
echo "v${NEW_VERSION}" > $FILE

exit 0
