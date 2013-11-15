#!/bin/bash
# Update Xmipp version number in the files that require it (see FILE variable below),
# and cannot update themselves

FILE=Doxyfile
if [ ! -f $FILE ]
then
	echo Please run this script from XMIPP_HOME
	exit 2
fi

if [ $# -eq 1 ]; then
    NEW_VERSION=$1
else
    NEW_VERSION=$(xmipp_python -c 'import protlib_filesystem; print protlib_filesystem.getXmippVersion()' )
fi

sed -i "s/PROJECT_NUMBER\(.*\)=\(.*\)[0-9]\.[0-9]/PROJECT_NUMBER\1=\2${NEW_VERSION}/g" $FILE

xmipp_python -c "import protlib_filesystem; protlib_filesystem.updateXmippVersion(\"${NEW_VERSION}\")"

# this files get the version directly from VERSION, so they no longer need to be updated from here
# It is interesting to keep them here anyway for reference purposes
FILE=protocols/compile_gui.py
FILE=protocols/about_gui.py

exit 0
