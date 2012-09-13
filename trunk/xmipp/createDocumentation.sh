#!/bin/sh
#edit file Doxyfile and update PROJECT_NUMBER

# Call local SCONS to compile
xmipp_python external/scons/scons.py mode=docs "$@"
