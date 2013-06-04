#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
#                Laura del Cano         (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
import sys
import os
from os.path import join, exists, isdir, relpath
from pyworkflow.utils.path import cleanPath, makePath
from pyworkflow.manager import Manager


if "SCIPION_HOME" not in os.environ:
    raise Exception("SCIPION_HOME is not defined as environment variable")

TESTS_HOME = join(os.environ['SCIPION_HOME'], 'tests')

def getInputPath(*filenames):
    """Return the path to the SCIPION_HOME/tests/input dir
    joined with filename"""
    return join(TESTS_HOME, "input", *filenames)

def getGoldPath(*filenames):
    """Return the path to the SCIPION_HOME/tests/gold dir
    joined with filename"""
    return join(TESTS_HOME, "gold", *filenames)

def getOutputPath(*filenames):
    """Return the path to the SCIPION_HOME/tests/output dir
    joined with filename"""
    return join(TESTS_HOME, "output", *filenames)

def getRelPath(filename):
    """Return the path relative to SCIPION_HOME/tests"""
    return relpath(filename, TESTS_HOME)

def setupOutput(test, outputDir):
    """ Define the output path for the calling test and 
    define a function to retrieve output path from this root. 
    """
    test.outputPath = getOutputPath(outputDir)
    cleanPath(test.outputPath)
    
def setupProject(testClass):
    """ Create and setup a project for this test. """
    testClass.projName = testClass.__name__
    testClass.proj = Manager().createProject(testClass.projName) # Now it will be loaded if exists
