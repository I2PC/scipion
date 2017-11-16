#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# * Authors:     roberto@cnb.csic.es
# *
# * [1] SciLifeLab, Stockholm University
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
# *  e-mail address 'delarosatrevin@scilifelab.se'
# *
# **************************************************************************

"""
This script will compute the statistics of the first SetOfCTFs
in a given project.
"""

# scipion python scripts scipionbox_report_statistics.py
# -p acquision.projname

import sys

from pyworkflow.manager import Manager
import pyworkflow.utils as pwutils
from pyworkflow.em.data import SetOfCTF
import json
import argparse


def getCtfStats(ctfSet):
    """ Compute some stats from a given set of ctf.
    so far resiolution related"""
    resolutionList = []
    defocusUList = []
    averageResolution = 0
    numberMovies = len(ctfSet)
    for ctf in ctfSet:
        resolution = ctf.getResolution()
        defocusU = ctf.getDefocusU()
        # remove inf and nan
        if resolution < 1000 and defocusU < 100000\
           and resolution > 0 and defocusU > -100000:
            resolutionList.append(resolution)
            defocusUList.append(defocusU)
            averageResolution += resolution

    averageResolution /= len(defocusUList)
    return resolutionList, defocusUList, averageResolution, numberMovies


def processCommandLine():
    '''This function parses and return arguments passed in'''
    # Assign description to the help doc
    parser = argparse.ArgumentParser(
         description='Script opens a project and '
                     'retrieve ctf resolution information')
    parser.add_argument(
        '-p', '--project', type=str, help='projectName to be processed',
        required=True, default=None)

    # Array for all arguments passed to script
    args = parser.parse_args()
    # Assign args to variables
    projectName = args.project
    return projectName


def getCTFData(projectName):
    """ get ctf related information"""
    manager = Manager()

    if not manager.hasProject(projectName):
        print("Unexistent project: %s" % pwutils.red(projectName))
        exit(0)

    project = manager.loadProject(projectName)

    for run in project.getRuns():
        for outputName, output in run.iterOutputAttributes(SetOfCTF):
            resolutionData, defocusData, averageResolution, numberMovies =
                getCtfStats(output)
            # return the first set of CTFs found
            return {"resolutionData": resolutionData,
                    "defocusData": defocusData,
                    "averageResolution": averageResolution,
                    "numberMovies": numberMovies}
    return {}


def main(projectName):
    # connect to scipion  and process the first setOfCTF
    # avaialble in project projecName
    return getCTFData(projectName)

if __name__ == '__main__':
    projectName = processCommandLine()
    resolutionDict = main(projectName)
    json.dump(resolutionDict, sys.stderr)
    exit(0)

