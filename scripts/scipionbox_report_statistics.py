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
This script will compute the statistics of the SetOfCTFs in a given project.
"""

import sys, os, argparse, datetime, time

from pyworkflow.manager import Manager
import pyworkflow.utils as pwutils
from pyworkflow.em.data import SetOfCTF

def getCtfStats(ctfSet):
    """ Compute some stats from a given set of ctf.
    so far resiolution related"""
    resolutionList   = []
    averageResolution = 0
    for ctf in ctfSet:
        resolution = ctf.getResolution()
        resolutionList.append(resolution)
        averageResolution += resolution

    averageResolution /= len(ctfSet)

    return resolution, averageResolution

def valid_date(s):
    try:
        return datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)

def processCommandLine():
    '''This function parses and return arguments passed in'''
    # Assign description to the help doc
    parser = argparse.ArgumentParser(
    description='Script send processing statistics either to a webservice or to a local database',
        epilog='--server and --dataBase are mutually exclusive')
    parser.add_argument(
        '-p', '--project', type=str, help='projectName to be processed',
                 required=True, default=None)
    parser.add_argument(
        '-u', '--username', type=str, help='user name',
                 required=True, default="noNameUser")
    parser.add_argument("-s",
                        "--startdate",
                        help="The Start Date - format YYYY-MM-DD",
                        required=True,
                        type=valid_date)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-s', '--server', type=str, help='Webservice waiting for the data.',
        required=False)
    group.add_argument(
        '-d', '--database', type=str, help='DataBase used to store the data.',
        required=False)

    # Array for all arguments passed to script
    args = parser.parse_args()
    # Assign args to variables
    server = args.server
    project = args.project
    dataBase = args.database
    userName = args.username
    startDate = args.startdate


    # Return all variable values
    return project, userName, startDate, server, dataBase

def getData():
    # Create a new project
    manager = Manager()

    if not manager.hasProject(projName):
        print("Unexistent project: %s" % pwutils.red(projName))

    project = manager.loadProject(projName)

    for run in project.getRuns():
        for outputName, output in run.iterOutputAttributes(SetOfCTF):
            print run.getRunName(), '-', outputName
            resolution, averageResolution = getCtfStats(output)
            return resolution, averageResolution
            #print "  defocus: (%0.4f - %0.4f)" % (minDefocus, maxDefocus)


def initializeDataBase(project, userName, startDate):
    #create tables is do not exist: project, magnitude,dataset

    #store userName project microscope session

def saveDataToDataBase():
    pass


def main(project, userName, startDate, server, dataBase):
    #store user_name, projectName, date
    #TODO: until end
    # get resolution data
    resolution, averageResolution = getData()
    #save to database
    if dataBase is not None:
        initializeDataBase(project, userName, startDate)


if __name__ == '__main__':
    project, userName, startDate, server, dataBase = processCommandLine()
    main(project, userName, startDate, server, dataBase)


    exit(0)

