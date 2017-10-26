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

import argparse
import datetime
import os
import sys
import sqlite3 as lite
import uuid

from pyworkflow.manager import Manager
import pyworkflow.utils as pwutils
from pyworkflow.em.data import SetOfCTF

PROJECTTABLE="project"
MAGNITUDETABLE="magnitude"
DATASETTABLE="dataset"

def getCtfStats(ctfSet):
    """ Compute some stats from a given set of ctf.
    so far resiolution related"""
    resolutionList = []
    averageResolution = 0
    for ctf in ctfSet:
        resolution = ctf.getResolution()
        resolutionList.append(resolution)
        averageResolution += resolution

    averageResolution /= len(ctfSet)

    return resolution, averageResolution


def valid_date(s):
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)


def processCommandLine():
    '''This function parses and return arguments passed in'''
    # Assign description to the help doc
    parser =\
        argparse.ArgumentParser(description='Script send processing '
                                            'statistics either to a '
                                            'webservice or to a local '
                                            'database',
                                epilog='--server and --dataBase are '
                                       'mutually exclusive')
    parser.add_argument('-p', '--project', type=str,
                        help='projectName to be processed',
                        required=True, default=None)
    parser.add_argument('-u', '--username', type=str, help='user name',
                        required=True, default="noNameUser")
    # parser.add_argument("-s",
    #                     "--startdate",
    #                     help="The Start Date - format YYYY-MM-DD",
    #                     type=valid_date, default="0001-01-01")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-w', '--webserviceurl', type=str, help='Webservice waiting for the data.',
        required=False)
    group.add_argument(
        '-d', '--database', type=str, help='DataBase used to store the data.',
        required=False)

    # Array for all arguments passed to script
    args = parser.parse_args()

    # Assign args to variables
    webserviceURL = args.webserviceurl
    projectName = args.project
    dataBaseName = args.database
    userName = args.username
    #startDate = args.startdate

    # Return all variable values
    return projectName, userName, webserviceURL, dataBaseName

def getProject(projectName):
    # Load a  project and return pointer to it
    manager = Manager()

    if not manager.hasProject(projectName):
        print("Unexistent project: %s" % pwutils.red(projectName))
        sys.exit(1)

    return manager.loadProject(projectName)  # return project

def _getUuidFileName(project):
    """project identifier is stored in this file"""
    return project.getLogPath("uuidStatistics.log")


def _getUuid(project):
    # Load (or create if not exits) a file
    # in the project Logs folder to store an unique
    # project identifier
    uuidFn = _getUuidFileName(project)
    try:
        with open(uuidFn) as f:
            uuidValue = f.readline()
    except IOError:
        uuidValue = str(uuid.uuid4())
        with open(uuidFn,'w') as f:
             f.write(uuidValue)

    return uuidValue



def getDataBaseFileName(project, dataBaseName='statistics.sqlite'):
    """proejct statistics are stored in this file (if saved on a file)"""
    return project.getLogPath(dataBaseName)


def getData(project):

    for run in project.getRuns():
        for outputName, output in run.iterOutputAttributes(SetOfCTF):
            print run.getRunName(), '-', outputName
            resolution, averageResolution = getCtfStats(output)
            return resolution, averageResolution
            # print "  defocus: (%0.4f - %0.4f)" % (minDefocus, maxDefocus)


def initializeDataBase(projectName, userName, projectUuid, conn):
    cur = conn.cursor()
    # create table project if it does not exist
    cur.execute("""CREATE TABLE IF NOT EXISTS  %s(
                timestamp DATE DEFAULT (datetime('now','localtime')),
                userName char(64),
                projectName char(128),
                projectUuid char(45) PRIMARY KEY-- project identifier, PK
                )
                """ % PROJECTTABLE)
    # insert project if new
    cur.execute("""INSERT OR IGNORE INTO %s
                          (userName, projectName, projectUuid)
                          VALUES ('%s','%s','%s')
                """% (PROJECTTABLE, userName, projectName, projectUuid))
    # create table magnitude
    cur.execute("""CREATE TABLE IF NOT EXISTS  %s(
                          magnitudeName char(128) PRIMARY KEY
                )
                """ % MAGNITUDETABLE)
    cur.execute("""INSERT OR IGNORE INTO %s
                          (magnitudeName)
                          VALUES ('%s')
                """% (MAGNITUDETABLE, 'ctfMaxResolution'))  #so far this is the only
                                                            # ctfMaxResolution but can be easily
                                                            # configured by the user
    # create table dataset
    cur.execute("""CREATE TABLE IF NOT EXISTS  %s(
                projectUuid char(45) REFERENCES %s(projectUuid),
                magnitudeName char(128) REFERENCES %s(magnitudeName),
                precomputedHistogram BLOB,
                dataArray BLOB,
                average FLOAT,
                nDataPoints INT
                )
                """ % (DATASETTABLE, PROJECTTABLE,  MAGNITUDETABLE ) )

# import sqlite, cPickle
# class Blob(object):
#     ''' automatic converter for binary strings '''
#     def _ _init_ _(self, s): self.s = s
#     def _quote(self): return "'%s'" % sqlite.encode(self.s)
#
# # make a test database in memory, get a cursor on it, and make a table
# connection = sqlite.connect(':memory:')
# cursor = connection.cursor( )
# cursor.execute("CREATE TABLE justatest (name TEXT, ablob BLOB)")
# # Prepare some BLOBs to insert in the table
# names = 'aramis', 'athos', 'porthos'
# data = {  }
# for name in names:
#     datum = list(name)
#     datum.sort( )
#     data[name] = cPickle.dumps(datum, 2)
# # Perform the insertions
# sql = 'INSERT INTO justatest VALUES(%s, %s)'
# for name in names:
#     cursor.execute(sql, (name, Blob(data[name])) )
# # Recover the data so you can check back
# sql = 'SELECT name, ablob FROM justatest ORDER BY name'
# cursor.execute(sql)
# for name, blob in cursor.fetchall( ):
#     print name, cPickle.loads(blob), cPickle.loads(data[name])
# # Done, close the connection (would be no big deal if you didn't, but...)
# connection.close( )

def closeDatabase(conn):
    conn.commit()
    conn.close()


def connectToDataBase(dataBaseName):
    conn = lite.connect(dataBaseName)
    return conn


def saveDataToDataBase():
    pass


def main(projectName, userName, webserviceURL, dataBaseName):
    # load the project changes the default directory
    if not os.path.isabs(dataBaseName):
        dataBaseName = os.path.join(os.getcwd(),dataBaseName)

    project = getProject(projectName)  # load project
    projectUuid = _getUuid(project)  # get project unique identifier
    if dataBaseName is not None:  #if output goes to database
        coon = connectToDataBase(dataBaseName)  #open cursor in data base
        # create tables -if needed- and insert project and category
        initializeDataBase(projectName, userName, projectUuid, coon) #create
    else:
        pass  # nothing to do here, database is already created
              # may be in the future allow category creation
              # initializeWebserver(projectName, userName, startDate, serverUrl)
    data, histogram = getData(project)#<<<<<<<<<<<<<<<<<
    average =
    histogram =

    #convert to string and send or store
    #create minimalistic test histogram
    # store user_name, projectName, date
    # TODO: until end
    # get resolution data
    #resolution, averageResolution = getData()
    # save to database


if __name__ == '__main__':
    projectName, userName, webserviceURL, dataBaseName = processCommandLine()
    main(projectName, userName, webserviceURL, dataBaseName)
    sys.exit(0)
