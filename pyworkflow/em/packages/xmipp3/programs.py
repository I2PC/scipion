#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from glob import glob

import pyworkflow.utils as pwutils
from xmipp3 import getXmippPath


# Create some short-cut functions to colors
failStr = pwutils.red
greenStr = pwutils.green


#------------- FUNCTION TO WORK WITH PROGRAMS META-INFORMATION -----------------
class LabelData():
    def __init__(self):
        pass
    
def getXmippLabels():
    ''' Parse the labels definition from the 'libraries/data/metadata_label.h' file '''
    labelHeader = getXmippPath('libraries', 'data', 'metadata_label.h')
    f = open(labelHeader)
    labels = []
    comments = {}
    labelsPrefixes = ['MDL_', 'RLN_', 'BSOFT_']
    
    for line in f:
        line = line.strip()
        for prefix in labelsPrefixes:
            if line.startswith(prefix) and '///<' in line:
                parts = line.split('///<')
                mdl = parts[0].strip()[:-1] # remove last comma
                comm = parts[1].strip()
                comments[mdl] = comm
            if line.startswith('MDL::addLabel(' + prefix):
                l = line.find('(')
                r = line.find(')')
                parts = line[l + 1:r].split(',')
                label = {}
                label['name'] = parts[2].replace('"', '').strip()
                label['type'] = parts[1].strip()
                label['enum'] = parts[0].strip()
                label['comment'] = comments.get(label['enum'], "")
                labels.append(label)
    return labels

def getXmippLabelsName():
    ''' Parse the labels name from the 'libraries/data/metadata_label.h' file '''
    labelHeader = getXmippPath('libraries', 'data', 'metadata_label.h')
    f = open(labelHeader)
    labels = []
    comments = {}
    
    for line in f:
        line = line.strip()
        if line.startswith('MDL_') and '///<' in line:
            parts = line.split('///<')
            mdl = parts[0].strip()[:-1] # remove last comma
            comm = parts[1].strip()
            comments[mdl] = comm
        if line.startswith('MDL::addLabel(MDL_'):
            l = line.find('(')
            r = line.find(')')
            parts = line[l + 1:r].split(',')
            labels.append(parts[2].replace('"', '').strip())
    return labels

def getXmippPrograms():
    '''Return the list of Xmipp's programs, taken from from bin/ folder'''     
    programs = [os.path.basename(p) for p in glob(getXmippPath('bin', 'xmipp_*'))]
    programs.sort()
    return programs

#FIXME: this is only while development
def skipProgram(programName):
    if programName in ['xmipp_sqlite3', 'xmipp_mpi_steps_runner',
                       'xmipp_angular_commonline', 'python',
                       'xmipp_transform_threshold', 'xmipp_mpi_write_test', 'xmipp_chimera_client',
                       'xmipp_imagej','xmipp_mpi_image_common_lines', 'xmipp_mpi_classify_CLTomo', 'xmipp_classify_CLTomo']:
        return True
    for p in ['xmipp_test', 'xmipp_template']:
        if programName.find(p) != -1:
            return True
    return False

def getProgramsDbName():
    return getXmippPath('.xmipp_programs.sqlite')

#Some helper functions
def createProgramsDb(dbName=None):
    if dbName is None:
        dbName = getProgramsDbName()
    db = ProgramDb(dbName)
    print 'Created db with name: %(dbName)s' % locals()
    db.create()
    #Create categories dictionary to classify programs
    #looking in program name and category prefixes
    categories = db.selectCategories()
    categoryDict = {}
    for c in categories:
        prefixes = c['prefixes'].split()
        for p in prefixes:
            categoryDict[p] = c
            
    programs = getXmippPrograms()
    for p in programs:
        p = os.path.basename(p)
        try:
            print greenStr(p), skipProgram(p)
            
            if not skipProgram(p):
                cmd = [p, "--xmipp_write_definition"]
                if p.find('_mpi') != -1:                    
                    cmd = ['mpirun', '-np', '1'] + cmd
                print ' '.join(cmd)
                from subprocess import Popen, PIPE
                ps = Popen(cmd, stdout=PIPE, stderr=PIPE)
                stderrdata = ps.communicate()[1]
                if stderrdata != '':
                    raise Exception(stderrdata)
                for prefix, category in categoryDict.iteritems():
                    if prefix in p:
                        db.updateProgramCategory(p, category)
                        break
        except Exception, e:
            print failStr("PROGRAM: " + p)
            print failStr("ERROR: " + str(e))
    labels = getXmippLabels()
    for l in labels:
        db.insertLabel(l)
    db.commit()
    return db


def createProgramsAutocomplete(script='.xmipp_programs.autocomplete'):
    programs = getXmippPrograms()
    
    if os.path.exists(script):
        os.remove(script)
        
    for p in programs:
        p = os.path.basename(p)
        try:
            print greenStr(p), skipProgram(p)
            
            if not skipProgram(p):
                cmd = [p, "--xmipp_write_autocomplete", script]
                if '_mpi' in p:                    
                    cmd = ['mpirun', '-np', '1'] + cmd
                print ' '.join(cmd)
                from subprocess import Popen, PIPE
                ps = Popen(cmd, stdout=PIPE, stderr=PIPE)
                stderrdata = ps.communicate()[1]
                if stderrdata != '':
                    raise Exception(stderrdata)                
        except Exception, e:
            print failStr("PROGRAM: " + p)
            print failStr("ERROR: " + str(e))

class ProgramKeywordsRank():
    def __init__(self, keywords=None):
        self.keywords = keywords
        self.weights = {'name': 5, 'keywords': 3, 'usage': 1}
        
    def getRank(self, program):
        if not self.keywords:
            return 1 
        rank = 0
        for k in self.keywords:
            for wkey, wvalue in self.weights.iteritems():
                if program[wkey].find(k) != -1:
                    rank += wvalue
        return rank
    

def escapeStr(text):
    return "'%s'" % text.replace("'", "''") 



class ProgramDb():
    ''' Class to handle Programs DB in the Xmipp installation folder '''
    DB_TIMEOUT = 1000 # ms for timeout waiting
    
    def __init__(self, dbName=None):
        if dbName is None:
            dbName = getProgramsDbName()
        self.dbName = dbName
        from sqlite3 import dbapi2 as sqlite
        self.connection = sqlite.Connection(dbName, timeout=self.DB_TIMEOUT)
        self.connection.row_factory = sqlite.Row
        self.cursor = self.connection.cursor()
        self.cursor.execute('pragma foreign_keys=ON')
            
    def create(self):
        self.createTables()
        
    def commit(self):
        self.connection.commit()
        
    def createTables(self):
            sqlCommand = """DROP TABLE IF EXISTS Category;
                            CREATE TABLE Category (
                               id INTEGER PRIMARY KEY ASC AUTOINCREMENT, 
                               name TEXT UNIQUE, 
                               desc TEXT, 
                               prefixes TEXT);
                            INSERT INTO Category VALUES(NULL, 'Classification', NULL, 'classify_ ml_ mlf_');
                            INSERT INTO Category VALUES(NULL, 'CTF', NULL, 'ctf_');
                            INSERT INTO Category VALUES(NULL, 'Images', NULL, 'image_ micrograph_');
                            INSERT INTO Category VALUES(NULL, 'Metadatas', NULL, 'metadata_');
                            INSERT INTO Category VALUES(NULL, 'Phantoms', NULL, 'phantom_ pdb_');
                            INSERT INTO Category VALUES(NULL, 'Angular assignment', NULL, 'angular_');
                            INSERT INTO Category VALUES(NULL, 'Tomography', NULL, 'tomo_ xray_');
                            INSERT INTO Category VALUES(NULL, 'Transformations', NULL, 'transform_');
                            INSERT INTO Category VALUES(NULL, 'Volumes', NULL, 'volume_ reconstruct_ resolution_');                             
                            
                            DROP TABLE IF EXISTS Program;
                            CREATE TABLE Program (
                               id INTEGER PRIMARY KEY ASC AUTOINCREMENT,
                               category_id INTEGER, 
                               name TEXT UNIQUE,
                               usage TEXT,
                               examples TEXT,
                               keywords TEXT);
                               
                            DROP TABLE IF EXISTS Label;
                            CREATE TABLE Label (
                               id INTEGER PRIMARY KEY ASC AUTOINCREMENT,
                               name TEXT UNIQUE,
                               type TEXT,
                               enum TEXT UNIQUE,
                               comment TEXT);
                         """
            self.cursor.executescript(sqlCommand)
            self.connection.commit()
            
    def insertProgram(self, program):
        program['desc'] = escapeStr(program['desc'])
        sqlCommand = """INSERT INTO Program VALUES (
                           NULL, %(category_id)d, %(name)s, %(desc)s, %(keywords)s);
                     """ % program
        self.cursor.execute(sqlCommand)
                     
    def selectPrograms(self, category=None):
        categoryWhere = ""
        if category:
            categoryWhere = "WHERE category_id=%d" % category['id']
        sqlCommand = "SELECT * FROM Program %s ORDER BY name;""" % categoryWhere
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchall()
    
    def selectProgram(self, program_name):
        sqlCommand = "SELECT * FROM Program WHERE name='%s';""" % program_name
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchone()
        
    def selectCategories(self):
        sqlCommand = "SELECT * FROM Category;"
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchall()
    
    def updateProgramCategory(self, program_name, category):
        sqlCommand = "UPDATE Program SET category_id = %d WHERE name='%s';" % (category['id'], program_name)
        self.cursor.execute(sqlCommand)
        self.connection.commit()
        
    def selectLabels(self):
        sqlCommand = "SELECT * FROM Label;"
        self.cursor.execute(sqlCommand)
        return self.cursor.fetchall()
    
    def insertLabel(self, labelData):
        labelData['comment'] = escapeStr(labelData['comment'])
        sqlCommand = """INSERT INTO Label VALUES (
                           NULL, '%(name)s', '%(type)s', '%(enum)s', %(comment)s);
                     """ % labelData
        self.cursor.execute(sqlCommand)


