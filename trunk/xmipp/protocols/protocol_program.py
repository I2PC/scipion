#!/usr/bin/env xmipp_python
'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

#------------------------------------------------------------------------------------------------
# Generic protocol for all Xmipp programs


from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
from protlib_xmipp import greenStr, redStr, blueStr
from protlib_sql import SqliteDb
from os.path import basename

class ProtXmippProgram(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.xmipp.name, scriptname, project)
        self.Import = 'from protocol_program import *'

    def validate(self):
        return []
        #return ["Protocol not implemented yet..."]
    
    def summary(self):
        return ["This is a test summary",
                "Need a real summary here..."
                ]
        
    def defineSteps(self):
        program = self.ProgramName
        params = buildCommandLine(self)
        self.NumberOfMpi = 1
        self.NumberOfThreads = 1
        #self.Db.insertStep('printCommandLine', programname=program, params=params)
        self.Db.insertStep('runJob', execution_mode=SqliteDb.EXEC_ALWAYS,
                             programname=program, 
                             params=params,
                             NumberOfMpi = self.NumberOfMpi,
                             NumberOfThreads = self.NumberOfThreads)
        
    def visualize(self):
        '''The generic visualization of programs just open the working dir'''
        from protlib_gui_ext import openLink
        openLink(self.WorkingDir)

# Some helper functions
def getParamName(paramId):
    if len(paramId) > 1: 
        return '--' + paramId
    return '-' + paramId

def printCommandLine(log, programname, params):
    print "Running program: ", programname
    print "         params: ", params
    
def buildCommandLine(prot):
    paramsLine = ""
    lastListValue = 'none'
    lastSubParamValue = None
    for k in dir(prot):
        if k.startswith('_K_'):
            from protlib_xmipp import redStr
            value = getattr(prot, k)
            # Check special cases of -o or --oroot 
            # and append the working dir to the value
            
            if ('_P_oroot_A_' in k or '_P_o_A_' in k) and value != "" and value == basename(value):
                value = prot.workingDirPath(value)
            if '_P_' in k: # Param variable (contains _P_)
                if '_L_' in k:  # This param is an option inside a list (OR group, contains _P_L_)
                    sep = '_P_L_'
                else:  # This is a common param, just '_P_'
                    sep = '_P_'
                myLine = ""
                key, suffix = k.split(sep)
                if '_A_' in suffix: # Arguments present (contains _A_)
                    args = suffix.split('_A_')
                    paramName = getParamName(args[0])
                    if len(args) > 2:
                        allowArg = lastSubParamValue == args[2]
                    else:
                        allowArg = True
                        lastSubParamValue = value
                        
                    if  allowArg and value != "":
                        if paramName not in paramsLine:
                            myLine = paramName
                        if 'arglist' not in args[1]:
                            value = escapeStr(value)
                        myLine += ' ' + value
                        
                else: #Param without Arguments, True or False value
                    paramName = getParamName(suffix)
                    if value:
                        myLine = paramName
                        
                if sep == '_P_': # normal param, not inside OR group
                    if lastListValue != 'none':
                        if lastListValue not in paramsLine:
                            paramsLine += " " + lastListValue # this is necessary for params inside OR group with no extra arguments
                        lastListValue = 'none'
                    paramsLine += " " + myLine
                elif paramName == lastListValue:
                    paramsLine += " " + myLine
            elif '_L_' in k: #exclusive list of some params
                lastListValue = value
    #Remove several white spaces            
    return ' '.join(paramsLine.split())

def escapeStr(value):
    '''If the value contains spaces and not start with ", add it'''
    if ' ' in value and not value.startswith('"'):
        value = '"%s"' % value
    return value
