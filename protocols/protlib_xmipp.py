#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
 '''

import os
import xmipp
from xmipp import Program, FileName
from protlib_utils import runImageJPlugin, runJavaIJapp
from protlib_filesystem import getXmippPath, xmippExists


class XmippScript():
    ''' This class will serve as wrapper around the XmippProgram class
    to have same facilities from Python scripts'''
    def __init__(self, runWithoutArgs=False):
        self._prog = Program(runWithoutArgs)
        
    def defineParams(self):
        ''' This function should be overwrited by subclasses for 
        define its own parameters'''
        pass
    
    def readParams(self):
        ''' This function should be overwrited by subclasses for 
        and take desired params from command line'''
        pass
    
    def checkParam(self, param):
        return self._prog.checkParam(param)
    
    def getParam(self, param, index=0):
        return self._prog.getParam(param, index)
    
    def getIntParam(self, param, index=0):
        return int(self._prog.getParam(param, index))
    
    def getDoubleParam(self, param, index=0):
        return float(self._prog.getParam(param, index))
    
    def getListParam(self, param):
        return self._prog.getListParam(param)
    
    def addUsageLine(self, line, verbatim=False):
        self._prog.addUsageLine(line, verbatim)

    def addExampleLine(self, line, verbatim=True):
        self._prog.addExampleLine(line, verbatim)
        
    def addParamsLine(self, line):
        self._prog.addParamsLine(line)
    
    def run(self):
        ''' This function should be overwrited by subclasses and
        it the main body of the script'''   
        pass
     
    def tryRun(self):
        ''' This function should be overwrited by subclasses and
        it the main body of the script'''
        try:
            import sys
            self.defineParams()
            doRun = self._prog.read(sys.argv)
            if doRun:
                self.readParams()
                self.run()
        except Exception:
            import traceback
            traceback.print_exc(file=sys.stderr)
            
class ScriptIJBase(XmippScript):
    def __init__(self, name):
        XmippScript.__init__(self)
        self.name = name
    
    def defineOtherParams(self):
        pass
    
    def readOtherParams(self):
        pass
    
    def defineParams(self):
        self.addParamsLine('  --input <...>                 : Input files to show');
        self.addParamsLine('         alias -i;');
        self.addParamsLine('  [--memory <mem="2g">]              : Memory ammount for JVM');
        self.addParamsLine('         alias -m;');
        self.defineOtherParams()
    
    def readInputFiles(self):
        self.inputFiles = self.getListParam('-i')
        
    def readParams(self):
        self.readInputFiles()
        self.memory = self.getParam('--memory')
        files = []
        missingFiles = []
        for f in self.inputFiles:
            if xmippExists(f):
                files.append(f)
            else:
                missingFiles.append(f)
        self.inputFiles = files
        if len(missingFiles):
            print "Missing files: \n %s" % '  \n'.join(missingFiles) 
        self.args = "-i %s" % ' '.join(self.inputFiles)
        self.readOtherParams()
 
class ScriptPluginIJ(ScriptIJBase):
    def __init__(self, macro):
        ScriptIJBase.__init__(self, macro)
                  
    def run(self):
        runImageJPlugin(self.memory, self.name, self.args)
     
class ScriptAppIJ(ScriptIJBase):
    def __init__(self, name):
        ScriptIJBase.__init__(self, name)
                  
    def run(self):
        if len(self.inputFiles) > 0:
            runJavaIJapp(self.memory, self.name, self.args, batchMode=False)
        else:
            print "No input files. Exiting..."
        
#------------- FUNCTION TO WORK WITH PROGRAMS META-INFORMATION -----------------    
def getImageData(img):
    ''' Function to get a matrix from an Image'''
#    xdim, ydim, z, n = img.getDimensions()
#    from numpy import zeros
#    Z = zeros((ydim, xdim), float)
#    for y in range(ydim):
#        for x in range(xdim):
#    #TODO: improve by avoiding use of getPixel
#            Z[y, x] = img.getPixel(y, x)
    Z = img.getData()
    return Z
    
def getFirstImage(mdFn):
    ''' Read the first image from Metadata'''
    md = xmipp.MetaData(mdFn)
    return md.getValue(xmipp.MDL_IMAGE, md.firstObject())

#------------- FUNCTION TO WORK WITH PROGRAMS META-INFORMATION -----------------
class LabelData():
    def __init__(self):
        pass
    
def getXmippLabels():
    ''' Parse the labels definition from the 'libraries/data/metadata_label.h' file '''
    labelHeader = getXmippPath(os.path.join('libraries', 'data', 'metadata_label.h'))
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
            label = {}
            label['name'] = parts[2].replace('"', '').strip()
            label['type'] = parts[1].strip()
            label['enum'] = parts[0].strip()
            label['comment'] = comments.get(label['enum'], "")
            labels.append(label)
    return labels

def getXmippPrograms():
    '''Return the list of Xmipp's programs, taken from from bin/ folder'''     
    from glob import glob
    programs = [os.path.basename(p) for p in glob(os.path.join(getXmippPath(), 'bin', 'xmipp_*'))]
    programs.sort()
    return programs

#FIXME: this is only while development
def skipProgram(programName):
    if programName in ['xmipp_sqlite3', 'xmipp_mpi_steps_runner',
                       'xmipp_angular_commonline', 'xmipp_python',
                       'xmipp_transform_threshold', 'xmipp_mpi_write_test']:
        return True
    for p in ['xmipp_test', 'xmipp_template']:
        if programName.find(p) != -1:
            return True
    return False

def getProgramsDbName():
    return os.path.join(getXmippPath(), '.xmipp_programs.sqlite')

#Some helper functions
def createProgramsDb(dbName=None):
    from protlib_sql import ProgramDb
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

#---------------------------------------------------------------------------
# Metadata stuff
#--------------------------------------------------------------------------- 
#create a metadata file with original image name, and two other 
#lines with variation over the original name
def intercalate_union_3(inFileName, outFileName, src1, targ1, src2, targ2):
    mD = xmipp.MetaData(inFileName)
    mDout = xmipp.MetaData()   
    for id in mD:       
        idOut = mDout.addObject()
        sIn = mD.getValue(xmipp.MDL_IMAGE, id)
        mDout.setValue(xmipp.MDL_IMAGE, sIn, idOut)
        enabled = mD.containsLabel(xmipp.MDL_ENABLED)
       
        if  (enabled):       
            i = int(mD.getValue(xmipp.MDL_ENABLED, id))
            mDout.setValue(xmipp.MDL_ENABLED, i, idOut)
       
        idOut = mDout.addObject()

        ss = sIn.replace(src1, targ1)
        mDout.setValue(xmipp.MDL_IMAGE, ss, idOut)
        
        if  (enabled):
            mDout.setValue(xmipp.MDL_ENABLED, i, idOut)
            
        idOut = mDout.addObject()
       
        ss = sIn.replace(src2, targ2)
        mDout.setValue(xmipp.MDL_IMAGE, ss, idOut)
        
        if  (enabled):
            mDout.setValue(xmipp.MDL_ENABLED, i, idOut)
       
    mDout.write(outFileName)

#set rot and tilt between -180,180 and -90,90
def check_angle_range(inFileName, outFileName):
    mD = xmipp.MetaData(inFileName)
    doWrite = False
    
    for id in mD: 
        doWrite2 = False
        rot = mD.getValue(xmipp.MDL_ANGLE_ROT, id)
        tilt = mD.getValue(xmipp.MDL_ANGLE_TILT, id)
        if tilt > 90.: 
            tilt = -(int(tilt) - 180)
            rot += 180.
            doWrite = True
            doWrite2 = True
        if tilt < -90.: 
            tilt = -(int(tilt) + 180)
            rot -= 180. 
            doWrite = True
            doWrite2 = True
        if (doWrite2):
            mD.setValue(xmipp.MDL_ANGLE_ROT , rot, id)
            mD.setValue(xmipp.MDL_ANGLE_TILT, tilt, id)
        
    if(doWrite or inFileName != outFileName):
        mD.write(outFileName)


#compute histogram
def compute_histogram(mD, bin, col, min, max):
    allMD = xmipp.MetaData()
    outMD = xmipp.MetaData()   
    _bin = (max - min) / bin
   
    for h in range(0, bin):
        outMD.removeObjects(xmipp.MDQuery("*"))
        if (h == 0):
            outMD.importObjects(mD, xmipp.MDValueRange(col, float(min), float(_bin * (h + 1) + min)))
        if (h > 0 and h < (bin - 1)):
            outMD.importObjects(mD, xmipp.MDValueRange(col, float(_bin * h + min), float(_bin * (h + 1) + min)))
        if (h == (bin - 1)):
            outMD.importObjects(mD, xmipp.MDValueRange(col, float(_bin * h + min), float(max)))
       
        _sum = float(outMD.aggregateSingle(xmipp.AGGR_SUM, xmipp.MDL_WEIGHT))
        outMD.addLabel(xmipp.MDL_COUNT)
        outMD.setValueCol(xmipp.MDL_COUNT, int(_sum + 0.1))
        allMD.unionAll(outMD)
       
    return allMD

def estimateFilenamesListMemory(input):
    n = min(len(input), 100) # 100 looks like a good value
    #for file in input:
    memory = n * estimateFileNameSize(input[0])
    return min(max(memory, 536870912), 2147483648) # minimum 512m, maximum 2gb

def estimateFileNameSize(input):
    from xmipp import FileName
    if '@' in input:
        trueFn = input.split('@')[1]
    else:
        trueFn = input
    fn = FileName(trueFn);
    if fn.isMetaData():
        return estimateMDSize(fn)
    elif fn.isImage():
        return estimateImageSize(fn)
    else:
        return 0

def estimateMDSize(input):
    from xmipp import FileName
    memory = 0
    fn = FileName(input);
    
    if fn.exists():
        MD = xmipp.MetaData(input)
        if(MD.containsLabel(xmipp.MDL_IMAGE)):
            for id in MD:
                fnImg = MD.getValue(xmipp.MDL_IMAGE, id)
                idMemory = estimateImageSize(fnImg)
                memory = max(memory, idMemory)
    
    return memory

def estimateImageSize(input):
    memory = 0
    input = FileName(input)
    
    if input.exists() and input.isImage():
        (Xdim, Ydim, Zdim, Ndim) = xmipp.SingleImgSize(input)
        memory = Xdim * Ydim * Zdim * Ndim * 8
    
    return memory

def convertBytes(bytes):
    from math import ceil
    
    bytes = float(bytes)
    if bytes >= 1099511627776:
        terabytes = bytes / 1099511627776
        size = '%dt' % ceil(terabytes)
    elif bytes >= 1073741824:
        gigabytes = bytes / 1073741824
        size = '%dg' % ceil(gigabytes)
    elif bytes >= 1048576:
        megabytes = bytes / 1048576
        size = '%dm' % ceil(megabytes)
    elif bytes >= 1024:
        kilobytes = bytes / 1024
        size = '%dk' % ceil(kilobytes)
    else:
        size = '%db' % ceil(bytes)
    return size

def estimateMemory(input):
    from math import log, ceil
    MD = xmipp.MetaData(input)
    memory = 0
    for id in MD:
        if MD.containsLabel(xmipp.MDL_IMAGE):
            fnImg = MD.getValue(xmipp.MDL_IMAGE, id)
        elif  MD.containsLabel(xmipp.MDL_MICROGRAPH):
            fnImg = MD.getValue(xmipp.MDL_MICROGRAPH, id)
        (Xdim, Ydim, Zdim, Ndim) = xmipp.SingleImgSize(fnImg)
        memory = max(memory, Xdim * Ydim * 8)
    memoryMb = int((2 ** ceil(log(memory, 2))) / (2 ** 20)); # Memory size in megabytes
    return memoryMb

#---------------------------------------------------------------------------
# Colors from Xmipp binding
#--------------------------------------------------------------------------- 
from xmipp import XMIPP_MAGENTA, XMIPP_BLUE, XMIPP_GREEN, XMIPP_RED, XMIPP_YELLOW, XMIPP_CYAN, colorStr

colorMap = {'red': XMIPP_RED, 'blue': XMIPP_BLUE,
                'green': XMIPP_GREEN, 'magenta': XMIPP_MAGENTA,
                'yellow': XMIPP_YELLOW, 'cyan': XMIPP_CYAN}


blueStr = lambda s: colorStr(XMIPP_BLUE, s)
greenStr = lambda s: colorStr(XMIPP_GREEN, s)
greenLowStr = lambda s: colorStr(XMIPP_GREEN, s, 0)
failStr = redStr = lambda s: colorStr(XMIPP_RED, s)
headerStr = magentaStr = lambda s: colorStr(XMIPP_MAGENTA, s)
yellowStr = lambda s: colorStr(XMIPP_YELLOW, s)
cyanStr = warnStr = cyanStr = lambda s: colorStr(XMIPP_CYAN, s)


def findColor(color):
    '''This function will search if there are color characters present
    on string and return the color and positions on string'''
    for k, v in colorMap.iteritems():
        x, y = colorStr(v, "_..._").split("_..._")
        fx = color.find(x)
        fy = color.find(y)
        if fx != -1 and fy != -1:
            color = color.replace(x, '').replace(y, '')
            return (k, fx, fy, color)
    return None
        
def validateSameSize(fileList, errors, errorPrefix='References'):
    '''Validate if a list of images(or volumes) have
    the same dimensions. 
    The dimensions tuple is returned'''
    from xmipp import SingleImgSize
    (xdim, ydim, zdim, ndim) = SingleImgSize(fileList[0])
    for filename in fileList[1:]:
        (xdim2, ydim2, zdim2, ndim2) = SingleImgSize(filename)
        if (xdim2, ydim2, zdim2, ndim2) != (xdim, ydim, zdim, ndim):
            errors.append("%s: %s and %s have not the same size" % \
                           (errorPrefix, fileList[0], filename)) 
    return (xdim, ydim, zdim, ndim)

def validateInputSize(references, images, md, errors):
    '''This function will validate that all references
    have the same size and also the input images
    The images metadata file should contains the MDL_IMAGE label
    '''
    from xmipp import MetaData, MDL_IMAGE, MetaDataInfo
    # Check reference size
    (xdim, ydim, zdim, ndim) = validateSameSize(references, errors)    
    # Check that volume and images have the same size
    if md.containsLabel(MDL_IMAGE):
        (xdimImg,ydimImg,_,_,_)=MetaDataInfo(md)    
        if (xdimImg, ydimImg) != (xdim, ydim):
            errors.append("References and images have not the same size")
    else:
        errors.append("Input metadata <%s> does not contain image column" % images)
    