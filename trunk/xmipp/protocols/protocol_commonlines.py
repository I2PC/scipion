#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for building initial references by common lines
# This protocol is based on EMAN 1
#
# Example use:
# ./xmipp_protocol_commonlines.py
#
# Author:Carlos Oscar Sorzano, January 2011
#
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob
from protlib_filesystem import changeDir, deleteFiles

class ProtCommonLines(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.commonlines.name, scriptname, project)
        self.Import = 'from protocol_commonlines import *'
    
    def defineSteps(self):
        fnOut=os.path.join(self.WorkingDir,"inputImages.hed")
        self.Db.insertStep('convertImages', [fnOut], Selfile=self.InSelFile, OutStack=fnOut)
        fnOut=os.path.join(self.WorkingDir,"ali.hed")
        self.Db.insertStep('centerImages', [fnOut], WorkingDir=self.WorkingDir, Radius=self.Radius)
        fnOut=os.path.join(self.WorkingDir,"threed.0a.mrc")
        self.Db.insertStep('commonlines', [fnOut], WorkingDir=self.WorkingDir, Radius=self.Radius, NumberOfMpi=self.NumberOfMpi,
                           Symmetry=self.Symmetry)

    def validate(self):
        from protlib_utils import which
        errors = []
        startAny=which('startAny')
        if startAny=='':
            errors.append("EMAN is not accesible")
        return errors
    
    def summary(self):
        message = []
        message.append("Initial volume by common lines constructed from " + self.InSelFile)
        message.append("Symmetry=" + self.Symmetry+" Radius="+str(self.Radius))
        return message

def convertImages(log,Selfile,OutStack):
    runJob(log,"xmipp_image_convert",' -i '+Selfile+' -o '+OutStack);

def centerImages(log,WorkingDir,Radius):
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    runJob(log,"cenalignint",'inputImages.hed mask='+str(Radius))
    deleteFiles(log,['avg.hed','avg.img'],True)
    changeDir(log,currentDir)

def commonlines(log,WorkingDir,Radius,NumberOfMpi,Symmetry):
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    params='ali.hed mask='+str(Radius)+" rounds=5"
    if NumberOfMpi>1:
        params+=" proc="+str(NumberOfMpi)
    if Symmetry!="c1":
        params+=" sym="+Symmetry
    runJob(log,"startAny",params)
    deleteFiles(log,['avg.hed','avg.img'],True)
    changeDir(log,currentDir)
