#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for performing operations over metadatas 
#
# Example use:
# ./xmipp_metadata_utilities.py
#
# Author: Javier Vargas  Sept 2013 

from xmipp import MetaData
from protlib_base import *
from protlib_utils import getListFromRangeString, runJob, runShowJ
from protlib_filesystem import copyFile, deleteFile, removeFilenamePrefix

class ProtMetadataUtils(XmippProtocol):
    
    def __init__(self, scriptname, project):

        XmippProtocol.__init__(self, protDict.metadata_utilities.name, scriptname, project)
        self.Import = 'from protocol_metadata_utilities import *'
        self.fnOut = self.extraPath("images.xmd")        
        
    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        
        if (self.OperationType == "set"):
            self.insertStep("setOption",OperationSetType=self.OperationSetType,md=self.MetadataInput,md2=self.md2,Label=self.Label,Label2=self.Label2,fnOut=self.fnOut,ModeType=self.ModeType)
        elif (self.OperationType == "operate"):
            self.insertStep("operateOption",Operation_TypeOp=self.Operation_TypeOp,md=self.MetadataInput,md2=self.md2,fnOut=self.fnOut, LabelOp=self.LabelOp, LabelSize=self.LabelSize,Label_KC=self.Label_KC,Label_DC=self.Label_DC,Label_RC=self.Label_RC,Expression=self.Expression,ModeType=self.ModeType)
        elif (self.OperationType == "query"):
            self.insertStep("queryOption",Query_Operation=self.Query_Operation,md=self.MetadataInput,md2=self.md2,fnOut=self.fnOut, SelectExpression=self.SelectExpression, LabelCount=self.LabelCount,LabelSum1=self.LabelSum1,LabelSum2=self.LabelSum2,ModeType=self.ModeType)
        elif (self.OperationType == "fill"):
            self.insertStep("fillOption",Fill_Operation=self.Fill_Operation,md=self.MetadataInput,fnOut=self.fnOut,LabelFill=self.LabelFill,ValueCont=self.ValueCont,InitValueLinealCont1=self.InitValueLinealCont1,InitValueLinealCont2=self.InitValueLinealCont2,InitValueRandUnifCont1=self.InitValueRandUnifCont1,InitValueRandUnifCont2=self.InitValueRandUnifCont2,InitValueRandGaussCont1=self.InitValueRandGaussCont1,InitValueRandGaussCont2=self.InitValueRandGaussCont2,InitValueRandStCont1=self.InitValueRandStCont1,InitValueRandStCont2=self.InitValueRandStCont2,InitValueRandStCont3=self.InitValueRandStCont3,ModeType=self.ModeType)

    def summary(self):        
        message=[]
        message.append("Input metadata: [%s]"%self.MetadataInput)
        message.append("Output metadata: [%s]"%self.fnOut)        
        return message
    
    def validate(self):
        errors = []            
        return errors

    def visualize(self):        
        runShowJ(self.fnOut)
        
def setOption(log,OperationSetType, md, md2, Label, Label2, fnOut,ModeType):
    if OperationSetType=="inner_join":
        runJob(log,'xmipp_metadata_utilities',"-i %s --set %s %s %s %s -o %s --mode %s"
               %(md,OperationSetType,md2,Label,Label2,fnOut,ModeType))
    else:
        runJob(log,'xmipp_metadata_utilities',"-i %s --set %s %s %s %s -o %s --mode %s"
               %(md,OperationSetType,md2,Label,Label2,fnOut,ModeType))

def operateOption(log,Operation_TypeOp,md,md2,fnOut,LabelOp,LabelSize,Label_KC,Label_DC,Label_RC,Expression,ModeType):
    
    if (Operation_TypeOp=="sort"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --operate %s %s -o %s --mode %s"
               %(md,Operation_TypeOp,LabelOp,fnOut,ModeType))
        
    elif (Operation_TypeOp=="random_subset"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --operate %s %s -o %s --mode %s"
               %(md,Operation_TypeOp,LabelSize,fnOut,ModeType))        

    elif (Operation_TypeOp=="bootstrap"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --operate %s -o %s --mode %s"
               %(md,Operation_TypeOp,fnOut,ModeType))        
    
    elif (Operation_TypeOp=="randomize"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --operate %s -o %s --mode %s"
               %(md,Operation_TypeOp,fnOut,ModeType))        

    elif (Operation_TypeOp=="keep_column"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --operate %s %s -o %s --mode %s"
               %(md,Operation_TypeOp,Label_KC,fnOut,ModeType))    

    elif (Operation_TypeOp=="drop_column"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --operate %s %s -o %s --mode %s"
               %(md,Operation_TypeOp,Label_DC,fnOut,ModeType))    

    elif (Operation_TypeOp=="rename_column"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --operate %s %s -o %s --mode %s"
               %(md,Operation_TypeOp,Label_RC,fnOut,ModeType))    

    elif (Operation_TypeOp=="modify_values"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --operate %s %s -o %s --mode %s"
               %(md,Operation_TypeOp,Expression,fnOut,ModeType))    
    

def queryOption(log,Query_Operation,md,md2,fnOut, SelectExpression, LabelCount,LabelSum1,LabelSum2,ModeType):

    if (Query_Operation=="select"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --query %s %s -o %s --mode %s"
               %(md,Query_Operation,SelectExpression,fnOut,ModeType))
        
    elif (Query_Operation=="count"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --query %s %s -o %s --mode %s"
               %(md,Query_Operation,LabelCount,fnOut,ModeType))
    
    elif (Query_Operation=="sum"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --query %s %s %s -o %s --mode %s"
               %(md,Query_Operation,LabelSum1,LabelSum2,fnOut,ModeType))  
  

def fillOption(log,Fill_Operation,md,fnOut,LabelFill,ValueCont,InitValueLinealCont1,InitValueLinealCont2,InitValueRandUnifCont1,InitValueRandUnifCont2,InitValueRandGaussCont1,InitValueRandGaussCont2,InitValueRandStCont1,InitValueRandStCont2,InitValueRandStCont3,ModeType):

    if (Fill_Operation=="constant"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --fill  %s %s %s -o %s --mode %s"
               %(md,LabelFill,Fill_Operation,ValueCont,fnOut,ModeType))

    elif (Fill_Operation=="lineal"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --fill  %s %s %s %s -o %s --mode %s"
               %(md,LabelFill,Fill_Operation,InitValueLinealCont1,InitValueLinealCont2,fnOut,ModeType))

    elif (Fill_Operation=="rand_uniform"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --fill  %s %s %s %s -o %s --mode %s"
               %(md,LabelFill,Fill_Operation,InitValueRandUnifCont1,InitValueRandUnifCont2,fnOut,ModeType))

    elif (Fill_Operation=="rand_gaussian"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --fill  %s %s %s %s -o %s --mode %s"
               %(md,LabelFill,Fill_Operation,InitValueRandGaussCont1,InitValueRandGaussCont2,fnOut,ModeType))

    elif (Fill_Operation=="rand_student"):
        runJob(log,'xmipp_metadata_utilities',"-i %s --fill  %s %s %s %s %s -o %s --mode %s"
               %(md,LabelFill,Fill_Operation,InitValueRandStCont1,InitValueRandStCont2,InitValueRandStCont3,fnOut,ModeType))        

    