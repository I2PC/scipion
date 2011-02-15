#make path abdolute in metadata file .value()
import os,sys
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
sys.path.append(scriptdir) # add default search path
#import XmippData
from xmipp import *
 

#change filename  by absolutepath/filename
def absolutePath(inFileName,outFileName):
   mD = MetaData(inFileName)
   mD.makeAbsPath(MDL_IMAGE)
   mD.write(outFileName)

#replace string in metada file
def replaceString(inFileName,outFileName,instring,outString):
    print '[replaceString] Need to be implemented'
#   mD = XmippData.MetaData(XmippData.FileName(inFileName))
#   sIn =XmippData.stringP() 
#   id=mD.firstObject()
#   while(id != -1 ):
#    XmippData.getValueString(mD, XmippData.MDL_IMAGE, sIn)
#    ss=sIn.value()
#        sIn.assign(ss.replace(instring,outString))
#    XmippData.setValueString(mD, XmippData.MDL_IMAGE, sIn)
#    id=mD.nextObject()
#   mD.write(XmippData.FileName(outFileName))

#replace string in medatada object
def replaceString(mD,instring,outString):
    print '[replaceString] Need to be implemented'
#   sIn =XmippData.stringP() 
#   id=mD.firstObject()
#   while(id != -1 ):
#    XmippData.getValueString(mD, XmippData.MDL_IMAGE, sIn)
#    ss=sIn.value()
#        sIn.assign(ss.replace(instring,outString))
#    XmippData.setValueString(mD, XmippData.MDL_IMAGE, sIn)
#    id=mD.nextObject()

#set enabled to -1 if exits
def deactivate_all_images(mD):
    print '[deactivate_all_images] Need to be implemented'
#   i =XmippData.intP() 
#   i.assign(-1)
#   id=mD.firstObject()
#   while(id != -1 ):
#    XmippData.setValueInt(mD, XmippData.MDL_ENABLED, i)
#    id=mD.nextObject()

#create a metadata file with original image name, and two other 
#lines with variation over the original name
def intercalate_union_3(inFileName,outFileName, src1,targ1,src2,targ2):
   
   mD = MetaData(inFileName)
   mDout = MetaData()
   
   for id in mD:       
       mDout.addObject()
       sIn = mD.getValue(MDL_IMAGE,id)
       mDout.setValue(MDL_IMAGE, sIn, id)
       enabled= mD.containsLabel(MDL_ENABLED)
       if  (enabled):
            i = int(mD.getValue(MDL_ENABLED,id))
            mDout.setValue(MDL_ENABLED, i, id)
       
       mDout.addObject()

       ss = sIn.replace(src1,targ1)
       mDout.setValue(MDL_IMAGE, ss, id)
       if  (enabled):
           mDout.setValue(MDL_ENABLED, i, id)
           
       mDout.addObject()
       
       ss = sIn.replace(src2,targ2)
       mDout.setValue(MDL_IMAGE, ss, id)
       if  (enabled):
           mDout.setValue(MDL_ENABLED, i, id)
       
   mDout.write(outFileName)

#set rot and tilt between -180,180 and -90,90
def check_angle_range(inFileName,outFileName):

    mD    = MetaData(inFileName)
    doWrite=False
    
    for id in mD: 
        doWrite2=False
        rot = mD.getValue(MDL_ANGLEROT,id)
        tilt = mD.getValue(MDL_ANGLETILT,id)
        if tilt > 90.: 
            tilt = -(int(tilt)-180)
            rot  += 180.
            doWrite=True
            doWrite2=True
        if tilt < -90.: 
            tilt = -(int(tilt)+180)
            rot  -= 180. 
            doWrite=True
            doWrite2=True
        if (doWrite2):
            mD.setValue(MDL_ANGLEROT , rot, id)
            mD.setValue(MDL_ANGLETILT, tilt, id)
        
    if(doWrite or inFileName != outFileName):
        mD.write(outFileName)


#compute histogram
def compute_histogram(mD,bin,col,min,max):
    
   allMD = MetaData()
   outMD = MetaData()   
   _bin = (max-min)/bin
   
   for h in range(0,bin):
       outMD.removeObjects(MDQuery("*"))
       if (h==0):
           outMD.importObjects(mD, MDValueRange(col, float(min), float(_bin*(h + 1)+min)))
       if (h>0 and h<(bin-1)):
           outMD.importObjects(mD, MDValueRange(col, float(_bin * h + min), float(_bin*(h + 1)+min)))
       if (h==(bin-1)):
           outMD.importObjects(mD, MDValueRange(col, float(_bin * h + min), float(max)))
       
       _sum=float(outMD.aggregateSingle(AGGR_SUM,MDL_WEIGHT))
       outMD.addLabel(MDL_COUNT)
       outMD.setValueCol(MDL_COUNT, int(_sum+0.1))
       allMD.unionAll(outMD)
       
   return allMD
