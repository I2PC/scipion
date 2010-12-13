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
       sIn = mD.getValue(MDL_IMAGE)
       mDout.setValue(MDL_IMAGE, sIn)
       enabled= mD.containsLabel(MDL_ENABLED)
       if  (enabled):
            i = int(mD.getValue(MDL_ENABLED))
            mDout.setValue(MDL_ENABLED, i)
       
       mDout.addObject()

       ss = sIn.replace(src1,targ1)
       mDout.setValue(MDL_IMAGE, ss)
       if  (enabled):
           mDout.setValue(MDL_ENABLED, i)
           
       mDout.addObject()
       
       ss = sIn.replace(src2,targ2)
       mDout.setValue(MDL_IMAGE, ss)
       if  (enabled):
           mDout.setValue(MDL_ENABLED, i)
       
   mDout.write(outFileName)

#set rot and tilt between -180,180 and -90,90
def check_angle_range(inFileName,outFileName):
    print '[check_angle_range] Need to be implemented' 
#
#   mD    = XmippData.MetaData(XmippData.FileName(inFileName))
#   rotP  = XmippData.doubleP() 
#   tiltP = XmippData.doubleP() 
#   id=mD.firstObject()
#   doWrite=False
#   while(id != -1 ):
#    doWrite2=False
#        XmippData.getValueDouble(mD, XmippData.MDL_ANGLEROT,   rotP)
#    rot  = rotP.value()
#    XmippData.getValueDouble(mD, XmippData.MDL_ANGLETILT, tiltP)
#    tilt = tiltP.value()
#    if tilt > 90.: 
#           tilt = -(tilt-180)
#           rot  += 180.
#       doWrite=True
#       doWrite2=True
#    if tilt < -90.: 
#           tilt = -(tilt+180)
#           rot  -= 180. 
#       doWrite=True
#       doWrite2=True
#    if (doWrite2):   
#       rotP.assign(rot)
#       tiltP.assign(tilt)
#       XmippData.setValueDouble(mD, XmippData.MDL_ANGLEROT , rotP)
#       XmippData.setValueDouble(mD, XmippData.MDL_ANGLETILT, tiltP)
#     id=mD.nextObject()
#   if(doWrite or inFileName != outFileName):
#     mD.write(XmippData.FileName(outFileName))


#compute histogram
def write_several(mD,bin,col,min,max):
    print '[write_several] Need to be implemented' 
#   #add column to metadata MD
#   outMD = XmippData.MetaData()
#   allMD = XmippData.MetaData()
#   a = XmippData.MDL()
#   iP=XmippData.intP()
#   _col =  a.label2Str(col)
#   _bin=(max-min)/bin
#   for h in range(0,bin):           
#       #delete next
#       outMD.clear()
#       expr = " "
#       if(h!=0):
#          expr = _col +" >= " + str(_bin * h + min ) 
#       if(h!=0 and h!=(bin-1) ):
#          expr += ' AND '
#       if(h!=(bin-1)):
#          expr +=  _col + " < " + str(_bin*(h + 1)+min) 
#       outMD.importObjects(mD, XmippData.MDExpression(expr));
#       id=outMD.firstObject()
#       _sum=outMD.aggregateSingle(XmippData.AGGR_SUM,XmippData.MDL_WEIGHT)
#       iP.assign(int(_sum+0.1))
#       XmippData.setValueColInt(outMD, XmippData.MDL_COUNT,iP)
#       allMD.unionAll(outMD)
#       #outMD.write(XmippData.FileName("/dev/stderr"))
#   allMD.setFilename(mD.getFilename())
#   return allMD

