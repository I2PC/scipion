#!/usr/bin/env python
"""/***************************************************************************
 *
 * Authors:     Roberto Marabini
 *
 * Universidad Autonoma de Madrid
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
"""
import os, glob, sys
scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
sys.path.append(scriptdir) # add default search path
import XmippData

from copy import copy
from optparse import Option, OptionValueError
class Code:
    def __init__(self):
       self.unknown=0
       self.union=1
       self.inter=2
       self.subs=3
       self.product=4
       self.copy=5
       self.move=6
       self.delete=7
       self.select=8
   
    def encode(self, _string):
       if (  _string=='-u' or _string=='--union'):
            return self.union
       elif (_string=='-n' or _string=='--inter' or _string=='--intersection'):
            return self.inter
       elif (_string=="-s" or _string == "--subs" or _string == "--substraction"):
            return self.subs
       elif (_string=="-p" or _string ==  "--product"or _string == "--naturalJoin"):
            return self.product
       elif (_string == '-c' or _string=='--copy'):
            return self.copy
       elif (_string == '-m' or _string=='--move'):
            return self.move
       elif (_string == '-d' or _string=='--delete'):
            return self.delete
       elif (_string == '-s' or _string=='--Select'):
            return self.select
       else:
            return self.unknown
        
    def decode(self, _operation):
       if (_operation==self.union):
           return 'Union'
       elif (_operation==self.inter):
            return 'Intersection'
       elif (_operation==self.subs):
            return 'Substraction'
       elif (_operation==self.product):
            return 'Natural Join (product)'
       elif (_operation==self.copy):
            return 'Copy'
       elif (_operation==self.move):
            return 'Move'
       elif (_operation==self.delete):
            return 'Delete'
       elif (_operation==self.select):
            return 'Select'
       else:
            return 'Unknown Operation'

_myCode=Code()
operation=_myCode.unknown
operationCounter=0

def check_operation(option, opt_str, value, parser):
    global operation
    global operationCounter
    global _myCode
    if (operationCounter > 0):
        raise OptionValueError("you cannot select two operations: "+ operation + ", " +_operation)
    operationCounter += 1    
    operation=_myCode.encode(opt_str)   
    if(operation!=_myCode.unknown):
        print "Operation: ", _myCode.decode(operation) 
    if(value!=None):    
         setattr(parser.values, option.dest, value)
         
def check_label(option, opt_str, value, parser):
    global operation
    if(operation==_myCode.unknown):
        raise OptionValueError("select operation before %s" % opt_str)
    if (operation == _myCode.union or 
        operation == _myCode.inter or
        operation == _myCode.subs or
        operation == _myCode.product or
        operation == _myCode.select):
       setattr(parser.values, option.dest, value)
    else:
       raise OptionValueError("selected option " + opt_str + " is not compatible with operation "+ _myCode.decode(operation))
#        
def check_value(option, opt_str, value, parser):
    global operation
    if(operation==_myCode.unknown):
        raise OptionValueError("select operation before %s" % opt_str)
    if (operation == _myCode.select):
       setattr(parser.values, option.dest, value)
    else:
       raise OptionValueError("selected option " + opt_str + " is not compatible with operation "+ _myCode.decode(operation))

def command_line_options():
      """ add command line options here"""
      import optparse
      _usage = """usage: %prog [options]
Example:
   operate_metadata  -i i1.doc --union i2.doc -o result.doc" 
   operate_metadata  -i i1.doc -o result.doc --copy /home/roberto/kk" 
   operate_metadata  -i i1.doc -o result.doc --naturalJoin i2.doc -l angleROT" 
   operate_metadata  -i i1.doc -o result.doc --select angleROT --minValue 0 --maxValue 23" 
   """
      global operation
      parser = optparse.OptionParser(_usage)
      parser.add_option("-i", "--inMetaDataFile", dest="inMetaDataFile",
			default="", type="string",
			help="MetaData input file name")	    
      parser.add_option("-o", "--outMetaDataFile", dest="outMetaDataFile",
			default="", type="string",
			help="MetaData output file name, by default overwrites inputfile")	    
      
      parser.add_option("-u","--union", dest="inMetaDataFile2",
			default="", type="string",action="callback", callback=check_operation,
			help="add two metadata files (repetitions are allowed)")
      parser.add_option("-n","--inter", "--intersection", dest="inMetaDataFile2",
            default="", type="string",action="callback", callback=check_operation,
            help="intersect two metadata files")
      parser.add_option("-s","--subs", "--substraction", dest="inMetaDataFile2",
            default="", type="string",action="callback", callback=check_operation,
            help="substract two metadata files) (result = entries in first metadata not in the second)")
      parser.add_option("-p", "--product", "--naturalJoin", dest="inMetaDataFile2",
            default="", type="string",action="callback", callback=check_operation,
            help=" natural join using label given in -l option, this label should be unique")
      parser.add_option("-l", "--label", dest="cLabel",
            default="", type="string",action="callback", callback=check_label,
            help="label used for natural join or select")

      parser.add_option("-c","--copy", dest="cPath",
            default="", type="string",action="callback", callback=check_operation,
            help="copy files in metadata file to new directory")
      parser.add_option("-m","--move", dest="cPath",
            default="", type="string",action="callback", callback=check_operation,
            help="move files in metadata file to new directory")
      parser.add_option("-d","--delete",
            action="callback", callback=check_operation,
            help="delete files in metadata file")

      parser.add_option("-S", "--Select", dest="cLabel",
            default="", type="string",action="callback", callback=check_operation,
            help="select entries with CLabel in the range given by --minvalue --maxvalue (only works with DOUBLES)")
      parser.add_option("-v", "--minValue", dest="minValue",
            default=0., type="float",action="callback", callback=check_value,
            help="copy files in metadata file to new directory")
      parser.add_option("-V", "--maxValue", dest="maxValue",
            default=0., type="float",action="callback", callback=check_value,
            help="copy files in metadata file to new directory")

            
      (options, args) = parser.parse_args()
      if(len(options.inMetaDataFile) < 1):
          parser.print_help()
          exit()
      #overwrite file	  
      if(len(options.outMetaDataFile) < 1):
             options.outMetaDataFile = options.inMetaDataFile    
    
      return(options.inMetaDataFile,
             options.outMetaDataFile,
             options.inMetaDataFile2,
             options.cLabel,
             options.cPath,
             options.minValue,
             options.maxValue)
	     
def process_union(inMetaDataFileS,inMetaDataFile2S,outMetaDataFileS):
    inMetaDataFile   = XmippData.FileName(inMetaDataFileS)
    inMetaDataFile2  = XmippData.FileName(inMetaDataFile2S)
    outMetaDataFile  = XmippData.FileName(outMetaDataFileS)
    MD1=XmippData.MetaData(inMetaDataFile)
    MD2=XmippData.MetaData(inMetaDataFile2)
    MD1.union_(MD2)
    MD1.write(outMetaDataFile)
    
def process_subs(inMetaDataFileS,inMetaDataFile2S,outMetaDataFileS,cLabel):
    #l must be given
    global operation
    if(len(cLabel)<1):
        print " Metadata label '-l' is required"
        sys.exit()
    inMetaDataFile   = XmippData.FileName(inMetaDataFileS)
    inMetaDataFile2  = XmippData.FileName(inMetaDataFile2S)
    outMetaDataFile  = XmippData.FileName(outMetaDataFileS)
    
    mdl=XmippData.MetaDataContainer()

    MD1=XmippData.MetaData(inMetaDataFile)
    MD2=XmippData.MetaData(inMetaDataFile2)
    MD3=XmippData.MetaData()
    if(operation==_myCode.subs):
        MD3.substraction(MD1,MD2,mdl.codifyLabel(cLabel))
        MD3.write(outMetaDataFile)
    elif(operation==_myCode.inter):
        MD3.intersection(MD1,MD2,mdl.codifyLabel(cLabel))
        MD3.write(outMetaDataFile)
    elif(operation==_myCode.product):
        MD1.combine(MD2,mdl.codifyLabel(cLabel))
        MD1.write(outMetaDataFile)

def process_copy(inMetaDataFileS,cPath):
    global operation
    import shutil
    if(len(cPath)<1):
        print "Product requires a Metadata label '-l'"
        sys.exit()
    _char= cPath[0:1]
    if(_char=='-'):
        print "Invalid directory name ", cPath 
        sys.exit()
    #create directory
    if not os.path.exists(cPath):
         os.makedirs(cPath)
    #loop
    inMetaDataFile   = XmippData.FileName(inMetaDataFileS)
    MD1=XmippData.MetaData(inMetaDataFile)
    id=MD1.firstObject()
    ss=XmippData.stringP()
    messsage=""
    if (operation==_myCode.copy):
        message="copying "
    elif (operation==_myCode.move):
        message="moving "
    elif (operation==_myCode.delete):
        message="delete "
        
    while(id!=XmippData.MetaData.NO_MORE_OBJECTS):
         XmippData.getValueString(MD1, XmippData.MDL_IMAGE, ss,id)
         if (operation==_myCode.copy):
             shutil.copy(ss.value(),cPath)
         elif (operation==_myCode.move):
             shutil.move(ss.value(),cPath)
         elif (operation==_myCode.delete):
             os.remove(ss.value())     
         print message + ss.value() 
         id=MD1.nextObject()
def process_select(inMetaDataFileS,outMetaDataFileS,minValue,maxValue,cLabel):
     print "NOT IMPLEMENTED (YET)"
#    inMetaDataFile   = XmippData.FileName(inMetaDataFileS)
#    outMetaDataFile  = XmippData.FileName(outMetaDataFileS)
#    MD1=XmippData.MetaData(inMetaDataFile)
#    MD2=XmippData.MetaData()
#    mdl=XmippData.MetaDataContainer()
#    XmippData.addObjectsInRangeDouble(MD1,MD2,mdl.codifyLabel(cLabel),minValue,
#                                    maxValue)
#    MD2.write(outMetaDataFile)
    
##########
##MAIN
##########

#process command line
(inMetaDataFile,
outMetaDataFile,
inMetaDataFile2,
cLabel,
cPath,
minValue,
maxValue
) = command_line_options()


#operate
if(operation==_myCode.union):
    process_union(inMetaDataFile,inMetaDataFile2,outMetaDataFile)
elif(operation==_myCode.inter or 
     operation==_myCode.subs  or
     operation==_myCode.product):
    process_subs(inMetaDataFile,inMetaDataFile2,outMetaDataFile,cLabel)

elif(operation==_myCode.copy or
     operation==_myCode.delete or
     operation==_myCode.move):
    process_copy(inMetaDataFile)
    
elif(operation==_myCode.select):
    process_select(inMetaDataFile,outMetaDataFile,minValue,maxValue,cLabel)
else:
    print "Error: Wrong operation:", operation
    

