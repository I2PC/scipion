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
import os,glob,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/lib'
sys.path.append(scriptdir) # add default search path
import XmippData
from math import *

"""valid math functions ['__builtins__', '__doc__', '__file__', '__name__', 'acos', 'asin ', 'atan', 'atan2', 'ceil',
'cos', 'cosh', 'degrees', 'e' , 'exp', 'fabs', 'floor', 'fmod', 'frexp', 'hidden_value', 'hypot',
'ldexp', 'log', 'log10', 'modf', 'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt', 'tan', 'tanh', 'user_func', 'x']  
"""

def command_line_options():
      """ add command line options here"""
      import optparse
      _usage="""usage: %prog [options]
Example:
   operate_metadata  -i a.doc -o b.doc -e  "angleRot=sin(angleRot*pi/180.)"
   operate_metadata  -i a.doc -o b.doc -e  "image=image.replace('xmp','spi')" """
      parser = optparse.OptionParser(_usage)
      parser.add_option("-i", "--inMetaDataFile",  dest="inMetaDataFile",
			default="", type="string",
			help="MetaData input file name")	    
      parser.add_option("-o", "--outMetaDataFile", dest="outMetaDataFile",
			default="", type="string",
			help="MetaData output file name, by default overwrites inputfile")	    
      parser.add_option("-e", "--expression", dest="expression",
			default="", type="string",
			help="any valid operation in python (no spaces allowed)")
			
      (options, args)      = parser.parse_args()
      if(len(options.inMetaDataFile)<1):
	  parser.print_help()
	  exit()
      #overwrite file	  
      if(len(options.outMetaDataFile)<1):
             options.outMetaDataFile = options.inMetaDataFile
      
      print  '**'
      print 'Input Metadata File  ', options.inMetaDataFile
      print 'Output Metadata File ', options.outMetaDataFile
      print 'Expression           ', options.expression
      print  '**'
    
      return(options.inMetaDataFile,
             options.outMetaDataFile,
             options.expression)
	     
#fill array with labels in expression
def fill_array(expression):
    _array=[]
    for i in range(XmippData.MDL_FIRST_LABEL,XmippData.MDL_LAST_LABEL):
        _label=XmippData.MetaDataContainer.decodeLabel(i)
        if(expression.find(_label) != -1):
	    _array.append(_label)
    return _array
##########
##MAIN
##########

#process command line
(inMetaDataFile,
 outMetaDataFile,
 expression)=command_line_options()

#create MetaData
inputMetaDataFile  = XmippData.FileName(inMetaDataFile)
outputMetaDataFile = XmippData.FileName(outMetaDataFile)
mD = XmippData.MetaData(inputMetaDataFile)

#parse expression
splitting=expression.partition('=')
metaDataLabelName=splitting[0]
expression=splitting[2]
labelsInExpression=fill_array(expression)

#loop through object
id=mD.firstObject()
_auxString=XmippData.stringP()
result=""
while(id!=XmippData.MetaData.NO_MORE_OBJECTS):
   #all operation but replace
   if (expression.find('replace')==-1):
       auxExpression=expression
       for i in labelsInExpression:
	   mD.writeValueToString(_auxString,i)
	   auxExpression=auxExpression.replace(i,_auxString.value())
       result=str(eval(auxExpression))    
   #replace
   else:
       auxExpression=metaDataLabelName
       for i in labelsInExpression:
	   mD.writeValueToString(_auxString,i)
	   auxExpression=auxExpression.replace(i,_auxString.value())
       replaceString="result=auxExpression.replace"+expression.split('replace')[1]
       exec replaceString   
   mD.setValue(metaDataLabelName,result) 
   
   id=mD.nextObject()
   
#write output   
mD.write(outMetaDataFile)
