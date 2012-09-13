#!/usr/bin/env python
FILENAMENUMBERLENTGH=6
#---------------------------------------------------------------------------
# utils_xmipp.composeFileName (Xmipp-like)
#---------------------------------------------------------------------------
def composeFileName(rootname,number,extension):
   if (number != -1):
      output = rootname + str(number).zfill(FILENAMENUMBERLENTGH)
   else:
      output = rootname

   if (extension != ''):
      output += '.' + extension 

   return output

def composeWildcardFileName(rootname,extension):
   output = rootname
   for i in range(FILENAMENUMBERLENTGH):
      output += '?'

   if (extension != ''):
      output += '.' + extension 

   return output

def getCommaSeparatedIntegerList(inputstring):
   import string
   lista=string.split(inputstring,",")
   output=[]
   for i in range(len(lista)):
      interval=string.split(lista[i],'-')
      if len(interval)==1:
         if not interval[0]=='':
            output+=[int(interval[0])]
      else:
         output+=range(int(interval[0]),
                       int(interval[1])+1)
   return output
