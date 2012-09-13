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

