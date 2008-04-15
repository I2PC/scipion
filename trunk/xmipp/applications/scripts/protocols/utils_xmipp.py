#!/usr/bin/env python
FILENAMENUMBERLENTGH=6
#---------------------------------------------------------------------------
# composeFileName (Xmipp-like)
#---------------------------------------------------------------------------
def composeFileName(rootname,number,extension):
   if (number != -1):
      output = rootname + str(number).zfill(FILENAMENUMBERLENTGH)
   else:
      output = rootname

   if (extension != ''):
      output += '.' + extension 

   return output

