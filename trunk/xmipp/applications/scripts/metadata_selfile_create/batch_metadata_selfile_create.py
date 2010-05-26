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
outFile = XmippData.FileName()
inFile = XmippData.FileName()

if len(sys.argv) == 1:
   print 'Usage: selfile_create "pattern"  metadataFile'
   sys.exit()
elif len(sys.argv) == 2:
   outFile = '/dev/stdout'
elif len(sys.argv) == 3:
   outFile = sys.argv[2]
else:
   print 'Usage   : xmipp_selfile_create "pattern"  metadataFile'
   print 'Example1: xmipp_selfile_create "Images/*xmp"    all_images.sel'
   print 'Example2: xmipp_selfile_create "Images/*xmp"  > all_images.sel'
   sys.exit()


mD = XmippData.MetaData()
sIn = XmippData.stringP()
sOut = XmippData.stringP()
ii = XmippData.intP()
files = glob.glob(sys.argv[1])
files.sort()
ii.assign(1)
x = XmippData.intP()
n = XmippData.intP()
nSize = 100

for file in files:
    sIn.assign(file)
    counter = 0
    inFile.compose(-1,sIn)
    XmippData.ImgSize(inFile, x, x, x, n)
    nSise = n.value()
    if nSize != 1:
        for jj in range (n.value()):
            mD.addObject()
            inFile.compose(counter, sIn)
            XmippData.setValueString(mD, XmippData.MDL_IMAGE, inFile, -1)
            XmippData.setValueInt(mD, XmippData.MDL_ENABLED, ii, -1)
            counter = counter + 1
    else:
        mD.addObject()
        XmippData.setValueString(mD, XmippData.MDL_IMAGE, sIn, -1)
        XmippData.setValueInt(mD, XmippData.MDL_ENABLED, ii, -1)
mD.write(outFile)
