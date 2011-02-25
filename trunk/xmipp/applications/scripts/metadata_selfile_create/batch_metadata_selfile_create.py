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
import os, glob, sys, optparse
   
def command_line_options():
      """ add command line options here"""
      _usage="""usage: %prog [options]
Example:
   xmipp_metadata_selfile_create -p "Images/*xmp" -o all_images.sel -isstack"""
      parser = optparse.OptionParser(_usage)
      parser.add_option("-p", "--pattern",  dest="pattern", type="string",
            help="The pattern to match")        
      parser.add_option("-o", "--outMetaDataFile", dest="outMetaDataFile",
            default="", type="string",
            help="MetaData output file name")
      parser.add_option("-s", "--isstack",
                  action="store_true", dest="isStack", default=False,
                  help="Check if the images are stacks")
      parser.add_option("-q", "--quiet",
                  action="store_true", dest="quiet", default=False,
                  help="Do not show any messages")
            
      (options, args)      = parser.parse_args()
      if(len(options.outMetaDataFile)<1):
          parser.print_help()
          exit()
      
      if not options.quiet:
          print 'Files pattern:        ', options.pattern
          print 'Output Metadata File: ', options.outMetaDataFile
          print 'Is stack:             ', options.isStack
    
      return(options.pattern, 
             options.outMetaDataFile,
             options.isStack,
             options.quiet)



scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/lib'
sys.path.append(scriptdir) # add default search path

from xmipp import *

pattern, outFile, isStack, quiet = command_line_options();

files = glob.glob(pattern)
files.sort()

mD = MetaData()
inFile = FileName()

nSize = 1
for file in files:
    fileAux=file
    if isStack:
        if file.endswith(".mrc"):
            fileAux=file+":mrcs"
        x, x, x, nSize = SingleImgSize(fileAux)
    if nSize != 1:
        counter = 1
        for jj in range(nSize):
            inFile.compose(counter, fileAux)
            objId = mD.addObject()
            mD.setValue(MDL_IMAGE, inFile, objId)
            mD.setValue(MDL_ENABLED, 1, objId)
            counter = counter + 1
    else:
        objId=mD.addObject()
        mD.setValue(MDL_IMAGE, fileAux, objId)
        mD.setValue(MDL_ENABLED, 1, objId)

mD.write(outFile)
