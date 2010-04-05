%{
#include "../metadata.h"
#include "../metadata_container.h"
%}
//Python does not have float/double variable, just reals that will be cast to double
%ignore MetaData::setValue( MetaDataLabel name, float value, long int objectID = -1 );
%ignore MetaDataContainer::addValue( MetaDataLabel name, float value );
%ignore MetaDataContainer::pairExists( MetaDataLabel name, float value );

%include ../metadata.h
%include ../metadata_container.h

/*
==================
First example
==================
python
import XmippData
#create metadata object
outFile=XmippData.FileName()
outFile='kk'
mD=XmippData.MetaData();
#create metadatacontainer (one line of metadata)
id=mD.addObject();
valor=10.3
#fill an object data
result=mD.setImage("image0001.xmp")
result=mD.setAngleRot(valor)
result=mD.setAnglePsi(123.2)

id=mD.addObject();
result=mD.setImage("image0002.xmp")
result=mD.setAngleRot(valor*2.)
result=mD.setAnglePsi(123.3)

mD.write(outFile)

#modify first added object
result=mD.setAngleRot(valor*3)

mD.write("prueba2.doc")

mD2=XmippData.MetaData();

mD2.read(outFile)

==================
ANOTHER EXAMPLE (xmipp_selfile_create)
==================

#!/usr/bin/env python
import os,glob,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/lib'
sys.path.append(scriptdir) # add default search path

import XmippData # import interface with xmipp

#some functions in xmipp require filenames, this is not a native 
#data type in python. So it need to be explicitely defined as follows
outFile=XmippData.FileName()
if len(sys.argv) == 1:
   print 'Usage: selfile_create "pattern"  metadataFile'
   sys.exit()
elif len(sys.argv) == 2:
   outFile='/dev/stdout'
elif len(sys.argv) == 3:
   outFile=sys.argv[2]
else:
   print 'Usage   : xmipp_selfile_create "pattern"  metadataFile'
   print 'Example1: xmipp_selfile_create "Images/*xmp"    all_images.sel'
   print 'Example2: xmipp_selfile_create "Images/*xmp"  > all_images.sel'
   sys.exit()

#create metadata object
mD=XmippData.MetaData()

#loop through the images
for file in glob.glob(sys.argv[1]):
    mD.addObject();
    mD.setImage(file)
    mD.setEnabled(1)
#write sel file    
mD.write(outFile)

*/
