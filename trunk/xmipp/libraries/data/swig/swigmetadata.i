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
outFile=XmippData.FileName('kk')
mD=XmippData.MetaData();
#create metadatacontainer (one line of metadata)
id=mD.addObject();
valor=10.3
#fill an object data
mD.setValue( XmippData.MDL_IMAGE, "image0001.xmp")
mD.setValue( XmippData.MDL_ANGLEPSI, valor)
mD.setValue( XmippData.MDL_ANGLEROT, 123.2)

id=mD.addObject();
mD.setValue( XmippData.MDL_IMAGE, "image0002.xmp")
mD.setValue( XmippData.MDL_ANGLEPSI, valor*2.)
mD.setValue( XmippData.MDL_ANGLEROT, 124.2)

mD.write(outFile)

#modify first added object
mD.setValue( XmippData.MDL_ANGLEROT, 444.2)

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
import XmippData
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


mD=XmippData.MetaData()

files = glob.glob(sys.argv[1])
files.sort()
for file in files:
    mD.addObject();
    mD.setValue( XmippData.MDL_IMAGE, file)
    mD.setValue( XmippData.MDL_ENABLED, 1)
mD.write(outFile)

================================
#Do not use de asigment operator rather use the copy constructor
mD2=XmippData.MetaData(mD);
mD2.write("two")
mD2.setValue( XmippData.MDL_ANGLEROT, 111.2)
mD2.write("three")
