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
namespace std {
   %template(vectorm) vector<MetaDataLabel>;
};

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
ss=XmippData.stringP()
dd=XmippData.doubleP()
mD2.getValue( XmippData.MDL_IMAGE, ss)
ss.value()
mD2.getValue( XmippData.MDL_ANGLEPSI, dd)
dd.value()



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

=================
import XmippData
outFile=XmippData.FileName('kk')
mD=XmippData.MetaData();
id=mD.addObject();
mD.setValue( XmippData.MDL_ANGLEPSI, 10.3)
b=XmippData.vectord()
b.push_back(12);
b.push_back(45);
mD.setValue(XmippData.MDL_NMA,b)
mD.setValue(XmippData.MDL_ANGLEROT, 11.0)
mD.write(outFile)

import XmippData
inFile=XmippData.FileName('kk')
mD=XmippData.MetaData();
mD.read(inFile)
b=XmippData.vectord();
mD.getValue(XmippData.MDL_NMA,b)
b[0]
b[1]
dd=XmippData.doubleP()
mD.getValue(XmippData.MDL_ANGLEROT,dd) 
print dd.value()

================

import os,glob,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/lib'
sys.path.append(scriptdir) # add default search path
import XmippData
outFile=XmippData.FileName("yy.doc")
dbFile=XmippData.FileName("kkk.db")
mD=XmippData.MetaData();
mD.read(outFile)
mD.toDataBase( dbFile, "prueba" );
vv=XmippData.vectorm()
vv.append(XmippData.MDL_IMAGE)
vv.append(XmippData.MDL_ANGLEPSI)

print "here"
mD.toDataBase( dbFile, "prueba",vv );
mD.toDataBase( dbFile, "prueba" );

====================

 XmippData.isBool(XmippData.MDL_ENABLED)
 
 ====================
import os,glob,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/lib'
sys.path.append(scriptdir) # add default search path
import XmippData
inputMetaDataFile  = XmippData.FileName('kk.ctfparam')
mD = XmippData.MetaData(inputMetaDataFile)
voltage=XmippData.doubleP()
defocusu=XmippData.doubleP()
mD.getValue(XmippData.MDL_DEFOCUSU,defocusu)
mD.getValue(XmippData.MDL_VOLTAGE,voltage)
print defocusu.value(),voltage.value()
mD.setValue(XmippData.MDL_DEFOCUSU,3333.)
mD.write('ww')
==================================
#combinewithfiles
import os,glob,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/lib'
sys.path.append(scriptdir) # add default search path
import XmippData
inputMetaDataFile   = XmippData.FileName('t.doc')
mD = XmippData.MetaData(inputMetaDataFile)


mD.combineWithFiles(XmippData.MDL_CTFMODEL)
mD.write('ww')

*/



