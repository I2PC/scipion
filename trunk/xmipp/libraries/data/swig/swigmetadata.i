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
*/
