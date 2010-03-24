%{
#include "../metadata.h"
#include "../metadata_container.h"
%}
//Python does not have float/double variable, just reals that will be cast to double
%ignore MetaData::setValue( MetaDataLabel name, float value, long int objectID = -1 );

%include ../metadata.h
%include ../metadata_container.h


/*

python
import XmippData
#create metadata object
mD=XmippData.MetaData();
#create metadatacontainer (one line of metadata)
id=mD.addObject();
valor=10.3
#fill an object data
result=mD.setValue(XmippData.MDL_ANGLEROT,valor)
result=mD.setValue(XmippData.MDL_ANGLETILT,-30.1)
result=mD.setValue(XmippData.MDL_ANGLEPSI,170.1)
result=mD.setValue(XmippData.MDL_IMAGE,"/home/jrbcast/test.xmp")

id=mD.addObject();
result=mD.setValue(XmippData.MDL_IMAGE,"/home/jrbcast/test2.xmp")
result=mD.setValue(XmippData.MDL_ANGLEPSI,1)
result=mD.setValue(XmippData.MDL_ANGLEROT,valor-1)
result=mD.setValue(XmippData.MDL_ANGLETILT,-35.1)
result=mD.setValue(XmippData.MDL_ANGLETILT,0.0)

mD.save("prueba1.doc")

#modify first added object
result=mD.setValue(XmippData.MDL_ANGLETILT,0.0,0)

mD.save("prueba2.doc")
*/
