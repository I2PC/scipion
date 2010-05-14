%{
#include "../image.h"
%}

/* ignore funcs for overload*/
//%ignore ImageImagicT;

//class Image {};
%include "../image.h"
PRINT(Image)
%template(Imaged) Image<double>;
%template(Imagei) Image<int>;
/*
%extend ImageXmippT {
   %template(set_eulerAngles)  set_eulerAngles<double>;
   %template(set_eulerAngles1) set_eulerAngles1<double>;
   %template(set_eulerAngles2) set_eulerAngles2<double>;
   %template(get_eulerAngles)  get_eulerAngles<double>;
   %template(get_eulerAngles1) get_eulerAngles1<double>;
   %template(get_eulerAngles2) get_eulerAngles2<double>;
};
*/
/*
%pythoncode
%{
    ImageT=ImageTd
    ImageT=ImageTi
    ImageXmippT=ImageXmippTd
    ImageXmippT=ImageXmippTi
%}
*/
/*
python
import os,glob,sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/lib'
sys.path.append(scriptdir) # add default search path
import XmippData
mD=XmippData.MetaData();
selFile=XmippData.FileName('sum.doc')
ss=XmippData.stringP()
mD.read(selFile)
XmippData.getValueString( mD,XmippData.MDL_IMAGE, ss)
inFile=XmippData.FileName(ss)
_image=XmippData.Imaged()
_image.read(inFile)
dd=XmippData.doubleP()

#mdc=XmippData.MetaDataContainer()
#mdc = mD.getObject()




vv=XmippData.vectorm()
vv.append(XmippData.MDL_IMAGE)
vv.append(XmippData.MDL_ANGLEROT)
vv.append(XmippData.MDL_ANGLETILT)

##################################vv.append(XmippData.MDL_ANGLEPSI)

_image.read(inFile,True,-1,True,False,mD,vv)
print ss.value()
XmippData.getValueDouble( _image.MD,XmippData.MDL_ANGLETILT,dd)
outFile=XmippData.FileName('kk.xmp')
_image.write(outFile)
*/
