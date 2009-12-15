%{
#include "../image.h"
%}

/* ignore funcs for overload*/
%ignore ImageImagicT;
%ignore ImageImagicinfo;
%ignore IMAGIC_TAG;
%ignore ImageOver;
%ignore ImageT::read(FILE*& fh, float fIform, int Ydim, int Xdim, bool reversed,
              Image_Type image_type);
%ignore ImageT::LoadImage;
%ignore ImageT::write(FILE*& fh, bool reversed, Image_Type image_type);
%ignore ImageXmippT::readImageHeaderAndContent;
%ignore ImageXmippT::readImageContent;
%ignore ImageXmippT::write(FILE *fp, bool force_reversed);

class Image {};
%include "../image.h"
PRINT(ImageT)
%template(ImageTd) ImageT<double>;
%template(ImageTi) ImageT<int>;

%extend ImageXmippT {
   %template(set_eulerAngles)  set_eulerAngles<double>;
   %template(set_eulerAngles1) set_eulerAngles1<double>;
   %template(set_eulerAngles2) set_eulerAngles2<double>;
   %template(get_eulerAngles)  get_eulerAngles<double>;
   %template(get_eulerAngles1) get_eulerAngles1<double>;
   %template(get_eulerAngles2) get_eulerAngles2<double>;
};

PRINT(ImageXmippT)
%template(ImageXmippTd) ImageXmippT<double>;
%template(ImageXmippTi) ImageXmippT<int>;

%pythoncode
%{
    ImageT=ImageTd
    ImageT=ImageTi
    ImageXmippT=ImageXmippTd
    ImageXmippT=ImageXmippTi
%}

/*
python
import XmippData


*/
