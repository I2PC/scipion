// Name of the output module
%module XmippData

// General includes
%{
#include <sstream>
#include <complex>
#include "../error.h"
//#include <string.h>
%}

// C and C++ string definitions
%include cstring.i
%include std_string.i

// Pointers to base classes
%include cpointer.i
%pointer_class(int,intP);
%pointer_class(char,charP);
%pointer_class(double,doubleP);
%pointer_class(float,floatP);

// Redefine all assignment operators of classes as the function .assign()
%rename(assign) *::operator=;

// The C++ insertion operator cannot be ported to Python. Instead, use the
// SWIG macro PRINT(type) to add printing capabilities to a class.
%ignore operator<<;
%define PRINT(type)
%extend type
{
    std::string __str__()
    {
        stringstream s;
        s << *self;
        return s.str();
    }
}
%enddef

// Exception handling
%exception {
    try {
        $action
    }
    catch (Xmipp_error XE) {
        PyErr_SetString(PyExc_RuntimeError,XE.msg.c_str());
        return NULL;
    }
}

// All interfaces being ported
%include swigdenoise.i
%include swigdocfile.i
%include swigfuncs.i

// rm libraries/data/swig/libXmippDataSwig.so libraries/data/swig/XmippData.py libraries/data/swig/swigXmippData_wrap.os libraries/data/swig/swigXmippData_wrap.cc ; ./scons.compile

