// Name of the output module
%module XmippData

// C and C++ string definitions
%include stl.i
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

// General includes
%{
#include <sstream>
#include <complex>
#include <string>
using namespace std;
%}

// Map reference types and pointers
%typemap(out) int &
{
$result=PyInt_FromLong((long)*$1);
}

%typemap(out) double &
{
$result=PyInt_FromLong((long)*$1);
}

%typemap(out) bool &
{
$result=PyInt_FromLong((long)*$1);
}

%typemap(out) float &
{
$result=PyInt_FromLong((long)*$1);
}

%typemap(out) double *
{
$result=PyInt_FromLong((long)*$1);
}

%typemap(out) bool *
{
$result=PyInt_FromLong((long)*$1);
}

%typemap(out) int *
{
$result=PyInt_FromLong((long)*$1);
}

%typemap(out) float *
{
$result=PyInt_FromLong((long)*$1);
}

// The C++ insertion operator cannot be ported to Python. Instead, use the
// SWIG macro PRINT(type) to add printing capabilities to a class.
%ignore operator<<;
%define PRINT(type)
%extend type
{
    std::string __str__()
    {
        std::stringstream s;
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
%include swigfuncs.i
%include swigdenoise.i
%include swigdocfile.i
%include swigargs.i

// rm libraries/data/swig/libXmippDataSwig.so libraries/data/swig/XmippData.py libraries/data/swig/swigXmippData_wrap.os libraries/data/swig/swigXmippData_wrap.cc ; ./scons.compile
