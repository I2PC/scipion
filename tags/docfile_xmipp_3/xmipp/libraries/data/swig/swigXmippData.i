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
%ignore operator>>;
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

// Arithmetic operators
%rename (add)   operator+;
%rename (sub)   operator-;
%rename (div)   operator/;
%rename (mul)   operator*;

// All interfaces being ported
%include swigfuncs.i
%include swigdocfile.i
%include swigargs.i
%include swigsqllite.i
%include swigmetadata.i
//%include swigmultidimensional_array.i
//%include swigmatrix1d.i
//%include swigmatrix2d.i
//%include swigmatrix3d.i
//%include swigimage.i
//%include swiggeometry.i
//%include swigmicrograph.i

// rm libraries/data/swig/_XmippData.so libraries/data/swig/XmippData.py libraries/data/swig/XmippData_wrap.os libraries/data/swig/XmippData_wrap.cc lib/_XmippData.so lib/XmippData.py ; ./scons.compile
