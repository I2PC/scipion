// Name of the output module
%module XmippData

// C and C++ string definitions
%include stl.i
%include std_string.i
%include "std_vector.i"

// Pointers to base classes
%include cpointer.i
%pointer_class(int, intP);
%pointer_class(char,charP);
%pointer_class(double,doubleP);
%pointer_class(float,floatP);
%pointer_class(std::string,stringP);

// General includes
%{
#include <sstream>
#include <complex>
#include <string>
using namespace std;
%}
%ignore operator<<;
%ignore operator>>;
%ignore *::operator();
%ignore *::operator=;
%ignore *::operator+=;
%ignore *::operator-=;
%ignore *::operator*=;
%ignore *::operator/=;
%rename (add)   operator+;
%rename (sub)   operator-;
%rename (div)   operator/;
%rename (mul)   operator*;

// rm libraries/data/swig/swigXmippData_wrap_java.cpp libraries/data/swig/*.java ; ./scons.compile

%{
#include "../funcs.h"
%}
%include "../funcs.h"

// Matrix3D
%ignore coreArrayByArray;
%ignore coreArrayByScalar;
%ignore coreScalarByArray;
%ignore applyGeometry<>;
%ignore applyGeometryBSpline<>;
%{
#include "../multidim_array.h"
%}
%template (MultidimArrayd) MultidimArray<double>;

%{
#include "../image.h"
%}
%include "../image.h"
%template (Image) Image<double>;

%ignore project_SimpleGrid;
%ignore project_GridVolume;
%ignore Projection::read;
%{
#include "../projection.h"
%}
%include "../projection.h"

%{
#include "../geometry.h"
%}
%include "../geometry.h"


