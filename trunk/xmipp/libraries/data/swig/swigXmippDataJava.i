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
%ignore FileName::FileName(char const *);
%{
#include "../funcs.h"
%}
class std::string {};

%include "../funcs.h"

// Matrix3D
%ignore MultidimArray::setXdim;
%ignore MultidimArray::setYdim;
%ignore MultidimArray::setZdim;

%ignore coreArrayByArray<>( const MultidimArray<T>& op1,
                            const MultidimArray<T>& op2,
                            MultidimArray<T>& result,
                            char operation);
%ignore coreArrayByScalar<>(const MultidimArray<T>& op1,
                            const T& op2,
                            MultidimArray<T>& result,
                            char operation);
%ignore coreScalarByArray<>(const T& op1,
                            const MultidimArray<T>& op2,
                            MultidimArray<T>& result,
                            char operation);
%ignore coreArrayByArray;
%ignore coreArrayByScalar;
%ignore coreScalarByArray;

%ignore applyGeometry<>;
%ignore applyGeometryBSpline<>;
%ignore getShifts;
%{
#include "../multidim_array.h"
%}
%include "../multidim_array.h"
%template (MultidimArrayd) MultidimArray<double>;

%{
#include "../image.h"
%}
%include "../image.h"
%template (ImageDouble) Image<double>;

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


