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
%ignore *::operator==;
%rename (add)   operator+;
%rename (sub)   operator-;
%rename (div)   operator/;
%rename (mul)   operator*;

// rm libraries/data/swig/swigXmippData_wrap_java.cpp libraries/data/swig/*.java ; ./scons.compile
%ignore FileName::FileName(char const *);
%{
#include "../filename.h"
%}
class std::string {};
%include "../filename.h"

%{
#include "../funcs.h"
%}
%include "../funcs.h"

// Matrix3D
%ignore MultidimArray::setXdim;
%ignore MultidimArray::setYdim;
%ignore MultidimArray::setZdim;
%ignore MultidimArray::setNdim;

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
%ignore getImage;
%{
#include "../multidim_array.h"
%}
%include "../multidim_array.h"
%template (MultidimArrayd) MultidimArray<double>;
//%template (vectorMDObject) std::vector<MDObject *>

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

%{
#include "../error.h"
%}
%include "../error.h"

// handle non-const string& parameters
%typemap(jni) std::string *INOUT, std::string &INOUT %{jobjectArray%}
%typemap(jtype) std::string *INOUT, std::string &INOUT "java.lang.String[]"
%typemap(jstype) std::string *INOUT, std::string &INOUT "java.lang.String[]"
%typemap(javain) std::string *INOUT, std::string &INOUT "$javainput"

%typemap(in) std::string *INOUT (std::string strTemp ), std::string &INOUT (std::string strTemp ) {
  if (!$input) {
    SWIG_JavaThrowException(jenv, SWIG_JavaNullPointerException, "array null");
    return $null;
  }
  if (JCALL1(GetArrayLength, jenv, $input) == 0) {
    SWIG_JavaThrowException(jenv, SWIG_JavaIndexOutOfBoundsException, "Array must contain at least 1 element");
    return $null;
  }

  jobject oInput = JCALL2(GetObjectArrayElement, jenv, $input, 0); 
  if ( NULL != oInput ) {
    jstring sInput = static_cast<jstring>( oInput );

    const char * $1_pstr = (const char *)jenv->GetStringUTFChars(sInput, 0); 
    if (!$1_pstr) return $null;
    strTemp.assign( $1_pstr );
    jenv->ReleaseStringUTFChars( sInput, $1_pstr);  
  }

  $1 = &strTemp;
}

%typemap(freearg) std::string *INOUT, std::string &INOUT ""

%typemap(argout) std::string *INOUT, std::string &INOUT
{ 
  jstring jStrTemp = jenv->NewStringUTF( strTemp$argnum.c_str() );
  JCALL3(SetObjectArrayElement, jenv, $input, 0, jStrTemp ); 
}

%apply std::string &INOUT { std::string & strOut };
%{
#include "../metadata.h"
%}
%include "../metadata.h"

// %template(getValueStr) MetaData::getValue<std::string>;

%ignore MDObject::toString;
%{
#include "../metadata_label.h"
%}
%include "../metadata_label.h"

%{
#include "../java_wrapper.h"
%}
%include "../java_wrapper.h"

%{
#include "../fft.h"
%}
%include "../fft.h"

%{
#include "../xmipp_java_interface.h"
%}
%include "../xmipp_java_interface.h"
