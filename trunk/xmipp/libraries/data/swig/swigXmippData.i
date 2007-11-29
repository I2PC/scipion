%module XmippData
%include cstring.i
%include std_string.i
%include cpointer.i
%pointer_class(int,intP);
%pointer_class(char,charP);
%pointer_class(double,doubleP);
%pointer_class(float,floatP);

%{
#include <string.h>
%}
