%{
#include "../funcs.h"
%}

// Many functions are ignored through the #ifndef SWIG mechanism.
// Open the ../funcs.h file to see which ones
class string {};
%include "../funcs.h"

PRINT(FileName)
