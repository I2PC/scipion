%{
#include "../filename.h"
%}

// Many functions are ignored through the #ifndef SWIG mechanism.
// Open the ../funcs.h file to see which ones
class std::string {};
%ignore FILENAMENUMBERLENGTH;
%include "../filename.h"
%ignore to_lowercase;


%pointer_class(FileName, FileNameP);

PRINT(FileName)
// %typemap(in) const FileName& {
//     $1 = &(FileName(PyString_AsString($input)));
// }