%{
#include "../filename.h"
%}

// Many functions are ignored through the #ifndef SWIG mechanism.
// Open the ../funcs.h file to see which ones
class std::string {};
%ignore FILENAMENUMBERLENGTH;
%include "../filename.h"
%ignore to_lowercase;
//python has its own create directory routines
%ignore mkpath(const FileName &, mode_t);
%ignore do_mkdir( char const *, mode_t);

%pointer_class(FileName, FileNameP);

PRINT(FileName)
// %typemap(in) const FileName& {
//     $1 = &(FileName(PyString_AsString($input)));
// }