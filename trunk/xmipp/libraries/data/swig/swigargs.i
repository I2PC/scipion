%{
#include "../args.h"
%}

%ignore  textToFloat(const std::string str, int _errno, std::string errmsg,
            int exit);
%ignore textToInteger(const std::string str, int _errno, std::string errmsg,
            int exit);
%ignore firstToken(const std::string& str);		

%include "../args.h"

