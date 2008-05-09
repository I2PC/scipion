%{
#include "../docfile.h"
%}

%ignore DocLine::operator[](int);
%ignore DocFile::FirstKey();
%ignore DocFile::operator=;
%ignore DocFile::show_line;
%rename("get") DocLine::operator[](int) const;
%include "../docfile.h"
PRINT(DocLine)
PRINT(DocFile)

/* Test code:
python
import XmippData
df=XmippData.DocFile()
df.read(XmippData.FileName("docfile.txt"))
print df
df.go_beginning()
line=0
while not df.eof():
    print df.get_current_line()
    if line>0:
        print df(1)
    df.next()
    line=line+1
*/
