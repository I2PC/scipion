%{
#include "../docfile.h"
%}

%ignore DocLine::operator[](int);
%ignore DocFile::show_line;
%rename("get") DocLine::operator[](int) const;
%include "../docfile.h"
PRINT(DocLine)
PRINT(DocFile)

/* Test code:
python
import XmippData
df=XmippData.DocFile()
df.read(XmippData.FileName("g0ta_movements.txt"))
print df
df.go_beginning()
line=0
while not df.eof():
    print df.get_current_line()
    if line>0:
        print df(1)
    df.next()
    line=line+1
df.go_beginning()
df.next()
print df.get_current_line()
print df.get_current_line().get(6)
*/
