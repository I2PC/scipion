%{
#include "../cppsqlite3.h"
%}

%ignore CppSQLite3Exception;
%ignore CppSQLite3Buffer;
%ignore CppSQLite3Binary;
%ignore CppSQLite3Query;
%ignore CppSQLite3Statement;
%include "../cppsqlite3.h"

/* Test code:
python
import XmippData
db=XmippData.CppSQLite3DB()
retval=db.open("mibase.db")
db.execDML("create table emp(empno int, empname char(20));")
db.execDML("insert into emp values (1,'Pepe')")
db.execDML("insert into emp values (2,'José')")
print db.execScalar("select count(*) from emp;")
t = db.getTable("select * from emp order by 1;");
for fld in range(t.numFields()):
    print t.fieldName(fld);
for row in range(t.numRows()):
    t.setRow(row);
    for fld in range(t.numFields()):
        if not t.fieldIsNull(fld):
            print t.fieldValue(fld)+"|",
        else:
            print "NULL|",
        print
db.close()
*/
