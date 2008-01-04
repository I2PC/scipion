%{
#include "../funcs.h"
%}

// Many functions are ignored through the #ifndef SWIG mechanism.
// Open the ../funcs.h file to see which ones
class string {};
%include "../funcs.h"

PRINT(FileName)

/* Test code
python
import XmippData

# A normal function with returned values and default arguments
print XmippData.gaussian1D(0,1);
print XmippData.gaussian1D(3,1);
print XmippData.gaussian1D(3,1,3);

# A normal function with returned valued and reference arguments
x1=doubleP()
x2=doubleP()
result=XmippData.solve_2nd_degree_eq(1.0,-1.0,-6.0,x1,x2)
print x1.value()
print x2.value()

# A normal function
XmippData.randomize_random_generator()
print XmippData.rnd_unif()

# FileNames with string parameters
fn=XmippData.FileName("g0ta00001.xmp")
print fn.get_root()
print fn.add_extension("ext")

# A function with a FileName as argument
print XmippData.exists(XmippData.FileName("docfile.txt"))
# La forma de construir un FileName es muy incómoda

*/
