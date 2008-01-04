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
x1=0.0
x2=0.0
result=XmippData.solve_2nd_degree_eq(1.0,-1.0,-6.0,x1,x2,1e-6)
# *** NO FUNCIONA 

print XmippData.gaussian1D(0,1);
print XmippData.gaussian1D(3,1);
print XmippData.gaussian1D(3,1,3);

XmippData.randomize_random_generator()
print XmippData.rnd_unif()

fn=XmippData.FileName("g0ta00001.xmp")
print fn.get_root()
print fn.add_extension("ext")

print XmippData.exists(XmippData.FileName("docfile.txt"))
# La forma de construir un FileName es muy incómoda

*/
