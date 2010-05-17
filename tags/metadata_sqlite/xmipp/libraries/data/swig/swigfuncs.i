%{
#include "../funcs.h"
%}

// Many functions are ignored through the #ifndef SWIG mechanism.
// Open the ../funcs.h file to see which ones
class std::string {};
%ignore FILENAMENUMBERLENGTH;
%ignore print(std::ostream& o, const bool b);
%ignore printb;
%ignore RealImag2Complex;
%ignore AmplPhase2Complex;
%ignore Complex2RealImag;
%ignore Complex2AmplPhase;
%ignore FileNameComparison;
%ignore exit_if_not_exists;
%ignore time_config;
%ignore annotate_time;
%ignore acum_time;
%ignore elapsed_time;
%ignore print_elapsed_time;
%ignore time_to_go;
%ignore xmippBaseListener;
%ignore xmippTextualListener;
%ignore FREAD;
%ignore FWRITE;
%ignore ByteSwap;
%ignore Marsaglia;
%include "../funcs.h"

PRINT(FileName)
//%typemap(in) const FileName& {
//    $1 = &(FileName(PyString_AsString($input)));
//}

/* Test code
python
import XmippData

# A normal function with returned values and default arguments
print XmippData.gaussian1D(0,1);
print XmippData.gaussian1D(3,1);
print XmippData.gaussian1D(3,1,3);

# A normal function with returned valued and reference arguments
x1=XmippData.doubleP()
x2=XmippData.doubleP()
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
fn=XmippData.FileName()
print fn
fn2=XmippData.FileName(fn)
print fn2
fn=XmippData.FileName("g0ta",2,"xmp")
print fn
fn=XmippData.FileName("new")
print fn
s="hello"
saux="txt"
fn.compose(s,0,saux)
print fn
fn.init_random(15)
print fn

# A function with a FileName as argument
print XmippData.exists(XmippData.FileName("docfile.txt"))
# La forma de construir un FileName es muy incómoda

# Tabsinc
t=XmippData.Tabsinc(0.01,10)
print t(0)

# Kaiser Bessel
t=XmippData.KaiserBessel(10,2,2,0.01,200)
t.sinhwin(32)

*/
