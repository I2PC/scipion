%{
#include "../matrix1d.h"
%}

// Open the matrix1d.h to see which are the ignored functions
%ignore adaptForNumericalRecipes;
%ignore killAdaptationForNumericalRecipes;
%ignore powellOptimizer;
%ignore coreArrayByArray;
%ignore coreArrayByScalar;
%ignore coreScalarByArray;
%ignore arrayByArray;
%ignore toPhysical;
%ignore toLogical;

%import swigmultidimensional_array.i
%include "../matrix1d.h"
PRINT(Matrix1D)

%template(Matrix1Dd) Matrix1D<double>;
%template(Matrix1Di) Matrix1D<int>; 

%pythoncode
%{
    Matrix1D=Matrix1Dd
    Matrix1D=Matrix1Di
%}

%template (vectorProduct)   vectorProduct<double>;
%template (vectorProduct)   vectorProduct<int>;
%template (dotProduct)      dotProduct<double>;
%template (dotProduct)      dotProduct<int>;
%template (cuttocommonsize) cutToCommonSize<double>;
%template (sorttwoVectors)  sortTwoVectors<double>;

/*
python
import XmippData

v=XmippData.Matrix1D(3,True)
print v
v=XmippData.Matrix1D(3,False)
print v
v=XmippData.Matrix1D(3)
print v
v2=XmippData.Matrix1D(v)
print v2

v.clear()
v.initLinear(0,10,6,"steps")
print v
v.printShape()

v.resize(3)
print v

v.setXmippOrigin()
v.xinit
v.printShape()
v.isCol()
v.isRow()

v(-1)
v(1)
v.setElement(-1,2)
print v
print v.X()

v2=v.transpose();
print v2
v.selfTranspose()
print v

v=XmippData.Matrix1D(3,True)
v.clear()
v.initLinear(0,10,6,"steps")
print v
i=v.indexSort();

v.window(-3,1)
print v

x=XmippData.intP()
v.maxIndex(x)
x.value()

vsplines=XmippData.Matrix1Dd()
v.produceSplineCoefficients(vsplines)
v.interpolatedElementBSpline(0.5)

v2d=XmippData.vectorR2(1,2)
print v2d

v.initLinear(0,4,2)
XmippData.dotProduct(v,v)

x=XmippData.vectorR3(1,0,0)
y=XmippData.vectorR3(0,1,0)
print XmippData.vectorProduct(x,y)

v=XmippData.vectorR3(0,0,1.0)
x=XmippData.vectorR3(1,1,1.0)
z=v+x
z=v*x
z=v-x
z=v/x

z=v+5
# COSS
# z=5.0+v

v.sameShape(v)
v.finishingX()
v.printStats()
v.computeMax()

minval=XmippData.doubleP()
maxval=XmippData.doubleP()
v.computeDoubleMinMax(minval,maxval)
minval.value()
maxval.value()

avg=XmippData.doubleP()
stddev=XmippData.doubleP()
minval=XmippData.doubleP()
maxval=XmippData.doubleP()
v.computeStats(avg,stddev,minval,maxval)
avg.value()
stddev.value()
minval.value()
maxval.value()

v.addNoise(0,1)
print v

v.threshold("abs_above",0,0)
print v
v.substitute(0,1)
print v

XmippData.MINnD(v,v,v)
print v

*/
