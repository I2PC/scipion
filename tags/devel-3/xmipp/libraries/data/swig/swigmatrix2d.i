%{
#include "../matrix2d.h"
%}

%ignore arrayByArray;
%ignore ludcmp<>;
%ignore lubksb<>;
%ignore svdcmp<>;
%ignore ludcmp;
%ignore lubksb;
%ignore svdcmp;

%rename(algebraicMultiplication) multiplyMatrix<>;
%rename(applyAffineBspline) applyGeometryBSpline<>;
%rename(applyAffine) applyGeometry<>(
    Matrix2D<T>& M2, Matrix2D< double > A, const Matrix2D<T>& M1, bool inv,
    bool wrap, T outside);
%rename(solveLinear) solve<>;

%import swigmultidimensional_array.i

%include "../matrix2d.h"
PRINT(Matrix2D)

%template(Matrix2Dd) Matrix2D<double>;
%template(Matrix2Di) Matrix2D<int>;

%pythoncode
%{
    Matrix2D=Matrix2Dd
    Matrix2D=Matrix2Di
%}

%template(multiplyMatrix) multiplyMatrix<double>;
%template(multiplyMatrix) multiplyMatrix<int>;
%template(solve) solve<double>;
%template(solveBySVD) solveBySVD<double>;
%template(ludcmp) ludcmp<double>;
%template(lubksb) lubksb<double>;
%template(svdcmp) svdcmp<double>;
%template(applyGeometry) applyGeometry<double>;
%template(applyGeometry) applyGeometry<int>;
%template(applyGeometryBSpline) applyGeometryBSpline<double>;
%template(applyGeometryBSpline) applyGeometryBSpline<int>;
%template(cutToCommonSize) cutToCommonSize<double>;
%template(cutToCommonSize) cutToCommonSize<int>;
%template(radialAverage) radialAverage<double>;
%template(radialAverage) radialAverage<int>;

/*
python
import XmippData

M=XmippData.Matrix2D()
print M
M=XmippData.Matrix2D(6,3)
print M
Maux=XmippData.Matrix2D(M)
print Maux

M=XmippData.Matrix2D(3,3)
M.initIdentity()
print M
M.initIdentity(4)
print M
M.initIdentity(2,3)
print M

M=XmippData.Matrix2D(3,3)
M.initZeros(6,3)
print M

M1=XmippData.Matrix1Di()
M1.initLinear(0,2,1)
print M1
M.fromVector(M1)
print M
M1.initZeros(3)
print M1
M.toVector(M1)
print M1

M.outside(1,1)
M.printShape()

M=XmippData.Matrix2D(3,3)
M.setXmippOrigin()
M.printShape()

Maux.sameShape(M)

M=XmippData.Matrix2Dd(3,3)
M.setXmippOrigin()
M(0,0)
M.setPixel(0,0,1.5)
print M
M.interpolatedElement(0.5,0.5)

M1=XmippData.Matrix1Dd()
M.profile(-1,-1,1,1,7,M1)
M.getRow(0,M1)
print M1
M1=M.Row(0)

M=XmippData.Matrix2D(3,3)
M1=XmippData.Matrix1Di()
M1.initLinear(0,2,1)
print M1
M.initIdentity(3)
print M
M.setCol(1,M1)
print M

Maux=M.transpose()
print Maux
M.det()
M.inv(Maux)
print Maux
Maux = M.inv()
print Maux

M=XmippData.Matrix2Dd(3,3)
M.initIdentity(3)
print M
M.selfRotate(45)
print M
M.initIdentity(3)
M.selfRotateBSpline(3,45)
print M

M1=XmippData.Matrix1Dd(3)
M1.initConstant(1)
Maux=XmippData.scale3DMatrix(M1)
print Maux

A=XmippData.Matrix2Dd(3,3)
b=XmippData.Matrix1Dd(3)
x=XmippData.Matrix1Dd(3)
A.initIdentity(3)
b.initLinear(1,3,1)
print A
print b
XmippData.solveNonNegative(A,b,x)
print x

v=XmippData.Matrix2Dd(3,3)
v.initRandom(0,1);
x=XmippData.Matrix2Dd(3,3)
x.initRandom(0,1);
x=-x

z=v+x
z=v*x
z=v-x
z=v/x

print x
print A*x
XmippData.multiplyElements(A,x,z)


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

A=XmippData.Matrix2Dd(3,3)
A.initIdentity()
center=XmippData.Matrix1D(2)
center.initZeros()
mean=XmippData.Matrix1Dd()
count=XmippData.Matrix1Di()
XmippData.radialAverage(A,center,mean,count);

A=XmippData.Matrix2Dd(3,3)
b=XmippData.Matrix1Dd(3)
x=XmippData.Matrix1Dd(3)
A.initIdentity(3)
b.initLinear(1,3,1)
print A
print b
XmippData.solveLinear(A,b,x)
print x
XmippData.solveBySVD(A,b,x,1e-6)
print x

*/
