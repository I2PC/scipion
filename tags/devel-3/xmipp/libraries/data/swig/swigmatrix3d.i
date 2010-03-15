%{
#include "../matrix3d.h"
%}

%ignore arrayByArray;

%rename(applyAffineBspline) applyGeometryBSpline<>;
%rename(applyAffine) applyGeometry<>;

%import swigmultidimensional_array.i
%include "../matrix3d.h"
PRINT(Matrix3D)

%template (applyGeometry) applyGeometry<double>;
%template (applyGeometry) applyGeometry<int>;
%template (applyGeometryBSpline) applyGeometryBSpline<double>;
%template (applyGeometryBSpline) applyGeometryBSpline<int>;
%template (Matrix3Dd) Matrix3D<double>;
%template (Matrix3Di) Matrix3D<int>;

%pythoncode
%{
    Matrix3D=Matrix3Dd
    Matrix3D=Matrix3Di
%}

%template (cutToCommonSize) cutToCommonSize<double>;
%template (cutToCommonSize) cutToCommonSize<int>;
%template (radialAverage) radialAverage<double>;
%template (radialAverage) radialAverage<int>;

/*
python
import XmippData

M=XmippData.Matrix3D(3,3,3)
print M
Maux=XmippData.Matrix3D(M)
print Maux
M.xinit
M.yinit
M.zinit
M.setXmippOrigin()
M.finishingX()

x=XmippData.intP() 
x.assign(4)
x.value()	   
y=XmippData.intP() 
z=XmippData.intP() 
M.xdim
M.ydim
M.zdim
M.getDimension(y,x,z)
print x.value(),y.value(),z.value()

M.selfReverseX()

M=XmippData.Matrix3Dd(3,3,3)
M.setVoxel(1,1,1,1)
M1=XmippData.vectorR3(0,1.0,0)
M.selfTranslate(M1)
print M

Maux=XmippData.Matrix3Dd()
M.pyramidExpand(Maux)
print Maux

M.initZeros(3,3,3)
M.setVoxel(1,1,1,1)
M.setVoxel(0,0,0,1)
center=XmippData.Matrix1Dd()
M.centerOfMass(center)
print center
mask=XmippData.Matrix3Di(3,3,3)
mask.setVoxel(1,0,0,1)
mask.setVoxel(1,0,1,1)
mask.setVoxel(1,0,2,1)
mask.setVoxel(1,1,0,1)
mask.setVoxel(1,1,1,1)
mask.setVoxel(1,1,2,1)
mask.setVoxel(1,2,0,1)
mask.setVoxel(1,2,1,1)
mask.setVoxel(1,2,2,1)
M.centerOfMass(center,mask)
print center
*/
