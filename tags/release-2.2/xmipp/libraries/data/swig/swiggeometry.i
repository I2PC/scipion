%{
#include "../geometry.h"
%}
%include "../geometry.h"

/*
python
import XmippData

point=XmippData.vectorR3(2.0,3.0,1.0)
direction=XmippData.vectorR3(1.0,1.0,1.0)
result=XmippData.Matrix1Dd()
XmippData.Uproject_to_plane(point, direction, 0.0, result)

A=XmippData.Matrix2Dd()
XmippData.Euler_angles2matrix(90,0,90,A)
print A


newrot=XmippData.doubleP()
newtilt=XmippData.doubleP()
newpsi=XmippData.doubleP()
XmippData.Euler_mirrorX(90,90,90,newrot,newtilt,newpsi)
print newrot.value(),newtilt.value(),newpsi.value()
*/
