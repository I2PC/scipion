%{
#include "../args.h"
%}
extern void checkAngle(const std::string& str);

/*
python
import XmippData

XmippData.checkAngle('rot')
XmippData.checkAngle('roT')
*/
