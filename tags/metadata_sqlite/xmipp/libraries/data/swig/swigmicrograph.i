%{
#include "../micrograph.h"
%}
%include "../micrograph.h"

/*
python
import XmippData

M=XmippData.Micrograph()
M.open_micrograph(XmippData.FileName("7773.raw"))
print M.depth()
# COSS
# listParticles=M.Particles()
# listParticles.size()

coord=XmippData.Particle_coords()
coord.label=1
coord.X=1000
coord.Y=1000
coord.valid=True

# COSS
# M.scissor(coord, result, 0, 255)

M(1000,1000)
M.set_val(1000,1000,0)
M(1000,1000)

Dmin=XmippData.doubleP()
Dmax=XmippData.doubleP()
M.computeDoubleMinMax(Dmin,Dmax)
print Dmin.value(),Dmax.value()

M.close_micrograph()

# COSS
# Normalize
*/
