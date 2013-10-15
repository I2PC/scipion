#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sanchez Sorzano, September 2013
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Volume parameters
#------------------------------------------------------------------------------------------------
# {file}(*.vol){validate}(PathExists) Input volume:
""" The input volume will be aligned to the reference volume"""
InputVolume=''

# Cylinder radius
CylinderRadius = -1

# Look for dihedrical symmetry
Dihedrical=False

#------------------------------------------------------------------------------------------------
# {section} Search space
#------------------------------------------------------------------------------------------------
# {expert} Minimum rotational angle
""" Helical symmetry is defined as V(r,rot,z)=V(r,rot+k*DeltaRot,z+k*Deltaz). This parameter defines the minimum value of DeltaRot. """
Rot0=0

# {expert} Maximum rotational angle
""" Helical symmetry is defined as V(r,rot,z)=V(r,rot+k*DeltaRot,z+k*Deltaz). This parameter defines the maximum value of DeltaRot. """
RotF=357

# {expert} Step rotational angle
""" Helical symmetry is defined as V(r,rot,z)=V(r,rot+k*DeltaRot,z+k*Deltaz). This parameter defines DeltaRot step. """
RotStep=3

# {expert} Minimum shiftZ
""" Helical symmetry is defined as V(r,rot,z)=V(r,rot+k*DeltaRot,z+k*Deltaz). This parameter defines the minimum value of DeltaZ. """
Z0=1

# {expert} Maximum shiftZ
""" Helical symmetry is defined as V(r,rot,z)=V(r,rot+k*DeltaRot,z+k*Deltaz). This parameter defines the maximum value of DeltaZ. """
ZF=10

# {expert} Step shiftZ
""" Helical symmetry is defined as V(r,rot,z)=V(r,rot+k*DeltaRot,z+k*Deltaz). This parameter defines DeltaZ step. """
ZStep=0.5

# {eval} expandParallel(mpi=0, threads=8, hours=4)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_helical_params import *

if __name__ == '__main__':
    protocolMain(ProtHelicalParams)
