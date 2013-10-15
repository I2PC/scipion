#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Author: Carlos Oscar Sanchez Sorzano, September 2013
#

# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Volume parameters
#------------------------------------------------------------------------------------------------
# {file}(*.vol){validate}(PathExists) Reference volume:
""" The input volume will be aligned to the reference volume"""
ReferenceVolume=''

# {file}(*.vol){validate}(PathExists) Input volume:
""" The input volume will be aligned to the reference volume"""
InputVolume=''

#-----------------------------------------------------------------------------
# {section}{has_question} Mask
#-----------------------------------------------------------------------------
# Apply mask?
DoMask = False

# {condition}(DoMask){list_combo}(circular, binary_file) Mask type
MaskType = "circular"

# {condition}(DoMask and MaskType=="circular"){wizard}(wizardSetMaskRadiusAlign) Mask radius
MaskRadius = -1

# {file}{condition}(DoMask and MaskType=="binary_file") Mask file
MaskFile = ""

#------------------------------------------------------------------------------------------------
# {section} Search method and space
#------------------------------------------------------------------------------------------------
# {list_combo}(Exhaustive, Local, Exhaustive+Local, Fast Fourier)Alignment algorithm
""" Exhaustive searches all possible combinations within a search space.
    Local searches around a given position.
    Be aware that the Fast Fourier algorithm requires a special compilation
    of Xmipp (DO_FRM=True in install.sh). It performs the same job as the exhaustive method but much faster.
"""
AlignmentMethod='Exhaustive'

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Minimum rotational angle
Rot0=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Maximum rotational angle
RotF=360

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Step rotational angle
RotStep=5

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Minimum tilt angle
Tilt0=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Maximum tilt angle
TiltF=180

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Step tilt angle
TiltStep=5

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Minimum in-plane angle
Psi0=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Maximum in-plane angle
PsiF=360

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Step in-plane angle
PsiStep=5

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Minimum shiftX
X0=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Maximum shiftX
XF=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Step shiftX
XStep=1

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Minimum shiftY
Y0=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Maximum shiftY
YF=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Step shiftY
YStep=1

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Minimum shiftZ
Z0=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Maximum shiftZ
ZF=0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local' or AlignmentMethod=='Fast Fourier') Step shiftZ
ZStep=1

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local') Minimum scale
Scale0=1.0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local') Maximum scale
ScaleF=1.0

# {expert}{condition}(AlignmentMethod=='Exhaustive' or AlignmentMethod=='Exhaustive+Local') Step scale
ScaleStep=0.005

# {expert}{condition}(AlignmentMethod=='Local') Initial rotational angle
RotCurrent=0.

# {expert}{condition}(AlignmentMethod=='Local') Initial tilt angle
TiltCurrent=0.

# {expert}{condition}(AlignmentMethod=='Local') Initial in-plane angle
PsiCurrent=0.

# {expert}{condition}(AlignmentMethod=='Local') Initial shiftX
ShiftXCurrent=0.

# {expert}{condition}(AlignmentMethod=='Local') Initial shiftY
ShiftYCurrent=0.

# {expert}{condition}(AlignmentMethod=='Local') Initial shiftZ
ShiftZCurrent=0.

# {expert}{condition}(AlignmentMethod=='Local') Initial scale
ScaleCurrent=1.

# {eval} expandParallel(mpi=0, threads=0, hours=10)

#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#		
# Main
#     
from protocol_align_volume import *

if __name__ == '__main__':
    protocolMain(ProtAlignVolume)
