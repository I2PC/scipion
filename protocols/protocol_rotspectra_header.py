#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for rotational spectra classification
# using self-organizing maps
#
# Example use:
# ./xmipp_rotational_spectra.py
#
# Author:Roberto Marabini, March 2007
#        Carlos Oscar Sorzano, January 2011
#
# required files: log.py, spider_header.py, Sel_Files.py
# by defaults searchs for python files in ../python directory
#
# {begin_of_header}

# {eval} expandCommentRun()

#------------------------------------------------------------------------------------------------
# {section} Rotational spectra parameters
#------------------------------------------------------------------------------------------------
# {file}(images*.xmd){validate}(PathExists) Selfile with the input images:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# {list_combo}(Use the middle of the image,Minimize first harmonic)How to find the center of rotation?
HowCenter='Use the middle of the image'

# Inner radius for rotational harmonics (%):
""" A percentage of the image radius """
SpectraInnerRadius = 15

# Outer radius for rotational harmonics (%):
""" A percentage of the image radius """
SpectraOuterRadius = 80

# {expert} Lowest harmonic to calculate
SpectraLowHarmonic = 1
# {expert} Highest harmonic to calculate
SpectraHighHarmonic = 15

#-----------------------------------------------------------------------------
# {section} Classification: classify_kerdensom 
#-----------------------------------------------------------------------------
# X-dimension of the self-organizing map:
SomXdim = 7
# Y-dimension of the self-organizing map:
SomYdim = 7
# {expert} Initial regularization factor:
""" 
The kerdenSOM algorithm anneals from an initial high regularization factor
to a final lower one, in a user-defined number of steps.
If the output map is too smooth, lower the regularization factors
If the output map is not organized, higher the regularization factors
See [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/KerDenSOM]
"""
SomReg0 = 1000

# {expert} Final regularization factor:
SomReg1 = 200

# {expert} Steps to lower the regularization factor:
SomSteps = 5
# {expert} Additional kerdenSOM parameters:
""" 
For a complete description 
See [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/KerDenSOM]
"""
KerdensomExtraCommand = ''

# {eval} expandParallel(threads=0, mpi=0, hours=96)

#-----------------------------------------------------------------------------
# {end_of_header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
from protocol_rotspectra import *
#
# main
#     
if __name__ == '__main__':
    protocolMain(ProtRotSpectra)
