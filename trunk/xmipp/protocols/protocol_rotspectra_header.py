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
#------------------------------------------------------------------------------------------
# {section}{has_question} Comment
#------------------------------------------------------------------------------------------
# Display comment
DisplayComment = False

# {text} Write a comment:
""" 
Describe your run here...
"""
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# Run name:
""" This will identify your protocol run. It need to be unique for each protocol. You could have run1, run2 for protocol X, but not two
run1 for it. This name together with the protocol output folder will determine the working dir for this run.
"""
RunName = "run_001"

# Delete working directory?
""" If TRUE the working directory will be deleted before run.
Set this option to TRUE if you want to start from scratch the same run
with previous parameters
"""
DoDeleteWorkingDir = False

""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName='all_images.sel'

#-----------------------------------------------------------------------------
# {section} Rotational spectra calculation
#-----------------------------------------------------------------------------
# {list}(Use the middle of the image,Minimize first harmonic, How to find the center of rotation?)
HowCenter='Use the middle of the image'
# Inner radius for rotational harmonics calculation:
""" These values are in pixels from the image center
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Makespectra
"""
SpectraInnerRadius=7
# Outer radius for rotational harmonics calculation:
SpectraOuterRadius=32
# {expert} Lowest harmonic to calculate
SpectraLowHarmonic=1
# {expert} Highest harmonic to calculate
SpectraHighHarmonic=15

#-----------------------------------------------------------------------------
# {section} Classification: classify_kerdensom 
#-----------------------------------------------------------------------------
# X-dimension of the self-organizing map:
SomXdim=7
# Y-dimension of the self-organizing map:
SomYdim=7
# {expert} Initial regularization factor:
""" The kerdenSOM algorithm anneals from an initial high regularization factor
    to a final lower one, in a user-defined number of steps.
    If the output map is too smooth, lower the regularization factors
    If the output map is not organized, higher the regularization factors
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
SomReg0=1000
# {expert} Final regularization factor:
SomReg1=200
# {expert} Number of steps to lower the regularization factor:
SomSteps=5
# {expert} Additional kerdenSOM parameters:
""" For a complete description 
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraCommand=''
#------------------------------------------------------------------------------------------------
# {hidden} Analysis of results
AnalysisScript='visualize_rotspectra.py'
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# {end_of_header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
from protocol_rotspectra import *
#
# main
#     
if __name__ == '__main__':
   protocolMain(ProtRotSpectra)