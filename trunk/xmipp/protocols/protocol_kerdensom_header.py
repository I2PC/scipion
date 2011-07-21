#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for image classification using self-organizing maps
#
# Example use:
# ./xmipp_protocol_kerdensom.py
#
# Author:Carlos Oscar Sorzano, January 2011
#
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

# {file} Selfile or stack with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName='CL2D/classes4/class_core_classes/class_core_00000.stk'

#------------------------------------------------------------------------------------------------
# {section} Mask parameters
#------------------------------------------------------------------------------------------------
#{list}(graphical,file) Choose mask type:
"""
Choose the type of mask you want to use. You can design it graphically
or use an existing mask from a file.
The first option will launch a graphical program to design your own mask.
Be careful NOT to submit your job via a queueing system!
"""
MaskType='graphical'
# {file}{condition}(MaskType=file) Mask file 
MaskFileName=''

#-----------------------------------------------------------------------------
# {section} Classification: classify_kerdensom 
#-----------------------------------------------------------------------------
# X-dimension of the map:
SomXdim=7
# Y-dimension of the map:
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
# {expert} Regularization steps
"""Number of steps to lower the regularization factor"""
SomSteps=5
# {expert} Additional kerdenSOM parameters:
""" For a complete description 
    See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraCommand=''

#------------------------------------------------------------------------------------------------
# {end_of_header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------

#
# Main
#     
from protocol_kerdensom import *

if __name__ == '__main__':
   protocolMain(ProtKerdensom)
