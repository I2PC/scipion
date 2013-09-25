#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of volumes 
# Author: Carlos Oscar, August 2013
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input
#-----------------------------------------------------------------------------
# {file}(*.vol){validate}(PathExists) Input volume:
"""This volume will be compared to the reference volume"""
InputVol = ''

# Calculate FSC and DPR
DoFSC=True

# {file}(*.vol){validate}(PathExists){condition}(DoFSC) Reference volume:
ReferenceVol = ''

# Calculate Structure Factor
DoStructureFactor=True

# Calculate Spectral SNR
""" The input volume for this option must be the raw output after reconstruction (reconstruction_Ref3D_???.vol) of
any of the iterations in a projection matching"""
DoSSNR=True

# Calculate Volumetric Spectral SNR
""" The input volume for this option must be the raw output after reconstruction (reconstruction_Ref3D_???.vol) of
any of the iterations in a projection matching"""
DoVSSNR=False

# {eval} expandParallel(jobsize=1)

#------------------------------------------------------------------------------------------------
# {section}{visualize} Visualization
#------------------------------------------------------------------------------------------------

# {condition}(DoFSC) Display Fourier Shell Correlation
DisplayFSC=True

# {condition}(DoFSC) Display Differential Phase Residual
DisplayDPR=True

# {condition}(DoStructureFactor) Display Structure factor
DisplayStructureFactor=True

# {condition}(DoStructureFactor) Display Guinier plot
DisplayGuinier=True

#{expert}{condition}(DisplayGuinier and DoStructureFactor) Use Matlab for Guinier
UseMatlab=False

# {condition}(DoSSNR) Display Spectral SNR
DisplaySSNR=True

# {condition}(DoVSSNR) Display Volumetric Spectral SNR
DisplayVSSNR=True

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_resolution3D import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtResolution3D)
