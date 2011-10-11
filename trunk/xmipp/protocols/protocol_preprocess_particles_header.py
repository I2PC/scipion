#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
# Author: Carlos Oscar, August 2011
#
# {begin_of_header}

# {eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Filter parameters
#-----------------------------------------------------------------------------
# {file}{validate}(PathExists) Selfile with the input images:
""" This selfile points to the stack or metadata containing your images 
"""
InSelFile=''

# Scale
""" Change the scale of the input images """
DoScale=False

#{condition}(DoScale) New image size
NewSize=0

# Fourier bandpass filter
""" You may do a lowpass filter by setting Freq_low to 0. You may do a high pass filter
    by setting Freq_high to 0.5."""
DoFourier=False 

#{condition}(DoFourier){wizard}(wizardChooseBandPassFilter) Freq_low (0<f<0.5)
""" Set it to 0 for low pass filters """
Freq_low=0.02

#{condition}(DoFourier) Freq_high (0<f<0.5)
""" Set it to 0.5 for high pass filters """
Freq_high=0.35

#{condition}(DoFourier){expert} Freq_decay (0<f<0.5)
""" It describes the length of the amplitude decay in a raised cosine """
Freq_decay=0.02

# Fourier Gaussian
""" Gaussian filter defined in Fourier space"""
DoGaussian=False

#{condition}(DoGaussian){wizard}(wizardChooseGaussianFilter) Frequency sigma
""" Remind that the Fourier frequency is normalized between 0 and 0.5"""
Freq_sigma=0.04

# Remove dust
""" Sets pixels with unusually large values to random values from a Gaussian with zero-mean and unity-standard deviation.
    It requires a previous normalization, i.e., Normalization must be set to Yes.
"""
DoRemoveDust=False

# {condition}(DoRemoveDust){wizard}(wizardChooseBadPixelsFilter)  Threshold for dust removal:
""" Pixels with a signal higher or lower than this value times the standard deviation of the image will be affected. For cryo, 3.5 is a good value.
    For high-contrast negative stain, the signal itself may be affected so that a higher value may be preferable.
"""
DustRemovalThreshold=3.5

# Normalize
""" It subtract a ramp in the gray values and normalizes so that in the background
    there is 0 mean and standard deviation 1 """
DoNorm=False

# {list_combo}(OldXmipp,NewXmipp,Ramp){condition}(DoNorm) Normalization type
""" OldXmipp (mean(Image)=0, stddev(Image)=1). NewXmipp (mean(background)=0, stddev(background)=1), Ramp (subtract background+NewXmipp)"""
NormType="Ramp"

# {condition}(DoNorm and NormType!="OldXmipp") Background radius
"""Pixels outside this circle are assumed to be noise and their stddev is set to 1.
   Radius for background circle definition (in pix.).
   If this value is 0, then the same as the particle radius is used. """
BackGroundRadius=0

# Apply mask
""" Apply mask from file """
DoMask=False

# {condition}(DoMask){wizard}(wizardDesignMask) Mask file
MaskFile=""

# {condition}(DoMask) Substitute value
"""Valid values are a number, min, max and avg """
Substitute="0"

# {eval} expandParallel(mpi=8, threads=0, hours=6)
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_preprocess_particles import *
#        
# Main
#     
 
if __name__ == '__main__':
    # create preprocess_particles_class object
    protocolMain(ProtPreprocessParticles)
