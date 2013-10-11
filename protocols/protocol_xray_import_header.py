#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

#-----------------------------------------------------------------------------
# {section} Input tomograms
#-----------------------------------------------------------------------------

# {list_combo}(Mistral, Bessy) Tomogram provider
"""Tomo series from Mistral microscope (Alba synchrotron) are expected to be
 in nexus format in hdf5 files, while from U41-TXM line (Bessy) you only set 
 the data folder and the indexes for initial and final images for both tomogram 
 and flatfields."""
TomogramProvider = "Mistral"

# {condition}(TomogramProvider=="Mistral") {list}(single file,folder) Import from 
ImportFrom = "single file"

# {condition}(TomogramProvider=="Mistral" and ImportFrom == "folder"){dir} Tomograms folder
"""Directory name from where to process all tomograms"""
DirTomograms = ''

# {condition}(TomogramProvider=="Mistral" and ImportFrom == "single file"){file} Tomogram file
"""Directory name from where to process all tomograms"""
Tomogram = ''

# {condition}(TomogramProvider=="Bessy"){dir} Raw data folder
"""Directory name from where to process all tomograms"""
DirBessyData = ''

# {condition}(TomogramProvider=="Bessy"){validate}(IsInt) Initial tomo index
TIni = ''
# {condition}(TomogramProvider=="Bessy"){validate}(IsInt) Final tomo index
TEnd = ''
# {condition}(TomogramProvider=="Bessy"){validate}(IsInt) Initial flatfield index
FIni = ''
# {condition}(TomogramProvider=="Bessy"){validate}(IsInt) Initial flatfield index
FEnd = ''

#------------------------------------------------------------------------------------------------
# {condition}{has_question}{section} Preprocess
#------------------------------------------------------------------------------------------------
# Preprocess tomograms?
DoPreprocess = True

# Log correction?
"""Because the intrinsic self-attenuation, pixel values of X-ray projections are 
not proportional to volume coefficients and, therefore, they are not optimal 
to be used with 3DEM reconstruction algorithm. This is fixed by applying: 
Ilog = log10(Inormalized)"""
DoLog = False

# {condition}(DoLog) Contrast inversion correction?
"""In addition to log correction, it applies a contrast inversion that allows 3DEM
 reconstruction algorithm returning volume coefficients close to the real expected
 absorption values. Icorrected = -log10(Inormalized)"""
DoCorrect = False 

# Crop borders?
""" Crop a given amount of pixels from each border. """
DoCrop = False

# {condition}(DoCrop) Pixels to crop
""" Amount of pixels you want to crop from borders """
Crop = 10

# Remove bad pixels from mask?
""" Nonnegative pixel values in the mask will be substituted by the local median."""
DoBadPixelsMask = False

# {condition}(DoBadPixelsMask){file} Mask file
BadPixelsMask = ""

# {eval} expandParallel(threads=0,hours=6)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
from protocol_xray_import import *
if __name__ == '__main__':
    protocolMain(ProtXrayImport)
