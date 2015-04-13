'''
#/***************************************************************************
# * Authors:     Joaquin Oton (joton@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
 '''
#------------------------------------------------------------------------------------------------
# Declaration of related X-ray protocols headers to be used either in independent protocols 
# and fast unattended protocols
# Author: Joaquin Oton, Oct 2013
 
def expandXrayImport():
    linesStr = '''
# {list_combo}(Alba, Bessy, Diamond) Synchrotron
"""Tomo series from Mistral microscope (Alba synchrotron) are expected to be
 in nexus format in hdf5 files, while from U41-TXM line (Bessy) you only set 
 the data folder and the indexes for initial and final images for both tomogram 
 and flatfields."""
Synchrotron = "Alba"

# {condition}(Synchrotron=="Alba") {list}(single file,folder) Import from 
ImportFrom = "single file"

# {condition}(Synchrotron=="Alba" and ImportFrom == "folder"){dir} Tomograms folder
"""Directory name from where to process all tomograms"""
DirTomograms = ''

# {condition}(Synchrotron=="Alba" and ImportFrom == "single file"){file}(*.hdf5){validate}(PathExists) Tomogram file
"""Tomogram Filename  to process"""
Tomogram = ''

# {condition}(Synchrotron=="Bessy"){dir} Raw data folder
"""Directory name from where to process all tomograms"""
DirBessyData = ''

# {condition}(Synchrotron=="Bessy"){validate}(IsInt) Initial tomo index
TIni = ''
# {condition}(Synchrotron=="Bessy"){validate}(IsInt) Final tomo index
TEnd = ''
# {condition}(Synchrotron=="Bessy"){validate}(IsInt) Initial flatfield index
FIni = ''
# {condition}(Synchrotron=="Bessy"){validate}(IsInt) Initial flatfield index
FEnd = ''

#------------------------------------------------------------------------------------------------
# {expert}{condition}{has_question}{section} Preprocess
#------------------------------------------------------------------------------------------------
# Preprocess tomograms?
DoPreprocess = True

# Log correction?
"""Because the intrinsic self-attenuation, pixel values of X-ray projections are 
not proportional to volume coefficients and, therefore, they are not optimal 
to be used with 3DEM reconstruction algorithm. This is fixed by applying: 
Ilog = log10(Inormalized)"""
DoLog = True

# {condition}(DoLog) Contrast inversion correction?
"""In addition to log correction, it applies a contrast inversion that allows 3DEM
 reconstruction algorithm returning volume coefficients close to the real expected
 absorption values. Icorrected = -log(Inormalized)"""
DoCorrect = False 

# Crop borders?
""" Crop a given amount of pixels from each border. """
DoCrop = False

# {condition}(DoCrop) Pixels to crop
""" Amount of pixels you want to crop from borders """
Crop = 10

# {list_combo}(No,Auto mask, From file)Remove bad pixels?
""" Nonnegative pixel values in the mask will be substituted by the local median.
 Auto Mask generates a mask from the flatfields average image being the threshold 
 mean - factor*std.
 From file read a mask from an image file."""
BadPixelMode = "Auto mask"

# {condition}(BadPixelMode=="Auto mask"){validate}(IsFloat) Factor threshold
BadPixelFactor = "3."    


# {condition}(BadPixelMode=="From file"){file} Mask file
BadPixelsMask = ""    
    '''
    return linesStr % locals()
 
 
 
def expandXrayFastAlign():
    linesStr = '''
#-----------------------------------------------------------------
# {expert}{condition}{has_question}{section} tiltxcorr
#-----------------------------------------------------------------
# Params
DoTiltXcorr = True
# rotation
rotation = 0.0
# sigma1
sigma1 = 0.03
# radius2
radius2 = 0.25
# sigma2
sigma2 = 0.05
# border
border = '49,49'
# size
size = '300,300'
# overlap
overlap = '0.33,0.33'

#-----------------------------------------------------------------
# {expert}{condition}{has_question}{section} tiltalign
#-----------------------------------------------------------------
# Params
DoTiltAlign = True
# RotationAngle
RotationAngle = 0.0
# AngleOffset
AngleOffset = 0.0 
# RotOption
RotOption = -1
# RotDefaultGrouping
RotDefaultGrouping = 5 
#TiltOption
TiltOption = 0 
# MagReferenceView
MagReferenceView = 1 
# MagOption
MagOption = 0 
# MagDefaultGrouping
MagDefaultGrouping = 4 
# XStretchOption
XStretchOption = 0 
# XStretchDefaultGrouping
XStretchDefaultGrouping = 7 
# SkewOption
SkewOption = 0 
# SkewDefaultGrouping
SkewDefaultGrouping = 11 
# ResidualReportCriterion
ResidualReportCriterion = 3.0  
# SurfacesToAnalyze
SurfacesToAnalyze = 1 
# MetroFactor
MetroFactor = 0.25 
# MaximumCycles
MaximumCycles = 1000 
# AxisZShift
AxisZShift = 0.0 
# LocalAlignments
LocalAlignments = 0 
# MinFidsTotalAndEachSurface
MinFidsTotalAndEachSurface = '8,3' 
# LocalOutputOptions
LocalOutputOptions = '0,0,0' 
# LocalRotOption
LocalRotOption = 3 
# LocalRotDefaultGrouping
LocalRotDefaultGrouping = 6 
# LocalTiltOption
LocalTiltOption = 5 
# LocalTiltDefaultGrouping
LocalTiltDefaultGrouping = 6 
# LocalMagReferenceView
LocalMagReferenceView = 1 
# LocalMagOption
LocalMagOption = 3 
# LocalMagDefaultGrouping
LocalMagDefaultGrouping = 7 
# LocalXStretchOption
LocalXStretchOption = 0 
# LocalXStretchDefaultGrouping
LocalXStretchDefaultGrouping = 7 
# LocalSkewOption
LocalSkewOption = 0
# LocalSkewDefaultGrouping
LocalSkewDefaultGrouping = 11 
# BeamTiltOption
BeamTiltOption = 0    
'''
    return linesStr % locals()   
 
 
 
 
def expandXrayReconstruct():
    linesStr = '''
# {list_combo}(imod, tomo3d) Reconstruction program
"""Tomo series from Mistral microscope (Alba synchrotron) are expected to be
 in nexus format in hdf5 files, while from U41-TXM line (Bessy) you only set 
 the data folder and the indexes for initial and final images for both tomogram 
 and flatfields."""
recProgram = "tomo3d"

# {validate}(IsInt) Thickness
thickness = ''

# {expert} Contrast inversion
"""Before reconstruct invert the contrast of the images. It is mandatory to 
   recover real absorption coefficient volumes."""
DoInvertContrast = True    
    '''
    return linesStr % locals()    
 
 
