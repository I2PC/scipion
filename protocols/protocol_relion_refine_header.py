#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Relion-based 3D classification
#
# Author: Roberto Marabini            roberto@cnb.csic.es    May 2013
#         J. M. de la Rosa Trevin     jmdelarosa@cnb.csic.es
#
#------------------------------------------------------------------------------------------------
# {begin_of_header}

# {eval} expandCommentRun()
# {cite}
CiteRelion3D = """
"A Bayesian view on cryo-EM structure determination" 
by Sjors H.W. Scheres (DOI: 10.1016/j.jmb.2011.11.010)
"""

#------------------------------------------------------------------------------------------
# {section} Input
#------------------------------------------------------------------------------------------
# {file}(images*.xmd){validate}(PathExists) Input images:
""" 
Provide a list of images from a stack <(Spider/MRC)> or metadata file that make up your data set.
The filenames should be relative to the <ProjectDir> where you are running the <Protocols>
If you have images outside the <ProjectDir> you should import them first.
Note that in the case of a stack, no metadata can be included and thus no CTF correction
can be performed.Nor will it be possible to perform noise spectra estimation or intensity
scale corrections in image groups.
"""
ImgMd = ""

# {hidden}Continue from here:
""" 
Select the *_optimiser.star file for the iteration from which you want to continue a previous run. 
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. 
If they are the same, the program will automatically add a '_ctX' to the output rootname, 
with X being the iteration from which one continues the previous run.Provide a list of images 
from a stack <(Spider/MRC)> or metadata file that make up your data set.
"""
ContinueFrom = ""

# {hidden}Number of classes
"""The number of classes (K) for a multi-reference refinement. These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.
"""
NumberOfClasses = 1

# {file}(*.vol, *.mrc){validate}(PathExists) Initial 3D reference volume:
"""
A 3D map in MRC/Spider format. Make sure this map has the same dimensions and 
the same pixel size as your input images.
"""
Ref3D = ""

# Ref. map is on absolute greyscale?
""" {expert}
Set this option to False unless you know what you are doing.
The probabilities are based on squared differences, so that the absolute grey scale is important.
Probabilities are calculated based on a Gaussian noise model, 
which contains a squared difference term between the reference and the experimental image. 
This has a consequence that the reference needs to be on the same absolute intensity 
grey-scale as the experimental images. RELION and XMIPP reconstruct maps at their absolute
intensity grey-scale. Other packages may perform internal normalisations of the reference 
density, which will result in incorrect grey-scales. Therefore: if the map was reconstructed
in RELION or in XMIPP, set this option to Yes, otherwise set it to No. If set to No, RELION 
will use a (grey-scale invariant) cross-correlation criterion in the first iteration, and 
prior to the second iteration the map will be filtered again using the initial low-pass filter.
This procedure is relatively quick and typically does not negatively affect the outcome of the
subsequent MAP refinement. Therefore, if in doubt it is recommended to set this option to No."""
IsMapAbsoluteGreyScale = False

#normalize input images
""" 
Normalize input images ?
average background value must be 0 and a stddev value must be 1. 
Note that the average and stddev values for the background are
 calculated outside a circle with the particle diameter 
"""
DoNormalizeInputImage = False

# {expert} Padding factor:
"""
The padding factor used for oversampling of the Fourier transform. The default is 3x padding, 
which is combined with nearest-neighbour interpolation. However, for large 3D maps, storing the 
entire 3x oversampled Fourier transform (as doubles) plus the real-space padded map in memory may 
be too much. Therefore, for large maps or in cases of many 3D references, in order to fit into memory 
one may need to use a smaller padding factor: e.g. 2, or (not recommended) 1. For padding factors smaller
than 3, (slower) linear interpolation will be used.
 
The approximate amount of memory (in Gb) required to store K maps of (size x size x size) voxels and a 
padding factor (pad) may be calculated as: 
K*2*8*(size*pad)^3/1024/1024/1024 
 
<Note>: also consider the use of threads if memory is an issue.
"""
PaddingFactor = 3.0

# Symmetry group
""" 
See [http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry]
for a description of the symmetry groups format
If no symmetry is present, give c1
"""
SymmetryGroup = 'c1'

#------------------------------------------------------------------------------------------------
# {section}{has_question} CTF
#------------------------------------------------------------------------------------------------

# Use CTF-amplitude correction?
"""
If set to Yes, CTFs will be corrected inside the MAP refinement. 
The resulting algorithm intrinsically implements the optimal linear, 
or Wiener filter. Note that CTF parameters for all images need to
be given in the input STAR file. 
"""
DoCTFCorrection = True

# Has reference been CTF-corrected?
"""
Set this option to Yes if the reference map represents CTF-unaffected density, 
e.g. it was created using Wiener filtering inside RELION or from a PDB. If set to No, 
then in the first iteration, the Fourier transforms of the reference projections 
are not multiplied by the CTFs.
"""
HasReferenceCTFCorrected = False

# {expert} Only flip phases?
"""
Set this to Yes to switch CTF-amplitude correction off. 
This option is NOT generally recommended.
"""
OnlyFlipPhases = False

# Have data been phase-flipped?
"""
Set this to Yes if the images have been ctf-phase corrected during the pre-processing steps. 
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, 
as this can be done inside the internal CTF-correction. However, if the phases do have been flipped, 
one should tell the program about it using this option.
"""
HaveDataPhaseFlipped = False

# Ignore CTFs until first peak?
"""
If set to Yes, then CTF-amplitude correction will only be performed from the first peak 
of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. 
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results.
Therefore, this option is not generally recommended.
"""
IgnoreCTFUntilFirstPeak = False

# {expert} Do intensity correction?
"""
An internal correction for differences in the intensity (grey-scale) of the signal between 
distinct micrographs is applied. This is useful if micrographs have very different 
signal-to-noise ratios, e.g. due to different ice thickness or contamination. 
Because one typically normalises the noise, this leads to distinct signal intensities in the data, 
and this procedure corrects for this. It is quite robust and therefore recommended for the general case.
"""
DoIntensityCorrection = False



#------------------------------------------------------------------------------------------------
# {section} Optimisation
#------------------------------------------------------------------------------------------------

# Initial low-pass filter (A): 
"""
It is recommended to strongly low-pass filter your initial reference map. 
If it has not yet been low-pass filtered, it may be done internally using this option. 
If set to 0, no low-pass filter will be applied to the initial reference(s).
"""
InitialLowPassFilterA = 60

# {hidden} Number of iterations:
"""
Number of iterations to be performed. Note that the current implementation does NOT comprise a convergence criterium. Therefore, the calculations will need to be stopped by the user if further iterations do not yield improvements in resolution or classes.
"""
NumberOfIterations = 1

# Particles mask diameter (A):
"""
The experimental images will be masked with a soft circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! If set to a value larger than the image size no masking will be performed.
"""
MaskDiameterA = 200

# Mask references structures?
"""
If set to yes, a mask will also be applied to the reconstructed references. This is useful to set the solvent region of your reconstruction to 0. Either a soft spherical mask (based on the diameter of the experimental image mask given above) or a user-provided mask (next option) may be used. The user-provided mask should have values between 0 and 1 only. Solvent flattening is recommended, but make sure not to mask any signal away.
"""
DoMaskParticles = True

# {file} Reference mask
"""
A Spider/mrc map containing a (soft) mask with the same dimensions as the reference(s), and values between 0 and 1, with 1 being 100% protein and 0 being 100% solvent. The reconstructed reference map will be multiplied by this mask.If no mask is given, a soft spherical mask based on the diameter of the mask for the experimental images will be applied.  
 
In some cases, for example for non-empty icosahedral viruses, it is also useful to use a second mask. For all white (value 1) pixels in this second mask the corresponding pixels in the reconstructed map are set to the average value of these pixels. Thereby, for example, the higher density inside the virion may be set to a constant. Note that this second mask should have one-values inside the virion and zero-values in the capsid and the solvent areas. To use a second mask, use the additional option --solvent_mask2, which may given in the Additional arguments line (in the Sampling tab).
"""
ReferenceMask = ""



#-----------------------------------------------------------------------------
# {section} Sampling: initial sampling rates will be auto-incremented
#-----------------------------------------------------------------------------

# 
# {list_combo}(30,15,7.5,3.7,1.8,0.9,0.5,0.2,0.1) Angular sampling interval (deg):
"""There are only a few discrete angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. The samplings are approximate numbers and vary slightly over the sphere.
"""
AngularSamplingDeg = '7.5'

# Offset search range (pix):
"""Probabilities will be calculated only for translations in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation for each image in the previous iteration.
"""
OffsetSearchRangePix = 5

# Offset search step (pix):
"""Translations will be sampled with this step-size (in pixels).Translational sampling is also done using the adaptive approach. Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.
"""
OffsetSearchStepPix = 1

# {list_combo}(30,15,7.5,3.7,1.8,0.9,0.5,0.2,0.1)  Local angular search range
"""
"""
LocalAngularSearchRange = 5.0

# {expert} Additional arguments
"""In this box command-line arguments may be provided that are not generated by the GUI. This may be useful for testing developmental options and/or expert use of the program. The command 'relion_refine' will print a list of possible options.
"""
AdditionalArguments = ""


# {eval} expandParallel()

#------------------------------------------------------------------------------------------------
# {section}{visualize} Preparation
#------------------------------------------------------------------------------------------------

# {view} Show the grey-scale corrected ref. map?
VisualizeCRVolume = True

# {view} Show the low-pass filtered ref. map?
VisualizeFRVolume = True

# {view} Show the generated seeds volumes?
VisualizeGSVolume = True

#------------------------------------------------------------------------------------------------
# {section}{visualize} Results per Iteration and Ref3D
#------------------------------------------------------------------------------------------------
# {hidden}{list_combo}( all, selection) Which ref3D you want to visualize?
""" 
   If you want to see the reference volume 2 and 5 write
   2 5. In relaion first reference is 1
"""
DisplayRef3DNo='all'

# {condition}(DisplayRef3DNo=='selection') Selected references 3D
""" Which iteration do you want to visualize 
If you want two see iterations 2 and 5 write
   2 5. In relion first iteration is 0"""
SelectedRef3DNo = ''

# {list_combo}(last, all, selection) Which iteration you want to visualize?
VisualizeIter = 'last'

# {condition}(VisualizeIter=='selection') Selected iterations
""" Which iteration do you want to visualize 
If you want two see iterations 2 and 5 write
   2 5. In relion first iteration is 0"""
SelectedIters = ''

# {list_combo}(x, y, z, surface) Display 3D maps along 
""" x -> Visualize volumes in slices along x
    y -> Visualize volumes in slices along y
    z -> Visualize volumes in slices along z
    surface: surface rendering, you need chimera installed!
"""
DisplayVolumeSlicesAlong='z'

# {view} Display reconstructed volume
""" Volume as given by the reconstruction algorithm.
    For each iteration two volumes are presented (each one reconstructed with
    the 50% of the images) plus a single final colume reconstructed using 
    all the data)
"""
DisplayReconstruction=False

#------------------------------------------------------------------------------------------------
# {section}{visualize} Overall Results
#------------------------------------------------------------------------------------------------
# {view} Display resolution plots (SSNR)
DisplayResolutionPlotsSSNR=True

###############################
# {view} Display resolution plots (FSC)
DisplayResolutionPlotsFSC=False

# {expert} Display a threshold in resolution plots (FSC)
ResolutionThreshold=0.5

# {view} Display angular distribution?
DisplayAngularDistribution=True
# {condition}(DisplayAngularDistribution) {list} (2D, 3D) Display Angular distribution in
""" 2D option uses matplotlib while 3D uses chimera
"""
DisplayAngularDistributionWith='2D'

# {expert} Scale RedSpheres
""" when using chimera for displaying red angular
distribution set radius of maximum sphere"""
SpheresMaxradius=-1.

# {hidden}{view} Display resolution plots (FSC)
DisplayResolutionPlotsFSC=False

# {hidden} No. Images assigned to class
""" Images assigned to each class per iteration"""
TableImagesPerClass=True

# {view} Changes Offset, Ang, No Part
""" changes in orientation, offset. number images assigned to each class"""
TableChange=True

# {view} LikeliHood Per Image
""" Max likelihood per image may be used to delete images with smaller value. 
The higher, the better. Considere remove particles with low values"""
Likelihood=True

# {view} AveragePmax
""" Average (per class) of the maximum value of normalized probability function """
AvgPMAX=False
#------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#------------------------------------------------------------------------------------------------

from protocol_relion_refine import *

if __name__ == '__main__':
    protocolMain(ProtRelionRefinner)

###############################
# {view} Display resolution plots (FSC)
#DisplayResolutionPlotsFSC=True

# {expert} Display a threshold in resolution plots (FSC)
#ResolutionThreshold=0.5
