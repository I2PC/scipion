# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This modules contains constants related to Relion protocols
"""

#------------------ Constants values --------------------------------------

import pyworkflow.em.metadata as md


MASK_FILL_ZERO = 0
MASK_FILL_NOISE = 1

ANGULAR_SAMPLING_LIST = ['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5', '0.2', '0.1']


# Map from Xmipp labels to Relion labels names
XMIPP_RELION_LABELS = {
                        md.MDL_ANGLE_ROT:         'rlnAngleRot'
                       ,md.MDL_ANGLE_TILT:        'rlnAngleTilt'
                       ,md.MDL_ANGLE_PSI:         'rlnAnglePsi'
                       ,md.MDL_AVG_CHANGES_ORIENTATIONS:'rlnChangesOptimalOrientations'
                       ,md.MDL_AVG_CHANGES_OFFSETS:     'rlnChangesOptimalOffsets'
                       ,md.MDL_AVG_CHANGES_CLASSES:     'rlnChangesOptimalClasses'
                       ,md.MDL_AVGPMAX:           'rlnAveragePmax'
                       ,md.MDL_CLASS_PERCENTAGE:  'rlnClassDistribution'
                       ,md.MDL_CTF_CA:            'rlnChromaticAberration'
                       ,md.MDL_CTF_CONVERGENCE_CONE: 'rlnConvergenceCone'
                       ,md.MDL_CTF_CS:            'rlnSphericalAberration'
                       ,md.MDL_CTF_DEFOCUSU:      'rlnDefocusU'
                       ,md.MDL_CTF_DEFOCUSV:      'rlnDefocusV'
                       ,md.MDL_CTF_DEFOCUS_ANGLE: 'rlnDefocusAngle'
                       ,md.MDL_CTF_ENERGY_LOSS:   'rlnEnergyLoss'
                       ,md.MDL_CTF_LENS_STABILITY:'rlnLensStability'
                       ,md.MDL_CTF_LONGITUDINAL_DISPLACEMENT:'rlnLongitudinalDisplacement'
                       ,md.MDL_CTF_TRANSVERSAL_DISPLACEMENT: 'rlnTransversalDisplacement'                       
                       ,md.MDL_CTF_Q0:            'rlnAmplitudeContrast'
                       ,md.MDL_CTF_SAMPLING_RATE: 'rlnDetectorPixelSize'
                       ,md.MDL_CTF_VOLTAGE:       'rlnVoltage'
                       ,md.MDL_DATATYPE:          'rlnDataType'
                       ,md.MDL_DEFGROUP:          'rlnGroupNumber'
                       ,md.MDL_ENABLED:           'rlnEnabled'
                       ,md.MDL_IMAGE:             'rlnImageName'
                       ,md.MDL_IMAGE_REF:         'rlnReferenceImage'
                       ,md.MDL_ITEM_ID:           'rlnImageId'
                       ,md.MDL_MAXCC:             'rlnLogLikelihood'
                       ,md.MDL_PSD:               'rlnCtfImage' # relion 1.3
                       ,md.MDL_LL:                'rlnLogLikeliContribution'
                       ,md.MDL_MAGNIFICATION:     'rlnMagnification'
                       ,md.MDL_MICROGRAPH:        'rlnMicrographName'
                       ,md.MDL_MICROGRAPH_ID:     'rlnMicrographId'
                       ,md.MDL_REF:               'rlnClassNumber'
                       ,md.MDL_RESOLUTION_FREQREAL:'rlnAngstromResolution'
                       ,md.MDL_RESOLUTION_FRC:     'rlnGoldStandardFsc'
                       ,md.MDL_RESOLUTION_FREQ:    'rlnResolution'
                       ,md.MDL_RESOLUTION_SSNR:    'rlnSsnrMap'
                       ,md.MDL_RANDOMSEED:         'rlnRandomSubset'
                       ,md.MDL_SAMPLINGRATE:       'rlnPixelSize'
                       #,md.MDL_SAMPLINGRATE_ORIGINAL: 'rlnPixelSize'
                       ,md.MDL_SCALE:              'rlnMagnificationCorrection'
                       ,md.MDL_SHIFT_X:           'rlnOriginX'
                       ,md.MDL_SHIFT_Y:           'rlnOriginY'
                       ,md.MDL_SHIFT_Z:           'rlnOriginZ'
                       ,md.MDL_PMAX:               'rlnMaxValueProbDistribution'
                       ,md.MDL_SAMPLINGRATE_X:     'rlnSamplingRateX'
                       ,md.MDL_SAMPLINGRATE_Y:     'rlnSamplingRateY'
                       ,md.MDL_SAMPLINGRATE_Z:     'rlnSamplingRateZ'
                       ,md.MDL_XCOOR:              'rlnCoordinateX'
                       ,md.MDL_XSIZE:              'rlnImageSizeX'
                       ,md.MDL_YCOOR:              'rlnCoordinateY'
                       ,md.MDL_YSIZE:              'rlnImageSizeY'
                       ,md.MDL_WEIGHT:             'rlnNrOfSignificantSamples'
                       ,md.MDL_ZSIZE:              'rlnImageSizeZ'
                       # relion 1.3
                       ,md.MDL_IMAGE2: 'rlnParticleName'
                       ,md.MDL_IMAGE_ORIGINAL: 'rlnOriginalParticleName'
                       ,md.MDL_SERIE: 'rlnGroupName'
                       }

XMIPP_RELION_LABELS_EXTRA = {
                       # Following labels have no correct matching, just to 
                       # pick one with the same datatype
                       md.MDL_ANGLE_Y2 : 'rlnOrientationalPriorMode' #int
                       ,md.MDL_ANGLE_ROT2 : 'rlnSigmaPriorRotAngle' #double
                       ,md.MDL_ANGLE_TILT2 : 'rlnSigmaPriorTiltAngle' #double
                       ,md.MDL_ANGLE_PSI2 : 'rlnSigmaPriorPsiAngle' #double
                       ,md.MDL_BLOCK_NUMBER:       'rlnGroupNumber' # just one
                       ,md.MDL_COUNT:             'rlnGroupNrParticles' # just one
                       ,md.MDL_CTF_CRIT_FITTINGSCORE:   'rlnCtfFigureOfMerit' #just one
                       ,md.MDL_CTF_CRIT_NORMALITY:   'rlnNormCorrection' #just one
                       ,md.MDL_CTF_CRIT_PSDVARIANCE: 'rlnCtfValue'         #just one
                       ,md.MDL_CTF_CRIT_PSDCORRELATION90: 'rlnCtfBfactor'  #just one
                       ,md.MDL_CRYSTAL_CELLX : 'rlnReferenceDimensionality'
                       ,md.MDL_CRYSTAL_CELLY : 'rlnOriginalImageSize'
                       ,md.MDL_CRYSTAL_DISAPPEAR_THRE : 'rlnCurrentResolution'
                       ,md.MDL_PICKING_PARTICLE_SIZE : 'rlnCurrentImageSize' #int
                       ,md.MDL_PICKING_AUTOPICKPERCENT : 'rlnPaddingFactor' #int
                       ,md.MDL_PICKING_TEMPLATES : 'rlnFourierSpaceInterpolator' #int
                       ,md.MDL_COLOR : 'rlnMinRadiusNnInterpolation' #int
                       ,md.MDL_DM3_NUMBER_TYPE : 'rlnNrClasses' #int
                       ,md.MDL_DM3_SIZE : 'rlnNrGroups' #int
                       ,md.MDL_NMA_COLLECTIVITY : 'rlnTau2FudgeFactor' #double
                       ,md.MDL_NMA_SCORE : 'rlnNormCorrectionAverage' #double
                       ,md.MDL_SIGMAOFFSET : 'rlnSigmaOffsets' #double
                       
                       ,md.MDL_MLF_WIENER: 'rlnOrientationDistribution' #double
                       ,md.MDL_IDX: 'rlnSpectralIndex' #int
                       ,md.MDL_MLF_NOISE: 'rlnSigma2Noise' #double
                       ,md.MDL_DM3_TAGNAME: 'rlnGroupName'  #string
                       ,md.MDL_MLF_SIGNAL: 'rlnGroupScaleCorrection' #double
                       
                       ,md.MDL_FOM: 'rlnAutopickFigureOfMerit'
                       ,md.MDL_ZSCORE_SHAPE1: 'rlnAccuracyRotations'
                       ,md.MDL_ZSCORE_SHAPE2: 'rlnAccuracyTranslations'
                       ,md.MDL_ZSCORE_SNR1: 'rlnReferenceSigma2'
                       ,md.MDL_ZSCORE_SNR2: 'rlnReferenceTau2'
                       ,md.MDL_ZSCORE_RESCOV: 'rlnSpectralOrientabilityContribution'

                       # Keep relion 1.3 new labels at the end
                       # just in case we want to keep 1.2 compatibility
                       ,md.MDL_ANGLE_ROT2:         'rlnAngleRotPrior'
                       ,md.MDL_ANGLE_TILT2:        'rlnAngleTiltPrior'
                       ,md.MDL_ANGLE_PSI2:         'rlnAnglePsiPrior'                       
                       ,md.MDL_SHIFT_X2 : 'rlnOriginXPrior' 
                       ,md.MDL_SHIFT_Y2 : 'rlnOriginYPrior' 
                       ,md.MDL_ZSCORE: 'rlnParticleSelectZScore' 
                       ,md.MDL_COUNT2: 'rlnNrOfFrames'
                       # Not the best labels, but just to grab some
                       ,md.MDL_CLASSIFICATION_DPR_05: 'rlnClassPriorOffsetX'
                       ,md.MDL_CLASSIFICATION_FRC_05: 'rlnClassPriorOffsetY'
                       #,md.MDL_CLASSIFICATION_INTRACLASS_DISTANCE: ''
                       }

# Supported versions:
V1_3 = '1.3'
V1_4 = '1.4'
V2_0 = '2.0'

