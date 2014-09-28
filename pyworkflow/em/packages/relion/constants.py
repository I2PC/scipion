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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This modules contains constants related to Relion protocols
"""

#------------------ Constants values --------------------------------------

import xmipp


MASK_FILL_ZERO = 0
MASK_FILL_NOISE = 1

ANGULAR_SAMPLING_LIST = ['30', '15', '7.5', '3.7', '1.8', '0.9', '0.5', '0.2', '0.1']


# Map from Xmipp labels to Relion labels names
XMIPP_RELION_LABELS = {
                        xmipp.MDL_ANGLE_ROT:         'rlnAngleRot'
                       ,xmipp.MDL_ANGLE_TILT:        'rlnAngleTilt'
                       ,xmipp.MDL_ANGLE_PSI:         'rlnAnglePsi'
                       ,xmipp.MDL_AVG_CHANGES_ORIENTATIONS:'rlnChangesOptimalOrientations'
                       ,xmipp.MDL_AVG_CHANGES_OFFSETS:     'rlnChangesOptimalOffsets'
                       ,xmipp.MDL_AVG_CHANGES_CLASSES:     'rlnChangesOptimalClasses'
                       ,xmipp.MDL_AVGPMAX:           'rlnAveragePmax'
                       ,xmipp.MDL_CLASS_PERCENTAGE:  'rlnClassDistribution'
                       ,xmipp.MDL_CTF_CA:            'rlnChromaticAberration'
                       ,xmipp.MDL_CTF_CONVERGENCE_CONE: 'rlnConvergenceCone'
                       ,xmipp.MDL_CTF_CS:            'rlnSphericalAberration'
                       ,xmipp.MDL_CTF_DEFOCUSU:      'rlnDefocusU'
                       ,xmipp.MDL_CTF_DEFOCUSV:      'rlnDefocusV'
                       ,xmipp.MDL_CTF_DEFOCUS_ANGLE: 'rlnDefocusAngle'
                       ,xmipp.MDL_CTF_ENERGY_LOSS:   'rlnEnergyLoss'
                       ,xmipp.MDL_CTF_LENS_STABILITY:'rlnLensStability'
                       ,xmipp.MDL_CTF_LONGITUDINAL_DISPLACEMENT:'rlnLongitudinalDisplacement'
                       ,xmipp.MDL_CTF_TRANSVERSAL_DISPLACEMENT: 'rlnTransversalDisplacement'                       
                       ,xmipp.MDL_CTF_Q0:            'rlnAmplitudeContrast'
                       ,xmipp.MDL_CTF_SAMPLING_RATE: 'rlnDetectorPixelSize'
                       ,xmipp.MDL_CTF_VOLTAGE:       'rlnVoltage'
                       ,xmipp.MDL_DATATYPE:          'rlnDataType'
                       ,xmipp.MDL_DEFGROUP:          'rlnGroupNumber'
                       ,xmipp.MDL_ENABLED:           'rlnEnabled'
                       ,xmipp.MDL_IMAGE:             'rlnImageName'
                       ,xmipp.MDL_IMAGE_REF:         'rlnReferenceImage'
                       ,xmipp.MDL_ITEM_ID:           'rlnImageId'
                       ,xmipp.MDL_MAXCC:             'rlnLogLikelihood'
                       ,xmipp.MDL_PSD:               'rlnCtfImage' # relion 1.3
                       ,xmipp.MDL_LL:                'rlnLogLikeliContribution'
                       ,xmipp.MDL_MAGNIFICATION:     'rlnMagnification'
                       ,xmipp.MDL_MICROGRAPH:        'rlnMicrographName'
                       ,xmipp.MDL_REF:               'rlnClassNumber'
                       ,xmipp.MDL_RESOLUTION_FREQREAL:'rlnAngstromResolution'
                       ,xmipp.MDL_RESOLUTION_FRC:     'rlnGoldStandardFsc'
                       ,xmipp.MDL_RESOLUTION_FREQ:    'rlnResolution'
                       ,xmipp.MDL_RESOLUTION_SSNR:    'rlnSsnrMap'
                       ,xmipp.MDL_RANDOMSEED:         'rlnRandomSubset'
                       ,xmipp.MDL_SAMPLINGRATE:       'rlnPixelSize'
                       #,xmipp.MDL_SAMPLINGRATE_ORIGINAL: 'rlnPixelSize'
                       ,xmipp.MDL_SCALE:              'rlnMagnificationCorrection'
                       ,xmipp.MDL_SHIFT_X:           'rlnOriginX'
                       ,xmipp.MDL_SHIFT_Y:           'rlnOriginY'
                       ,xmipp.MDL_SHIFT_Z:           'rlnOriginZ'
                       ,xmipp.MDL_PMAX:               'rlnMaxValueProbDistribution'
                       ,xmipp.MDL_SAMPLINGRATE_X:     'rlnSamplingRateX'
                       ,xmipp.MDL_SAMPLINGRATE_Y:     'rlnSamplingRateY'
                       ,xmipp.MDL_SAMPLINGRATE_Z:     'rlnSamplingRateZ'
                       ,xmipp.MDL_XCOOR:              'rlnCoordinateX'
                       ,xmipp.MDL_XSIZE:              'rlnImageSizeX'
                       ,xmipp.MDL_YCOOR:              'rlnCoordinateY'
                       ,xmipp.MDL_YSIZE:              'rlnImageSizeY'
                       ,xmipp.MDL_WEIGHT:             'rlnNrOfSignificantSamples'
                       ,xmipp.MDL_ZSIZE:              'rlnImageSizeZ'
                       # relion 1.3
                       ,xmipp.MDL_IMAGE2: 'rlnParticleName'
                       ,xmipp.MDL_IMAGE_ORIGINAL: 'rlnOriginalParticleName'
                       }

XMIPP_RELION_LABELS_EXTRA = {
                       # Following labels have no correct matching, just to 
                       # pick one with the same datatype
                       xmipp.MDL_ANGLE_Y2 : 'rlnOrientationalPriorMode' #int
                       ,xmipp.MDL_ANGLE_ROT2 : 'rlnSigmaPriorRotAngle' #double
                       ,xmipp.MDL_ANGLE_TILT2 : 'rlnSigmaPriorTiltAngle' #double
                       ,xmipp.MDL_ANGLE_PSI2 : 'rlnSigmaPriorPsiAngle' #double
                       ,xmipp.MDL_BLOCK_NUMBER:       'rlnGroupNumber' # just one
                       ,xmipp.MDL_COUNT:             'rlnGroupNrParticles' # just one
                       ,xmipp.MDL_CTF_CRIT_FITTINGSCORE:   'rlnCtfFigureOfMerit' #just one
                       ,xmipp.MDL_CTF_CRIT_NORMALITY:   'rlnNormCorrection' #just one
                       ,xmipp.MDL_CTF_CRIT_PSDVARIANCE: 'rlnCtfValue'         #just one
                       ,xmipp.MDL_CTF_CRIT_PSDCORRELATION90: 'rlnCtfBfactor'  #just one
                       ,xmipp.MDL_CRYSTAL_CELLX : 'rlnReferenceDimensionality'
                       ,xmipp.MDL_CRYSTAL_CELLY : 'rlnOriginalImageSize'
                       ,xmipp.MDL_CRYSTAL_DISAPPEAR_THRE : 'rlnCurrentResolution'
                       ,xmipp.MDL_PICKING_PARTICLE_SIZE : 'rlnCurrentImageSize' #int
                       ,xmipp.MDL_PICKING_AUTOPICKPERCENT : 'rlnPaddingFactor' #int
                       ,xmipp.MDL_PICKING_TEMPLATES : 'rlnFourierSpaceInterpolator' #int
                       ,xmipp.MDL_COLOR : 'rlnMinRadiusNnInterpolation' #int
                       ,xmipp.MDL_DM3_NUMBER_TYPE : 'rlnNrClasses' #int
                       ,xmipp.MDL_DM3_SIZE : 'rlnNrGroups' #int
                       ,xmipp.MDL_NMA_COLLECTIVITY : 'rlnTau2FudgeFactor' #double
                       ,xmipp.MDL_NMA_SCORE : 'rlnNormCorrectionAverage' #double
                       ,xmipp.MDL_SIGMAOFFSET : 'rlnSigmaOffsets' #double
                       
                       ,xmipp.MDL_MLF_WIENER: 'rlnOrientationDistribution' #double
                       ,xmipp.MDL_IDX: 'rlnSpectralIndex' #int
                       ,xmipp.MDL_MLF_NOISE: 'rlnSigma2Noise' #double
                       ,xmipp.MDL_DM3_TAGNAME: 'rlnGroupName'  #string
                       ,xmipp.MDL_MLF_SIGNAL: 'rlnGroupScaleCorrection' #double
                       
                       ,xmipp.MDL_FOM: 'rlnAutopickFigureOfMerit'
                       ,xmipp.MDL_ZSCORE_SHAPE1: 'rlnAccuracyRotations'
                       ,xmipp.MDL_ZSCORE_SHAPE2: 'rlnAccuracyTranslations'
                       ,xmipp.MDL_ZSCORE_SNR1: 'rlnReferenceSigma2'
                       ,xmipp.MDL_ZSCORE_SNR2: 'rlnReferenceTau2'
                       ,xmipp.MDL_ZSCORE_RESCOV: 'rlnSpectralOrientabilityContribution'

                       # Keep relion 1.3 new labels at the end
                       # just in case we want to keep 1.2 compatibility
                       ,xmipp.MDL_ANGLE_ROT2:         'rlnAngleRotPrior'
                       ,xmipp.MDL_ANGLE_TILT2:        'rlnAngleTiltPrior'
                       ,xmipp.MDL_ANGLE_PSI2:         'rlnAnglePsiPrior'                       
                       ,xmipp.MDL_SHIFT_X2 : 'rlnOriginXPrior' 
                       ,xmipp.MDL_SHIFT_Y2 : 'rlnOriginYPrior' 
                       ,xmipp.MDL_ZSCORE: 'rlnParticleSelectZScore' 
                       ,xmipp.MDL_COUNT2: 'rlnNrOfFrames'
                       # Not the best labels, but just to grab some
                       ,xmipp.MDL_CLASSIFICATION_DPR_05: 'rlnClassPriorOffsetX'
                       ,xmipp.MDL_CLASSIFICATION_FRC_05: 'rlnClassPriorOffsetY'
                       #,xmipp.MDL_CLASSIFICATION_INTRACLASS_DISTANCE: ''
                       }