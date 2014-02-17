#!/usr/bin/env xmipp_python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from xmipp import *
import os
from os.path import  exists

CTF_BASIC_LABELS = [
                     MDL_CTF_CA
                    ,MDL_CTF_CS
                    ,MDL_CTF_DEFOCUS_ANGLE 
                    ,MDL_CTF_DEFOCUSU 
                    ,MDL_CTF_DEFOCUSV 
                    ,MDL_CTF_K
                    ,MDL_CTF_Q0 
                    ,MDL_CTF_SAMPLING_RATE
                    ,MDL_CTF_VOLTAGE
                    ]

# Map from Xmipp labels to Relion labels names
XMIPP_RELION_LABELS = {
                        MDL_ANGLE_ROT:         'rlnAngleRot'
                       ,MDL_ANGLE_TILT:        'rlnAngleTilt'
                       ,MDL_AVG_CHANGES_ORIENTATIONS:'rlnChangesOptimalOrientations'
                       ,MDL_AVG_CHANGES_OFFSETS:     'rlnChangesOptimalOffsets'
                       ,MDL_AVG_CHANGES_CLASSES:     'rlnChangesOptimalClasses'
                       ,MDL_ANGLE_PSI:         'rlnAnglePsi'
                       ,MDL_CLASS_PERCENTAGE:  'rlnClassDistribution'
                       ,MDL_CTF_CA:            'rlnChromaticAberration'
                       ,MDL_CTF_CONVERGENCE_CONE: 'rlnConvergenceCone'
                       ,MDL_CTF_CS:            'rlnSphericalAberration'
                       ,MDL_CTF_DEFOCUSU:      'rlnDefocusU'
                       ,MDL_CTF_DEFOCUSV:      'rlnDefocusV'
                       ,MDL_CTF_DEFOCUS_ANGLE: 'rlnDefocusAngle'
                       ,MDL_CTF_ENERGY_LOSS:   'rlnEnergyLoss'
                       ,MDL_CTF_LENS_STABILITY:'rlnLensStability'
                       ,MDL_CTF_LONGITUDINAL_DISPLACEMENT:'rlnLongitudinalDisplacement'
                       ,MDL_CTF_TRANSVERSAL_DISPLACEMENT: 'rlnTransversalDisplacement'                       
                       ,MDL_CTF_Q0:            'rlnAmplitudeContrast'
                       ,MDL_CTF_SAMPLING_RATE: 'rlnDetectorPixelSize'
                       ,MDL_CTF_VOLTAGE:       'rlnVoltage'
                       ,MDL_DATATYPE:          'rlnDataType'
                       ,MDL_DEFGROUP:          'rlnGroupNumber'
                       ,MDL_ENABLED:           'rlnEnabled'
                       ,MDL_IMAGE:             'rlnImageName'
                       ,MDL_IMAGE_REF:         'rlnReferenceImage'
                       ,MDL_ITEM_ID:           'rlnImageId'
                       ,MDL_MAXCC:             'rlnLogLikelihood'
                       ,MDL_LL:                'rlnLogLikeliContribution'
                       ,MDL_MAGNIFICATION:     'rlnMagnification'
                       ,MDL_MICROGRAPH:        'rlnMicrographName'
                       ,MDL_AVGPMAX:           'rlnAveragePmax'
                       ,MDL_REF3D:             'rlnClassNumber'
                       ,MDL_RESOLUTION_FREQREAL:'rlnAngstromResolution'
                       ,MDL_RESOLUTION_FRC:     'rlnGoldStandardFsc'
                       ,MDL_RESOLUTION_FREQ:    'rlnResolution'
                       ,MDL_RESOLUTION_SSNR:    'rlnSsnrMap'
                       ,MDL_RANDOMSEED:         'rlnRandomSubset'
                       ,MDL_SAMPLINGRATE:       'rlnPixelSize'
                       #,MDL_SAMPLINGRATE_ORIGINAL: 'rlnPixelSize'
                       ,MDL_SCALE:              'rlnMagnificationCorrection'
                       ,MDL_ORIGIN_X:           'rlnOriginX'
                       ,MDL_ORIGIN_Y:           'rlnOriginY'
                       ,MDL_ORIGIN_Z:           'rlnOriginZ'
                       ,MDL_PMAX:               'rlnMaxValueProbDistribution'
                       ,MDL_SAMPLINGRATE_X:           'rlnSamplingRateX'
                       ,MDL_SAMPLINGRATE_Y:           'rlnSamplingRateY'
                       ,MDL_SAMPLINGRATE_Z:           'rlnSamplingRateZ'
                       ,MDL_XCOOR:              'rlnCoordinateX'
                       ,MDL_XSIZE:              'rlnImageSizeX'
                       ,MDL_YCOOR:              'rlnCoordinateY'
                       ,MDL_YSIZE:              'rlnImageSizeY'
                       ,MDL_WEIGHT:             'rlnNrOfSignificantSamples'
                       ,MDL_ZSIZE:              'rlnImageSizeZ'
                       }
XMIPP_RELION_LABELS_EXTRA = {
                       # Following labels have no correct matching, just to 
                       # pick one with the same datatype
                       MDL_BLOCK_NUMBER:       'rlnGroupNumber' # just one
                       ,MDL_COUNT:             'rlnGroupNrParticles' # just one
                       ,MDL_CTF_CRIT_FITTINGSCORE:   'rlnCtfFigureOfMerit' #just one
                       ,MDL_CTF_CRIT_NORMALITY:   'rlnNormCorrection' #just one
                       ,MDL_CTF_CRIT_PSDVARIANCE: 'rlnCtfValue'         #just one
                       ,MDL_CTF_CRIT_PSDCORRELATION90: 'rlnCtfBfactor'  #just one
                       ,MDL_CRYSTAL_CELLX : 'rlnReferenceDimensionality'
                       ,MDL_CRYSTAL_CELLY : 'rlnOriginalImageSize'
                       ,MDL_CRYSTAL_DISAPPEAR_THRE : 'rlnCurrentResolution'
                       ,MDL_PICKING_PARTICLE_SIZE : 'rlnCurrentImageSize' #int
                       ,MDL_PICKING_AUTOPICKPERCENT : 'rlnPaddingFactor' #int
                       ,MDL_PICKING_TEMPLATES : 'rlnFourierSpaceInterpolator' #int
                       ,MDL_COLOR : 'rlnMinRadiusNnInterpolation' #int
                       ,MDL_DM3_NUMBER_TYPE : 'rlnNrClasses' #int
                       ,MDL_DM3_SIZE : 'rlnNrGroups' #int
                       ,MDL_NMA_COLLECTIVITY : 'rlnTau2FudgeFactor' #double
                       ,MDL_NMA_SCORE : 'rlnNormCorrectionAverage' #double
                       ,MDL_SIGMAOFFSET : 'rlnSigmaOffsets' #double
                       ,MDL_ANGLE_Y2 : 'rlnOrientationalPriorMode' #int
                       ,MDL_ANGLE_ROT2 : 'rlnSigmaPriorRotAngle' #double
                       ,MDL_ANGLE_TILT2 : 'rlnSigmaPriorTiltAngle' #double
                       ,MDL_ANGLE_PSI2 : 'rlnSigmaPriorPsiAngle' #double

                       ,MDL_MLF_WIENER: 'rlnOrientationDistribution' #double
                       ,MDL_IDX: 'rlnSpectralIndex' #int
                       ,MDL_MLF_NOISE: 'rlnSigma2Noise' #double
                       ,MDL_DM3_TAGNAME: 'rlnGroupName'  #string
                       ,MDL_MLF_SIGNAL: 'rlnGroupScaleCorrection' #double
                       
                       ,MDL_ZSCORE_SHAPE1: 'rlnAccuracyRotations'
                       ,MDL_ZSCORE_SHAPE2: 'rlnAccuracyTranslations'
                       ,MDL_ZSCORE_SNR1: 'rlnReferenceSigma2'
                       ,MDL_ZSCORE_SNR2: 'rlnReferenceTau2'
                       ,MDL_ZSCORE_RESCOV: 'rlnSpectralOrientabilityContribution'


                       }
# from data.star
#WARNING: Ignoring unknown column: rlnMaxValueProbDistribution

def convertCtfparam(oldCtf):
    '''Convert the old format (Xmipp2.4) of the CTF 
    and return a new MetaData'''
    oldLabelsName = ['sampling_rate', 'voltage', 'defocusU', 'defocusV', 'azimuthal_angle', 
                     'spherical_aberration', 'chromatic_aberration', 'Q0', 'K']                   
    conversionDict = dict(zip(oldLabelsName, CTF_BASIC_LABELS))
    
    f = open(oldCtf)
    md = MetaData()
    md.setColumnFormat(False)
    objId = md.addObject()
    
    for line in f:
        line = line.strip()
        if len(line) and not line.startswith('#'):
            parts = line.split('=')
            old_key = parts[0].strip()
            value = float(parts[1].strip())
            label = conversionDict[old_key]
            md.setValue(label, value, objId)

    f.close()
  
    #Change sign of defocusU, defocusV and Q0
    labels = [label2Str(l) for l in [MDL_CTF_DEFOCUSU, MDL_CTF_DEFOCUSV, MDL_CTF_Q0]]
    exps = ["%s=%s*-1" % (l, l) for l in labels]
    expression = ','.join(exps)
    
    md.operate(expression)
    return md
    
def convertBox(boxFile, posFile, ysize, family='DefaultFamily', particleSize=None):
    ''' Convert the .box files with EMAN coordinates for one micrographs
    to the .pos metadata file on Xmipp, change the coordinate from left-upper 
    corner to center. As parameter the family could be provided and the 
    particle size(if different from EMAN. If ysize is provide, the y-coordinate is inverted'''
    f = open(boxFile)
    
    md = MetaData()
    md.readPlain(boxFile, 'xcoor ycoor')
    eman1 = False
    
    if particleSize is None: #Find the EMAN particle size
        for line in f:
            if len(line.strip()):
                parts = line.strip().split()
                if len(parts) == 5:
                    eman1 = True # the y should be inverted
                particleSize = int(parts[2])                    
                break
                
    half = particleSize / 2
    # Move coordinates to the center and also invert the y-coordinate
    if eman1: #invert y
        md.operate('xcoor=xcoor+%(half)d,ycoor=%(ysize)d-(ycoor+%(half)d)' % locals())
    else:
        md.operate('xcoor=xcoor+%(half)d,ycoor=ycoor+%(half)d' % locals())
    
    mdFamily = MetaData()
    objId = mdFamily.addObject()
    mdFamily.setValue(MDL_PICKING_FAMILY, family, objId)
    mdFamily.setValue(MDL_PICKING_MICROGRAPH_FAMILY_STATE, 'Manual', objId)
    mdFamily.write('families@%s' % posFile)
    md.write("%s@%s" % (family, posFile), MD_APPEND)
    
def renameMdLabels(inputMd, outputMd, labelsDict):
    '''Change the labels' name on inputMd and write as outputMd
    The change will be made using the labelsDict, changing
    key by value 
    If dictString is False, the keys in the dictionary are Xmipp MDL_* labels
    '''
    fIn = open(inputMd)
    fOut = open(outputMd, 'w+')
    for l in fIn:        
        label = l.split()[0].strip()[1:] # remove starting _ character
        if label in labelsDict:
            l = l.replace(label, labelsDict[label])
        fOut.write(l)
    fIn.close()
    fOut.close()
    
def removeMdHeader(inputMd, outputFile):
    ''' Remove header from md and writing only data lines'''
    fIn = open(inputMd)
    fOut = open(outputFile, 'w+')
    for l in fIn:
        if not (l.startswith('data') or l.startswith('loop') or l.startswith('_')):
            fOut.write(l)
    fIn.close()
    fOut.close()
    
def exportSpiderParticles(inputMd, outputFile):
    ''' Read an Xmipp3.0 .pos metadata containing particles
    from one micrograph and export as expected in Spider'''
    md = MetaData(inputMd)
    fOut = open(outputFile, 'w+')
    
    i = 0
    for objId in md:
        i = i + 1
        x = float(md.getValue(MDL_XCOOR, objId))
        y = float(md.getValue(MDL_YCOOR, objId))
        print >> fOut, " %(i)04d 6        %(i)04d  %(x)f  %(y)f  %(x)f  %(y)f           1" % locals()
    fOut.close()
    
def readMdFromSpider(inputDoc, labelsStr):
    """ Read an Spider docfile as a metadata.
    Params:
        inputDoc: Spider docfile path
        labelsStr: an string containing the names of labels to use.
    """
    labels = [str2Label(l) for l in labelsStr.split()]
    fDoc = open(inputDoc, 'r')
    md = MetaData()
    
    for line in fDoc:
        # Read only non comment lines
        if not line.strip().startswith(';'):
            values = line.split()
            objId = md.addObject()
            md.setValue(MDL_ITEM_ID, long(values[0]), objId) # Read Spider key as item_id
            for l, v in zip(labels, values[2:]): # Exclude key and number of colums
                md.setValue(l, float(v), objId)
    fDoc.close()
    
    return md
    
def exportEman2Boxes(inputMd, outputFile, dim):
    ''' Read an Xmipp3.0 .pos metadata containing particles
    from one micrograph and export as expected in Eman2'''
    md = MetaData(inputMd)
    fOut = open(outputFile, 'w+')
    half = dim / 2
    
    for objId in md:
        x = md.getValue(MDL_XCOOR, objId) - half
        y = md.getValue(MDL_YCOOR, objId) - half
        print >> fOut, "%(x)d  %(y)d  %(dim)d %(dim)d" % locals()
    fOut.close()    
        
def exportMdToRelion(md, outputRelion):
    """ This function will receive a Xmipp metadata and will
    convert to the one expected by Relion.    
    All labels not recognized by Relion will be dropped.
    Params:
     md: input xmipp metadata.
     outputRelion: output filename to store the Relion star file.
    """
    for label in md.getActiveLabels():
        if not label in XMIPP_RELION_LABELS:
            md.removeLabel(label)
    tmpFile = outputRelion + '.tmp'
    md.write(tmpFile)
    # Create a dict with the names
    d = {}
    for k, v in XMIPP_RELION_LABELS.iteritems():
        d[label2Str(k)] = v
        
    #print "dict: ", d
        
    renameMdLabels(tmpFile, outputRelion, d)
    from protlib_filesystem import deleteFile
    deleteFile(None, tmpFile)
    
def addRelionLabels(replace=False, extended=False):
    '''Add relion labels as aliases
    '''
    from xmipp import addLabelAlias
    for k, v in XMIPP_RELION_LABELS.iteritems():
        addLabelAlias(k, v, replace)
    if extended:
        for k, v in XMIPP_RELION_LABELS_EXTRA.iteritems():    
            addLabelAlias(k, v, replace)        
        
def relionLabelString():
    """ create an string that can be used for XMIPP_EXTRA_ALIASES
    for adding the labels of Relion.
    """
    from xmipp import label2Str
    pairs = []
    for k, v in XMIPP_RELION_LABELS.iteritems():
        pairs.append('%s=%s' % (label2Str(k), v))
    for k, v in XMIPP_RELION_LABELS_EXTRA.iteritems():
        pairs.append('%s=%s' % (label2Str(k), v))        
    return ';'.join(pairs)
        
def exportReliontoMetadataFile(inputRelion,outputXmipp):
    """ This function will receive a relion file and will
    convert to the one xmipp metadata file.    
    """
    if str2Label('rlnImageName')== -1:
       addRelionLabels()
    #add xmipp header
    line1="""# XMIPP_STAR_1 * 
# 
"""
    tmpFile = inputRelion + '.tmp'
    fOut = open(tmpFile,'w')
    fOut.write(line1)
    #copy and paste relion file
    fIn = open(inputRelion,'r')
    fOut.write(fIn.read())
    fOut.flush()
    fOut.close()
    #addRelionLabels MUST BE CALLED FROM CODE 
    blocklist = getBlocksInMetaDataFile(tmpFile)
    if exists (outputXmipp):
        os.remove(outputXmipp)
    for block in blocklist:
        md = MetaData(block + '@'+tmpFile)
        if len(md.getActiveLabels())!=0:
            md.write(block +'@'+ outputXmipp,MD_APPEND)
        
    from protlib_filesystem import deleteFile
    deleteFile(None, tmpFile,False)
    
    
