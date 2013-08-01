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
                       ,MDL_CTF_DEFOCUSU:      'rlnDefocusU'
                       ,MDL_CTF_DEFOCUSV:      'rlnDefocusV'
                       ,MDL_CTF_DEFOCUS_ANGLE: 'rlnDefocusAngle'
                       ,MDL_CTF_VOLTAGE:       'rlnVoltage'
                       ,MDL_CTF_CS:            'rlnSphericalAberration'
                       ,MDL_CTF_Q0:            'rlnAmplitudeContrast'
                       ,MDL_IMAGE:             'rlnImageName'
                       ,MDL_LL:                'rlnLogLikeliContribution'
                       ,MDL_MICROGRAPH:        'rlnMicrographName'
                       ,MDL_AVGPMAX:           'rlnAveragePmax'
                       ,MDL_REF3D:             'rlnClassNumber'
                       ,MDL_RESOLUTION_FREQREAL:'rlnAngstromResolution'
                       ,MDL_RESOLUTION_FRC:     'rlnGoldStandardFsc'
                       ,MDL_RESOLUTION_FREQ:    'rlnResolution'
                       ,MDL_RESOLUTION_SSNR:    'rlnSsnrMap'
                       ,MDL_SCALE:              'rlnMagnificationCorrection'
                       ,MDL_SHIFT_X:            'rlnOriginX'
                       ,MDL_SHIFT_Y:            'rlnOriginY'
                       ,MDL_WEIGHT:             'rlnNrOfSignificantSamples'
                       }

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
    
def addRelionLabels():
    '''Add relion labels as aliases
    '''
    from xmipp import addLabelAlias
    for k, v in XMIPP_RELION_LABELS.iteritems():
        addLabelAlias(k,v)
        
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
    
    
