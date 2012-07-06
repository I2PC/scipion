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

from xmipp import MDL_CTF_SAMPLING_RATE, MDL_CTF_VOLTAGE, MDL_CTF_DEFOCUSU, MDL_CTF_DEFOCUSV, \
MDL_CTF_DEFOCUS_ANGLE, MDL_CTF_CS, MDL_CTF_CA, MDL_CTF_Q0, MDL_CTF_K, label2Str, MetaData,\
MDL_XCOOR, MDL_YCOOR, MDL_PICKING_FAMILY, MDL_PICKING_MICROGRAPH_FAMILY_STATE

CTF_BASIC_LABELS = [MDL_CTF_SAMPLING_RATE, 
                    MDL_CTF_VOLTAGE, 
                    MDL_CTF_DEFOCUSU, 
                    MDL_CTF_DEFOCUSV, 
                    MDL_CTF_DEFOCUS_ANGLE, 
                    MDL_CTF_CS,
                    MDL_CTF_CA,
                    MDL_CTF_Q0, 
                    MDL_CTF_K]

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
    md.readPlain(boxFile, 'Xcoor Ycoor')
    
    if particleSize is None: #Find the EMAN particle size
        for line in f:
            if len(line.strip()):
                parts = line.strip().split()
                particleSize = int(parts[2])                    
                break
                
    half = particleSize / 2
    # Move coordinates to the center and also invert the y-coordinate
    md.operate('Xcoor=Xcoor+%(half)d,Ycoor=Ycoor+%(half)d' % locals())
    
    mdFamily = MetaData()
    objId = mdFamily.addObject()
    mdFamily.setValue(MDL_PICKING_FAMILY, family, objId)
    mdFamily.setValue(MDL_PICKING_MICROGRAPH_FAMILY_STATE, 'Manual', objId)
    mdFamily.write('families@%s' % posFile)
    md.writeBlock(posFile, family)
    
def renameMdLabels(inputMd, outputMd, labelsDict):
    '''Change the labels' name on inputMd and write as outputMd
    The change will be made using the labelsDict, changing
    key by value '''
    fIn = open(inputMd)
    fOut = open(outputMd, 'w+')
    for l in fIn:
        if l.strip() in labelsDict:
            l = l.strip()
            l = labelsDict.get(l, l) + '\n'
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
    