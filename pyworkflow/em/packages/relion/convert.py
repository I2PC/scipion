# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This module contains converter functions that will serve to:
1. Write from base classes to Relion specific files
2. Read from Relion files to base classes
"""

import os
from pyworkflow.em.constants import NO_INDEX
from pyworkflow.object import String
#Following imports should not be used when a more generic library is used to deal with metadatas
import xmipp
from pyworkflow.em.packages.xmipp3 import XmippMdRow, getLabelPythonType, LABEL_TYPES

XMIPP_RELION_LABELS = {
                        xmipp.MDL_ANGLE_ROT:         'rlnAngleRot'
                       ,xmipp.MDL_ANGLE_TILT:        'rlnAngleTilt'
                       ,xmipp.MDL_AVG_CHANGES_ORIENTATIONS:'rlnChangesOptimalOrientations'
                       ,xmipp.MDL_AVG_CHANGES_OFFSETS:     'rlnChangesOptimalOffsets'
                       ,xmipp.MDL_AVG_CHANGES_CLASSES:     'rlnChangesOptimalClasses'
                       ,xmipp.MDL_ANGLE_PSI:         'rlnAnglePsi'
                       ,xmipp.MDL_CTF_DEFOCUSU:      'rlnDefocusU'
                       ,xmipp.MDL_CTF_DEFOCUSV:      'rlnDefocusV'
                       ,xmipp.MDL_CTF_DEFOCUS_ANGLE: 'rlnDefocusAngle'
                       ,xmipp.MDL_CTF_VOLTAGE:       'rlnVoltage'
                       ,xmipp.MDL_CTF_CS:            'rlnSphericalAberration'
                       ,xmipp.MDL_CTF_Q0:            'rlnAmplitudeContrast'
                       ,xmipp.MDL_IMAGE:             'rlnImageName'
                       ,xmipp.MDL_LL:                'rlnLogLikeliContribution'
                       ,xmipp.MDL_MICROGRAPH:        'rlnMicrographName'
                       ,xmipp.MDL_AVGPMAX:           'rlnAveragePmax'
                       ,xmipp.MDL_REF3D:             'rlnClassNumber'
                       ,xmipp.MDL_RESOLUTION_FREQREAL:'rlnAngstromResolution'
                       ,xmipp.MDL_RESOLUTION_FRC:     'rlnGoldStandardFsc'
                       ,xmipp.MDL_RESOLUTION_FREQ:    'rlnResolution'
                       ,xmipp.MDL_RESOLUTION_SSNR:    'rlnSsnrMap'
                       ,xmipp.MDL_SCALE:              'rlnMagnificationCorrection'
                       ,xmipp.MDL_SHIFT_X:            'rlnOriginX'
                       ,xmipp.MDL_SHIFT_Y:            'rlnOriginY'
                       ,xmipp.MDL_WEIGHT:             'rlnNrOfSignificantSamples'
                       }

#def getLabelPythonType(label):
#    """ From xmipp label to python variable type """
#    for k,v in XMIPP_RELION_LABELS.iteritems():
#        if v == label:
#            labelType = xmipp.labelType(k)
#            return LABEL_TYPES.get(labelType, str)

def objectToRow(obj, row, attrDict):
    """ This function will convert an EMObject into a RelionMdRow.
    Params:
        obj: the EMObject instance (input)
        row: the RelionMdRow instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Relion (values).
    """
    
    for attr, label in attrDict.iteritems():
        if hasattr(obj, attr):
            print attr
            valueType = getLabelPythonType(label)
            row.setValue(label, valueType(getattr(obj, attr).get()))
            
def _rowToObject(row, obj, attrDict):
    """ This function will convert from a XmippMdRow to an EMObject.
    Params:
        row: the RelionMdRow instance (input)
        obj: the EMObject instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).
    """
    for attr, label in attrDict.iteritems():
        if not hasattr(obj, attr):
            setattr(obj, attr, String()) #TODO: change string for the type of label
            
        getattr(obj, attr).set(row.getValue(label))
    
def rowToObject(md, objId, obj, attrDict):
    """ Same as rowToObject, but creating the row from md and objId. """
    row = XmippMdRow()
    row.readFromMd(md, objId)
    _rowToObject(row, obj, attrDict)
    
def imageToRow(img, imgRow):
    imgDict = { 
               "_filename": xmipp.MDL_IMAGE,
               }
    index, filename = img.getLocation()
    fn = locationToRelion(index, filename)
    #TODO: que pasa con los pars relativos al acquisition?
    objectToRow(img, imgRow, imgDict)
        
def ctfToRow(ctf, imgRow):
    ctfDict = {
               "defocusU": xmipp.MDL_CTF_DEFOCUSU,
               "defocusV": xmipp.MDL_CTF_DEFOCUSV,
               "defocusAngle": xmipp.MDL_CTF_DEFOCUS_ANGLE
               }
    objectToRow(ctf, imgRow, ctfDict)
    
def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    #TODO: Maybe we need to add more logic dependent of the format
    if index != NO_INDEX:
        return "%d@%s" % (index, filename)

def writeSetOfImages(imgSet, filename):
    """ This function will write a SetOfImages as Relion metadata.
    Params:
        imgSet: the SetOfImages instance.
        filename: the filename where to write the metadata.
    """   
    print "in writeSetOfImages"
    md = xmipp.MetaData()
    
    for img in imgSet:
        objId = md.addObject()
        imgRow = XmippMdRow()
        imageToRow(img, imgRow)
        if img.hasCTF():
            imgRow.setValue(xmipp.MDL_MICROGRAPH, img.getCTF().micFile)
            ctfToRow(img.getCTF(), imgRow)  
        imgRow.writeToMd(md, objId)
    print "antes de escribir el star"
    tmpFile = filename + '.tmp'
    md.write(tmpFile)
    # Create a dict with the names
    d = {}
    for k, v in XMIPP_RELION_LABELS.iteritems():
        d[xmipp.label2Str(k)] = v
    #Rename labels on Metadata to use relion labels
    renameMdLabels(tmpFile, filename, d)
    imgSet._relionStar = String(filename)
    
def createRelionInputImages(self, imgSet, imagesFn=None): 
    imgsStar = getattr(imgSet, '_relionStar', None)
    if imgsStar is None:
        imgsFn = self._getPath(imagesFn or 'input_images.star')
        writeSetOfImages(imgSet, imgsFn)
    else:
        imgsFn = imgsStar.get()
    return imgsFn

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