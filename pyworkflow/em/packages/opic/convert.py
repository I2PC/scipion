# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains converter functions that will serve to:
1. Write from base classes to Relion specific files
2. Read from Relion files to base classes
"""

import os
import sys
import re
from os.path import join, basename
from itertools import izip
import numpy
from collections import OrderedDict

from pyworkflow.object import ObjectWrap, String, Integer
from pyworkflow.utils import Environ
from pyworkflow.utils.path import (createLink, cleanPath, copyFile,
                                   replaceBaseExt, getExt, removeExt)
import pyworkflow.em as em
import pyworkflow.em.metadata as md

# This dictionary will be used to map
# between CTFModel properties and Xmipp labels
ACQUISITION_DICT = OrderedDict([ 
       ("_amplitudeContrast", md.RLN_CTF_Q0),
       ("_sphericalAberration", md.RLN_CTF_CS),
       ("_voltage", md.RLN_CTF_VOLTAGE),
        ("_magnification", md.RLN_CTF_MAGNIFICATION)
       ])

COOR_DICT = OrderedDict([
             ("_x", md.RLN_IMAGE_COORD_X),
             ("_y", md.RLN_IMAGE_COORD_Y)
             ])

CTF_DICT = OrderedDict([
       ("_defocusU", md.RLN_CTF_DEFOCUSU),
       ("_defocusV", md.RLN_CTF_DEFOCUSV),
       ("_defocusAngle", md.RLN_CTF_DEFOCUS_ANGLE)
       ])

CTF_PSD_DICT = OrderedDict([
       ("_psdFile", md.RLN_CTF_IMAGE)
       ])

CTF_EXTRA_LABELS = [   
    md.RLN_CTF_FOM,
    # In Relion the ctf also contains acquisition information
    md.RLN_CTF_Q0,
    md.RLN_CTF_CS,
    md.RLN_CTF_VOLTAGE,
    md.RLN_CTF_MAGNIFICATION,
    md.RLN_CTF_DETECTOR_PIXEL_SIZE
    ]

# Some extra labels to take into account the zscore
IMAGE_EXTRA_LABELS = [
    md.RLN_SELECT_PARTICLES_ZSCORE,
    md.RLN_IMAGE_FRAME_NR,
    ]

# ANGLES_DICT = OrderedDict([
#        ("_angleY", md.RLN_ANGLE_Y),
#        ("_angleY2", md.RLN_ANGLE_Y2),
#        ("_angleTilt", md.RLN_ANGLE_TILT)
#        ])
 
ALIGNMENT_DICT = OrderedDict([ 
       ("_rlnOriginX", md.RLN_ORIENT_ORIGIN_X),
       ("_rlnOriginY", md.RLN_ORIENT_ORIGIN_Y),
       ("_rlnOriginZ", md.RLN_ORIENT_ORIGIN_Z),
       ("_rlnAngleRot", md.RLN_ORIENT_ROT),
       ("_rlnAngleTilt", md.RLN_ORIENT_TILT),
       ("_rlnAnglePsi", md.RLN_ORIENT_PSI),
       ])


def getRelionEnviron():
    """ Setup the environment variables needed to launch Relion. """
    environ = Environ(os.environ)
    environ.update({
            'PATH': join(os.environ['RELION_HOME'], 'bin'),
            'LD_LIBRARY_PATH': join(os.environ['RELION_HOME'], 'lib') + ":" + join(os.environ['RELION_HOME'], 'lib64'),
            'SCIPION_MPI_FLAGS': os.environ.get('RELION_MPI_FLAGS', ''),
            }, position=Environ.BEGIN)
    return environ

def setEnviron():
    """ Setup the environment variables needed to import localrec classes. """
    os.environ.update({'PYTHONPATH': join(os.environ['LOCALREC_HOME'], 'lib')
                       })
    os.environ.update(getRelionEnviron())
    sys.path.append(os.path.join(os.environ["LOCALREC_HOME"], "lib"))


def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    if index != em.NO_INDEX:
        return "%06d@%s" % (index, filename)
    
    return filename


def objectToRow(obj, row, attrDict, extraLabels={}):
    """ This function will convert an EMObject into a XmippMdRow.
    Params:
        obj: the EMObject instance (input)
        row: the XmippMdRow instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).        
        extraLabels: a list with extra labels that could be included
            as _xmipp_labelName
    """
    if obj.isEnabled():
        enabled = True
    else:
        enabled = False
    row.setValue(md.RLN_IMAGE_ENABLED, enabled)
    
    for attr, label in attrDict.iteritems():
        if hasattr(obj, attr):
            valueType = md.label2Python(label)
            row.setValue(label, valueType(getattr(obj, attr).get()))

    attrLabels = attrDict.values()
    
    for label in extraLabels:
        attrName = '_' + md.label2Str(label)
        if label not in attrLabels and hasattr(obj, attrName):
            value = obj.getAttributeValue(attrName) 
            row.setValue(label, value)


def rowToObject(row, obj, attrDict, extraLabels={}):
    """ This function will convert from a XmippMdRow to an EMObject.
    Params:
        row: the XmippMdRow instance (input)
        obj: the EMObject instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as properties with the label name such as: _rlnSomeThing
    """
    obj.setEnabled(row.getValue(md.RLN_IMAGE_ENABLED, 1) > 0)
    
    for attr, label in attrDict.iteritems():
        value = row.getValue(label)
        if not hasattr(obj, attr):
            setattr(obj, attr, ObjectWrap(value))
        else:
            getattr(obj, attr).set(value)
        
    attrLabels = attrDict.values()
    
    for label in extraLabels:
        if label not in attrLabels and row.hasLabel(label):
            labelStr = md.label2Str(label)
            setattr(obj, '_' + labelStr, row.getValueAsObject(label))


def setRowId(mdRow, obj, label=md.RLN_IMAGE_ID):
    mdRow.setValue(label, long(obj.getObjId()))


def acquisitionToRow(acquisition, ctfRow):
    """ Set labels values from acquisition to md row. """
    objectToRow(acquisition, ctfRow, ACQUISITION_DICT)


def ctfModelToRow(ctfModel, ctfRow):
    """ Set labels values from ctfModel to md row. """
    objectToRow(ctfModel, ctfRow, CTF_DICT, extraLabels=CTF_EXTRA_LABELS)
    

def geometryFromMatrix(matrix, inverseTransform):
    from pyworkflow.em.transformations import translation_from_matrix, euler_from_matrix

    if inverseTransform:
        from numpy.linalg import inv
        matrix = inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    angles = -numpy.rad2deg(euler_from_matrix(matrix, axes='szyz'))
    return shifts, angles


def matrixFromGeometry(shifts, angles, inverseTransform):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """
    from pyworkflow.em.transformations import euler_matrix
    from numpy import deg2rad
    radAngles = -deg2rad(angles)
    
    M = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if inverseTransform:
        from numpy.linalg import inv
        M[:3, 3] = -shifts[:3]
        M = inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M


def alignmentToRow(alignment, alignmentRow):
    matrix = alignment.getMatrix()
    shifts, angles = geometryFromMatrix(matrix, True)

    alignmentRow.setValue(md.RLN_ORIENT_ORIGIN_X, shifts[0])
    alignmentRow.setValue(md.RLN_ORIENT_ORIGIN_Y, shifts[1])
    alignmentRow.setValue(md.RLN_ORIENT_ORIGIN_Z, shifts[2])
    alignmentRow.setValue(md.RLN_ORIENT_ROT,  angles[0])
    alignmentRow.setValue(md.RLN_ORIENT_TILT, angles[1])
    alignmentRow.setValue(md.RLN_ORIENT_PSI,  angles[2])


def rowToAlignment(subpart):
    inverseTransform = True
    
    alignment = em.Transform()
    angles = numpy.zeros(3)
    shifts = numpy.zeros(3)
    shifts[0] = subpart.rlnOriginX
    shifts[1] = subpart.rlnOriginY
    shifts[2] = subpart.rlnOriginZ or 0
    angles[0] = subpart.rlnAngleRot
    angles[1] = subpart.rlnAngleTilt
    angles[2] = subpart.rlnAnglePsi

    M = matrixFromGeometry(shifts, angles, inverseTransform)
    alignment.setMatrix(M)
    
    return alignment


def coordinateToRow(coord, coordRow, label, copyId=True):
    """ Set labels values from Coordinate coord to md row. """
    if copyId:
        setRowId(coordRow, coord)
    objectToRow(coord, coordRow, COOR_DICT)
    if coord.getMicId():
        coordRow.setValue(label, long(coord.getMicId()))


def imageToRow(img, imgRow, imgLabel, **kwargs):
    # Provide a hook to be used if something is needed to be 
    # done for special cases before converting image to row
    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, imgRow)
        
    setRowId(imgRow, img) # Set the id in the metadata as MDL_ITEM_ID
    index, fn = img.getLocation()
    # check if the is a file mapping
    filesDict = kwargs.get('filesDict', {})
    filename = filesDict.get(fn, fn)
     
    imgRow.setValue(imgLabel, locationToRelion(index, filename))

    if kwargs.get('writeCtf', True) and img.hasCTF():
        ctfModelToRow(img.getCTF(), imgRow)
        
    # alignment is mandatory at this point, it shoud be check
    # and detected defaults if not passed at readSetOf.. level
    
    alignmentToRow(img.getTransform(), imgRow)
                
    if kwargs.get('writeAcquisition', True) and img.hasAcquisition():
        acquisitionToRow(img.getAcquisition(), imgRow)
    
    # Write all extra labels to the row    
    objectToRow(img, imgRow, {}, extraLabels=IMAGE_EXTRA_LABELS)

    # Provide a hook to be used if something is needed to be 
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, imgRow)


def particleToRow(part, partItem, **kwargs):
    """ Set labels values from Particle to md row. """
    
    partRow = md.Row()
    coord = part.getCoordinate()
    if coord is not None:
        coordinateToRow(coord, partRow, md.RLN_MICROGRAPH_ID, copyId=False)
    if part.hasMicId():
        partRow.setValue(md.RLN_MICROGRAPH_ID, long(part.getMicId()))
        # If the row does not contains the micrgraphs name
        # use a fake micrograph name using id to relion
        # could at least group for CTF using that
        if not partRow.hasLabel(md.RLN_MICROGRAPH_NAME):
            partRow.setValue(md.RLN_MICROGRAPH_NAME, 'fake_micrograph_%06d.mrc' % part.getMicId())
    if part.hasAttribute('_rlnParticleId'):
        partRow.setValue(md.RLN_PARTICLE_ID, long(part._rlnParticleId.get()))
    imageToRow(part, partRow, md.RLN_IMAGE_NAME, **kwargs)
    
    for label, value in partRow._labelDict.iteritems():
        labelStr = md.label2Str(label)
        setattr(partItem, labelStr, value)


def rowToSubcoordinate(subpart, coord, particle):
    coord.setMicId(particle.getObjId())

    if particle.hasCTF():
        ctf = particle.getCTF()
        ctf.setDefocusU(subpart.rlnDefocusU)
        ctf.setDefocusV(subpart.rlnDefocusV)
    
    
    particle.setTransform(rowToAlignment(subpart))
    
    
    coord.setX(subpart.rlnCoordinateX)
    coord.setY(subpart.rlnCoordinateY)
    coord._subparticle = particle