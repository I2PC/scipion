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
In this module are protocol base classes related to EM.
Them should be sub-classes in the different sub-packages from
each EM-software packages.
"""


from pyworkflow.protocol import Protocol
from pyworkflow.object import Set
from pyworkflow.em.data import (SetOfMicrographs, SetOfCoordinates, SetOfParticles,
                                SetOfClasses2D, SetOfClasses3D, SetOfClassesVol,
                                SetOfVolumes, SetOfCTF, SetOfMovies,
                                SetOfMovieParticles, SetOfAverages, SetOfNormalModes)
from pyworkflow.em.constants import RELATION_SOURCE, RELATION_TRANSFORM, RELATION_CTF
from pyworkflow.em.data_tiltpairs import SetOfAngles, CoordinatesTiltPair
from pyworkflow.utils.path import cleanPath
from pyworkflow.mapper.sqlite_db import SqliteDb



class EMProtocol(Protocol):
    """ Base class to all EM protocols.
    It will contains some common functionalities. 
    """
    _base = True
    
    def __init__(self, **kwargs):
        Protocol.__init__(self, **kwargs)
        
    def __createSet(self, SetClass, template, suffix):
        """ Create a set and set the filename using the suffix. 
        If the file exists, it will be delete. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        cleanPath(setFn)
        
        SqliteDb.closeConnection(setFn)        
        setObj = SetClass(filename=setFn)
        return setObj
    
    def _createSetOfMicrographs(self, suffix=''):
        return self.__createSet(SetOfMicrographs, 'micrographs%s.sqlite', suffix)
    
    def _createSetOfCoordinates(self, micSet, suffix=''):
        coordSet = self.__createSet(SetOfCoordinates, 'coordinates%s.sqlite', suffix)
        coordSet.setMicrographs(micSet)       
        return coordSet
    
    def _createSetFromName(self, className, suffix=''):
        """ Create a new set from the given className. """
        _createSetFunc = getattr(self, '_create%s' % className)
        return _createSetFunc(suffix)
        
    def _createSetOfParticles(self, suffix=''):
        return self.__createSet(SetOfParticles, 'particles%s.sqlite', suffix)

    def _createSetOfAverages(self, suffix=''):
        return self.__createSet(SetOfAverages, 'averages%s.sqlite', suffix)
        
    def _createSetOfMovieParticles(self, suffix=''):
        return self.__createSet(SetOfMovieParticles, 'movie_particles%s.sqlite', suffix)
    
    def _createSetOfClasses2D(self, imgSet, suffix=''):
        classes = self.__createSet(SetOfClasses2D, 'classes2D%s.sqlite', suffix)
        classes.setImages(imgSet)
        return classes
    
    def _createSetOfClasses3D(self, imgSet, suffix=''):
        classes =  self.__createSet(SetOfClasses3D, 'classes3D%s.sqlite', suffix)
        classes.setImages(imgSet)
        return classes
    
    def _createSetOfClassesVol(self, suffix=''):
        return self.__createSet(SetOfClassesVol, 'classesVol%s.sqlite', suffix)
    
    def _createSetOfVolumes(self, suffix=''):
        return self.__createSet(SetOfVolumes, 'volumes%s.sqlite', suffix)
    
    def _createSetOfCTF(self, suffix=''):
        return self.__createSet(SetOfCTF, 'ctfs%s.sqlite', suffix)
    
    def _createSetOfNormalModes(self, suffix=''):
        return self.__createSet(SetOfNormalModes, 'modes%s.sqlite', suffix)
    
    def _createSetOfMovies(self, suffix=''):
        return self.__createSet(SetOfMovies, 'movies%s.sqlite', suffix)
    
    def _createSetOfAngles(self, suffix=''):
        return self.__createSet(SetOfAngles, 'tiltpairs_angles%s.sqlite', suffix)
       
    def _defineSourceRelation(self, srcObj, dstObj):
        """ Add a DATASOURCE relation between srcObj and dstObj """
        self._defineRelation(RELATION_SOURCE, srcObj, dstObj)
    
    def _defineTransformRelation(self, srcObj, dstObj):
        self._defineRelation(RELATION_TRANSFORM, srcObj, dstObj)
        # A transform relation allways implies a source relation
        self._defineSourceRelation(srcObj, dstObj)
        
    def _defineCtfRelation(self, micsObj, ctfsObj):
        self._defineRelation(RELATION_CTF, micsObj, ctfsObj)
        # A ctf relation allways implies a source relation
        self._defineSourceRelation(micsObj, ctfsObj)
    
    def _insertChild(self, key, child):
        if isinstance(child, Set):
            child.write()
        Protocol._insertChild(self, key, child)
        
    def _validateDim(self, obj1, obj2, errors, label1='Input 1', label2='Input 2'):
        """ Validate that obj1 and obj2 has the same dimensions.
        Params:
            obj1, obj2: input objects that can be Images or SetOfImages subclasses.
            errors: an error list where to append the error if not same dimensions
            label1, label2: labels for input objects used for the error message.
        """
        if obj1 is not None and obj2 is not None:
            d1 = obj1.getXDim()
            d2 = obj2.getXDim()

            if d1 is None:
                errors.append("Can not get dimensions from %s." % label1)
            elif d2 is None:
                errors.append("Can not get dimensions from %s." % label2)
            elif d1 != d2:
                msg = '%s and %s have not the same dimensions, \n' % (label1, label2)
                msg += 'which are %d and %d, respectively' % (d1, d2)
                errors.append(msg)
            
    def __str__(self):
        return self.getObjLabel()
    
    def allowsDelete(self, obj):
        if (isinstance(obj, SetOfCoordinates) or
            isinstance(obj, CoordinatesTiltPair)):
            return True
        return False
        
    
