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

import time
from itertools import izip

from pyworkflow.protocol import Protocol
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from pyworkflow.em.data import (SetOfMicrographs, SetOfCoordinates,
                                SetOfParticles, SetOfImages,
                                SetOfClasses2D, SetOfClasses3D, SetOfClassesVol,
                                SetOfVolumes, SetOfCTF, SetOfMovies, SetOfFSCs,
                                SetOfMovieParticles, SetOfAverages,
                                SetOfNormalModes, SetOfPDBs)
from pyworkflow.em.constants import (RELATION_SOURCE, RELATION_TRANSFORM,
                                     RELATION_CTF)
from pyworkflow.em.data_tiltpairs import (SetOfAngles, CoordinatesTiltPair,
                                          TiltPair)
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
        return self.__createSet(SetOfMicrographs,
                                'micrographs%s.sqlite', suffix)
    
    def _createSetOfCoordinates(self, micSet, suffix=''):
        coordSet = self.__createSet(SetOfCoordinates,
                                    'coordinates%s.sqlite', suffix)
        coordSet.setMicrographs(micSet)       
        return coordSet

    def _createCoordinatesTiltPair(self, micTiltPairs, uCoords, tCoords,
                                   angles, suffix):
        coordTiltPairs = self.__createSet(CoordinatesTiltPair,
                                          'coordinates_pairs%s.sqlite', suffix)
        coordTiltPairs.setUntilted(uCoords)
        coordTiltPairs.setTilted(tCoords)
        coordTiltPairs.setAngles(angles)
        coordTiltPairs.setMicsPair(micTiltPairs)

        for coordU, coordT in izip(uCoords, tCoords):
            coordTiltPairs.append(TiltPair(coordU, coordT))

        return coordTiltPairs
    
    def _createSetFromName(self, className, suffix=''):
        """ Create a new set from the given className. """
        _createSetFunc = getattr(self, '_create%s' % className)
        return _createSetFunc(suffix)
        
    def _createSetOfImages(self, suffix=''):
        return self.__createSet(SetOfImages, 'images%s.sqlite', suffix)

    def _createSetOfParticles(self, suffix=''):
        return self.__createSet(SetOfParticles, 'particles%s.sqlite', suffix)

    def _createSetOfAverages(self, suffix=''):
        return self.__createSet(SetOfAverages, 'averages%s.sqlite', suffix)
        
    def _createSetOfMovieParticles(self, suffix=''):
        return self.__createSet(SetOfMovieParticles,
                                'movie_particles%s.sqlite', suffix)
    
    def _createSetOfClasses2D(self, imgSet, suffix=''):
        classes = self.__createSet(SetOfClasses2D,
                                   'classes2D%s.sqlite', suffix)
        classes.setImages(imgSet)
        return classes
    
    def _createSetOfClasses3D(self, imgSet, suffix=''):
        classes =  self.__createSet(SetOfClasses3D,
                                    'classes3D%s.sqlite', suffix)
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
        return self.__createSet(SetOfAngles,
                                'tiltpairs_angles%s.sqlite', suffix)

    def _createSetOfFSCs(self, suffix=''):
        return self.__createSet(SetOfFSCs, 'fscs%s.sqlite', suffix)
       
    def _createSetOfPDBs(self, suffix=''):
        return self.__createSet(SetOfPDBs, 'pdbs%s.sqlite', suffix)

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
        
    def _validateDim(self, obj1, obj2, errors,
                     label1='Input 1', label2='Input 2'):
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
        
    # ------ Methods for Streaming picking --------------

    def _defineStreamingParams(self, form):
        """ This function can be called during the _defineParams method
        of some protocols that support stream processing.
        It will add an Streaming section together with the following
        params:
            streamingSleepOnWait: Some streaming protocols are quite fast,
                so, checking input/ouput updates creates an IO overhead.
                This params allows them to sleep (without consuming resources)
                to wait for new work to be done.
            streamingBatchSize: For some programs it is more efficient to process
                many items at once and not one by one. So this parameter will
                allow to group a number of items to be processed in the same
                protocol step. This can also reduce some IO overhead and spawing
                new OS processes.
        """
        form.addSection("Streaming")
        form.addParam("streamingWarning", params.LabelParam, important=True,
                      label="The following params are related to how "
                            "streaming is done in Scipion.")
        form.addParam("streamingSleepOnWait", params.IntParam, default=0,
                      label="Sleep when waiting (secs)",
                      help="If you specify a value greater than zero, "
                           "it will be the number of seconds that the "
                           "protocol will sleep when waiting for new "
                           "input data in streaming mode. ")
        form.addParam("streamingBatchSize", params.IntParam, default=1,
                      label="Batch size",
                      help="This value allows to group several items to be "
                           "processed inside the same protocol step. You can "
                           "use the following values: \n"
                           "*1*    The default behavior, the items will be "
                           "processed one by one.\n"
                           "*0*    Put in the same step all the items "
                           "available. If the sleep time is short, it could be "
                           "practically the same of one by one. If not, you "
                           "could have steps with more items. If the steps will "
                           "be executed in parallel, it is better not to use "
                           "this option.\n"
                           "*>1*   The number of items that will be grouped into "
                           "a step.")

    def _getStreamingSleepOnWait(self):
        return self.getAttributeValue('streamingSleepOnWait', 0)

    def _getStreamingBatchSize(self):
        return self.getAttributeValue('streamingBatchSize', 1)

    def _streamingSleepOnWait(self):
        """ This method should be used by protocols that want to sleep
        when there is not more work to do.
        """
        sleepOnWait = self._getStreamingSleepOnWait()
        if sleepOnWait > 0:
            self.info("Not much work to do now, sleeping %s seconds."
                      % sleepOnWait)
            time.sleep(sleepOnWait)


    def _insertNewMics(self, inputMics, getMicKeyFunc,
                       insertStepFunc, insertStepListFunc, *args):
        """ Insert steps of new micrographs taking into account the batch size.
        It is assumed that a self.micDict exists mapping between micKey and mic.
        It is also assumed that self.streamClosed is defined...with True value
        if the input stream is closed.
        This function can be used from several base protocols that support
        streaming and batch:

        - Prot
        - ProtParticlePickingAuto
        - ProtExtractParticles
        Params:
            inputMics: the input micrographs to be inserted into steps
            getMicKeyFunc: function to get the key of a micrograph
                (usually mic.getMicName()
            insertStepFunc: function used to insert a single step
            insertStepListFunc: function used to insert many steps.
            *args: argument list to be passed to step functions
        Returns:
            The list of step Ids that can be used as dependencies.
        """
        deps = []
        insertedMics = inputMics

        # Despite this function only should insert new micrographs
        # let's double check that they are not inserted already
        micList = [mic for mic in inputMics
                   if getMicKeyFunc(mic) not in self.micDict]

        def _insertSubset(micSubset):
            stepId = insertStepListFunc(micSubset, self.initialIds, *args)
            deps.append(stepId)

        # Now handle the steps depending on the streaming batch size
        batchSize = self._getStreamingBatchSize()

        if batchSize == 1: # This is one by one, as before the batch size
            for mic in micList:
                stepId = insertStepFunc(mic, self.initialIds, *args)
                deps.append(stepId)
        elif batchSize == 0:  # Greedy, take all available ones
            _insertSubset(micList)
        else:  # batchSize > 0, insert only batches of this size
            n = len(inputMics)
            d = n / batchSize # number of batches to insert
            nd = d * batchSize
            for i in range(d):
                _insertSubset(micList[i*batchSize:(i+1)*batchSize])

            if n > nd and self.streamClosed:  # insert last ones
                _insertSubset(micList[nd:])
            else:
                insertedMics = micList[:nd]

        for mic in insertedMics:
            self.micDict[getMicKeyFunc(mic)] = mic

        return deps
