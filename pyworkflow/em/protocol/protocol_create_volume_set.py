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

import pyworkflow.protocol.params as params
from protocol_3d import Prot3D
from pyworkflow.em import Volume, SetOfVolumes



class ProtCreateVolumeSet(Prot3D):
    """ 
    Group output volumes into a new set.
    
    Some protocols receive as input a set of volumes.
    So this protocols is useful to create a new set
    from output volumes from other run.
    Input can be either volumes or set of volumes.
     """
    _label = 'create volume set'
    
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolumes', params.MultiPointerParam, 
                      pointerClass='SetOfVolumes,Volume',  
                      label="Input volume(s)", important=True, 
                      help='Select one or more volumes (Volume or SetOfVolumes)\n'
                           'to be grouped into a new SetOfVolumes')
        
    #--------------------------- INSERT steps functions --------------------------------------------    

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions --------------------------------------------
    
    def createOutputStep(self):
        volSet = self._createSetOfVolumes()
        oldSampling = None
        for vol in self._iterInputVolumes():
            samplingRate = vol.getSamplingRate()
            if (oldSampling is not None and
                oldSampling != samplingRate):
                self.warning('All the volumes in the output set ')
            # Only copy basic attributes 
            outVol = Volume()
            outVol.setObjId(None)
            outVol.setLocation(vol.getLocation())
            outVol.setObjLabel(vol.getObjLabel())
            outVol.setObjComment(vol.getObjComment())
            volSet.append(outVol)
            oldSampling = samplingRate
            
        volSet.setSamplingRate(samplingRate)

        self._defineOutputs(outputVolumes=volSet)
        for pointer in self.inputVolumes:
            self._defineSourceRelation(pointer, volSet)

    #--------------------------- INFO functions --------------------------------------------
    
    def _validate(self):
        errors = []
        for pointer in self.inputVolumes:
            if pointer.pointsNone():
                errors.append('Invalid input, pointer: %s' % pointer.getObjValue())
                errors.append('              extended: %s' % pointer.getExtended())
        return errors    
    
    def _summary(self):
        summary = []
        nVols = self._getNumberOfInputs()
            
        if nVols > 0:
            summary.append("Volumes grouped: *%d* " % nVols)
        else:
            summary.append("No volumes grouped.")
                
        return summary
        
    #--------------------------- UTILS functions --------------------------------------------
    def _iterInputVolumes(self):
        """ Iterate over all the input volumes. """
        for pointer in self.inputVolumes:
            item = pointer.get()
            if item is None:
                break
            itemId = item.getObjId()
            if isinstance(item, Volume):
                item.outputName = self._getExtraPath('output_vol%06d.vol' % itemId)
                yield item
            elif isinstance(item, SetOfVolumes):
                for vol in item:
                    vol.outputName = self._getExtraPath('output_vol%06d_%03d.vol' % (itemId, vol.getObjId()))
                    yield vol
                    
    def _getNumberOfInputs(self):
        """ Return the total number of input volumes. """
        nVols = 0
        for _ in self._iterInputVolumes():
            nVols += 1    
            
        return nVols    
