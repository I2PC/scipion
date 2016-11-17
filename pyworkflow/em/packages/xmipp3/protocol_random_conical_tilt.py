# **************************************************************************
# *
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
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

from os.path import exists

import xmipp
from pyworkflow.object import String
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtInitialVolume
from pyworkflow.em.data import Volume, SetOfParticles
from pyworkflow.em.constants import ALIGN_2D

from convert import getImageLocation, alignmentToRow
from xmipp3 import XmippMdRow


class XmippProtRCT(ProtInitialVolume):
    """Creates initial volumes by using a set of projections/classes
    from a tilted-pair picking process and using RCT algorithm. """
    _label = 'random conical tilt'
    
    def __init__(self, **args):
        ProtInitialVolume.__init__(self, **args)

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        # TODO: Input can be either a SetOfParticles or
        # a SetOfClasses2D
        form.addParam('inputParticlesTiltPair', params.PointerParam,
                      label="Input particles tilt pair",
                      important=True,
                      pointerClass='ParticlesTiltPair',
                      help='Select the input particles tilt pair file '
                           'that will be used.  file. This file is used '
                           'to associate each micrograph with its tilted '
                           'equivalent.')
        
        form.addParam('inputParticles', params.PointerParam,
                      label="Input classes", important=True,
                      pointerClass='SetOfParticles,SetOfClasses',
                      help='Select the input images or classes from '
                           'the project.')
        
        form.addSection(label='Alignment parameters')
        form.addParam('thinObject', params.BooleanParam,
                      default=False, label='Thin Object',
                      help='If the object is thin, then the tilted '
                           'projections can be stretched to match the '
                           'untilted projections')
                       
        form.addParam('maxShift', params.IntParam,
                      default="10", expertLevel=LEVEL_ADVANCED,
                      label="Maximum allowed shift for tilted particles (pixels)", 
                      help='Particles that shift more will be discarded. '
                           'A value larger than the image size will not '
                           'discard any particle.')
        
        form.addParam('skipTranslation', params.BooleanParam,
                      default=False, expertLevel=LEVEL_ADVANCED,
                      label='Skip tilted translation alignment', 
                      help='If the tilted image quality is very low, '
                           'then this alignment might result in poor '
                           'estimates.')

        form.addSection(label='Reconstruction')
        form.addParam('additionalParams', params.StringParam,
                      default="-n 5 -l 0.01", expertLevel=LEVEL_ADVANCED,
                      label='Additional reconstruction parameters', 
                      help='See: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Reconstruct_art_v31')
        
        form.addParam('doFilter', params.BooleanParam, default=True,
                      label='Filter reconstructed volumes?', 
                      help='Filtering may be useful to remove noise, '
                           'especially when few particles '
                           'contribute to the reconstruction.')
        
        form.addParam('resoLowPassFilter', params.FloatParam, default=0.2,
                      label='Resolution of the low-pass filter (dig.freq)',
                      help='Resolution of the low-pass filter (dig.freq)')        

        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- STEPS functions -------------------------------

    def _insertAllSteps(self):
        self.inputSet = self.inputParticles.get()
        self.rctClassesFn = self._getExtraPath('rct_classes.xmd')
        
        firstStepId = self._insertFunctionStep('createRctImagesStep')
        reconSteps = []
                
        if isinstance(self.inputSet, SetOfParticles):
            reconStep = self._reconstructImages(self.inputSet, firstStepId)
            reconSteps.append(reconStep)
        else:
            for class2D in self.inputSet:
                reconStep = self._reconstructImages(class2D, firstStepId)
                reconSteps.append(reconStep)
                
        self._insertFunctionStep('createOutputStep', prerequisites=reconSteps)      
        
    def createRctImagesStep(self):
        """ Function to create the Xmipp metadata needed to run the
        Xmipp protocol """
        if isinstance(self.inputSet, SetOfParticles):
            self._appendRctImages(self.inputSet)
        else:
            for class2D in self.inputSet:
                self._appendRctImages(class2D)

    def _appendRctImages(self, particles):
        blockMd = "class%06d_images@%s" % (particles.getObjId(),
                                           self.rctClassesFn)
        classMd = xmipp.MetaData()

        partPairs = self.inputParticlesTiltPair.get()
        uImages = partPairs.getUntilted()
        tImages = partPairs.getTilted()
        sangles = partPairs.getCoordsPair().getAngles()

        micPairs = partPairs.getCoordsPair().getMicsPair()
        uMics = micPairs.getUntilted()
        tMics = micPairs.getTilted()
        
        scaleFactor = uImages.getSamplingRate() / particles.getSamplingRate()
        
        for img in particles:
            imgId = img.getObjId()
                       
            uImg = uImages[imgId]
            tImg = tImages[imgId]
            
            if uImg is None or tImg is None:
                print (">>> Warning, for id %d, tilted or untilted particle "
                       "was not found. Ignored." % imgId)
            else:
                objId = classMd.addObject()
                pairRow = XmippMdRow()
                pairRow.setValue(xmipp.MDL_IMAGE, getImageLocation(uImg))
                uCoord = uImg.getCoordinate()
                micId = uCoord.getMicId()
                uMic = uMics[micId]
                angles = sangles[micId]
                pairRow.setValue(xmipp.MDL_MICROGRAPH, uMic.getFileName())
                pairRow.setValue(xmipp.MDL_XCOOR, uCoord.getX())
                pairRow.setValue(xmipp.MDL_YCOOR, uCoord.getY())
                pairRow.setValue(xmipp.MDL_ENABLED, 1)
                pairRow.setValue(xmipp.MDL_ITEM_ID, long(imgId))
                pairRow.setValue(xmipp.MDL_REF, 1)
    
                alignment = img.getTransform()
    
                # Scale alignment by scaleFactor
                alignment.scale(scaleFactor)
                alignmentToRow(alignment, pairRow, alignType=ALIGN_2D)
                                   
                pairRow.setValue(xmipp.MDL_IMAGE_TILTED, getImageLocation(tImg))
                tMic = tMics[micId]
                pairRow.setValue(xmipp.MDL_MICROGRAPH_TILTED, tMic.getFileName())
                (angleY, angleY2, angleTilt) = angles.getAngles()
                pairRow.setValue(xmipp.MDL_ANGLE_Y, float(angleY))
                pairRow.setValue(xmipp.MDL_ANGLE_Y2, float(angleY2))
                pairRow.setValue(xmipp.MDL_ANGLE_TILT, float(angleTilt))
                
                pairRow.writeToMd(classMd, objId)
        
        classMd.write(blockMd, xmipp.MD_APPEND)
            
    def _reconstructImages(self, particles, deps):
        """ Function to insert the step needed to reconstruct a
        class (or setOfParticles) """
        classNo = particles.getObjId()
        blockMd = "class%06d_images@%s" % (classNo, self.rctClassesFn)
        classNameIn = blockMd
        classNameOut = self._getExtraPath("rct_images_%06d.xmd" % classNo)
        classVolumeOut = self._getPath("rct_%06d.vol" % classNo)
        
        if particles.hasRepresentative():
            classImage = getImageLocation(particles.getRepresentative())
        else:
            classImage = None
        
        reconStep = self._insertFunctionStep('reconstructClass',
                                             classNameIn, classNameOut,
                                             classImage, classVolumeOut,
                                             prerequisites=[deps])
        return reconStep

    def reconstructClass(self, classIn, classOut, classImage, classVolumeOut):
        # If class image doesn't exists, generate it by averaging
        if classImage is None:
            classRootOut = classOut.replace(".xmd", "") + "_"
            statsFn = self._getExtraPath('stats.xmd')
            args = "-i %(classIn)s --save_image_stats %(classRootOut)s -o %(statsFn)s"
            self.runJob("xmipp_image_statistics", args % locals(),
                        numberOfMpi=1, numberOfThreads=1)

            classImage = classRootOut + "average.xmp"
            
        centerMaxShift = self.maxShift.get()

        args = "-i %(classIn)s -o %(classOut)s --ref %(classImage)s " % locals()
        args += " --max_shift %(centerMaxShift)d" % locals()
        
        if self.thinObject.get():
            args += " --do_stretch"
        
        if self.skipTranslation.get():
            args += " --do_not_align_tilted"
        
        self.runJob("xmipp_image_align_tilt_pairs", args,
                    numberOfMpi=1, numberOfThreads=1)
        
        reconstructAdditionalParams = self.additionalParams.get()

        args = "-i %(classOut)s -o %(classVolumeOut)s " % locals()
        args += " %(reconstructAdditionalParams)s" % locals()

        if self.numberOfMpi >1:
            # threading is not supported in mpi version
            self.runJob("xmipp_reconstruct_art", args, numberOfThreads=1)
        else:
            args += " --thr %d" % self.numberOfThreads.get()
            self.runJob("xmipp_reconstruct_art", args)

        if exists(classVolumeOut):
            mdFn = self._getPath('volumes.xmd')
            md = xmipp.MetaData()
            
            if exists(mdFn):
                md.read(mdFn)
            objId = md.addObject()
            md.setValue(xmipp.MDL_IMAGE, classVolumeOut, objId)
                        
            if self.doFilter.get():
                filteredVolume = classVolumeOut.replace('.vol', '_filtered.vol')
                lowPassFilter = self.resoLowPassFilter.get()
                args = "-i %(classVolumeOut)s -o %(filteredVolume)s " % locals()
                args += " --fourier low_pass %(lowPassFilter)f " % locals()
                args += " --thr %d" % self.numberOfThreads.get()
                self.runJob("xmipp_transform_filter", args, numberOfMpi=1)
                objId = md.addObject()
                md.setValue(xmipp.MDL_IMAGE, filteredVolume, objId)
            md.write(mdFn)
                    
    def createOutputStep(self):
        # TODO: Refactor following code if possible
        self.volumesSet = self._createSetOfVolumes()
        self.volumesSet.setStore(False)
        self.sampling = self.inputParticlesTiltPair.get().getUntilted().getSamplingRate()
        self.volumesSet.setSamplingRate(self.sampling)

        if self.doFilter.get():
            self.volumesFilterSet = self._createSetOfVolumes('filtered')
            self.volumesFilterSet.setStore(False)
            self.volumesFilterSet.setSamplingRate(self.sampling)

        if isinstance(self.inputSet, SetOfParticles):
            volumeOut = self._getPath("rct_%06d.vol" % self.inputSet.getObjId())
            self._appendOutputVolume(volumeOut)
        else:
            for class2D in self.inputSet:
                volumeOut = self._getPath("rct_%06d.vol" % class2D.getObjId())
                self._appendOutputVolume(volumeOut)

        self._defineOutputs(outputVolumes=self.volumesSet)

        if self.doFilter.get():
            self._defineOutputs(outputFilteredVolumes=self.volumesFilterSet)

        self._defineSourceRelation(self.inputParticlesTiltPair, self.volumesSet)

    def _appendOutputVolume(self, volumeOut):
        vol = Volume()
        vol.setFileName(volumeOut)
        vol.setSamplingRate(self.sampling)
        self.volumesSet.append(vol)

        if self.doFilter.get():
            volumeFilterOut = volumeOut.replace('.vol', '_filtered.vol')
            volf = Volume()
            volf.setFileName(volumeFilterOut)
            volf.setSamplingRate(self.sampling)
            self.volumesFilterSet.append(volf)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []

        if isinstance(self.inputParticles.get(), SetOfParticles):
            if not self.inputParticles.get().hasAlignment():
                errors.append("The input particles should have "
                              "alignment information.")
        else:
            for class2D in self.inputParticles.get():
                if not class2D.hasAlignment():
                    errors.append("The input classes should have "
                                  "alignment information.")

        return errors
    
    def _summary(self):
        summary = []

        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            if isinstance(self.inputParticles.get(), SetOfParticles):
                summary.append("Input particles: %d"
                               % self.inputParticles.get().getSize())
            else:
                summary.append("Input classes: %d"
                               % self.inputParticles.get().getSize())
            summary.append("Output volumes: %d"
                           % self.outputVolumes.getSize())

            if self.doFilter.get():
                summary.append("Output filtered volumes: %d"
                               % self.outputFilteredVolumes.getSize())

        return summary
        
    def _methods(self):
        methods = []

        if not hasattr(self, 'outputVolumes'):
            methods.append("Output volumes not ready yet.")
        else:
            if isinstance(self.inputParticles.get(), SetOfParticles):
                methods.append('Set of %d particles %s was employed to create '
                               'an initial volume using RCT method.'
                               % (len(self.inputParticles.get()),
                                  self.getObjectTag('inputParticles')))
            else:
                particlesArray = [len(s) for s in self.inputParticles.get()]
                particlesArrayString = String(particlesArray)
                methods.append('Set of %d classes %s was employed to create %d '
                               'initial volumes using RCT method. '
                               % (len(self.inputParticles.get()),
                                  self.getObjectTag('inputParticles'),
                                  len(self.inputParticles.get())))
                methods.append('For each initial volume were used respectively '
                               '%s particles' % particlesArrayString)

            methods.append("Output volumes: %s" %
                           self.getObjectTag('outputVolumes'))

            if self.doFilter.get():
                methods.append("Output filtered volumes: %s"
                               % self.getObjectTag('outputFilteredVolumes'))

        return methods
            
    def _citations(self):
        return ['Sorzano2015b']
