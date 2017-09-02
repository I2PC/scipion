# **************************************************************************
# *
# * Authors:     Grigory Sharov     (sharov@igbmc.fr)
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

import pyworkflow.object as pwobj
from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import (PointerParam, FloatParam, StringParam,
                                        BooleanParam, IntParam, LEVEL_ADVANCED)
from pyworkflow.utils import removeExt
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtParticles
from pyworkflow.em.data import SetOfClasses3D, SetOfParticles, SetOfClasses
from convert import (convertBinaryVol, readSetOfParticles,
                     writeSetOfParticles, writeReferences)
import pyworkflow.em.metadata as md


class ProtRelionSortParticles(ProtParticles):
    """
    Relion particle sorting protocol.
    It calculates difference images between particles and their aligned
    (and CTF-convoluted) references, and produces Z-score on the characteristics
    of these difference images (such as mean, standard deviation, skewness,
    excess kurtosis and rotational symmetry).

    """
    _label = 'sort particles'
    _lastUpdateVersion = VERSION_1_1

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfParticles,SetOfClasses2D,SetOfClasses3D',
                      label='Input particles', important=True,
                      help='Select a set of particles after Relion auto-picking'
                           ' or 3D refinement. You can also select Relion 2D/3D'
                           ' classes. Particles should have at least in-plane '
                           'alignment parameters and class assignment.')

        form.addParam('referenceAverages', PointerParam,
                      pointerClass="SetOfAverages",
                      condition='inputSet and isInputAutoPicking',
                      label='Reference 2D averages',
                      help='Select references 2D averages used for auto-picking')

        form.addParam('referenceVolume', PointerParam,
                      pointerClass="Volume",
                      condition='inputSet and isInputAutoRefine',
                      label='Reference volume',
                      help='Select reference volume 2D after 3D auto-refine')

        form.addParam('maskDiameterA', IntParam, default=-1,
                      label='Particle mask diameter (A)',
                      help='The experimental images will be masked with a '
                           'soft circular mask with this <diameter>. Make '
                           'sure this diameter is not set too small because '
                           'that may mask away part of the signal! If set to '
                           'a value larger than the image  size no masking '
                           'will be performed.\n\n'
                           'The same diameter will also be used for a '
                           'spherical mask of the reference structures if no '
                           'user-provided mask is specified.')
        form.addParam('doLowPass', IntParam, default=-1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Low pass filter references to (A):',
                      help='Lowpass filter in Angstroms for the references '
                           '(prevent Einstein-from-noise!)')
        form.addParam('doInvert', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Invert contrast of references?',
                      help='Density in particles is inverted compared to the '
                           'density in references')
        form.addParam('doCTF', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Do CTF-correction?',
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement.  The resulting algorithm '
                           'intrinsically implements the optimal linear,  '
                           'or Wiener filter. Note that input particles '
                           'should contains CTF parameters.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED, condition='doCTF',
                      label='Ignore CTFs until their first peak?',
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak of each CTF '
                           'onward. This can be useful if the CTF model is '
                           'inadequate at the lowest resolution. Still, in '
                           'general using higher amplitude contrast on the '
                           'CTFs (e.g. 10-20%) often yields better results. '
                           'Therefore, this option is not generally '
                           'recommended.')
        form.addParam('minZ', FloatParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label='Min Z-value?',
                      help='Minimum Z-value to count in the sorting of '
                           'outliers')
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED, label='Additional parameters',
                      help='In this box command-line arguments may be '
                           'provided that are not generated by the GUI. This '
                           'may be useful for testing  developmental options '
                           'and/or expert use of the program, e.g: \n'
                           '--verb 1\n')

        form.addParallelSection(threads=0, mpi=1)
            
    #--------------------------- INSERT steps functions ------------------------

    def isInputAutoPicking(self):
        inputSet = self.inputSet.get()
        return (isinstance(inputSet, SetOfParticles) and
                not inputSet.hasAlignmentProj())

    def isInputAutoRefine(self):
        inputSet = self.inputSet.get()
        return (isinstance(inputSet, SetOfParticles) and
                inputSet.hasAlignmentProj())

    def isInputClasses(self):
        return isinstance(self.inputSet.get(), SetOfClasses)

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        refAvg = self.referenceAverages.get()
        avgId = refAvg.getObjId() if refAvg is not None else None
        refVol = self.referenceVolume.get()
        volId = refVol.getObjId() if refVol is not None else None

        self._insertFunctionStep('convertInputStep',
                                 self.inputSet.get().getObjId(), avgId, volId)
        self._insertRelionStep()
        self._insertFunctionStep('createOutputStep')

    def _insertRelionStep(self):
        """ Prepare the command line arguments before calling Relion. """
        # Join in a single line all key, value pairs of the args dict
        args = {}
        self._setArgs(args)
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.iteritems()])

        if self.extraParams.hasValue():
            params += ' ' + self.extraParams.get()

        self._insertFunctionStep('runRelionStep', params)

    #--------------------------- STEPS functions -------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'input_particles': self._getExtraPath('input_particles.star'),
            'input_refs': self._getExtraPath('input_references.star'),
            'output_star': self._getExtraPath('input_particles_sorted.star'),
            'input_refvol': self._getTmpPath('input_vol.mrc')
        }
        self._updateFilenamesDict(myDict)

    def convertInputStep(self, inputId, avgId, volId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        Params:
            particlesId: use this parameters just to force redo of convert if
                the input particles are changed.
        """
        inputSet = self.inputSet.get()
        imgStar = self._getFileName('input_particles')
        refStar = self._getFileName('input_refs')
        # Pass stack file as None to avoid write the images files
        self.info("Converting set from '%s' into '%s'"
                  % (inputSet.getFileName(), imgStar))

        refSet = None
        # case refine3D
        if self.isInputClasses():
            refSet = self.inputSet.get()  # 2D or 3D classes

        else:
            if self.isInputAutoRefine():
                em.ImageHandler().convert(self.referenceVolume.get(),
                                       self._getFileName('input_refvol'))
            else: # Autopicking case
                refSet = self.referenceAverages.get()

        self.classDict = {}

        if refSet:
            self.info("Converting reference from '%s' into %s"
                      % (refSet.getFileName(), refStar))

            # Compute class mapping
            classList = [cls.getObjId() for cls in refSet]
            classList.sort()
            for i, c in enumerate(classList):
                self.classDict[c] = i + 1

            writeReferences(refSet, removeExt(refStar),
                            postprocessImageRow=self._updateClasses)

        # Write particles star file
        allParticles = self._allParticles(iterate=False)
        writeSetOfParticles(allParticles, imgStar, self._getPath(),
                            postprocessImageRow=self._postProcessImageRow)

    def runRelionStep(self, params):
        """ Execute relion steps with given params. """
        self.runJob(self._getProgram(), params)
        
    def createOutputStep(self):
        sortedImgSet = self._createSetOfParticles()

        particles = self._sampleParticles()
        sortedImgSet.copyInfo(particles)

        particlesIter = self._allParticles(iterate=False)

        sortedStar = self._getFileName('output_star')
        sortedImgSet.copyItems(particlesIter,
                               updateItemCallback=self._updateZScore,
                               itemDataIterator=md.iterRows(sortedStar))

        self._defineOutputs(outputParticles=sortedImgSet)
        self._defineSourceRelation(self.inputSet, sortedImgSet)


    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        #imgSet = self.inputParticles.get()
        #if not imgSet.hasCTF() and self.doCTF:
        #    errors.append('Input particles have no CTF information!')

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Input %s particles were sorted by Z-score: %s" %
                           (self.inputSet.get().getSize(),
                            self.getObjectTag('outputParticles')))
        return summary
    
    #--------------------------- UTILS functions -------------------------------
    def _allParticles(self, iterate=False):
        # A handler function to iterate over the particles
        inputSet = self.inputSet.get()

        if self.isInputClasses():
            iterParticles = inputSet.iterClassItems()
            if iterate:
                return iterParticles
            else:
                particles = SetOfParticles(filename=":memory:")
                particles.copyInfo(inputSet.getFirstItem())
                particles.copyItems(iterParticles)
                return particles
        else:
            if iterate:
                return inputSet.iterItems()
            else:
                return inputSet

    def _sampleParticles(self):
        """ Return the input particles set, or the first class
        in the case the input are classes.
        """
        if self.isInputClasses():
            return self.inputSet.get().getFirstItem()
        else:
            return self.inputSet.get()

    def _setArgs(self, args):
        from pyworkflow.em.packages.relion.convert import getVersion
        particles = self._sampleParticles()

        if self.maskDiameterA <= 0:
            maskDiameter = particles.getSamplingRate() * particles.getXDim()
        else:
            maskDiameter = self.maskDiameterA.get()

        args.update({'--i': self._getFileName('input_particles'),
                     '--particle_diameter': maskDiameter,
                     '--angpix': particles.getSamplingRate(),
                     '--min_z': self.minZ.get()
                     })
        
        if getVersion() == "2.0":
            args['--o'] = self._getFileName('output_star')
        else:
            args['--o'] = 'sorted'
            
        #if inputReferences is a volume, convert it to mrc here
        if self.isInputAutoRefine():
            args['--ref'] = self._getFileName('input_refvol')
        else:
            args['--ref'] = self._getFileName('input_refs')

        if self.doInvert:
            args['--invert'] = ''

        if self.ignoreCTFUntilFirstPeak:
            args['--ctf_intact_first_peak'] = ''

        if self.doCTF:
            args['--ctf'] = ''

    def _getProgram(self, program='relion_particle_sort'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    def _postProcessImageRow(self, img, imgRow):
        if self.classDict:
            classId = img.getClassId()

            if classId is not None:
                if not classId in self.classDict:
                    raise Exception("Class Id %s from particle %s is not found"
                                    % (classId, img.getObjId()))
                newClassId = self.classDict[classId]
            else:
                newClassId = imgRow.getValue(md.RLN_PARTICLE_CLASS)
        else:
            newClassId = 1

        imgRow.setValue(md.RLN_PARTICLE_CLASS, newClassId)

    def _updateClasses(self, img, imgRow):
        index, _ = img.getLocation()
        # renumber class numbers
        imgRow.setValue(md.RLN_PARTICLE_CLASS, int(index))

    def _updateZScore(self, img, imgRow):
        zscore = imgRow.getValue(md.RLN_SELECT_PARTICLES_ZSCORE)
        img._rlnSelectParticlesZscore = pwobj.Float(zscore)
