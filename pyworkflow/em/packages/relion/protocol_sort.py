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

from pyworkflow.protocol.params import (PointerParam, FloatParam, StringParam,
                                        BooleanParam, IntParam, LEVEL_ADVANCED)
from pyworkflow.utils import removeExt
from pyworkflow.em.protocol import ProtParticles
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

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label='Input particles',
                      help='Select a set of particles after Relion 2D/3D classification, refinement '
                           'or auto-picking. Particles should have at least in-plane alignment '
                           'parameters and class assignment.')
        form.addParam('inputReferences', PointerParam,
                      pointerClass="SetOfClasses2D, SetOfClasses3D, SetOfVolumes, Volume, SetOfAverages",
                      label='Input references',
                      help='Select references: a set of classes/averages or a 3D volume')
        form.addParam('maskDiameterA', IntParam, default=-1,
                      label='Particle mask diameter (A)',
                      help='The experimental images will be masked with a soft circular mask '
                           'with this <diameter>. '
                           'Make sure this diameter is not set too small because that may mask '
                           'away part of the signal! If set to a value larger than the image '
                           'size no masking will be performed.\n\n'
                           'The same diameter will also be used for a spherical mask of the '
                           'reference structures if no user-provided mask is specified.')
        form.addParam('doLowPass', IntParam, default=-1,
                      label='Low pass filter references to (A):',
                      help='Lowpass filter in Angstroms for the references (prevent Einstein-from-noise!)')
        form.addParam('doInvert', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Invert contrast of references?',
                      help='Density in particles is inverted compared to the density in references')
        form.addParam('doCTF', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Do CTF-correction?',
                      help='If set to Yes, CTFs will be corrected inside the MAP refinement. '
                           'The resulting algorithm intrinsically implements the optimal linear, '
                           'or Wiener filter. Note that input particles should contains CTF parameters.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED, condition='doCTF',
                      label='Ignore CTFs until their first peak?',
                      help='If set to Yes, then CTF-amplitude correction will only be performed from the '
                           'first peak of each CTF onward. This can be useful if the CTF model is '
                           'inadequate at the lowest resolution. Still, in general using higher '
                           'amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. '
                           'Therefore, this option is not generally recommended.')
        form.addParam('minZ', FloatParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label='Min Z-value?',
                      help='Minimum Z-value to count in the sorting of outliers')
        form.addParam('extraParams', StringParam, default='', expertLevel=LEVEL_ADVANCED,
                      label='Additional parameters',
                      help='In this box command-line arguments may be provided that '
                           'are not generated by the GUI. This may be useful for testing '
                           'developmental options and/or expert use of the program, e.g: \n'
                           '--verb 1\n')

        form.addParallelSection(threads=0, mpi=1)
            
    #--------------------------- INSERT steps functions ------------------------

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.get().getObjId(),
                                 self.inputReferences.get().getObjId())
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
            'output_star': self._getExtraPath('input_particles_sorted.star')}
        self._updateFilenamesDict(myDict)

    def convertInputStep(self, particlesId, referencesId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        Params:
            particlesId: use this parameters just to force redo of convert if
                the input particles are changed.
        """
        imgSet = self.inputParticles.get()
        refSet = self.inputReferences.get()
        imgStar = self._getFileName('input_particles')
        refStar = self._getFileName('input_refs')
        self.classesList = []

        for img in imgSet:
            # Take classId from img
            classId = img.getClassId()
            if not classId in self.classesList:
                self.classesList.append(classId)
        self.classesList.sort()
        # Pass stack file as None to avoid write the images files
        self.info("Converting set from '%s' into '%s'"
                  % (imgSet.getFileName(), imgStar))

        # case refine3D
        if self.classesList[0] == None and refSet.getClassName() == 'Volume':
            writeSetOfParticles(imgSet, imgStar, self._getPath(),
                                postprocessImageRow=self._resetClasses)
        # case auto-pick, no classIds
        elif self.classesList[0] == None:
            writeSetOfParticles(imgSet, imgStar, self._getPath())
        # case cl2d or cl3d
        else:
            writeSetOfParticles(imgSet, imgStar, self._getPath(),
                                postprocessImageRow=self._postProcessImageRow)

        if not refSet.getClassName() == 'Volume':
            self.info("Converting set from '%s' into %s"
                      % (refSet.getFileName(), refStar))
            writeReferences(refSet, removeExt(refStar),
                            postprocessImageRow=self._updateClasses)

    def runRelionStep(self, params):
        """ Execute relion steps with given params. """
        self.runJob(self._getProgram(), params)
        
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        sortedImgSet = self._createSetOfParticles()
        sortedImgSet.copyInfo(imgSet)
        readSetOfParticles(self._getFileName('output_star'),
                           sortedImgSet,
                           alignType=imgSet.getAlignment())

        self._defineOutputs(outputParticles=sortedImgSet)
        self._defineSourceRelation(imgSet, sortedImgSet)

    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        imgSet = self.inputParticles.get()
        if not imgSet.hasCTF() and self.doCTF:
            errors.append('Input particles have no CTF information!')

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append("Input %s particles: %s were sorted by Z-score" %
                           (self.inputParticles.get().getSize(),
                            self.inputParticles.get().getNameId()))
        return summary
    
    #--------------------------- UTILS functions -------------------------------
    def _setArgs(self, args):
        maskDiameter = self.maskDiameterA.get()
        if maskDiameter <= 0:
            x, _, _ = self.inputParticles.get().getDim()
            maskDiameter = self.inputParticles.get().getSamplingRate() * x

        args.update({'--i': self._getFileName('input_particles'),
                     '--particle_diameter': maskDiameter,
                     '--angpix': self.inputParticles.get().getSamplingRate(),
                     '--min_z': self.minZ.get(),
                     '--o': 'sorted'
                     })
        #if inputReferences is a volume, convert it to mrc here
        if self.inputReferences.get().getClassName() == 'Volume':
            args['--ref'] = convertBinaryVol(self.inputReferences.get(),
                                             self._getTmpPath())

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
        if img.hasClassId() and imgRow.hasLabel(md.RLN_PARTICLE_AUTOPICK_FOM):
            # classId from classification has higher priority then autopick classId
            imgRow.setValue(md.RLN_PARTICLE_CLASS, int(img.getClassId()))

        #now renumber classes
        classId = imgRow.getValue(md.RLN_PARTICLE_CLASS)
        newClassId = [i for i, x in enumerate(self.classesList) if x == classId]
        #print "oldcls=", classId, "newcls=", newClassId[0] + 1
        imgRow.setValue(md.RLN_PARTICLE_CLASS, newClassId[0] + 1)

    def _updateClasses(self, img, imgRow):
        index, _ = img.getLocation()
        # renumber class numbers
        imgRow.setValue(md.RLN_PARTICLE_CLASS, int(index))

    def _resetClasses(self, img, imgRow):
        # after auto-refine 3D we have just 1 class
        imgRow.setValue(md.RLN_PARTICLE_CLASS, 1)
