# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
# *           Qiyu Jin
# *           Slavica Jonic                (jonic@impmc.upmc.fr)
# * Ported to Scipion:
# *           J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Jan 2014
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


from os.path import basename

from pyworkflow.utils import isPower2, getListFromRangeString
from pyworkflow.utils.path import copyFile, cleanPath 
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.packages.xmipp3 import XmippMdRow
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles,\
     xmippToLocation, getImageLocation, setXmippAttributes
from pyworkflow.protocol.params import NumericRangeParam

from convert import modeToRow

import pyworkflow.em.metadata as md


NMA_ALIGNMENT_WAV = 0
NMA_ALIGNMENT_PROJ = 1    

        
class XmippProtAlignmentNMA(ProtAnalysis3D):
    """ Protocol for flexible angular alignment. """
    _label = 'nma alignment'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputModes', params.PointerParam, pointerClass='SetOfNormalModes',
                      label="Normal modes",                        
                      help='Set of normal modes to explore.')
        form.addParam('modeList', NumericRangeParam, expertLevel=params.LEVEL_ADVANCED,
                      label="Modes selection",
                      help='Select which modes do you want to use from all of them.\n'
                           'If you leave the field empty, all modes will be used.\n'
                           'You have several ways to specify selected modes.\n'
                           '   Examples:\n'
                           ' "7,8-10" -> [7,8,9,10]\n'
                           ' "8, 10, 12" -> [8,10,12]\n'
                           ' "8 9, 10-12" -> [8,9,10,11,12])\n')
        form.addParam('inputParticles', params.PointerParam, label="Input particles", 
                      pointerClass='SetOfParticles',
                      help='Select the set of particles that you want to use for flexible analysis.')  

        form.addParam('copyDeformations', params.PathParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Precomputed results (for development)',
                      help='Enter a metadata file with precomputed elastic  \n'
                           'and rigid-body alignment parameters and perform \n'
                           'all remaining steps using this file.')
        
        form.addSection(label='Angular assignment')
        form.addParam('trustRegionScale', params.IntParam, default=1,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Trust region scale',
                      help='For elastic alignment, this parameter scales the initial \n'
                           'value of the trust region radius for optimization purposes.\n'
                           'Use larger values for larger expected deformation amplitudes.')    
        form.addParam('alignmentMethod', params.EnumParam, default=NMA_ALIGNMENT_WAV,
                      choices=['wavelets & splines', 'projection matching'],
                      label='Alignment method',
                      help='For rigid-body alignment, use Projection Matching (faster) instead\n'
                           'of Wavelets and Splines (more accurate). In the case of Wavelets \n'
                           'and Splines, the size of images should be a power of 2.')
        form.addParam('discreteAngularSampling', params.FloatParam, default=10,
                      label="Discrete angular sampling (deg)", 
                      help='This parameter is used in Projection Matching and Wavelets methods\n'
                           'for a rough rigid-body alignment. It is the angular step (in degrees)\n'
                           'with which the library of reference projections is computed. This \n'
                           'alignment is refined with Splines method if Wavelets and Splines \n'
                           'alignment is chosen.')
                      
        form.addParallelSection(threads=0, mpi=8)    
    
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def getInputPdb(self):
        """ Return the Pdb object associated with the normal modes. """
        return self.inputModes.get().getPdb()
    
    def _insertAllSteps(self):
        atomsFn = self.getInputPdb().getFileName()
        # Define some outputs filenames
        self.imgsFn = self._getExtraPath('images.xmd') 
        self.atomsFn = self._getExtraPath(basename(atomsFn))
        self.modesFn = self._getExtraPath('modes.xmd')
        
        self._insertFunctionStep('convertInputStep', atomsFn) 
        
        if self.copyDeformations.empty(): #ONLY FOR DEBUGGING
            self._insertFunctionStep("performNmaStep", self.atomsFn, self.modesFn)
        else:   
            # TODO: for debugging and testing it will be useful to copy the deformations
            # metadata file, not just the deformation.txt file         
            self._insertFunctionStep('copyDeformationsStep', self.copyDeformations.get())
            
        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions --------------------------------------------   
    def convertInputStep(self, atomsFn):
        # Write the modes metadata taking into account the selection
        self.writeModesMetaData()
        # Write a metadata with the normal modes information
        # to launch the nma alignment programs
        writeSetOfParticles(self.inputParticles.get(), self.imgsFn)
        # Copy the atoms file to current working dir
        copyFile(atomsFn, self.atomsFn)
            
    def writeModesMetaData(self):
        """ Iterate over the input SetOfNormalModes and write
        the proper Xmipp metadata.
        Take into account a possible selection of modes (This option is 
        just a shortcut for testing. The recommended
        way is just create a subset from the GUI and use that as input)
        """
        modeSelection = []
        if self.modeList.empty():
            modeSelection = []
        else:
            modeSelection = getListFromRangeString(self.modeList.get())
            
        mdModes = md.MetaData()
        
        inputModes = self.inputModes.get()
        for mode in inputModes:
            # If there is a mode selection, only
            # take into account those selected
            if not modeSelection or mode.getObjId() in modeSelection:
                row = XmippMdRow()
                modeToRow(mode, row)
                row.writeToMd(mdModes, mdModes.addObject())
        mdModes.write(self.modesFn)
            
    def copyDeformationsStep(self, deformationMd):
        copyFile(deformationMd, self.imgsFn)
        # We need to update the image name with the good ones
        # and the same with the ids.
        inputSet = self.inputParticles.get()
        mdImgs = md.MetaData(self.imgsFn)
        for objId in mdImgs:
            imgPath = mdImgs.getValue(md.MDL_IMAGE, objId)
            index, fn = xmippToLocation(imgPath)
            # Conside the index is the id in the input set
            particle = inputSet[index]
            mdImgs.setValue(md.MDL_IMAGE, getImageLocation(particle), objId)
            mdImgs.setValue(md.MDL_ITEM_ID, long(particle.getObjId()), objId)
        mdImgs.write(self.imgsFn)
        
    def performNmaStep(self, atomsFn, modesFn):
        sampling = self.inputParticles.get().getSamplingRate()
        discreteAngularSampling = self.discreteAngularSampling.get()
        trustRegionScale = self.trustRegionScale.get()
        odir = self._getTmpPath()
        imgFn = self.imgsFn
        
        args = "-i %(imgFn)s --pdb %(atomsFn)s --modes %(modesFn)s --sampling_rate %(sampling)f "
        args += "--discrAngStep %(discreteAngularSampling)f --odir %(odir)s --centerPDB "
        args += "--trustradius_scale %(trustRegionScale)d --resume "
        
        if self.getInputPdb().getPseudoAtoms():
            args += "--fixed_Gaussian "
        
        if self.alignmentMethod == NMA_ALIGNMENT_PROJ:
            args += "--projMatch "
    
        self.runJob("xmipp_nma_alignment", args % locals())
        
        cleanPath(self._getPath('nmaTodo.xmd'))
    
    def createOutputStep(self):
        inputSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        pdbPointer = self.inputModes.get()._pdbPointer
        
        partSet.copyInfo(inputSet)
        partSet.copyItems(inputSet,
                          updateItemCallback=self._updateParticle,
                          itemDataIterator=md.iterRows(self.imgsFn, sortByLabel=md.MDL_ITEM_ID))
        
        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(pdbPointer, partSet)
        self._defineTransformRelation(self.inputParticles, partSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        errors = []
        xdim = self.inputParticles.get().getDim()[0]
        if not isPower2(xdim):
            errors.append("Image dimension (%s) is not a power of two, consider resize them" % xdim)
        return errors
    
    def _citations(self):
        return ['Jonic2005', 'Sorzano2004b']
    
    def _methods(self):
        pass
    
    #--------------------------- UTILS functions --------------------------------------------
    def _printWarnings(self, *lines):
        """ Print some warning lines to 'warnings.xmd', 
        the function should be called inside the working dir."""
        fWarn = open("warnings.xmd",'w')
        for l in lines:
            print >> fWarn, l
        fWarn.close()

    def _getLocalModesFn(self):
        modesFn = self.inputModes.get().getFileName()
        return self._getBasePath(modesFn)
    
    def _updateParticle(self, item, row):
        setXmippAttributes(item, row, md.MDL_NMA, md.MDL_COST)
