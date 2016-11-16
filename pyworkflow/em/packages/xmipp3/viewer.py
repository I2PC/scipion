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

import os

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO, CommandView
from pyworkflow.em.data import *
from pyworkflow.em.protocol import *

from xmipp3 import getXmippPath, getEnviron
from pyworkflow.em.data_tiltpairs import MicrographsTiltPair, ParticlesTiltPair, CoordinatesTiltPair
from convert import *
from os.path import dirname, join
from pyworkflow.utils import makePath, runJob, copyTree, cleanPath
import pyworkflow as pw
import xmipp
import pyworkflow.gui.dialog as dialog

from protocol_cl2d_align import XmippProtCL2DAlign
from protocol_cl2d import XmippProtCL2D
from protocol_compare_reprojections import XmippProtCompareReprojections
from protocol_ctf_discrepancy import XmippProtCTFDiscrepancy
from protocol_extract_particles import XmippProtExtractParticles
from protocol_extract_particles_pairs import XmippProtExtractParticlesPairs
from protocol_helical_parameters import XmippProtHelicalParameters
from protocol_kerdensom import XmippProtKerdensom
from protocol_particle_pick_automatic import XmippParticlePickingAutomatic
from protocol_particle_pick import XmippProtParticlePicking
from protocol_particle_pick_pairs import XmippProtParticlePickingPairs
from protocol_preprocess import XmippProtPreprocessVolumes
from protocol_preprocess_micrographs import XmippProtPreprocessMicrographs
from protocol_rotational_spectra import XmippProtRotSpectra
from protocol_screen_particles import XmippProtScreenParticles
from protocol_ctf_micrographs import XmippProtCTFMicrographs
from pyworkflow.em.showj import *
from protocol_validate_nontilt import XmippProtValidateNonTilt
from protocol_multireference_alignability import XmippProtMultiRefAlignability
from protocol_assignment_tilt_pair import XmippProtAssignmentTiltPair





class XmippViewer(Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [                  
                CoordinatesTiltPair, 
                Image, 
                MicrographsTiltPair, 
                ParticlesTiltPair, 
                ProtExtractParticles, 
                SetOfClasses2D, 
                SetOfClasses3D, 
                SetOfCoordinates, 
                SetOfCTF, 
                SetOfImages, 
                SetOfMovies, 
                SetOfNormalModes, 
                XmippProtCompareReprojections, 
                XmippParticlePickingAutomatic, 
                XmippProtExtractParticlesPairs, 
                XmippProtKerdensom, 
                ProtParticlePicking, 
                XmippProtParticlePickingPairs,
                XmippProtRotSpectra, 
                XmippProtScreenParticles,
                XmippProtCTFMicrographs, 
                XmippProtValidateNonTilt,
                XmippProtAssignmentTiltPair,
                XmippProtMultiRefAlignability
                ]
    
    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)
        self._views = []   

    # FIXME: JMRT: I think this function is not necessary, we should remove it
    def visualize(self, obj, **kwargs):
        self._visualize(obj, **kwargs)
        
        for v in self._views:
            v.show()

    def __createTemporaryCtfs(self, obj, setOfMics):
        pwutils.cleanPath(obj._getPath("ctfs_temporary.sqlite"))
        self.protocol._createFilenameTemplates()
        ctfSet = self.protocol._createSetOfCTF("_temporary")

        for mic in setOfMics:
            micDir = obj._getExtraPath(removeBaseExt(mic.getFileName()))
            ctfparam = self.protocol._getFileName('ctfparam', micDir=micDir)

            if exists(ctfparam) or exists('xmipp_default_ctf.ctfparam'):
                if not os.path.exists(ctfparam):
                    ctfparam = 'xmipp_default_ctf.ctfparam'
                ctfModel = readCTFModel(ctfparam, mic)
                self.protocol._setPsdFiles(ctfModel, micDir)
                ctfSet.append(ctfModel)

        if not ctfSet.isEmpty():
            ctfSet.write()
            ctfSet.close()

        return ctfSet

    def _visualize(self, obj, **kwargs):
        cls = type(obj)

        if issubclass(cls, Volume):
            fn = getImageLocation(obj)
            self._views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={RENDER: 'image',
                                                      SAMPLINGRATE: obj.getSamplingRate()}))
                 
        elif issubclass(cls, Image):
            fn = getImageLocation(obj)
            self._views.append(ObjectView(self._project, obj.strId(), fn))
            
        elif issubclass(cls, SetOfNormalModes):
            fn = obj.getFileName()
            from nma.viewer_nma import OBJCMD_NMA_PLOTDIST, OBJCMD_NMA_VMD
            objCommands = "'%s' '%s'" % (OBJCMD_NMA_PLOTDIST, OBJCMD_NMA_VMD)
            self._views.append(ObjectView(self._project, self.protocol.strId(),
                                          fn, obj.strId(),
                                          viewParams={OBJCMDS: objCommands},
                                          **kwargs))

        elif issubclass(cls, SetOfMovies):
            fn = obj.getFileName()
            # Enabled for the future has to be available
            labels = 'id _filename _samplingRate  '
            self._views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={ORDER: labels, 
                                                      VISIBLE: labels, 
                                                      MODE: MODE_MD, RENDER: "no"}))

        elif issubclass(cls, SetOfMicrographs):            
            fn = obj.getFileName()
            self._views.append(ObjectView(self._project, obj.strId(), fn, **kwargs))
            
        elif issubclass(cls, MicrographsTiltPair):
            labels = 'id enabled _untilted._filename _tilted._filename'
            self._views.append(ObjectView(self._project, obj.strId(), obj.getFileName(),
                                          viewParams={ORDER: labels, 
                                                      VISIBLE: labels, 
                                                      MODE: MODE_MD,
                                                      RENDER: '_untilted._filename _tilted._filename'}))
            
        elif issubclass(cls, ParticlesTiltPair):          
            labels = 'id enabled _untilted._filename _tilted._filename'
            self._views.append(ObjectView(self._project, obj.strId(), obj.getFileName(),
                                          viewParams={ORDER: labels, 
                                                      VISIBLE: labels,
                                                      RENDER:'_untilted._filename _tilted._filename',
                                                      MODE: MODE_MD}))

        elif issubclass(cls, SetOfCoordinates):
            micSet = obj.getMicrographs()  # accessing mics to provide metadata file
            if micSet is None:
                raise Exception('visualize: SetOfCoordinates has no micrographs set.')
            
            mdFn = getattr(micSet, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:  # happens if protocol is not an xmipp one
                fn = self._getTmpPath(micSet.getName() + '_micrographs.xmd')
                writeSetOfMicrographs(micSet, fn)
            tmpDir = self._getTmpPath(obj.getName())
            cleanPath(tmpDir)
            makePath(tmpDir)
            # FIXME: (JMRT) We are always writing the SetOfCoordinates and removing
            # the tmpDir, we need to take into account if the user have pick
            # some particles in the tmpDir and have not save them, that now
            # will loose all picked particles.
            # A possible solution could be to alert that changes have not been
            # written during modification of tmpDir or create a new Xmipp picking
            # protocol to continue picking later without loosing the coordinates.
            writeSetOfCoordinates(tmpDir, obj)
            self._views.append(CoordinatesObjectView(self._project, fn, tmpDir,
                                                     self.protocol,
                                                     inTmpFolder=True))

        elif issubclass(cls, SetOfParticles):
            fn = obj.getFileName()
            labels = 'id enabled _index _filename _xmipp_zScore _xmipp_cumulativeSSNR _sampling '
            labels += '_ctfModel._defocusU _ctfModel._defocusV _ctfModel._defocusAngle _transform._matrix'
            self._views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={ORDER: labels, 
                                                      VISIBLE: labels, 
                                                      'sortby': '_xmipp_zScore asc',
                                                      RENDER:'_filename'}))
               
        elif issubclass(cls, SetOfVolumes):
            fn = obj.getFileName()
            labels = 'id enabled comment _filename '
            self._views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={MODE: MODE_MD,
                                                      ORDER: labels,
                                                      VISIBLE: labels,
                                                      RENDER: '_filename'}))

        elif issubclass(cls, SetOfClasses2D):
            fn = obj.getFileName()
            self._views.append(ClassesView(self._project, obj.strId(), fn, **kwargs))
            
        elif issubclass(cls, SetOfClasses3D):
            fn = obj.getFileName()
            self._views.append(Classes3DView(self._project, obj.strId(), fn))

        elif issubclass(cls, SetOfImages):
            fn = obj.getFileName()
            self._views.append(
                ObjectView(self._project, obj.strId(), fn, **kwargs))
        
        if issubclass(cls, XmippProtCTFMicrographs):
            if obj.hasAttribute('outputCTF'):
                ctfSet = obj.outputCTF
            else:
                mics = obj.inputMicrographs.get()
                ctfSet = self.__createTemporaryCtfs(obj, mics)

            if ctfSet.isEmpty():
                self._views.append(self.infoMessage("No CTF estimation has finished yet"))
            else:
                self._views.append(CtfView(self._project, ctfSet))

        elif issubclass(cls, SetOfCTF):
            self._views.append(CtfView(self._project, obj))

        elif issubclass(cls, CoordinatesTiltPair):
            tmpDir = self._getTmpPath(obj.getName()) 
            makePath(tmpDir)

            mdFn = join(tmpDir, 'input_micrographs.xmd')
            writeSetOfMicrographsPairs(obj.getUntilted().getMicrographs(),
                                        obj.getTilted().getMicrographs(), 
                                        mdFn) 
            parentProtId = obj.getObjParentId()
            parentProt = self.getProject().mapper.selectById(parentProtId)
            extraDir = parentProt._getExtraPath()
            
 #           extraDir = parentProt._getExtraPath()
            #TODO: Review this if ever a non Xmipp CoordinatesTiltPair is available
            writeSetOfCoordinates(tmpDir, obj.getUntilted()) 
            writeSetOfCoordinates(tmpDir, obj.getTilted()) 
            launchTiltPairPickerGUI(mdFn, tmpDir, self.protocol)


        elif issubclass(cls, XmippProtExtractParticles) or issubclass(cls, XmippProtScreenParticles):
            particles = obj.outputParticles
            self._visualize(particles)
            
            fn = obj._getPath('images.xmd')
            md = xmipp.MetaData(fn) 
            # If Zscore on output images plot Zscore particle sorting
            if md.containsLabel(xmipp.MDL_ZSCORE):
                from plotter import XmippPlotter
                xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
                xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
                xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_ZSCORE)
                self._views.append(xplotter)
    
        elif issubclass(cls, XmippProtRotSpectra):
            self._visualize(obj.outputClasses,
                            viewParams={'columns': obj.SomXdim.get(),
                                        RENDER: ' spectraPlot._filename average._filename',
                                        ZOOM: 30,
                                        VISIBLE:  'enabled id _size average._filename spectraPlot._filename',
                                        'labels': 'id _size',
                                        SORT_BY: 'id'})
        
        elif issubclass(cls, XmippProtKerdensom):
            self._visualize(obj.outputClasses,
                            viewParams={'columns': obj.SomXdim.get(),
                                       'render': 'average._filename _representative._filename',
                                       'labels': '_size',
                                       'sortby': 'id'})
        

        
        elif issubclass(cls, XmippProtCompareReprojections):
                fn = obj.outputParticles.getFileName()
                labels = 'id enabled _index _xmipp_image._filename _xmipp_imageRef._filename _xmipp_imageResidual._filename _xmipp_imageCovariance._filename _xmipp_cost _xmipp_zScoreResCov _xmipp_zScoreResMean _xmipp_zScoreResVar _xmipp_continuousA _xmipp_continuousB _xmipp_continuousX _xmipp_continuousY'
                labelRender = "_xmipp_image._filename _xmipp_imageRef._filename _xmipp_imageResidual._filename _xmipp_imageCovariance._filename"
                self._views.append(ObjectView(self._project, obj.outputParticles.strId(), fn,
                                              viewParams={ORDER: labels, 
                                                      VISIBLE: labels, 
                                                      SORT_BY: '_xmipp_cost asc', RENDER:labelRender,
                                                      MODE: MODE_MD}))
            
        elif issubclass(cls, XmippParticlePickingAutomatic):
            micSet = obj.getInputMicrographs()
            mdFn = getattr(micSet, '_xmippMd', None)
            inTmpFolder = False
            if mdFn:
                micsfn = mdFn.get()
            else:  # happens if protocol is not an xmipp one
                micsfn = self._getTmpPath(micSet.getName() + '_micrographs.xmd')
                writeSetOfMicrographs(micSet, micsfn)
                inTmpFolder = True
                
            posDir = obj._getExtraPath()  
            memory = '%dg'%obj.memory.get(), 
            launchSupervisedPickerGUI(micsfn, posDir, obj, mode='review', memory=memory, inTmpFolder=inTmpFolder)

         # We need this case to happens before the ProtParticlePicking one
        elif issubclass(cls, XmippProtAssignmentTiltPair):
            if obj.getOutputsSize() >= 1:
                coordsSet = obj.getCoordsTiltPair()
                self._visualize(coordsSet)  
                
        elif issubclass(cls, ProtParticlePicking):
            if obj.getOutputsSize() >= 1:
                coordsSet = obj.getCoords()
                self._visualize(coordsSet)
            
        elif issubclass(cls, XmippProtValidateNonTilt):
            outputVols = obj.outputVolumes
            labels = 'id enabled comment _filename weight'
            self._views.append(ObjectView(self._project, outputVols.strId(), outputVols.getFileName(),
                                          viewParams={MODE: MODE_MD, VISIBLE:labels, ORDER: labels,
                                                      SORT_BY: 'weight desc', RENDER: '_filename'}))
        
        elif issubclass(cls, XmippProtMultiRefAlignability):
            outputVols = obj.outputVolumes
            labels = 'id enabled comment _filename weightAlignabilityPrecision weightAlignabilityAccuracy'
            self._views.append(ObjectView(self._project, outputVols.strId(), outputVols.getFileName(),
                                          viewParams={MODE: MODE_MD, VISIBLE:labels, ORDER: labels,
                                                      SORT_BY: 'weightAlignabilityAccuracy desc', RENDER: '_filename'}))
            
            fn = obj.outputParticles.getFileName()
            labels = 'id enabled _index _filename _xmipp_scoreAlignabilityAccuracy _xmipp_scoreAlignabilityPrecision'
            labelRender = "_filename"
            self._views.append(ObjectView(self._project, obj.outputParticles.strId(), fn,
                                            viewParams={ORDER: labels, 
                                                      VISIBLE: labels, 
                                                      SORT_BY: '_xmipp_scoreAlignabilityAccuracy desc', RENDER:labelRender,
                                                      MODE: MODE_MD}))
            
            fn = obj._getExtraPath('vol001_pruned_particles_alignability.xmd')
            md = xmipp.MetaData(fn)
            from plotter import XmippPlotter
            from pyworkflow.em.plotter import EmPlotter
            plotter = XmippPlotter()
            plotter.createSubPlot('Soft-alignment validation plot','Angular Precision', 'Angular Accuracy')
            plotter.plotMdFile(md, xmipp.MDL_SCORE_BY_ALIGNABILITY_PRECISION, xmipp.MDL_SCORE_BY_ALIGNABILITY_ACCURACY,
                               marker='.', markersize=.55, color='red', linestyle='')
            self._views.append(plotter)


        elif issubclass(cls, XmippProtExtractParticlesPairs):
            self._visualize(obj.outputParticlesTiltPair)

        return self._views
    
    

        
