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
This module implement the wrappers around xmipp_showj
visualization program.
"""

import os
from pyworkflow.viewer import Viewer, Wizard, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.data import *
from pyworkflow.em.protocol import *
from pyworkflow.utils.process import runJob
from xmipp3 import getXmippPath
from protocol_preprocess_micrographs import XmippProtPreprocessMicrographs
from protocol_ctf_micrographs import XmippProtCTFMicrographs
from protocol_extract_particles import XmippProtExtractParticles, ProtImportParticles
from protocol_cl2d_align import XmippProtCL2DAlign
from protocol_cl2d import XmippProtCL2D
from protocol_kerdensom import XmippProtKerdensom
from protocol_rotational_spectra import XmippProtRotSpectra
from protocol_create_mask import XmippProtCreateMask3D
from convert import writeSetOfMicrographs, writeSetOfParticles, writeSetOfClasses2D, writeSetOfCoordinates, writeSetOfCTFs, locationToXmipp, \
                    writeSetOfClasses3D
from os.path import dirname, join
from pyworkflow.utils.path import makePath


import xmipp


class XmippViewer(Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [Image, SetOfImages, SetOfCoordinates, SetOfClasses2D, SetOfClasses3D,
                ProtImportMicrographs, XmippProtPreprocessMicrographs, ProtCTFMicrographs,
                ProtParticlePicking, ProtImportParticles, XmippProtExtractParticles, ProtUserSubSet,
                ProtAlign, ProtProcessParticles, XmippProtKerdensom, XmippProtRotSpectra,  XmippProtCreateMask3D,
                SetOfClasses2D, SetOfCTF, NormalModes]
    
    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):
        cls = type(obj)
        
        if issubclass(cls, Image):
            print "visualizing Image"
            fn = locationToXmipp(*obj.getLocation())
            runShowJ(fn)
            
        elif issubclass(cls, NormalModes):
            runShowJ(obj.getFileName())         
              
        elif issubclass(cls, SetOfMicrographs):
            
            mdFn = getattr(obj, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:
                fn = self._getTmpPath(obj.getName() + '_micrographs.xmd')
                writeSetOfMicrographs(obj, fn)
                
            extra = ''
            if obj.hasCTF():
                extra = ' --mode metadata --render first'
            #runShowJ(fn, extraParams=extra)
            runScipionShowJ(fn, "Set Of Micrographs", obj)    
        
        elif issubclass(cls, SetOfCoordinates):
            obj_mics = obj.getMicrographs()
            mdFn = getattr(obj_mics, '_xmippMd', None)
            
            if mdFn:
                fn = mdFn.get()
            else:
                fn = self._getTmpPath(obj_mics.getName() + '_micrographs.xmd')
                writeSetOfMicrographs(obj_mics, fn)
                
            extraFn = getattr(obj, '_xmippMd', None)
            if extraFn:
                extraDir = dirname(extraFn.get())
            else:
                extraDir = self._getTmpPath(obj.getName()) # TODO: CHECK to create an extra for the coordinates obj
                makePath(extraDir)
                writeSetOfCoordinates(extraDir, obj)            
                
            runParticlePicker(fn, extraDir, extraParams='readonly')
        
        elif issubclass(cls, SetOfParticles) or issubclass(cls, SetOfVolumes):
            mdFn = getattr(obj, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:
                fn = self._getTmpPath(obj.getName() + '_images.xmd')
                #Set hasCTF to False to avoid problems
                writeSetOfParticles(obj, fn, self._getTmpPath())
            if issubclass(cls, SetOfParticles):
                runScipionShowJ(fn, "Set Of Particles", obj)  
            else:
                runShowJ(fn)
            md = xmipp.MetaData(fn) 
            #print "MD=%s" % obj.outputParticles.getFileName()
            if md.containsLabel(xmipp.MDL_ZSCORE):
                print "MD contains ZSCORE"
                from plotter import XmippPlotter
                xplotter = XmippPlotter(windowTitle="Zscore particles sorting")
                xplotter.createSubPlot("Particle sorting", "Particle number", "Zscore")
                xplotter.plotMd(md, False, mdLabelY=xmipp.MDL_ZSCORE)
                figFn = fn.replace('.xmd', '.png')
                xplotter.savefig(figFn)
                xplotter.show()
        
        elif issubclass(cls, SetOfClasses2D):
            mdFn = getattr(obj, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:
                fn = self._getTmpPath(obj.getName() + '_classes.xmd')
                writeSetOfClasses2D(obj, fn, self._getTmpPath())
            runShowJ(fn, extraParams=args.get('extraParams', ''))  
        elif issubclass(cls, SetOfClasses3D):
            mdFn = getattr(obj, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:
                fn = self._getTmpPath(obj.getName() + '_classes.xmd')
                writeSetOfClasses3D(obj, fn, self._getTmpPath())
            runShowJ("classes@"+fn, extraParams=args.get('extraParams', ''))  
        elif issubclass(cls, SetOfCTF):
            mdFn = getattr(obj, '_xmippMd', None)
            if mdFn:
                fn = mdFn.get()
            else:
                fn = self._getTmpPath(obj.getName() + '_ctfs.xmd')
                writeSetOfCTFs(obj, fn)
            runShowJ(fn, extraParams=' --mode metadata --render first')  
        
        elif (issubclass(cls, ProtImportMicrographs) or
              issubclass(cls, XmippProtPreprocessMicrographs)):
            self.visualize(obj.outputMicrographs)
        
        elif issubclass(cls, ProtParticlePicking):
            self.visualize(obj.outputCoordinates)
        
        elif (issubclass(cls, ProtImportParticles) or
              issubclass(cls, XmippProtExtractParticles)):
            self.visualize(obj.outputParticles)
            # If Zscore on output images plot Zscore particle sorting
        elif (issubclass(cls, ProtUserSubSet)):
            self.visualize(obj.outputset)
                          
        elif issubclass(cls, ProtAlign):
            self.visualize(obj.outputAverage)
            self.visualize(obj.outputParticles) 
        elif issubclass(cls, XmippProtRotSpectra):
            self.visualize(obj.outputClasses, extraParams='--mode rotspectra --columns %d' % obj.SomXdim.get())
        elif issubclass(cls, XmippProtCreateMask3D):
            self.visualize(obj.outputMask)
        elif issubclass(cls, XmippProtKerdensom):
            self.visualize(obj.outputClasses, extraParams='--columns %d' % obj.SomXdim.get())
        elif issubclass(cls, ProtCTFMicrographs):
            self.visualize(obj.outputCTF)
        elif issubclass(cls, ProtProcessParticles):
            self.visualize(obj.outputParticles)
        else:
            raise Exception('XmippViewer.visualize: can not visualize class: %s' % obj.getClassName())
        
    @classmethod
    def getView(self):
        return "viewerXmipp"
        
# ------------- Xmipp utilities function to launch Java applications ------------

def getArchitecture():
    import platform
    arch = platform.architecture()[0]
    for a in ['32', '64']:
        if a in arch:
            return a
    return 'NO_ARCH' 
    
def getJavaIJappArguments(memory, appName, appArgs):
    """ Build the command line arguments to launch 
    a Java application based on ImageJ. 
    memory: the amount of memory passed to the JVM.
    appName: the qualified name of the java application.
    appArgs: the arguments specific to the application.
    """ 
    if len(memory) == 0:
        memory = "1g"
        print "No memory size provided. Using default: " + memory
    imagej_home = getXmippPath("external", "imagej")
    lib = getXmippPath("lib")
    javaLib = getXmippPath('java', 'lib')
    plugins_dir = os.path.join(imagej_home, "plugins")
    arch = getArchitecture()
    args = "-Xmx%(memory)s -d%(arch)s -Djava.library.path=%(lib)s -Dplugins.dir=%(plugins_dir)s -cp %(imagej_home)s/*:%(javaLib)s/* %(appName)s %(appArgs)s" % locals()

    return args
    
def runJavaIJapp(memory, appName, args, batchMode=True):
    args = getJavaIJappArguments(memory, appName, args)
    runJob(None, "java", args, runInBackground=batchMode)
    
def runShowJ(inputFiles, memory="1g", extraParams=""):
    runJavaIJapp(memory, "'xmipp.viewer.Viewer'", "-i %s %s" % (inputFiles, extraParams), True)
    
def runScipionShowJ(inputFiles, set, obj, memory="1g", extraParams=""):
    SCIPION_PYTHON = os.environ["SCIPION_PYTHON"]
    PW_HOME = os.environ["PW_HOME"]
    
    script = "%s %s/pyworkflow/apps/%s" % (SCIPION_PYTHON, PW_HOME, "pw_create_image_subset.py") 

    runJavaIJapp(memory, "'xmipp.viewer.scipion.ScipionViewer'", "-i %s %s --command \"%s\" \"%s\" %s" % (inputFiles, extraParams, set, script, obj.strId()), True)

def runParticlePicker(inputMics, inputCoords, memory="1g", extraParams=""):
    runJavaIJapp(memory, "xmipp.viewer.particlepicker.training.Main", "%s %s %s" % (inputMics, inputCoords, extraParams), True)



