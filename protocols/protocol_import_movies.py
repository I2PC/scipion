#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of movies
# Author: Carlos Oscar Sorzano, July 2011

from glob import glob
from protlib_base import *
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_MOVIE, MDL_SAMPLINGRATE, MDL_CTF_VOLTAGE, \
    MDL_CTF_CS, MDL_CTF_SAMPLING_RATE, MDL_MAGNIFICATION, MDL_ENABLED, MDL_REF, checkImageFileSize, \
    MD_APPEND, MD_OVERWRITE, getImageSize
from protlib_filesystem import  replaceFilenameExt
from protlib_xmipp import RowMetaData


class ProtImportMovies(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.import_movies.name, scriptname, project)
        self.Import = "from protocol_import_movies import *"
        self.PatternMovies = join(self.DirMovies, self.ExtMovies)
        self.MicrographsMd = self.getFilename('micrographs')
            
    def defineSteps(self):
        # Create microscope
        self.insertCreateMicroscope()
        # Decide name after preprocessing
        movies = self.getMovies()
        if self.CopyMovies:
            self.insertStep("createDir",verifyfiles=[self.ExtraDir],path=self.ExtraDir)
            func = self.insertCopyMicrograph
            funcOutput = lambda i: os.path.join(self.ExtraDir,os.path.basename(i))
        else:
            func = lambda i, o: i #Do nothing
            funcOutput = lambda i: i 
        filenameDict = {}
        for m in movies:
            output = funcOutput(m)
            filenameDict[m] = output
            func(m, output)
        # Insert step for result metadatas creation     
        self.insertCreateResults(filenameDict)
            
    def validate(self):
        errors = []
        # Check that there are any micrograph to process
        micrographList = self.getMovies()
        if len(micrographList) == 0:
            errors.append("There are no movies to process in " + self.PatternMovies)
        else:
            for micrograph in micrographList:
                try:
                    if not checkImageFileSize(micrograph):
                        errors.append(micrograph+" seems to be corrupted")
                except Exception:
                    errors.append(micrograph+" seems to be corrupted")

        return errors

    def summary(self):
        message = []
        message.append("Import of <%d> movies from [%s]" % (len(self.getMovies()), self.DirMovies))
        
        fnAcquisition = self.getFilename('acquisition')
        if os.path.exists(fnAcquisition):
            md = MetaData(fnAcquisition)
            sampling = md.getValue(MDL_SAMPLINGRATE, md.firstObject())
            message.append("\nSampling rate = <%f> A/pix" % sampling)
        return message
    
    def getMovies(self):
        ''' Return a list with movies in WorkingDir'''
        return glob(self.PatternMovies)
    
    def visualize(self):
        summaryFile = self.MicrographsMd
        if os.path.exists(summaryFile):
            from protlib_utils import runShowJ
            runShowJ(summaryFile)

    def insertCreateMicroscope(self):    
        if self.SamplingRateMode == "From image":
            self.AngPix = float(self.SamplingRate)
        else:
            scannedPixelSize=float(self.ScannedPixelSize)
            self.AngPix = (10000. * scannedPixelSize) / self.Magnification
        fnOut = self.getFilename('microscope')
        self.insertStep("createMicroscope", verifyfiles=[fnOut], fnOut=fnOut, Voltage=self.Voltage,
                        SphericalAberration=self.SphericalAberration,SamplingRate=self.AngPix, SamplingRateMode=self.SamplingRateMode,
                        Magnification=self.Magnification)
            
    def insertCreateResults(self, filenameDict):
        micFn = self.MicrographsMd
        self.insertStep('createResults', 
                        verifyfiles=[micFn], 
                        WorkingDir=self.WorkingDir, 
                        FilenameDict=filenameDict, 
                        MicrographFn=micFn, 
                        PixelSize=self.AngPix)
        
    def insertCopyMicrograph(self, inputMic, outputMic):
        self.insertStep('copyFile', source=inputMic, dest=outputMic)
        if inputMic.endswith('.raw'):
            self.insertStep('copyFile', source=replaceFilenameExt(inputMic, '.inf'), 
                            dest=replaceFilenameExt(outputMic, '.inf'))
                

#def getWdFile(wd, key):
#    return getProtocolFilename(key, WorkingDir=wd)


def createMicroscope(log, fnOut, Voltage, SphericalAberration, SamplingRate, SamplingRateMode, Magnification):
    md = RowMetaData()
    md.setValue(MDL_CTF_VOLTAGE, float(Voltage))    
    md.setValue(MDL_CTF_CS, float(SphericalAberration))    
    md.setValue(MDL_CTF_SAMPLING_RATE, float(SamplingRate))
    if SamplingRateMode != "From image":
        md.setValue(MDL_MAGNIFICATION, float(Magnification))
    md.write(fnOut)    

def createResults(log, WorkingDir, FilenameDict, MicrographFn, PixelSize):
    ''' Create a metadata movies.xmd with all movies
    '''
    
    movies = FilenameDict.values()
    movies.sort()
    md = MetaData()
    mdClass = MetaData()
    [xdim,ydim,zdim,ndim]=getImageSize(movies[0])
     
    md.write("classes@"+MicrographFn,MD_OVERWRITE)
    for m in movies:
        id = md.addObject()
        md.setValue(MDL_REF, int(id), id)        
        md.setValue(MDL_MICROGRAPH, ('%05d@%s'%(1,m)), id)
        md.setValue(MDL_ENABLED, 1, id)
        md.setValue(MDL_MICROGRAPH_MOVIE, m, id)        
        
        mdClass.clear()
        for num in range(1,ndim+1):
            idClass = mdClass.addObject()
            mdClass.setValue(MDL_MICROGRAPH, ('%05d@%s'%(num,m)), idClass)
            mdClass.setValue(MDL_MICROGRAPH_MOVIE, m, idClass)
            mdClass.setValue(MDL_REF, int(id), idClass)
            mdClass.write('class%05d_movies@%s'%(id,MicrographFn),MD_APPEND)                 
        
    md.write("classes@"+MicrographFn,MD_APPEND)
    mdAcquisition = MetaData()
    mdAcquisition.setValue(MDL_SAMPLINGRATE,float(PixelSize),mdAcquisition.addObject())
    mdAcquisition.write(os.path.join(WorkingDir,"acquisition_info.xmd"))
    
