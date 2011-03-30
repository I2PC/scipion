#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs: 
#  - downsampling
#  - power spectral density (PSD) and CTF estimation on the micrograph
#
# For each micrograph given, this script will perform 
# the requested operations below.
# For each micrograph a subdirectory will be created
#
# Author: Sjors Scheres, March 2007
#         Roberto Marabini (mpi extension)
#         Carlos Oscar Sorzano, November 2010
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Working subdirectory:
WorkingDir='Preprocessing'
# {dir} Directory name from where to process all scanned micrographs
DirMicrographs='Micrographs'
# Which files in this directory to process
""" This is typically *.tif or *.ser, but may also be *.mrc, *.spi 
    (see the expert options)
    Note that any wildcard is possible, e.g. *3[1,2].tif
"""
ExtMicrographs='*.mrc'
# Rootname for these micrographs
""" Several files will be created called <rootname>_...
"""
RootName='all'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir ='/media/usbdisk/Experiments/TestProtocols/Protocol_Small'
#------------------------------------------------------------------------------------------------
# {section} Preprocess
#------------------------------------------------------------------------------------------------
# Do proceprocess
# Perform preprocessing? 
DoPreprocess=False
# Crop borders
""" Crop a given amount of pixels from each border.
    Set this option to -1 for not applying it."""
Crop=-1
# Remove bad pixels
""" Values will be thresholded to this multiple of standard deviations. Typical
    values are about 5, i.e., pixel values beyond 5 times the standard deviation will be
    substituted by the local median. Set this option to -1 for not applying it."""
Stddev=-1
# Downsampling factor 
""" Set to 1 for no downsampling. Non-integer downsample factors are possible with
    the Fourier kernel. """
Down=1
# {expert}{list}|Fourier|Rectangle|Sinc| Which method to use for downsampling?
""" Fourier is theoretically the best option, but it may take more memory than your machine
    can handle. Then, Rectangle is the fastest, but much less accurate. Sinc is reasonably
    accurate, but painfully slow...
"""
DownKernel='Fourier'
#------------------------------------------------------------------------------------------------
# {section} CTF estimation
#------------------------------------------------------------------------------------------------
# Perform CTF estimation?
DoCtfEstimate=True
# Microscope voltage (in kV)
Voltage=200
# Spherical aberration
SphericalAberration=2.26
# Magnification rate
Magnification=50000
# Scanned pixel size (in um)
ScannedPixelSize=7
# Amplitude Contrast
AmplitudeContrast=0.07
# {expert} Only perform power spectral density estimation?
""" Skip the CTF estimation part, and only estimate the PSD
"""
OnlyEstimatePSD=False
# {expert} Lowest resolution for CTF estimation
""" Give a value in digital frequency (i.e. between 0.0 and 0.5)
    This cut-off prevents the typically peak at the center of the PSD to interfere with CTF estimation.  
    The default value is 0.05, but for micrographs with a very fine sampling this may be lowered towards 0.0
"""
LowResolCutoff=0.05
# {expert} Highest resolution for CTF estimation
""" Give a value in digital frequency (i.e. between 0.0 and 0.5)
    This cut-off prevents high-resolution terms where only noise exists to interfere with CTF estimation.  
    The default value is 0.35, but it should be increased for micrographs with signals extending beyond this value
"""
HighResolCutoff=0.35
# {expert} Minimum defocus to search (in Ang.)
""" Minimum defocus value (in Angstrom) to include in defocus search
"""
MinFocus=5000
# {expert} Maximum defocus to search (in Ang.)
""" Maximum defocus value (in Angstrom) to include in defocus search
"""
MaxFocus=100000
# {expert} Window size for CTFFIND
WinSize=256
# {expert} Defocus step for CTFFIND (in Ang.)
""" Step size for defocus search (in Angstrom)
"""
StepFocus=500

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# distributed-memory parallelization (MPI)?
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
DoParallel=True

# Number of MPI processes to use:
NumberOfMpiProcesses=3

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. The following values are available: 
"""
SystemFlavour=''

#------------------------------------------------------------------------------------------------
# {hidden} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_preprocess_micrographs.py'

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
# Taken from Python 2.6
import posixpath
def relpath(path, start=posixpath.curdir):
    """Return a relative version of a path"""
    if not path:
        raise ValueError("no path specified")
    start_list=posixpath.abspath(start).split(posixpath.sep)
    path_list=posixpath.abspath(path).split(posixpath.sep)
    # Work out how much of the filepath is shared by start and path.
    i=len(posixpath.commonprefix([start_list, path_list]))
    rel_list=[posixpath.pardir] * (len(start_list) - i) + path_list[i:]
    if not rel_list:
        return curdir
    return posixpath.join(*rel_list)

def stepPerformed(step,filename):
    import re
    f = open(filename, 'r')
    lines=f.readlines()
    f.close()
    expr = re.compile(step)
    return len(filter(expr.search,lines))>0

import glob, os, shutil, sys
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
sys.path.append(scriptdir)

class preprocess_A_class:
    def saveAndCompareParameters(self, listOfParameters):
        fnOut=self.WorkingDir + "/protocolParameters.txt"
        linesNew=[];
        for prm in listOfParameters:
            eval("linesNew.append('"+prm +"='+str("+prm+")+'\\n')")
        if os.path.exists(fnOut):
            f = open(fnOut, 'r')
            linesOld=f.readlines()
            f.close()
            same=True;
            if len(linesOld)==len(linesNew):
                for i in range(len(linesNew)):
                    if not linesNew[i]==linesOld[i]:
                        same=False
                        break;
            else:
                same=False
            if not same:
                print("Deleting")
                self.log.info("Deleting working directory since it is run with different parameters")
                shutil.rmtree(self.WorkingDir)
                os.makedirs(self.WorkingDir)
        f = open(fnOut, 'w')
        f.writelines(linesNew)
        f.close()
    
    def __init__(self,
                 WorkingDir,
                 DirMicrographs,
                 ExtMicrographs,
                 RootName,
                 ProjectDir,
                 DoPreprocess,
                 Crop,
                 Stddev,
                 Down,
                 DownKernel,
                 DoCtfEstimate,
                 Voltage,
                 SphericalAberration,
                 Magnification,
                 ScannedPixelSize,
                 AmplitudeContrast,
                 OnlyEstimatePSD,
                 LowResolCutoff,
                 HighResolCutoff,
                 MinFocus,
                 MaxFocus,
                 WinSize,
                 StepFocus,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour
                 ):
        
        import log

        self.WorkingDir=os.path.abspath(WorkingDir)
        self.DirMicrographs=os.path.abspath(DirMicrographs)
        self.ExtMicrographs=ExtMicrographs.strip()
        self.ProjectDir=os.path.abspath(ProjectDir)
        self.LogDir="Logs"
        self.RootName=RootName.strip()
        self.DoPreprocess=DoPreprocess
        self.Crop=Crop
        self.Stddev=Stddev
        if (float(Down) == int(Down)):
            self.Down=int(Down)
        else:
            self.Down=float(Down)
        self.DownKernel=DownKernel
        self.DoCtfEstimate=DoCtfEstimate
        self.Voltage=Voltage
        self.SphericalAberration=SphericalAberration
        self.Magnification=Magnification
        self.ScannedPixelSize=ScannedPixelSize
        self.AmplitudeContrast=AmplitudeContrast
        self.OnlyEstimatePSD=OnlyEstimatePSD
        self.LowResolCutoff=LowResolCutoff
        self.HighResolCutoff=HighResolCutoff
        self.MinFocus=MinFocus
        self.MaxFocus=MaxFocus
        self.WinSize=WinSize
        self.StepFocus=StepFocus
        self._MySystemFlavour=SystemFlavour
        self._DoParallel=DoParallel
        self._MyNumberOfMpiProcesses=NumberOfMpiProcesses

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     self.LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)

        # Check ctffind executable
        import which
        self.CtffindExec=which.which('ctffind3.exe')
        self.DoCtffind=self.CtffindExec!=''

        # Delete working directory if exists, make a new one
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Save parameters and compare to possible previous runs
        self.saveAndCompareParameters([
                 "DirMicrographs",
                 "ExtMicrographs",
                 "RootName",
                 "DoPreprocess",
                 "Crop",
                 "Stddev",
                 "Down",
                 "DownKernel",
                 "DoCtfEstimate",
                 "Voltage",
                 "SphericalAberration",
                 "Magnification",
                 "ScannedPixelSize",
                 "AmplitudeContrast",
                 "OnlyEstimatePSD",
                 "LowResolCutoff",
                 "HighResolCutoff",
                 "MinFocus",
                 "MaxFocus",
                 "WinSize",
                 "StepFocus"]);
        
        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                   os.path.abspath(self.WorkingDir))
        xmpi_run_file=self.WorkingDir + "/preprocess_micrographs" 
        self.xmpi_run_file=os.path.abspath(xmpi_run_file)

        # Execute protocol in the working directory
        os.chdir(self.WorkingDir)
        self.process_all_micrographs(self.xmpi_run_file)

    def process_all_micrographs(self, xmpi_run_file):
        import time
        print '*********************************************************************'
        print '*  Processing the following micrographs: '
        for self.filename in glob.glob(self.DirMicrographs + '/' + self.ExtMicrographs):
            (self.filepath, self.name)=os.path.split(self.filename)
            print '*  ' + self.name

        self.SFinputparams=[]
        self.SFmicrograph=[]
        self.SFctffindmrc=[]
        self.SFquadrant=[]
        self.SFctffind=[]
        self.SFshort=[]
        self.SFhalf=[]
        self.SFctf=[]
        self.SFpsd=[]

        fh_mpi=os.open(xmpi_run_file + '.sh', os.O_WRONLY | os.O_TRUNC | os.O_CREAT, 0700)
        # Preprocessing
        for filename in glob.glob(self.DirMicrographs + '/' + self.ExtMicrographs):

            # Get the shortname and extension
            (filepath, micrographName)=os.path.split(filename)
            (shortname, extension)=os.path.splitext(micrographName)
            self.SFshort.append(shortname)
            
            # Create directory for this micrograph
            if not os.path.exists(shortname):
                os.makedirs(shortname)
                fh=open(shortname + "/status.txt", "w")
                fh.write("Process started at " + time.asctime() + "\n")
                fh.write("Step 0: Directory created at " + time.asctime() + "\n")
                fh.close()
            else:
                fh=open(shortname + "/status.txt", "a")
                fh.write("Process started at " + time.asctime() + "\n")
                fh.close()

            # Preprocess
            finalName=self.perform_preprocessing(filename, fh_mpi)
            self.SFmicrograph.append(finalName)

        # Stop here until preprocessing is done
        if self._DoParallel:
            os.write(fh_mpi, "MPI_Barrier\n");

        # Estimate CTF
        idx=0;
        for filename in self.SFmicrograph:
            if self.DoCtfEstimate:
                shortname=self.SFshort[idx]
                names=self.perform_ctfestimate_xmipp(shortname, filename, fh_mpi)
                self.SFpsd.append(names[0])
                self.SFinputparams.append(names[1])
                if not self.OnlyEstimatePSD:                    
                    self.SFctf.append(names[2])
                    self.SFhalf.append(names[3])
                    self.SFquadrant.append(names[4])
                if self.DoCtffind:
                    self.perform_ctfestimate_ctffind(shortname, filename, fh_mpi)
            idx += 1

        # Launch Preprocessing and calculation of the CTF
        os.close(fh_mpi)
        self.launchCommandFile(xmpi_run_file + '.sh')
        
        # Pickup results from CTFFIND
        if self.DoCtfEstimate and self.DoCtffind:
            idx=0;
            for filename in self.SFmicrograph:
                shortname=self.SFshort[idx]
                ctfparam=self.SFctf[idx]
                results=self.pickup_ctfestimate_ctffind(shortname, filename)
                self.SFctffind.append(results[0])
                self.SFctffindmrc.append(results[1])
                idx += 1
        
        # Write the different selfiles
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
        sys.path.append(scriptdir)
        import xmipp
        MD=xmipp.MetaData()
        idx=0
        for filename in self.SFmicrograph:
            fh=open(self.SFshort[idx] + "/status.txt", "a")
            fh.write("Step F: Finished on " + time.asctime() + "\n")
            fh.close()
            objId=MD.addObject()
            MD.setValue(xmipp.MDL_IMAGE, filename,objId)
            if self.DoCtfEstimate:
                MD.setValue(xmipp.MDL_PSD,   self.SFpsd[idx],objId)
            if len(self.SFinputparams) > 0:
                MD.setValue(xmipp.MDL_CTFINPUTPARAMS, self.SFinputparams[idx],objId)
            if len(self.SFctf) > 0:
                MD.setValue(xmipp.MDL_CTFMODEL,          self.SFctf[idx],objId)
                MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE1, self.SFhalf[idx],objId)
                MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE2, self.SFquadrant[idx],objId)
            if len(self.SFctffind) > 0:
                MD.setValue(xmipp.MDL_CTFMODEL2, self.SFctffind[idx],objId)
                MD.setValue(xmipp.MDL_ASSOCIATED_IMAGE3, self.SFctffindmrc[idx],objId)
            idx += 1
        MD.sort(xmipp.MDL_IMAGE);
        MD.write(self.WorkingDir + "/" + self.RootName + "_micrographs.sel")

        # CTF Quality control
        if self.DoCtfEstimate:
            command="xmipp_ctf_sort_psds -i " + self.RootName + "_micrographs.sel\n"
            self.log.info(command)     
            os.system(command)     
        
        message=" Done pre-processing of micrographs"
        print '* ', message
        print '*********************************************************************'
        self.log.info(message)

    def launchCommandFile(self, commandFile):
        import launch_job, log
        log.cat(self.log, commandFile)
        if self._DoParallel:
            command=' -i ' + commandFile
            launch_job.launch_job("xmipp_run", command, self.log, True,
                  self._MyNumberOfMpiProcesses, 1, self._MySystemFlavour)
        else:
            self.log.info(commandFile)     
            os.system(commandFile)     

    def perform_preprocessing(self, filename, fh_mpi):
        # Decide name after preprocessing
        filename=relpath(filename)
        (filepath, micrographName)=os.path.split(filename)
        (shortname, extension)=os.path.splitext(micrographName)
        finalname=shortname + '/down' + str(self.Down) + '_' + shortname
        if not self.Stddev == -1 or not self.Crop == -1 or not self.Down == 1:
            if not self.Down == 1:
                finalname += ".spi"
            else:
                finalname += ".raw"
        else:
            finalname += extension
            if not os.path.exists(finalname):
                command='ln -s ' + relpath(filename, shortname) + ' ' + finalname + "; "
                if micrographName.endswith(".raw"):
                    command+='ln -s ' + relpath(filename, shortname) + '.inf ' + finalname + ".inf; "
                command += "if [ -e " + finalname + ' ]; then ' + \
                            'echo "Step 1: Preprocessed image created " `date` >> ' + shortname + "/status.txt; " + \
                          "fi"
                command += "\n"
                os.write(fh_mpi, command)
                return finalname
        if not self.DoPreprocess:
            return finalname
        
        # Check if the preprocessing has already been done
        if stepPerformed("Step 1",shortname + "/status.txt"):
            return finalname
        
        # Crop
        iname=filename
        command="";
        if not self.Crop == -1:
            command += "xmipp_window -i " + iname + " -o " + finalname + " -crop " + str(self.Crop) + " ; "
            iname=finalname
        
        # Remove bad pixels
        if not self.Stddev == -1:
            command += "xmipp_filter -i " + iname + " -desvPixels " + str(self.Stddev)
            if not iname == finalname:
                command += " -o " + finalname
                iname=finalname
            command += " ; "
        
        # Downsample
        if not self.Down == 1:
            command += "xmipp_micrograph_downsample -i " + iname + " -o " + shortname + "/tmp.spi " + \
                     "-datatype float "
            if (self.DownKernel == 'Fourier'):
                scale=1. / self.Down
                command += ' -fourier ' + str(scale)
            elif (self.DownKernel == 'Sinc'):
                command += ' -Xstep ' + str(self.Down) + ' -kernel sinc 0.02 0.1'
            elif (self.DownKernel == 'Rectangle'):
                command += ' -Xstep ' + str(self.Down) + ' -kernel rectangle ' + str(self.Down) + ' ' + str(self.Down)
            command += " ; rm -f " + finalname + " ; mv -i " + shortname + "/tmp.spi " + finalname + " ; "
        
        # Postprocessing
        command += "if [ -e " + finalname + ' ]; then ' + \
                     'echo "Step 1: Preprocessed image created " `date` >> ' + shortname + "/status.txt; " + \
                   "fi"
        
        # Write the preprocessing command
        if not command == "":
            command += "\n"
            os.write(fh_mpi, command)
        return finalname

    def perform_ctfestimate_xmipp(self, shortname, filename, fh_mpi):
        (filepath, micrographName)=os.path.split(filename)
        (fnRoot, extension)=os.path.splitext(micrographName)
        pname=shortname + '/' + fnRoot + '_estimate_ctf_input.param'
        retval=[]
        
        retval.append(shortname + '/' + fnRoot + "_Periodogramavg.psd")
        retval.append(pname)

        # prepare parameter file
        paramlist=[]
        paramlist.append('image= ' + filename + '\n')
        paramlist.append('micrograph_averaging= yes \n')
        paramlist.append('N_horizontal= 512 \n')
        paramlist.append('periodogram= yes \n')
        if not self.OnlyEstimatePSD:
            AngPix=(10000. * self.ScannedPixelSize * self.Down) / self.Magnification
            paramlist.append('voltage= ' + str(self.Voltage) + '\n')
            paramlist.append('spherical_aberration= ' + str(self.SphericalAberration) + '\n')
            paramlist.append('sampling_rate= ' + str(AngPix) + '\n')
            paramlist.append('particle_horizontal= 256 \n')
            paramlist.append('Q0= -' + str(self.AmplitudeContrast) + '\n')
            paramlist.append('min_freq= ' + str(self.LowResolCutoff) + '\n')
            paramlist.append('max_freq= ' + str(self.HighResolCutoff) + '\n')
            retval.append(shortname + '/' + fnRoot + "_Periodogramavg.ctfparam")
            retval.append(shortname + '/' + fnRoot + "_Periodogramavg_ctfmodel_halfplane.xmp")
            retval.append(shortname + '/' + fnRoot + "_Periodogramavg_ctfmodel_quadrant.xmp")

        # Check if the preprocessing has already been done
        if stepPerformed("Step 2",shortname + "/status.txt"):
            return retval

        # Perform CTF estimation
        fh=open(pname, "w")
        fh.writelines(paramlist)
        fh.close()
        command='if grep -q "Step 1" ' + shortname + '/status.txt ; then ' + \
            'xmipp_ctf_estimate_from_micrograph -i ' + pname + ' ; ' + \
            'if [ -e ' + shortname + '/' + fnRoot + '_Periodogramavg.ctfparam ] ; then ' + \
               'echo "Step 2: CTF estimated with Xmipp " `date` >> ' + shortname + '/status.txt;' + \
            ' fi; fi'
        command += '\n'
        os.write(fh_mpi, command);
        return retval

    def perform_ctfestimate_ctffind(self, shortname, filename, fh_mpi):
        # Check if the preprocessing has already been done
        if stepPerformed("Step 3",shortname + "/status.txt"):
            return

        # The new line is different if we are running in parallel or not
        theNewLine='\n'
        if(self._DoParallel):
            theNewLine='MPI_NEWLINE'

        # Convert image to MRC
        command='if grep -q "Step 1" ' + shortname + '/status.txt' + theNewLine + \
            'then' + theNewLine + '(' + theNewLine
        command += 'xmipp_image_convert -i ' + filename + ' -o ' + shortname + '/tmp.mrc -v 0; '

        # Prepare parameters for CTFTILT
        AngPix=(10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        (filepath, micrographName)=os.path.split(filename)
        (fnRoot, extension)=os.path.splitext(micrographName)
        command += "export NATIVEMTZ=kk" + theNewLine
        command += self.CtffindExec + '  << eof > ' + shortname + '/' + fnRoot + '_ctffind.log' + theNewLine
        command += shortname + '/tmp.mrc' + theNewLine
        command += shortname + '/' + fnRoot + '_ctffind_spectrum.mrc' + theNewLine
        command += str(self.SphericalAberration) + ',' + \
                  str(self.Voltage) + ',' + \
                  str(self.AmplitudeContrast) + ',' + \
                  str(self.Magnification) + ',' + \
                  str(self.ScannedPixelSize * self.Down) + theNewLine
        command += str(self.WinSize) + ',' + \
                  str(AngPix / self.LowResolCutoff) + ',' + \
                  str(AngPix / self.HighResolCutoff) + ',' + \
                  str(self.MinFocus) + ',' + \
                  str(self.MaxFocus) + ',' + \
                  str(self.StepFocus) + theNewLine
        command += 'eof' + theNewLine + ')' + theNewLine + 'fi\n'
        os.write(fh_mpi, command)

    def pickup_ctfestimate_ctffind(self, shortname, filename):
        import xmipp, time
        (filepath, micrographName)=os.path.split(filename)
        (fnRoot, extension)=os.path.splitext(micrographName)
        fnOut=shortname + '/' + fnRoot + '_ctffind.ctfparam'

        # Check if the preprocessing has already been done
        if stepPerformed("Step 3",shortname + "/status.txt"):
            retval=[];
            retval.append(fnOut);
            retval.append(shortname + '/' + fnRoot + '_ctffind_spectrum.mrc');
            return retval

        # Remove temporary files
        if os.path.exists(shortname + '/tmp.mrc'):
            os.remove(shortname + '/tmp.mrc')

        # Pick values from ctffind
        if not os.path.exists(shortname + '/' + fnRoot + '_ctffind.log'):
            retval=[]
            retval.append('NA')
            retval.append('NA')
            return retval
        
        # Effectively pickup results
        fh=open(shortname + '/' + fnRoot + '_ctffind.log', 'r')
        lines=fh.readlines()
        fh.close()
        DF1=0.
        DF2=0.
        Angle=0.
        found=False
        for i in range(len(lines)):
            if not (lines[i].find('Final Values') == -1):
                words=lines[i].split()
                DF1=float(words[0])
                DF2=float(words[1])
                Angle=float(words[2])
                found=True
                break

        if not found:
            retval=[]
            retval.append('NA')
            retval.append('NA')
            return retval
        
        retval=[];
        retval.append(fnOut);
        retval.append(shortname + '/' + fnRoot + '_ctffind_spectrum.mrc');

        # Generate Xmipp .ctfparam file:
        MD=xmipp.MetaData()
        MD.setColumnFormat(False)
        AngPix=(10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        objId=MD.addObject()
        MD.setValue(xmipp.MDL_CTF_SAMPLING_RATE, AngPix, objId)
        MD.setValue(xmipp.MDL_CTF_VOLTAGE,       float(self.Voltage), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUSU,      float(-DF2), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUSV,      float(-DF1), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, float(Angle), objId)
        MD.setValue(xmipp.MDL_CTF_CS,            float(self.SphericalAberration), objId)
        MD.setValue(xmipp.MDL_CTF_Q0,            float(-self.AmplitudeContrast), objId)
        MD.setValue(xmipp.MDL_CTF_K,             1.0, objId)
        MD.write(fnOut)

        fh=open(shortname + "/status.txt", "a")
        fh.write("Step 3: CTF estimated with CTFFind " + time.asctime() + "\n")
        fh.close()

        return retval

# Preconditions
def preconditions(gui):
    retval=True
    # Check if there is workingdir
    if WorkingDir == "":
        message="No working directory given"
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    # Check that there are any micrograph to process
    listOfMicrographs=glob.glob(DirMicrographs + '/' + ExtMicrographs)
    if len(listOfMicrographs) == 0:
        message="There are no micrographs to process in " + DirMicrographs + '/' + ExtMicrographs
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    return retval

#        
# Main
#     
if __name__ == '__main__':
    # create preprocess_A_class object
    if not preconditions(False):
        sys.exit(1)
    preprocessA=preprocess_A_class(
                 WorkingDir,
                 DirMicrographs,
                 ExtMicrographs,
                 RootName,
                 ProjectDir,
                 DoPreprocess,
                 Crop,
                 Stddev,
                 Down,
                 DownKernel,
                 DoCtfEstimate,
                 Voltage,
                 SphericalAberration,
                 Magnification,
                 ScannedPixelSize,
                 AmplitudeContrast,
                 OnlyEstimatePSD,
                 LowResolCutoff,
                 HighResolCutoff,
                 MinFocus,
                 MaxFocus,
                 WinSize,
                 StepFocus,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour)
