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
#         Carlos Oscar
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Working subdirectory:
WorkingDir='Preprocessing'
# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=False
# {dir} Directory name from where to process all scanned micrographs
DirMicrographs='Micrographs'
# Which files in this directory to process
""" This is typically *.tif or *.ser, but may also be *.mrc, *.spi 
    (see the expert options)
    Note that any wildcard is possible, e.g. *3[1,2].tif
"""
ExtMicrographs='*.tif'
# Rootname for these micrographs
""" Several files will be created called <rootname>_...
"""
RootName='all'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir='/media/usbdisk/Experiments/TestProtocols'
# {expert} Directory name for logfiles:
LogDir='Logs'
#------------------------------------------------------------------------------------------------
# {section} Preprocess
#------------------------------------------------------------------------------------------------
# Do proceprocess
# Perform preprocessing? 
DoPreprocess=True
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
#------------------------------------------------------------------------------------------------
# {section} CTFFIND
#------------------------------------------------------------------------------------------------
# Use N. Grigorieffs CTFFIND beside Xmipp?
# {file} Location of the CTFFIND executable
""" If available, CTFFIND will be used to double check the validity of the defoci """
CtffindExec='/media/usbdisk/OtherPackages/ctffind/ctf/ctffind3.exe'
# {expert} Window size
WinSize=256
# {expert} Minimum resolution (in norm. freq.)
""" Lowest resolution to include in CTF estimation. Valid range 0-0.5
"""
MinResCTF=0.02
# {expert} Maximum resolution factor
""" highest resolution to include in CTF estimation. Valid range 0-0.5
"""
MaxResCTF=0.35
# {expert} Defocus step (in Ang.)
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
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_preprocess_micrographs.py'

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
import glob,os,shutil,sys
class preprocess_A_class:
    def __init__(self,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 DirMicrographs,
                 ExtMicrographs,
                 RootName,
                 ProjectDir,
                 LogDir,
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
                 CtffindExec,
                 WinSize,
                 MinResCTF,
                 MaxResCTF,
                 StepFocus,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour
                 ):
        
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir)
        import log

        self.WorkingDir=os.path.abspath(WorkingDir)
        self.DirMicrographs=os.path.abspath(DirMicrographs)
        self.ExtMicrographs=ExtMicrographs
        self.ProjectDir=os.path.abspath(ProjectDir)
        self.LogDir=LogDir
        self.RootName=RootName
        self.DoPreprocess=DoPreprocess
        self.Crop=Crop
        self.Stddev=Stddev
        if (float(Down)==int(Down)):
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
        self.DoCtffind=(not CtffindExec=="")
        self.WinSize=WinSize
        self.MinResCTF=MinResCTF
        self.MaxResCTF=MaxResCTF
        self.StepFocus=StepFocus
        self._MySystemFlavour=SystemFlavour
        self._DoParallel=DoParallel
        self._MyNumberOfMpiProcesses=NumberOfMpiProcesses

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)

        # Check ctffind executable
        if (self.DoCtffind):
            self.CtffindExec=os.path.abspath(CtffindExec)

        # Delete working directory if exists, make a new one
        if (DoDeleteWorkingDir): 
            if (self.WorkingDir==""):
                raise RuntimeError,"No working directory given"
            if os.path.exists(self.WorkingDir):
                shutil.rmtree(self.WorkingDir)
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                   os.path.abspath(self.WorkingDir))
        unique_number=os.getpid()
        xmpi_run_file=self.WorkingDir+"/"+'xmipp_preprocess_micrographs_' + str(unique_number) 
        self.xmpi_run_file=os.path.abspath(xmpi_run_file)

        # Execute protocol in the working directory
        os.chdir(self.WorkingDir)
        self.process_all_micrographs(self.xmpi_run_file)

    def process_all_micrographs(self,xmpi_run_file):
        print '*********************************************************************'
        print '*  Processing the following micrographs: '
        for self.filename in glob.glob(self.DirMicrographs+'/'+self.ExtMicrographs):
            (self.filepath, self.name) = os.path.split(self.filename)
            print '*  '+self.name

        self.SFinputparams = []
        self.SFmicrograph = []
        self.SFquadrant = []
        self.SFctffind = []
        self.SFshort = []
        self.SFhalf = []
        self.SFctf = []
        self.SFpsd = []

        fh_mpi  = os.open(self.xmpi_run_file+ '_1.sh',os.O_WRONLY|os.O_TRUNC|os.O_CREAT, 0700)
        # Preprocessing
        for filename in glob.glob(self.DirMicrographs+'/'+self.ExtMicrographs):

            # Get the shortname and extension
            (filepath, micrographName) = os.path.split(filename)
            (shortname,extension) = os.path.splitext(micrographName)
            self.SFshort.append(shortname)
            
            # Create directory for this micrograph
            if not os.path.exists(shortname):
                os.makedirs(shortname)

            # Preprocess
            finalName=self.perform_preprocessing(filename,fh_mpi)
            self.SFmicrograph.append(finalName)

        # Stop here until preprocessing is done
        if self._DoParallel:
            os.write(fh_mpi,"MPI_Barrier"+"\n");

        # Estimate CTF
        idx=0;
        for filename in self.SFmicrograph:
            if self.DoCtfEstimate:
                shortname=self.SFshort[idx]
                names=self.perform_ctfestimate_xmipp(shortname,filename,fh_mpi)
                self.SFpsd.append(names[0])
                self.SFinputparams.append(names[1])
                if not self.OnlyEstimatePSD:                    
                    self.SFctf.append(names[2])
                    self.SFhalf.append(names[3])
                    self.SFquadrant.append(names[4])
                if self.DoCtffind:
                    self.perform_ctfestimate_ctffind(shortname,filename,fh_mpi)
            idx+=1

        # Launch Preprocessing and calculation of the CTF
        os.close(fh_mpi)
        self.launchCommandFile(xmpi_run_file + '_1.sh')
        
        # Pickup results from CTFFIND
        fh_mpi  = os.open(self.xmpi_run_file+ '_2.sh',os.O_WRONLY|os.O_TRUNC|os.O_CREAT, 0700)
        idx=0;
        for filename in self.SFmicrograph:
            shortname=self.SFshort[idx]
            ctfparam=self.SFctf[idx]
            # Pickup result from CTFTILT
            if self.DoCtfEstimate and self.DoCtffind:
                fnCTFfind=self.pickup_ctfestimate_ctffind(shortname,filename,fh_mpi)
                self.SFctffind.append(fnCTFfind)

            os.write(fh_mpi,"\n");
            idx+=1
        
        # Write the different selfiles
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir)
        import XmippData
        MD=XmippData.MetaData()
        ptr=XmippData.stringP()
        idx=0
        for filename in self.SFmicrograph:
            MD.addObject()
            ptr.assign(filename)
            XmippData.setValueString(MD,XmippData.MDL_IMAGE,ptr)
            ptr.assign(self.SFpsd[idx])
            XmippData.setValueString(MD,XmippData.MDL_PSD,ptr)
            if len(self.SFinputparams)>0:
                ptr.assign(self.SFinputparams[idx])
                XmippData.setValueString(MD,XmippData.MDL_CTFINPUTPARAMS,ptr)
            if len(self.SFctf)>0:
                ptr.assign(self.SFctf[idx])
                XmippData.setValueString(MD,XmippData.MDL_CTFMODEL,ptr)
                ptr.assign(self.SFhalf[idx])
                XmippData.setValueString(MD,XmippData.MDL_ASSOCIATED_IMAGE1,ptr)
                ptr.assign(self.SFquadrant[idx])
                XmippData.setValueString(MD,XmippData.MDL_ASSOCIATED_IMAGE2,ptr)
            if len(self.SFctffind)>0:
                ptr.assign(self.SFctffind[idx])
                XmippData.setValueString(MD,XmippData.MDL_CTFMODEL2,ptr)
            idx+=1
        MD.write(XmippData.FileName(self.WorkingDir+"/"+self.RootName+"_micrographs.sel"))

        # CTF Quality control
        command="xmipp_psd_sort -i "+self.RootName+"_micrographs.sel\n"
        os.write(fh_mpi,"echo * " + command)
        os.write(fh_mpi,command)
        
        # Launch second shell
        os.close(fh_mpi)
        self.launchCommandFile(xmpi_run_file + '_2.sh')
        
        message=" Done pre-processing of micrographs"
        print '* ',message
        print '*********************************************************************'
        self.log.info(message)

    def launchCommandFile(self,commandFile):
        import launch_job, log
        log.cat(self.log,commandFile)
        if self._DoParallel:
            command =' -i '+ commandFile
            launch_job.launch_job("xmipp_run",command,self.log,True,
                  self._MyNumberOfMpiProcesses,1,self._MySystemFlavour)
        else:
            self.log.info(commandFile)     
            os.system(commandFile)     
        # COSS os.remove(commandFile)

    def perform_preprocessing(self,filename,fh_mpi):
        # Decide name after preprocessing
        (filepath, micrographName) = os.path.split(filename)
        (shortname,extension) = os.path.splitext(micrographName)
        finalname=shortname+'/down'+str(self.Down)+'_'+shortname
        if not self.Stddev==-1 or not self.Crop==-1 or not self.Down==1:
            finalname+=".spi"
        else:
            finalname+=extension
            if not os.path.exists(finalname):
                command='ln -s '+filename+' '+finalname + "\n"
                os.write(fh_mpi,"echo * " + command)
                os.write(fh_mpi,command)
        if not self.DoPreprocess:
            return finalname
        
        # Crop
        iname=filename
        command="";
        if not self.Crop==-1:
            command+="xmipp_window -i "+iname+" -o "+finalname+" -crop "+str(self.Crop)+" ; "
            iname=finalname
        
        # Remove bad pixels
        if not self.Stddev==-1:
            command+="xmipp_filter -i "+iname+" -stdfilter "+str(self.Stddev)
            if not iname==finalname:
                command+=" -o "+finalname
            command+=" ; "
        
        # Downsample
        if not self.Down==1:
            command+="xmipp_micrograph_downsample -i "+iname+" -o "+shortname+"/tmp.spi "+\
                     "--output_bits 32 "
            if (self.DownKernel=='Fourier'):
                scale = 1./self.Down
                command+=' -fourier '+str(scale)
            elif (self.DownKernel=='Sinc'):
                command+=' -Xstep '+str(self.Down)+' -kernel sinc 0.02 0.1'
            elif (self.DownKernel=='Rectangle'):
                command+=' -Xstep '+str(self.Down)+' -kernel rectangle '+str(self.Down)+' '+str(self.Down)
            command+=" ; rm -f "+finalname+" ; mv -i "+shortname+"/tmp.spi "+finalname
        
        # Write the preprocessing command
        if not command=="":
            command+="\n"
            os.write(fh_mpi,"echo * " + command)
            os.write(fh_mpi,command)
        return finalname

    def perform_ctfestimate_xmipp(self,shortname,filename,fh_mpi):
        (filepath, micrographName) = os.path.split(filename)
        (fnRoot,extension) = os.path.splitext(micrographName)
        pname=shortname+'/'+fnRoot+'_estimate_ctf_input.param'
        print '*********************************************************************'
        print '*  Estimate PSD/CTF for micrograph: '+filename
        retval=[]
        
        retval.append(shortname+'/'+fnRoot+"_Periodogramavg.psd")
        retval.append(pname)

        # prepare parameter file
        paramlist = []
        paramlist.append('image= '+filename+'\n')
        paramlist.append('micrograph_averaging= yes \n')
        paramlist.append('N_horizontal= 512 \n')
        paramlist.append('periodogram= yes \n')
        if not self.OnlyEstimatePSD:
            AngPix = (10000. * self.ScannedPixelSize * self.Down) / self.Magnification
            paramlist.append('voltage= '+str(self.Voltage)+'\n')
            paramlist.append('spherical_aberration= '+str(self.SphericalAberration)+'\n')
            paramlist.append('sampling_rate= '+str(AngPix)+'\n')
            paramlist.append('particle_horizontal= 256 \n')
            paramlist.append('Q0= -'+str(self.AmplitudeContrast)+'\n')
            paramlist.append('min_freq= '+str(self.LowResolCutoff)+'\n')
            paramlist.append('max_freq= '+str(self.HighResolCutoff)+'\n')
            retval.append(shortname+'/'+fnRoot+"_Periodogramavg.ctfparam")
            retval.append(shortname+'/'+fnRoot+"_Periodogramavg_ctfmodel_halfplane.xmp")
            retval.append(shortname+'/'+fnRoot+"_Periodogramavg_ctfmodel_quadrant.xmp")

        # Perform CTF estimation
        fh=open(pname,"w")
        fh.writelines(paramlist)
        fh.close()
        command='xmipp_ctf_estimate_from_micrograph -i '+pname
        os.write(fh_mpi,"echo * " + command+"\n")
        os.write(fh_mpi,command+"\n");
        return retval

    def perform_ctfestimate_ctffind(self,shortname,filename,fh_mpi):
        # Convert image to MRC
        command='xmipp_convert_image -i '+filename+' -o '+shortname+'/tmp.mrc ; '

        # The new line is different if we are running in parallel or not
        theNewLine='\n'
        if(self._DoParallel):
            theNewLine='MPI_NEWLINE'
        
        # Prepare parameters for CTFTILT
        AngPix = (10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        (filepath, micrographName) = os.path.split(filename)
        (fnRoot,extension) = os.path.splitext(micrographName)
        command+="export NATIVEMTZ=kk ; "
        command+= self.CtffindExec+'  << eof > '+shortname+'/'+fnRoot+'_ctffind.log'+theNewLine
        command+= shortname+'/tmp.mrc'+theNewLine
        command+= shortname+'/'+fnRoot+'_ctffind_spectrum.mrc'+theNewLine
        command+= str(self.SphericalAberration) + ',' + \
                  str(self.Voltage) + ',' + \
                  str(self.AmplitudeContrast) + ',' + \
                  str(self.Magnification) + ',' + \
                  str(self.ScannedPixelSize*self.Down) +theNewLine
        command+= str(self.WinSize) + ',' + \
                  str(AngPix/self.MinResCTF) + ',' + \
                  str(AngPix/self.MaxResCTF) + ',' + \
                  str(self.MinFocus) + ',' + \
                  str(self.MaxFocus) + ',' + \
                  str(self.StepFocus) +theNewLine
        command+= 'eof\n'
        os.write(fh_mpi,"echo * " + command)
        os.write(fh_mpi,command)

    def pickup_ctfestimate_ctffind(self,shortname,filename,fh_mpi):
        import XmippData
        (filepath, micrographName) = os.path.split(filename)
        (fnRoot,extension) = os.path.splitext(micrographName)

        # Pick values from ctffind
        fh=open(shortname+'/'+fnRoot+'_ctffind.log','r')
        lines=fh.readlines()
        fh.close()
        DF1=0.
        DF2=0.
        Angle=0.
        for i in range(len(lines)):
            if not (lines[i].find('Final Values')==-1):
                words=lines[i].split()
                DF1=float(words[0])
                DF2=float(words[1])
                Angle=float(words[2])
                break

        # Generate Xmipp .ctfparam file:
        MD=XmippData.MetaData()
        MD.setColumnFormat(False)
        ptr=XmippData.doubleP()
        
        AngPix = (10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        ptr.assign(AngPix)
        MD.addObject()
        XmippData.setValueDouble(MD,XmippData.MDL_CTF_SAMPLING_RATE,ptr)

        ptr.assign(self.Voltage)
        XmippData.setValueDouble(MD,XmippData.MDL_CTF_VOLTAGE,ptr)

        ptr.assign(-DF2)
        XmippData.setValueDouble(MD,XmippData.MDL_CTF_DEFOCUSU,ptr)
        
        ptr.assign(-DF1)
        XmippData.setValueDouble(MD,XmippData.MDL_CTF_DEFOCUSV,ptr)

        ptr.assign(Angle)
        XmippData.setValueDouble(MD,XmippData.MDL_CTF_DEFOCUS_ANGLE,ptr)

        ptr.assign(self.SphericalAberration)
        XmippData.setValueDouble(MD,XmippData.MDL_CTF_CS,ptr)

        ptr.assign(-self.AmplitudeContrast)
        XmippData.setValueDouble(MD,XmippData.MDL_CTF_Q0,ptr)

        ptr.assign(1.0)
        XmippData.setValueDouble(MD,XmippData.MDL_CTF_K,ptr)

        fnOut=shortname+'/'+fnRoot+'_ctffind.ctfparam'
        MD.write(XmippData.FileName(fnOut))

        # Remove temporary files
        command="rm " + shortname+'/tmp.mrc\n'
        os.write(fh_mpi,"echo * " + command)
        os.write(fh_mpi,command)
        return fnOut

# Preconditions
def preconditions(gui):
    retval=True
    # Check ctffind executable
    if not CtffindExec=="":
        if (not os.path.exists(CtffindExec)):
            message="Cannot find ctffind executable: " + CtffindExec
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
                 DoDeleteWorkingDir,
                 DirMicrographs,
                 ExtMicrographs,
                 RootName,
                 ProjectDir,
                 LogDir,
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
                 CtffindExec,
                 WinSize,
                 MinResCTF,
                 MaxResCTF,
                 StepFocus,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour)
