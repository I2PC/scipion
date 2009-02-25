#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based pre-processing of micrographs: 
#  - tif2raw conversion, 
#  - downsampling
#  - power spectral density (PSD) and CTF estimation on the micrograph
#  - CTF phase flipping on the micrograph
#
# For each file <micrograph>.tif given, this script will perform 
# the requested operations below.
# For each micrograph a subdirectory will be created
#
# Example use:
# ./xmipp_preprocess_micrographs.py *.tif
#   or if you wish to use several CPUs
# mpirun -np 8 ./xmipp_preprocess_micrographs.py *.tif
# NOTE1: you need the packahe mpi4py in order to use the mpi extension
# NOTE2: each process will create a different log file. The name is the same that
#        the one used by the secuencial program plus the CPU number     
# Author: Sjors Scheres, March 2007
#         Roberto Marabini (mpi extension)
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
""" This is typically *.tif, but may also be *.mrc or *.spi (see the expert options)
    Note that any wildcard is possible, e.g. *3[1,2].tif
"""
ExtMicrographs='*.tif'
# Name for the output micrograph selfile:
""" Be aware that this file will be overwritten if it already exists!
"""
MicrographSelfile='all_micrographs.sel'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir='/home/coss/temp/Far_from_focus'
# {expert} Directory name for logfiles:
LogDir='Logs'
#------------------------------------------------------------------------------------------------
# {section} Initial conversion to raw
#------------------------------------------------------------------------------------------------
# Perform micrograph conversion? 
DoConversion=True
# {list}|Tif2Raw|Mrc2Raw|Spi2Raw|Raw2Raw| Which conversion to perform?
""" Some TIF formats are not recognized. In that case, save your micrographs as spider, mrc or raw and try to convert those. Note that raw2raw assumes the raw files have an Xmipp-like raw.info file with the same rootname, and that in this case no conversion takes place, but only the required directory structure is made.
    
"""
ConversionTask='Tif2Raw'
#------------------------------------------------------------------------------------------------
# {section} Downsampling
#------------------------------------------------------------------------------------------------
# Perform downsampling?
DoDownSample=True
# Downsampling factor 
Down=2
# {expert}{list}|Fourier|Rectangle|Sinc| Which method to use for downsampling?
""" Fourier is theoretically the best option, but it may take more memory than your machine can handle. Then, Rectangle is the fastest, but least accurate. Since is reasonably accurate, but painfully slow...
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
AmplitudeContrast=0.1
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
#------------------------------------------------------------------------------------------------
# {section} CTFFIND
#------------------------------------------------------------------------------------------------
# Use N. Grigorieffs CTFFIND instead of Xmipp?
""" Some people prefer the faster CTFFIND program.
    Note however that this option will yield no information about the CTF envelope, and therefore this option cannot be used with the high-resolution refinement protocol.
"""
DoCtffind=True
# {file} Location of the CTFFIND executable
CtffindExec='../../bin/ctffind3.exe'
# {expert} Window size
WinSize=128
# {expert} Minimum resolution (in Ang.)
MinRes=200.0
""" Lowest resolution to include in CTF estimation
"""
# {expert} Maximum resolution factor
""" The highest resolution used in the CTF estimation will be:
    <this factor> x 10000 x Down x ScannedPixelSize / Magnification
    recommended value: 3
"""
MaxResFactor=3
# {expert} Minimum defocus to search (in Ang.)
""" Minimum defocus value (in Angstrom) to include in defocus search
"""
MinFocus=5000
# {expert} Maximum defocus to search (in Ang.)
""" Maximum defocus value (in Angstrom) to include in defocus search
"""
MaxFocus=100000
# {expert} Defocus step (in Ang.)
""" Step size for defocus search (in Angstrom)
"""
StepFocus=500
#------------------------------------------------------------------------------------------------
# {section} CTF phase flipping
#------------------------------------------------------------------------------------------------
# Perform CTF phase flipping on the micrographs?
""" The phase-flipped micrographs will be saved with a different format (spider) than the original raw-format.
"""
DoCtfPhaseFlipping=True
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_preprocess_micrographs.py'
#
#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# This script has been parallelized at python level
""" This python script can be run as a parallel job,
    In order to use several CPUs you will need to
    answer yes to the "Use a job queing system" pop-up
    window and after the queueing command add
    mpirun -np N -machinefile machine.txt 
    or whatever is adecuate in your computer
    (N is the number of CPUs).  Each process will create a
    different log file. The names of these files are the
    same as the ones used by the sequential program plus the process number
"""
UselessVariable='read help'
# {expert} Abort if mpi is not installed?
"""If set to True the script will abort if mpi4py is not installed and one executes this script through mpirun
"""
CheckMpi=False
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
#Next lines will import mpi4py package, if import fails a dummy 
#mpi class will be created. 
#In this way the same file may be use for mpi and
#serial jobs  
try:
    from mpi4py import MPI
    mpi = MPI.COMM_WORLD
    mpiActive=True
except ImportError:
    class DummyMPI:
        '''A dummy MPI class for running serial jobs.'''
        def Get_rank(dummy):
            return 0
        def Get_size(dummy):
            return 1
        def Get_processor_name(dummy):
            return 'localhost' 
        def Barrier(dummy):
            return 1
    MPI = DummyMPI()
    mpi = DummyMPI()
    mpiActive=False
#
class preprocess_A_class:
    def __init__(self,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 DirMicrographs,
                 ExtMicrographs,
                 MicrographSelfile,
                 ProjectDir,
                 LogDir,
                 DoConversion,
                 ConversionTask,
                 DoDownSample,
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
                 DoCtffind,
                 CtffindExec,
                 WinSize,
                 MinRes,
                 MaxResFactor,
                 MinFocus,
                 MaxFocus,
                 StepFocus,
                 DoCtfPhaseFlipping,
                 CheckMpi
                 ):
        
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir)
        import log

        self.WorkingDir=os.path.abspath(WorkingDir)
        self.DirMicrographs=os.path.abspath(DirMicrographs)
        self.ExtMicrographs=ExtMicrographs
        self.ProjectDir=os.path.abspath(ProjectDir)
        self.LogDir=LogDir
        self.DoConversion=DoConversion
        self.ConversionTask=ConversionTask
        self.DoDownSample=DoDownSample
        self.Down=Down
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
        self.DoCtffind=DoCtffind
        self.WinSize=WinSize
        self.MinRes=MinRes
        self.MaxResFactor=MaxResFactor
        self.MinFocus=MinFocus
        self.MaxFocus=MaxFocus
        self.StepFocus=StepFocus
        self.DoCtfPhaseFlipping=DoCtfPhaseFlipping
        self.MicrographSelfile=MicrographSelfile
        self.myrank = mpi.Get_rank()
        self.nprocs = mpi.Get_size()

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0]+"_"+str(self.myrank),
                                     self.WorkingDir)
        # Check ctffind executable
        if (self.DoCtffind):
            if (not os.path.exists(CtffindExec)):
                message = "Error: cannot find ctffind executable: ",CtffindExec
                print '*',message
                self.log.error(message)
                sys.exit()
            else:
                self.CtffindExec=os.path.abspath(CtffindExec)

        if (mpiActive):
            print "Python has mpi4py package, python parallel jobs are possible"
            self.log.info("Python has mpi4py package, python parallel jobs are possible")        
        else:
            print "Python has NOT mpi4py package, python parallel jobs are NOT possible"
            self.log.info("Python has NOT mpi4py package, python parallel jobs are NOT possible")
            if (CheckMpi):
                print "mpi is not installed and checkmpi flag is set to True"
                print "disable mpi checking flag in the script" 
                print "and launch it as a serial job" 
                print "that is, without mpirun" 
                exit(1)
        # Delete working directory if exists, make a new one
        if (self.myrank==0):
            if (DoDeleteWorkingDir): 
                if os.path.exists(self.WorkingDir):
                    shutil.rmtree(self.WorkingDir)
            if not os.path.exists(self.WorkingDir):
                os.makedirs(self.WorkingDir)

            # Backup script
            log.make_backup_of_script_file(sys.argv[0],
                                       os.path.abspath(self.WorkingDir))
        #wait until myrank=0 has arrived here
        mpi.Barrier() 
        # Execute protocol in the working directory
        os.chdir(self.WorkingDir)
        self.process_all_micrographs()

        # Return to parent dir
        os.chdir(os.pardir)

        
    def process_all_micrographs(self):
        import os
        import glob
        print '*********************************************************************'
        print '*  Processing the following micrographs: '
        for self.filename in glob.glob(self.DirMicrographs+'/'+self.ExtMicrographs):
            (self.filepath, self.name) = os.path.split(self.filename)
            print '*  '+self.name

        self.psdselfile = []
        self.ctfselfile = []
        self.inputselfile = []
        self.micselfile = []
	job_index=0;
        for self.filename in glob.glob(self.DirMicrographs+'/'+self.ExtMicrographs):
            (self.filepath, self.name) = os.path.split(self.filename)
            (self.shortname,extension) = os.path.splitext(self.name)
            self.downname='down'+str(self.Down)+'_'+self.shortname
            job_index = job_index + 1
		
	    #if I am not the right node skip this micrograph	
            #print job_index, self.nprocs, job_index % self.nprocs , self.myrank
	    if ((job_index % self.nprocs)==self.myrank):
                if not os.path.exists(self.shortname):
                    os.makedirs(self.shortname)

                if (self.DoConversion):
                    if (self.ConversionTask=='Tif2Raw'):
                        self.perform_tif2raw()
                    elif (self.ConversionTask=='Mrc2Raw'):
                        self.perform_mrc2raw()
                    elif (self.ConversionTask=='Spi2Raw'):
                        self.perform_spi2raw()
                    elif (self.ConversionTask=='Raw2Raw'):
                        self.perform_raw2raw()
                    else:
                        message="Unrecognized ConversionTask: choose from list options"
                        print '*',message
                        self.log.error(message)
                        sys.exit()

                if (self.DoDownSample):
                    self.perform_downsample()

            if not os.path.exists(self.shortname+'/'+self.downname+'.raw'):
                self.downname=self.shortname

	    if ((job_index % self.nprocs)==self.myrank):
                if (self.DoCtfEstimate):
                    if not (self.DoCtffind):
                        if (self.OnlyEstimatePSD):
                            self.perform_only_psdestimate()
                        else:
                            self.perform_ctfestimate_xmipp()
                    else:
                        self.perform_ctfestimate_ctffind()

                if (self.DoCtfPhaseFlipping):
                    self.perform_ctf_phase_flipping()

            #only one node should write sel file with micrographs
	    if (self.myrank==0):
	        self.append_micrograph_selfile()
        	if (self.DoCtfEstimate):
                    if not (self.DoCtffind):
                	if (self.OnlyEstimatePSD):
                            #self.perform_only_psdestimate()
        		    oname=self.shortname+'/'+self.downname+'_Periodogramavg.psd'
        		    self.psdselfile.append(oname+' 1\n')
			    fh = open("all_psds.sel","w")
        		    fh.writelines(self.psdselfile)
        		    fh.close()
                	else:
                            #self.perform_ctfestimate_xmipp()
        		    # Add entry to the ctfselfile (for visualization of all CTFs)
        		    oname=self.shortname+'/'+self.downname+'_Periodogramavg.ctfmodel_halfplane'
        		    self.ctfselfile.append(oname+' 1\n')
        		    fh=open('all_ctfs.sel','w')
        		    fh.writelines(self.ctfselfile)
        		    fh.close()
                            # Add entry to the inputparamselfile
        		    oname=self.shortname+'/'+self.shortname+'_input.param'
                            self.inputselfile.append(oname+' 1\n')
        		    fh=open('all_inputparams.sel','w')
        		    fh.writelines(self.inputselfile)
        		    fh.close()

        		    #Also add enrty to psdselfile
        		    oname=self.shortname+'/'+self.downname+'_Periodogramavg.psd'
        		    self.psdselfile.append(oname+' 1\n')
        		    fh = open("all_psds.sel","w")
        		    fh.writelines(self.psdselfile)
        		    fh.close()
                    else:
                	#self.perform_ctfestimate_ctffind()
        		ctfname = self.shortname + '/ctffind_' + self.downname + '_ctfmodel.xmp'
        		self.ctfselfile.append(ctfname+' 1 \n')
			fh=open('all_ctfs.sel','w')
        		fh.writelines(self.ctfselfile)
        		fh.close()

    def perform_tif2raw(self):
        import os
        import launch_job
        oname=self.shortname+'/'+self.shortname+'.raw'
        print '*********************************************************************'
        print '*  Generating RAW for micrograph: '+self.name
        command=' '+self.filename+' '+oname
        launch_job.launch_job("xmipp_convert_tiff2raw",
                              command,
                              self.log,
                              False,1,1,'')

    def perform_mrc2raw(self):
        import os
        import launch_job
        oname=self.shortname+'/'+self.shortname+'.raw'
        tname=self.shortname+'/'+self.shortname+'.spi'
        print '*********************************************************************'
        print '*  Generating RAW for micrograph: '+self.name
        command=' -i '+self.filename+' -o '+tname
        launch_job.launch_job("xmipp_convert_spi22ccp4",
                              command,
                              self.log,
                              False,1,1,'')
        command=' -generate_inf -f -i '+tname+' -o '+oname
        launch_job.launch_job("xmipp_convert_raw22spi",
                              command,
                              self.log,
                              False,1,1,'')
        os.remove(tname)
        
    def perform_spi2raw(self):
        import os
        import launch_job
        oname=self.shortname+'/'+self.shortname+'.raw'
        print '*********************************************************************'
        print '*  Generating RAW for micrograph: '+self.name
        command=' -generate_inf -f -i '+self.filename+' -o '+oname
        launch_job.launch_job("xmipp_convert_raw22spi",
                              command,
                              self.log,
                              False,1,1,'')

    def perform_raw2raw(self):
        import os
        oname=self.shortname+'/'+self.shortname+'.raw'
        print '*********************************************************************'
        print '*  Generating RAW for micrograph: '+self.name
        command='cp '+self.filename+' '+oname
        print '* ',command
        self.log.info(command)
        os.system(command)
        command='cp '+self.filename+'.inf '+oname+'.inf'
        print '* ',command
        self.log.info(command)
        os.system(command)
    
    def perform_downsample(self):
        import os
        import launch_job
        iname=self.shortname+'/'+self.shortname+'.raw'
        oname=self.shortname+'/'+self.downname+'.raw'
        print '*********************************************************************'
        print '*  Downsampling micrograph: '+iname
        if (self.DownKernel=='Fourier'):
            scale = 1./self.Down
            command=' -i '+iname+' -o '+oname+' -output_bits 32 -fourier '+str(scale)
        elif (self.DownKernel=='Sinc'):
            command=' -i '+iname+' -o '+oname+' -output_bits 32 -Xstep '+str(self.Down)+' -kernel sinc 0.02 0.1'
        elif (self.DownKernel=='Rectangle'):
            command=' -i '+iname+' -o '+oname+' -output_bits 32 -Xstep '+str(self.Down)+' -kernel rectangle '+str(self.Down)+' '+str(self.Down)
        else:
            message="Unrecognized DownKernel: choose from list options"
            print '*',message
            self.log.error(message)
            sys.exit()

        launch_job.launch_job("xmipp_micrograph_downsample",
                              command,
                              self.log,
                              False,1,1,'')

    def perform_only_psdestimate(self):
        import os
        import launch_job
        iname=self.shortname+'/'+self.downname+'.raw'
        pname=self.shortname+'/'+self.shortname+'_psd.param'
        print '*********************************************************************'
        print '*  Estimate PSD for micrograph: '+iname
        paramlist = []
        paramlist.append('image= '+iname+'\n')
        paramlist.append('micrograph_averaging= yes \n')
        paramlist.append('N_horizontal= 512 \n')
        paramlist.append('periodogram= yes \n')
        paramlist.append('dont_adjust_CTF= yes \n')
    	
        fh = open(pname,"w")
        fh.writelines(paramlist)
        fh.close()
        command=' -i '+pname
        launch_job.launch_job("xmipp_ctf_estimate_from_micrograph",
                              command,
                              self.log,
                              False,1,1,'')
        #oname=self.shortname+'/'+self.downname+'_Periodogramavg.psd'
        #self.psdselfile.append(oname+' 1\n')
	#fh = open("all_psds.sel","w")
        #fh.writelines(self.psdselfile)
        #fh.close()
    
    def perform_ctfestimate_xmipp(self):
        import os
        import launch_job
        iname=self.shortname+'/'+self.downname+'.raw'
        pname=self.shortname+'/'+self.shortname+'_input.param'
        print '*********************************************************************'
        print '*  Estimate CTF for micrograph: '+iname

        # prepare parameter file
        AngPix = (10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        paramlist = []
        paramlist.append('image= '+iname+'\n')
        paramlist.append('micrograph_averaging= yes \n')
        paramlist.append('voltage= '+str(self.Voltage)+'\n')
        paramlist.append('spherical_aberration= '+str(self.SphericalAberration)+'\n')
        paramlist.append('sampling_rate= '+str(AngPix)+'\n')
        paramlist.append('particle_horizontal= 128 \n')
        paramlist.append('Q0= -'+str(self.AmplitudeContrast)+'\n')
        paramlist.append('N_horizontal= 512 \n')
        paramlist.append('min_freq= '+str(self.LowResolCutoff)+'\n')
        paramlist.append('max_freq= '+str(self.HighResolCutoff)+'\n')
        paramlist.append('periodogram= yes \n')

        # Perform CTF estimation
        fh=open(pname,"w")
        fh.writelines(paramlist)
        fh.close()
        command=' -i '+pname
        launch_job.launch_job("xmipp_ctf_estimate_from_micrograph",
                              command,
                              self.log,
                              False,1,1,'')

        ## Add entry to the ctfselfile (for visualization of all CTFs)
        #oname=self.shortname+'/'+self.downname+'_Periodogramavg.ctfmodel_halfplane'
        #self.ctfselfile.append(oname+' 1\n')
        #fh=open('all_ctfs.sel','w')
        #fh.writelines(self.ctfselfile)
        #fh.close()

        #Also add enrty to psdselfile
        #oname=self.shortname+'/'+self.downname+'_Periodogramavg.psd'
        #self.psdselfile.append(oname+' 1\n')
        #fh = open("all_psds.sel","w")
        #fh.writelines(self.psdselfile)
        #fh.close()


    def perform_ctfestimate_ctffind(self):
        import os
        # Prepare stuff
        DStep = self.ScannedPixelSize * self.Down
        MaxRes = (self.MaxResFactor * 10000. * DStep) / self.Magnification
        self.convert_raw_to_mrc()

        # Execute CTFFIND
        #next line tell ctfind to skip endian checking
        #I could have modified the program spi22ccp4 but
        #I do not agree with cctfind interpretation of the flag
	os.putenv('NATIVEMTZ', "kk")
        command=  self.CtffindExec+'  << eof > '+self.shortname+'/ctffind_'+self.downname+'.log\n'
        command+= self.shortname+'/tmp.mrc\n'
        command+= self.shortname+'/spectrum.mrc\n'
        command+= str(self.SphericalAberration) + ',' + \
                  str(self.Voltage) + ',' + \
                  str(self.AmplitudeContrast) + ',' + \
                  str(self.Magnification) + ',' + \
                  str(DStep) + '\n'
        command+= str(self.WinSize) + ',' + \
                  str(self.MinRes) + ',' + \
                  str(MaxRes) + ',' + \
                  str(self.MinFocus) + ',' + \
                  str(self.MaxFocus) + ',' + \
                  str(self.StepFocus) + '\n'
        command+= 'eof\n'
        print '* ',command
        self.log.info(command)
        os.system(command )

        # Convert output to Xmipp ctfparam files
        self.convert_ctffind_output_to_xmipp_style()        

        # Remove temporary files
        os.remove(self.shortname+'/tmp.mrc')
        os.remove(self.shortname+'/tmp.spi')
        os.remove(self.shortname+'/spectrum.mrc')

    def convert_raw_to_mrc(self):
        import os
        import launch_job
        command= ' -i '+ self.shortname+'/'+self.downname+'.raw ' + \
                 ' -o '+ self.shortname+'/tmp.spi ' + \
                 ' -is_micrograph -f'
        launch_job.launch_job("xmipp_convert_raw22spi",
                              command,
                              self.log,
                              False,1,1,'')
        command= ' -i '+ self.shortname+'/tmp.spi ' + \
                 ' -o '+ self.shortname+'/tmp.mrc '
        launch_job.launch_job("xmipp_convert_spi22ccp4",
                              command,
                              self.log,
                              False,1,1,'')

    def convert_ctffind_output_to_xmipp_style(self):
        import os;
        import launch_job
        logfile=self.shortname+'/ctffind_'+self.downname+'.log'
        fh=open(logfile,'r')
        lines=fh.readlines()
        fh.close()
        newlines=[]
        DF1=0.
        DF2=0.
        Angle=0.
        for i in range(len(lines)):
            if not (lines[i].find('Final Values')==-1):
                words=lines[i].split()
                DF1=words[0]
                DF2=words[1]
                Angle=words[2]
            if (lines[i].find('CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]')>-1):
                newlines.append(str(DF1)+' '+str(DF2)+' '+str(Angle)+' '+lines[i+1])

        # Write CTFFIND .ctf file
        ctffile=self.shortname+'/ctffind_'+self.downname+'.ctf'
        fh=open(ctffile,'w')
        fh.writelines(newlines)
        fh.close()

        # Convert MRC ctf model to Xmipp image
        ctfname = self.shortname + '/ctffind_' + self.downname + '_ctfmodel.xmp'
        command= ' -i ' + self.shortname + '/spectrum.mrc ' + \
                 ' -o '+ ctfname
        launch_job.launch_job("xmipp_convert_spi22ccp4",
                              command,
                              self.log,
                              False,1,1,'')

        # Add entry to the ctfselfile (for visualization of all CTFs)
        #ctfname = self.shortname + '/ctffind_' + self.downname + '_ctfmodel.xmp'
        #self.ctfselfile.append(ctfname+' 1 \n')
	#fh=open('all_ctfs.sel','w')
        #fh.writelines(self.ctfselfile)
        #fh.close()
  
        # Generate Xmipp .ctfparam file:
        paramname=self.shortname+'/'+self.downname+'_Periodogramavg.ctfparam'
        if os.path.exists(paramname):
            os.remove(paramname)
        AngPix = (10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        paramlist = []
        paramlist.append('defocusU= -'+str(DF1)+'\n')
        paramlist.append('defocusV= -'+str(DF2)+'\n')
        paramlist.append('azimuthal_angle= '+str(Angle)+'\n')
        paramlist.append('sampling_rate= '+str(AngPix)+'\n')
        paramlist.append('voltage= '+str(self.Voltage)+'\n')
        paramlist.append('spherical_aberration= '+str(self.SphericalAberration)+'\n')
        paramlist.append('Q0= -'+str(self.AmplitudeContrast)+'\n')
        paramlist.append('K= 1.\n')
        fh=open(paramname,"w")
        fh.writelines(paramlist)
        fh.close()


    def perform_ctf_phase_flipping(self):
        import os
        import launch_job
        iname=self.shortname+'/'+self.downname+'.raw'
        oname=self.shortname+'/'+self.downname+'.spi'
        paramname=self.shortname+'/'+self.downname+'_Periodogramavg.ctfparam'
        command= ' -i   ' + iname + \
                 ' -o   ' + oname + \
                 ' -ctf ' + paramname
        launch_job.launch_job("xmipp_micrograph_phase_flipping",
                              command,
                              self.log,
                              False,1,1,'')

    def append_micrograph_selfile(self):
        if self.DoCtfPhaseFlipping:
            name=self.shortname+'/'+self.downname+'.spi'
        else:
            name=self.shortname+'/'+self.downname+'.raw'
        self.micselfile.append(name+' 1\n')
        fh = open(self.MicrographSelfile,'w')
        fh.writelines(self.micselfile)
        fh.close()


    def close(self):
        message=" Done pre-processing all"
        print '* ',message
        print '*********************************************************************'
        self.log.info(message)
#		
# Main
#     
if __name__ == '__main__':

    # create preprocess_A_class object
    
    preprocessA=preprocess_A_class(WorkingDir,
                                   DoDeleteWorkingDir,
                                   DirMicrographs,
                                   ExtMicrographs,
                                   MicrographSelfile,
                                   ProjectDir,
                                   LogDir,
                                   DoConversion,
                                   ConversionTask,
                                   DoDownSample,
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
                                   DoCtffind,
                                   CtffindExec,
                                   WinSize,
                                   MinRes,
                                   MaxResFactor,
                                   MinFocus,
                                   MaxFocus,
                                   StepFocus,
                                   DoCtfPhaseFlipping,
                                   CheckMpi)

    # close 
    preprocessA.close()

