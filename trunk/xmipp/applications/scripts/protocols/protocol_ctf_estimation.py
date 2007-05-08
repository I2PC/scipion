#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for CTF estimation of micrographs: 
#
# It is assumed that you have already ran the preprocess_micrographs protocol,
# You also need the xmipp_preprocess_micrographs.py file in the current directory
#
# Example use:
# ./protocol_ctf_estimation.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {dir} Working subdirectory:
""" Use the same directory where you executed protocol_preprocess_micrographs.py
"""
WorkingDir="Preprocessing"
# {file} Selfile with micrographs on which to perform processing
MicrographSelfile='Preprocessing/all_micrographs.sel'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir="/home2/bioinfo/scheres/work/protocols/G40P"
# {expert} Directory name for logfiles:
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Contrast Transfer Function estimation
#------------------------------------------------------------------------------------------------
# Microscope voltage (in kV)
Voltage=200
# Spherical aberration
SphericalAberration=2.26
# Magnification rate
Magnification=50000
# Scanned pixels size in micrometer
ScannedPixelSize=14
# Amplitude Contrast
AmplitudeContrast=0.1
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
# Use N. Grigorieffs CTFFIND instead of Xmipp to estimate CTF?
DoCtffind=False
# Location where CTFFIND is installed
CtffindExec='/home/cnbb13/indivi/bin/ctffind3b.exe'
# {expert} Window size
WinSize=128
# {expert} Minimum resolution (in Ang.)
MinRes=200.0
# {expert} Maximum resolution factor
""" The maximum resolution used in the CTF estimation will be:
    this factor x 1000 x Down x ScannedPixelSize / Magnification
    recommended value: 3
"""
MaxResFactor=3
# {expert} Minimum defocus to search
MinFocus=5000
# {expert} Maximum defocus to search
MaxFocus=100000
# {expert} Defocus step
StepFocus=500
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_ctf_estimation.py"
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#  {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#

class preprocess_particles_class:

    #init variables
    def __init__(self,
                 WorkingDir,
                 MicrographSelfile,
                 ProjectDir,
                 LogDir,
                 Voltage,
                 SphericalAberration,
                 Magnification,
                 ScannedPixelSize,
                 AmplitudeContrast,
                 LowResolCutoff,
                 HighResolCutoff,
                 DoCtffind,
                 CtffindExec,,
                 WinSize,
                 MinRes,
                 MaxResFactor,
                 MinFocus,
                 MaxFocus,
                 StepFocus
                 ):
	     
        import os,sys
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log,protocol_preprocess_micrographs
        
        self.WorkingDir=WorkingDir
        self.MicrographSelfile=os.path.abspath(MicrographSelfile)
        self.ProjectDir=ProjectDir
        self.LogDir=LogDir
        self.Voltage=Voltage
        self.SphericalAberration=SphericalAberration
        self.Magnification=Magnification
        self.ScannedPixelSize=ScannedPixelSize
        self.AmplitudeContrast=AmplitudeContrast
        self.DoCtffind=DoCtffind
        self.LowResolCutoff=LowResolCutoff
        self.HighResolCutoff=HighResolCutoff
        self.DoCtffind=DoCtffind
        self.CtffindExec=CtffindExec
        self.WinSize=WinSize
        self.MinRes=MinRes
        self.MaxResFactor=MaxResFactor
        self.MinFocus=MinFocus
        self.MaxFocus=MaxFocus
        self.StepFocus=StepFocus
        
        # Parameters set from outside
        self.Down=protocol_preprocess_micrographs.Down

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Make working directory if it does not exist yet
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                       os.path.abspath(self.WorkingDir))
    
        # Execute protocol in the working directory
        os.chdir(self.WorkingDir)        
        self.process_all_micrographs()

        # Return to parent dir
        os.chdir(os.pardir)
            
    def process_all_micrographs(self):
        import os
        import selfile
        print '*********************************************************************'
        print '*  Estimate CTF for micrographs in '+os.path.basename(self.MicrographSelfile)

        self.ctfselfile = []
        self.allselfile = []
        self.allctflibfile = []

        mysel=selfile.selfile()
        mysel.read(self.MicrographSelfile)
        for name,state in mysel.sellines:
            self.shortname,self.downname=os.path.split(name)
            self.downname=self.downname.replace('.raw','')
            if (state.find('-1') < 0):
                if (self.UseCtffind):
                    self.perform_estimate_ctffind()
                else:
                    self.perform_estimate_xmipp()
                    

    def perform_estimate_ctffind(self):
        print 'Error: CTFFIND CTF-estimation not implemented yet...'
        import sys
        sys.exit()
        
        
    def perform_estimate_xmipp(self):
        import os
        iname=self.shortname+'/'+self.downname+'.raw'
        pname=self.shortname+'/'+self.shortname+'_input.param'
        selname=self.shortname+'/'+self.downname+'.raw.sel' 
        ctfname=self.downname+'.raw.ctf.sel' 
        ctfname2=self.shortname+'/'+self.downname+'.raw.ctf.sel' 
        print '*********************************************************************'
        print '*  Estimate CTF for micrograph: '+iname

        # Calculate sampling in Angstrom/pixel
        AngPix = (10000. * self.ScannedPixelSize * self.Down) / self.Magnification
        
        # prepare parameter file
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
        paramlist.append('selfile= '+selname+'\n')
        paramlist.append('periodogram= yes \n')

        # Perform CTF estimation
        fh=open(pname,"w")
        fh.writelines(paramlist)
        fh.close()
        command='xmipp_ctf_estimate_from_micrograph -i '+pname
        print '* ',command
        self.log.info(command)
        os.system(command )

        # Add entry to the ctfselfile (for visualization of all CTFs)
        currdir=os.getcwd()
        oname=currdir+'/'+self.shortname+'/'+self.downname+'_Periodogramavg.ctfmodel'
        self.ctfselfile.append(oname+' 1 \n')

        # Move ctf.sel file into subdirectory
        if os.path.exists(ctfname):
            os.rename(ctfname,ctfname2)

        # Fill selfile with all CTFs
        fh=open(ctfname2,'r')
        text = fh.readlines()
        fh.close()
        outctfselname=self.ProjectDir+'/'+OutCTFSelFile
        fh=open(outctfselname,'a')
        fh.writelines(text)
        fh.close()

        # Update all_ctfs.sel
        fh=open('all_ctfs.sel','w')
        fh.writelines(self.ctfselfile)
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

   	# create preprocess_particles_class object

	preprocess_particles=preprocess_particles_class(WorkingDir,
                                                        MicrographSelfile,
                                                        ProjectDir,
                                                        LogDir,
                                                        UseCtffind,
                                                        Voltage,
                                                        SphericalAberration,
                                                        Magnification,
                                                        ScannedPixelSize,
                                                        AmplitudeContrast,
                                                        LowResolCutoff,
                                                        HighResolCutoff,
                                                        OutCTFSelFile)

	# close 
	preprocess_particles.close()

