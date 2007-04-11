#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
#  - CTF estimation
#  - extraction of particles
#  - phase flipping
#  - normalization
#  - sort_junk
#
# It is assumed that you have already ran the preprocess_micrographs protocol,
#  and that you have picked the particles for each micrograph
# You also need the xmipp_preprocess_micrographs.py file in the current directory
#
# Example use:
# ./xmipp_preprocess_particles.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Selfile with micrographs on which to perform processing
MicrographSelfile='all_micrographs.sel'
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Extract particles
#------------------------------------------------------------------------------------------------
# Extract particles from micrographs?
DoExtract=True
## These file are supposed to be in the corresponding subdirectories
# Extension for the picked coordinates file
PosFile="raw.Common.pos"
# Dimension of the particles to extract (in pix.)
Size=64
# {expert} Directory name for particle images:
ImagesDir="Images"
# {expert} Selfile with all particles (is placed in ProjectDir)
OutSelFile="all_images.sel"
#------------------------------------------------------------------------------------------------
# {section} Normalization
#------------------------------------------------------------------------------------------------
# Perform particle normalization?
DoNormalize=True
# Pixels outside this circle are assumed to be noise and their stddev is set to 1.
# Radius for background circle definition (in pix.)
BackGroundRadius=30
# Perform ramping background correction?
DoUseRamp=True
# Perform black dust particles removal?
DoRemoveBlackDust=False
# Perform white dust particles removal?
DoRemoveWhiteDust=False
#------------------------------------------------------------------------------------------------
# {section} Contrast Transfer Function estimation
#------------------------------------------------------------------------------------------------
# Perform CTF estimation?
DoEstimateCTF=True
# Microscope voltage (in kV)
Voltage=200
# Spherical aberration
SphericalAberration=2.26
# Pixel size in the particles (in Angstrom/pixel)
Sampling=2.8
# Amplitude Contrast
AmplitudeContrast=0.1
# Visualize CTF-estimation?
DoVisualizeCTF=True
# {expert} Selfile with CTFs for all particles (is placed in ProjectDir)
OutCTFSelFile="all_images_ctf.sel"
#------------------------------------------------------------------------------------------------
# {section} Phase correction
#------------------------------------------------------------------------------------------------
# Perform CTF-phase flipping?
DoPhaseFlipping=True
#------------------------------------------------------------------------------------------------
# {section} Particle sorting
#------------------------------------------------------------------------------------------------
# Perform particle sorting to identify junk particles?
DoSorting=True
# Visualize sorted particles?
DoVisualizeSorting=True
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
                     MicrographSelfile,
                     ProjectDir,
                     LogDir,
                     DoExtract,
                     PosFile,
                     Size,
                     ImagesDir,
                     OutSelFile,
                     DoNormalize,
                     BackGroundRadius,
                     DoUseRamp,
                     DoRemoveBlackDust,
                     DoRemoveWhiteDust,
		     DoEstimateCTF,
                     Voltage,
                     SphericalAberration,
                     Sampling,
                     AmplitudeContrast,
                     DoVisualizeCTF,
                     OutCTFSelFile,
                     DoPhaseFlipping,
                     DoSorting,
                     DoVisualizeSorting,
                     ):
	     
            import os,sys
            scriptdir=os.path.expanduser('~')+'/scripts/'
            sys.path.append(scriptdir) # add default search path
            import log,protocol_preprocess_micrographs
        
            self.MicrographSelfile=MicrographSelfile
            self.ProjectDir=ProjectDir
            self.LogDir=LogDir
            self.DoExtract=DoExtract
            self.PosFile=PosFile
            self.Size=Size
            self.ImagesDir=ImagesDir
            self.OutSelFile=OutSelFile
            self.DoNormalize=DoNormalize
            self.BackGroundRadius=BackGroundRadius
            self.DoUseRamp=DoUseRamp
            self.DoRemoveBlackDust=DoRemoveBlackDust
            self.DoRemoveWhiteDust=DoRemoveWhiteDust
            self.DoEstimateCTF=DoEstimateCTF
            self.Voltage=Voltage
            self.SphericalAberration=SphericalAberration
            self.Sampling=Sampling
            self.AmplitudeContrast=AmplitudeContrast
            self.DoVisualizeCTF=DoVisualizeCTF
            self.OutCTFSelFile=OutCTFSelFile
            self.DoPhaseFlipping=DoPhaseFlipping
            self.DoSorting=DoSorting
            self.DoVisualizeSorting=DoVisualizeSorting
            # Parameters set from outside
            self.Down=protocol_preprocess_micrographs.Down
            self.log=log.init_log_system(self.ProjectDir,
                                         self.LogDir,
                                         sys.argv[0],
                                         '.')
            
            # Execute program
            self.process_all_micrographs()

            
        def check_have_marked(self):
            import os
            posname=self.shortname+'/'+self.downname+'.'+PosFile 
            if os.path.exists(posname):
                return True
            else:
                return False

        def perform_extract(self,is_for_flipping):
            import os
            iname=self.shortname+'/'+self.downname+'.raw' 
            selname=self.downname+'.raw.sel' 
            selname2=self.shortname+'/'+self.downname+'.raw.sel' 
            posname=self.shortname+'/'+self.downname+'.'+PosFile 
            imgsubdir=self.ProjectDir+'/'+self.ImagesDir+'/'+self.shortname
            rootname=imgsubdir+'/'+self.shortname+'_'
            logname=self.shortname+'/scissor.log'

	    # Use double image size for phase flipping
            if (is_for_flipping):
                size=2*self.Size
            else:
                size=self.Size

	    # Make directory if necessary
	    if not os.path.exists(imgsubdir):
                os.makedirs(imgsubdir)

	    # perform scissor operation
            command='xmipp_micrograph_scissor -i '+iname+' -pos '+posname+' -root '+rootname+' -Xdim '+str(size)+'|grep "corresponding image is set to blank"> '+logname
            print '* ',command
            self.log.info(command)
            os.system(command)

            # Move output selfile inside the sub-directory:
            os.rename(selname,selname2)

	    # Remove particles near the border: delete empty images and rewrite selfile & posfile
            fh=open(logname,"r")
            text=fh.read()
	    fh.close()            
            num_lines=text.count("\n")
            print '*  Removed '+str(num_lines)+' particles from selfile because they were too near to the border'
            fh=open(posname,"r")
            positions=fh.readlines()
	    fh.close()
            fh=open(selname2,"r")
            text = fh.readlines()
	    fh.close()
            newtext=[]
            newpos=[]
            newpos.append(positions[0])
            for i in range(len(text)):
                args=text[i].split(' ')
                if args[1].find('-1') > -1):
                    os.remove(args[0])
                else:
                    newtext.append(text[i])
                    newpos.append(positions[i+1])
                    self.allselfile.append(text[i])
            fh=open(selname2, 'w')
	    fh.writelines(newtext)
	    fh.close()
            fh=open(posname, 'w')
	    fh.writelines(newpos)
	    fh.close()
            outselfname=self.ProjectDir+'/'+self.OutSelFile
            fh=open(outselfname, 'w')
	    fh.writelines(self.allselfile)
	    fh.close()


        def perform_ctfestimate(self):
            import os
            iname=self.shortname+'/'+self.downname+'.raw'
            pname=self.shortname+'/'+self.shortname+'_input.param'
            posname=self.shortname+'/'+self.downname+'.'+PosFile 
	    selname=self.shortname+'/'+self.downname+'.raw.sel' 
            ctfname=self.downname+'.raw.ctf.sel' 
            ctfname2=self.shortname+'/'+self.downname+'.raw.ctf.sel' 
            print '*********************************************************************'
            print '*  Estimate CTF for micrograph: '+iname

            # prepare parameter file
            paramlist = []
            paramlist.append('image= '+iname+'\n')
            paramlist.append('micrograph_averaging= yes \n')
            paramlist.append('voltage= '+str(self.Voltage)+'\n')
            paramlist.append('spherical_aberration= '+str(self.SphericalAberration)+'\n')
            paramlist.append('sampling_rate= '+str(self.Sampling)+'\n')
            paramlist.append('particle_horizontal= 128 \n')
	    paramlist.append('Q0= '+str(self.AmplitudeContrast)+'\n')
            paramlist.append('N_horizontal= 512 \n')
	    paramlist.append('selfile= '+selname+'\n')
            paramlist.append('picked='+posname+'\n')
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


        def perform_normalize(self):
            import os
            iname=self.shortname+'/'+self.downname+'.raw.sel'
            print '*********************************************************************'
            print '*  Normalize particles in: '+iname
            param=' -i ' +iname+' -background circle '+str(self.BackGroundRadius)
            if (self.DoUseRamp):
                param=param+' -method Ramp'
            if (self.DoRemoveBlackDust):
                param=param+'  -remove_black_dust'
            if (self.DoRemoveWhiteDust):
                param=param+'  -remove_white_dust'
            command='xmipp_normalize '+param
            print '* ',command
            self.log.info(command)
            os.system(command)


        def perform_phaseflip(self):
            import os
            import spider_header
            selname=self.shortname+'/'+self.downname+'.raw.sel' 
            ctfparam=self.shortname+'/'+self.downname+'_Periodogramavg.ctfparam'

            # Check that images are 2 times self.Size, if not: re-extract and re-normalize them
            fh=open(selname,"r")
            text = fh.readlines()
	    fh.close()
            args=text[0].split()
            myheader=spider_header.spiderheader(args[0])
            if not myheader.nx==2*self.Size:
                self.perform_extract(True)
                self.perform_normalize()
                
            # Perform phase flipping operation
            print '*********************************************************************'
            print '*  Flipping phases for images in: '+selname
            command='xmipp_ctf_correct_phase -i '+selname+' -ctf '+ctfparam
            print '* ',command
            self.log.info(command)
            os.system(command)

            # Re-window the images to their original size
            command='xmipp_window -i '+selname+' -size '+str(self.Size)
            print '* ',command
            self.log.info(command)
            os.system(command)
            
        def visualize_ctfs(self):
            import os
            print '*********************************************************************'
            print '*  Visualizing all CTFs: '
            command='xmipp_show -ctfsel all_ctfs.sel &'
            print '* ',command
            self.log.info(command)
            os.system(command )

        def perform_sort_junk(self):
            import os
            print '*********************************************************************'
            print '*  Sorting images to identify junk particles in: '+self.OutSelFile
            command='xmipp_sort_by_statistics -i '+self.OutSelFile
            print '* ',command
            self.log.info(command)
            os.system(command)
	    if (self.DoVisualizeSorting):
		    command='xmipp_show -sel sort_junk.sel & '
		    print '* ',command
		    self.log.info(command)
		    os.system(command)


        def process_all_micrographs(self):
            import os
            import glob
            import SelFiles
            print '*********************************************************************'
            print '*  Pre-processing micrographs in '+str(self.MicrographSelfile)

            dirname=self.ProjectDir+'/'+self.ImagesDir
            if not os.path.exists(dirname):
                os.makedirs(dirname)

	    self.ctfselfile = []
            self.allselfile = []
            self.allctflibfile = []
            mysel=SelFiles.selfile()
            mysel.read(self.MicrographSelfile)
            for line in mysel.sellines:
                name,state=line[:-1].split(" ")
                self.shortname,self.downname=os.path.split(name)
                self.downname=self.downname.replace('.raw','')

                if (state.find('-1') < 0):
                    if (self.check_have_marked()):
                        if (self.DoExtract):
                            self.perform_extract(self.DoPhaseFlipping)

                            # Normalize before phase flipping to remove dust particles
                        if (self.DoNormalize):
                            self.perform_normalize()
                                  
                        if (self.DoEstimateCTF):
                            self.perform_ctfestimate()

                        if (self.DoPhaseFlipping):
                            self.perform_phaseflip()

            if (self.DoVisualizeCTF):
                self.visualize_ctfs()

            if (self.DoSorting):
                self.perform_sort_junk()

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

	preprocess_particles=preprocess_particles_class(MicrographSelfile,
                                                        ProjectDir,
                                                        LogDir,
                                                        DoExtract,
                                                        PosFile,
                                                        Size,
                                                        ImagesDir,
                                                        OutSelFile,
                                                        DoNormalize,
                                                        BackGroundRadius,
                                                        DoUseRamp,
                                                        DoRemoveBlackDust,
                                                        DoRemoveWhiteDust,
                                                        DoEstimateCTF,
                                                        Voltage,
                                                        SphericalAberration,
                                                        Sampling,
                                                        AmplitudeContrast,
                                                        DoVisualizeCTF,
                                                        OutCTFSelFile,
                                                        DoPhaseFlipping,
                                                        DoSorting,
                                                        DoVisualizeSorting)

	# close 
	preprocess_particles.close()

