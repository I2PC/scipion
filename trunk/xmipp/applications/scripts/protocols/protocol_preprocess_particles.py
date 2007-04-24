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
# Is this selfile a list of untilted-tilted pairs?
""" True for RCT-processing. In that case, provide a 3-column selfile as follows:
    untilted_pair1.raw tilted_pair1.raw 1
    untilted_pair2.raw tilted_pair2.raw 1
    etc...

    Note that for pair lists, only extraction and normalization is performed (ie no CTF or sorting)
    Also note that the extension of the coordinates files has to be raw.Common.pos!
"""
IsPairList=False
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
""" For tilted pairs, extension .sel is replaced by _untilted.sel and _tilted.sel
"""
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
DoEstimateCTF=False
# Microscope voltage (in kV)
Voltage=200
# Spherical aberration
SphericalAberration=2.26
# Pixel size in the particles (in Angstrom/pixel)
Sampling=2.8
# Amplitude Contrast (typically negative!)
AmplitudeContrast=-0.1
# Visualize CTF-estimation?
DoVisualizeCTF=True
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
                 IsPairList,
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
                 LowResolCutoff,
                 HighResolCutoff,
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
        self.IsPairList=IsPairList
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
        self.LowResolCutoff=LowResolCutoff
        self.HighResolCutoff=HighResolCutoff
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
        
        # Check pairlist file
        if (not self.IsPairList):
            self.process_all_micrographs()
        else:
            if (self.DoEstimateCTF):
                message='Warning: CTFs will not be estimated for tilted pairs!'
                print '*',message
                self.log.error(message)
            if (self.DoPhaseFlipping):
                message='Warning: CTF phases will not be flipped for tilted pairs!'
                print '*',message
                self.log.error(message)
            if (self.DoSorting):
                message='Warning: images will not be sorted for tilted pairs!'
                print '*',message
                self.log.error(message)
            if (not self.PosFile=="raw.Common.pos"):
                message='Error: for tilted pairs The coordinate file extension has to be raw.Common.pos!'
                print '*',message
                self.log.error(message)
                sys.exit()

            self.process_all_pairs()

            
    def process_all_micrographs(self):
        import os
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
        for name,state in mysel.sellines:
            self.shortname,self.downname=os.path.split(name)
            self.downname=self.downname.replace('.raw','')
            if (state.find('-1') < 0):
                if (self.check_have_marked()):
                    if (self.DoExtract):
                        self.perform_extract(self.DoPhaseFlipping)

                    # Normalize before phase flipping to remove dust particles
                    if (self.DoNormalize):
                        self.perform_normalize(self.shortname+'/'+self.downname+'.raw.sel')
                              
                    if (self.DoEstimateCTF):
                        self.perform_ctfestimate()

                    if (self.DoPhaseFlipping):
                        self.perform_phaseflip()

        if (self.DoVisualizeCTF):
            self.visualize_ctfs()

        if (self.DoSorting):
            self.perform_sort_junk()

    def process_all_pairs(self):
        import os
        print '*********************************************************************'
        print '*  Pre-processing all micrograph pairs in '+str(self.MicrographSelfile)

        dirname=self.ProjectDir+'/'+self.ImagesDir
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        self.allselfile = []
        self.allselfile2 = []

        fh=open(self.MicrographSelfile,'r')
        self.pairlines=fh.readlines()
        fh.close()
        words=self.pairlines[0].split()
        if (len(words)<3):
            message='Error: Selfile is not a pairlist file!'
            print '*',message
            self.log.error(message)
            sys.exit()
        for line in self.pairlines:
            words=line.split()
            self.shortname,self.downname=os.path.split(words[0])
            self.shortname2,self.downname2=os.path.split(words[1])
            self.downname=self.downname.replace('.raw','')
            self.downname2=self.downname2.replace('.raw','')
            state=words[2]
            if (state.find('-1') < 0):
                if (self.check_have_marked()):
                    if (self.DoExtract):
                        self.perform_extract_pairs()

                    if (self.DoNormalize):
                        self.perform_normalize(self.shortname+'/'+self.downname+'.raw.sel')
                        self.perform_normalize(self.shortname2+'/'+self.downname2+'.raw.sel')


    def check_file_exists(self,name):
        import os
        if not os.path.exists(name):
            message='Error: File '+name+' does not exist, exiting...'
            print '*',message
            self.log.error(message)
            sys.exit()


    def check_have_marked(self):
        import os
        posname=self.shortname+'/'+self.downname+'.'+self.PosFile 
        if os.path.exists(posname):
            return True
        else:
            return False


    def perform_extract_pairs(self):
        import os
        
        iname=self.shortname+'/'+self.downname+'.raw' 
        iname2=self.shortname2+'/'+self.downname2+'.raw' 
        imgsubdir=self.ProjectDir+'/'+self.ImagesDir+'/'+self.shortname
        imgsubdir2=self.ProjectDir+'/'+self.ImagesDir+'/'+self.shortname2
        rootname=imgsubdir+'/'+self.shortname+'_' 
        rootname2=imgsubdir+'/'+self.shortname2+'_'
        selname=self.downname+'.raw.sel' 
        selname2=self.downname2+'.raw.sel' 
        selnameb=self.shortname+'/'+self.downname+'.raw.sel' 
        selnameb2=self.shortname2+'/'+self.downname2+'.raw.sel' 
        posname=self.shortname+'/'+self.downname+'.'+self.PosFile
        posname2=self.shortname2+'/'+self.downname2+'.'+self.PosFile
        angname=self.shortname+'/'+self.downname+'.ang'
        logname=self.shortname+'/scissor.log'

        # Make directories if necessary
        if not os.path.exists(imgsubdir):
            os.makedirs(imgsubdir)
        if not os.path.exists(imgsubdir2):
            os.makedirs(imgsubdir2)

        # Check existence of second pos and ang files
        self.check_file_exists(posname)
        self.check_file_exists(posname2)
        self.check_file_exists(angname)

#        command='xmipp_micrograph_scissor -i ' + iname + ' -root ' + rootname + \
        command='xmipp_scissor -i ' + iname + ' -root ' + rootname + \
                 ' -tilted ' + iname2 + ' -root_tilted ' + rootname2 + \
                 ' -Xdim ' + str(self.Size) + \
                 '|grep "corresponding image is set to blank"> ' + logname            
        print '* ',command
        self.log.info(command)
        os.system(command)

        # Move output selfiles inside the sub-directory:
        os.rename(selname,selnameb)
        os.rename(selname2,selnameb2)

        # Remove pairs with one image near the border
        self.remove_empty_images_pairs(posname,posname2,selnameb,selnameb2)
        

    def remove_empty_images_pairs(self,posname,posname2,selname,selname2):
        import os
        newsel=[]
        newpos=[]
        newsel2=[]
        newpos2=[]

        count = 0
        # read old selfile and posfile
        fh  = open(posname,"r")
        pos = fh.readlines()
        fh.close()
        fh   = open(posname2,"r")
        pos2 = fh.readlines()
        fh.close()
        fh  = open(selname,"r")
        sel = fh.readlines()
        fh.close()
        fh   = open(selname2,"r")
        sel2 = fh.readlines()
        fh.close()
        newpos.append(pos[0]) # append header line
        newpos2.append(pos2[0]) # append header line
        for i in range(len(sel)):
            args=sel[i].split(' ')
            args2=sel2[i].split(' ')
            if ((args[1].find('-1') > -1) or (args2[1].find('-1') > -1)):
                # remove empty images
                if (os.path.exists(args[0])):
                    os.remove(args[0])
                if (os.path.exists(args2[0])):
                    os.remove(args2[0])
                count += 1
            else:
                # or append to new selfile and posfile
                newsel.append(sel[i])
                newsel2.append(sel2[i])
                newpos.append(pos[i+1])
                newpos2.append(pos2[i+1])
                self.allselfile.append(sel[i])
                self.allselfile2.append(sel2[i])
        # write new selfiles and posfiles
        fh=open(selname, 'w')
        fh.writelines(newsel)
        fh.close()
        fh=open(selname2, 'w')
        fh.writelines(newsel2)
        fh.close()
        fh=open(posname, 'w')
        fh.writelines(newpos)
        fh.close()
        fh=open(posname2, 'w')
        fh.writelines(newpos2)
        fh.close()
        # Update allselfiles
        outselfname=self.ProjectDir+'/'+self.OutSelFile.replace('.sel','_untilted.sel')
        fh=open(outselfname, 'w')
        fh.writelines(self.allselfile)
        fh.close()
        outselfname=self.ProjectDir+'/'+self.OutSelFile.replace('.sel','_tilted.sel')
        fh=open(outselfname, 'w')
        fh.writelines(self.allselfile2)
        fh.close()
        # Output to screen
        message='Removed '+str(count)+' pairs from selfiles because at least one of the particles was too near to the border'
        print '* ',message
        self.log.info(message)


    def perform_extract(self,is_for_flipping):
        import os
        iname=self.shortname+'/'+self.downname+'.raw' 
        selname=self.downname+'.raw.sel' 
        selnameb=self.shortname+'/'+self.downname+'.raw.sel' 
        posname=self.shortname+'/'+self.downname+'.'+self.PosFile 
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

#        command='xmipp_micrograph_scissor -i ' + iname + ' -pos ' + posname + \
        command='xmipp_scissor -i ' + iname + ' -pos ' + posname + \
                 ' -root ' + rootname + ' -Xdim ' + str(size) + \
                 '|grep "corresponding image is set to blank"> ' + logname
        print '* ',command
        self.log.info(command)
        os.system(command)
        os.rename(selname,selnameb)

        # Remove particles near the border:
        self.remove_empty_images(posname,selnameb)

    def remove_empty_images(self,posname,selname):
        import os
        newsel=[]
        newpos=[]
        count = 0

        # read old selfile and posfile
        fh  = open(posname,"r")
        pos = fh.readlines()
        fh.close()
        fh  = open(selname,"r")
        sel = fh.readlines()
        fh.close()
        newpos.append(pos[0])
        for i in range(len(sel)):
            args=sel[i].split(' ')
            if (args[1].find('-1') > -1):
                # remove empty image
                os.remove(args[0])
                count += 1
            else:
                # or append to new selfile and posfile
                newsel.append(sel[i])
                newpos.append(pos[i+1])
                self.allselfile.append(sel[i])
        # write new selfile and posfile
        fh=open(selname, 'w')
        fh.writelines(newsel)
        fh.close()
        fh=open(posname, 'w')
        fh.writelines(newpos)
        fh.close()
        # Write updated allselfile
        outselfname=self.ProjectDir+'/'+self.OutSelFile
        fh=open(outselfname, 'w')
        fh.writelines(self.allselfile)
        fh.close()
        # Output to screen
        message='Removed '+str(count)+' particles from selfile because they were too near to the border'
        print '* ',message
        self.log.info(message)


    def perform_normalize(self,iname):
        import os
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


    def perform_ctfestimate(self):
        import os
        iname=self.shortname+'/'+self.downname+'.raw'
        pname=self.shortname+'/'+self.shortname+'_input.param'
        posname=self.shortname+'/'+self.downname+'.'+self.PosFile 
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
        paramlist.append('min_freq= '+str(self.LowResolCutoff)+'\n')
        paramlist.append('max_freq= '+str(self.HighResolCutoff)+'\n')
        paramlist.append('selfile= '+selname+'\n')
        paramlist.append('picked='+posname+'\n')
        paramlist.append('periodogram= yes \n')

        # Perform CTF estimation
        fh=open(pname,"w")
        fh.writelines(paramlist)
        fh.close()
#        command='xmipp_ctf_estimate_from_micrograph -i '+pname
        command='xmipp_assign_CTF -i '+pname
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
            self.perform_normalize(self.shortname+'/'+self.downname+'.raw.sel')
            
        # Perform phase flipping operation
        print '*********************************************************************'
        print '*  Flipping phases for images in: '+selname
#        command='xmipp_ctf_correct_phase -i '+selname+' -ctf '+ctfparam
        command='xmipp_correctphase -i '+selname+' -ctf '+ctfparam
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
#        command='xmipp_sort_by_statistics -i '+self.ProjectDir+'/'+self.OutSelFile
        command='xmipp_sort_junk -i '+self.ProjectDir+'/'+self.OutSelFile
        print '* ',command
        self.log.info(command)
        os.system(command)
        if (self.DoVisualizeSorting):
    	    command='xmipp_show -sel sort_junk.sel & '
    	    print '* ',command
    	    self.log.info(command)
    	    os.system(command)


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
                                                        IsPairList,
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
                                                        LowResolCutoff,
                                                        HighResolCutoff,
                                                        DoVisualizeCTF,
                                                        OutCTFSelFile,
                                                        DoPhaseFlipping,
                                                        DoSorting,
                                                        DoVisualizeSorting)

	# close 
	preprocess_particles.close()

