#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
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
# {dir} Working subdirectory:
""" Use the same directory where you executed protocol_preprocess_micrographs.py
"""
WorkingDir="Preprocessing"
# {file} Selfile with micrographs on which to perform processing
MicrographSelfile='Preprocessing/all_micrographs.sel'
# Is this selfile a list of untilted-tilted pairs?
""" True for RCT-processing. In that case, provide a 3-column selfile as follows:
    untilted_pair1.raw tilted_pair1.raw 1
    untilted_pair2.raw tilted_pair2.raw 1
    etc...

    Note that for pair lists, only extraction and normalization is performed (ie no CTF or sorting)
    Also note that the extension of the coordinates files has to be raw.Common.pos!
"""
IsPairList=False
# {expert} Name for the output selfile:
""" This name should have extension .sel
"""
OutSelFile='all_images.sel'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir="/home2/bioinfo/scheres/work/protocols/"
# {expert} Directory name for logfiles:
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Extract particles
#------------------------------------------------------------------------------------------------
# Extract particles from micrographs?
DoExtract=True
# Family name for the picked coordinates file
""" By default this is Common, and the posfiles are called <mic>.raw.Common.pos
    These files are supposed to be in the micrograph-related subdirectories
"""
PosFile="Common"
# Dimension of the particles to extract (in pix.)
Size=64
# {expert} Directory name for particle images:
""" This directory will be placed in the project directory
"""
ImagesDir="Images"
#------------------------------------------------------------------------------------------------
# {section} Normalization
#------------------------------------------------------------------------------------------------
# Perform particle normalization?
DoNormalize=True
# Pixels outside this circle are assumed to be noise and their stddev is set to 1.
# Radius for background circle definition (in pix.)
BackGroundRadius=30
# Perform ramping background correction?
""" Correct for inclined background densities by fitting a least-squares plane through the background pixels
"""
DoUseRamp=True
# Perform black dust particles removal?
""" Sets pixels with unusually low values (i.e. stddev<-3.5) to zero
"""
DoRemoveBlackDust=False
# Perform white dust particles removal?
""" Sets pixels with unusually high values (i.e. stddev>3.5) to zero
"""
DoRemoveWhiteDust=False
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
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_preprocess_particles.py"
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
                 IsPairList,
                 ProjectDir,
                 LogDir,
                 DoExtract,
                 PosFile,
                 Size,
                 ImagesDir,
                 DoNormalize,
                 BackGroundRadius,
                 DoUseRamp,
                 DoRemoveBlackDust,
                 DoRemoveWhiteDust,
                 DoPhaseFlipping,
                 DoSorting,
                 ):
	     
        import os,sys
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,protocol_preprocess_micrographs
        
        self.WorkingDir=WorkingDir
        self.MicrographSelfile=os.path.abspath(MicrographSelfile)
        self.IsPairList=IsPairList
        self.ProjectDir=ProjectDir
        self.LogDir=LogDir
        self.DoExtract=DoExtract
        self.PosFile=PosFile
        self.Size=Size
        self.ImagesDir=ImagesDir
        self.DoNormalize=DoNormalize
        self.BackGroundRadius=BackGroundRadius
        self.DoUseRamp=DoUseRamp
        self.DoRemoveBlackDust=DoRemoveBlackDust
        self.DoRemoveWhiteDust=DoRemoveWhiteDust
        self.DoPhaseFlipping=DoPhaseFlipping
        self.DoSorting=DoSorting
        self.OutSelFile=OutSelFile

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
        
        # Parameters set from outside
        self.Down=protocol_preprocess_micrographs.Down
        
        # Check pairlist file
        if (not self.IsPairList):
            self.process_all_micrographs()
        else:
            if (self.DoPhaseFlipping):
                message='Warning: CTF phases will not be flipped for tilted pairs!'
                print '*',message
                self.log.error(message)
            if (self.DoSorting):
                message='Warning: images will not be sorted for tilted pairs!'
                print '*',message
                self.log.error(message)
            if (not self.PosFile=="Common"):
                message='Error: for tilted pairs The coordinate family name has to be Common!'
                print '*',message
                self.log.error(message)
                sys.exit()

            self.process_all_pairs()

        # Return to parent dir
        os.chdir(os.pardir)
            
    def process_all_micrographs(self):
        import os
        import selfile
        print '*********************************************************************'
        print '*  Pre-processing micrographs in '+os.path.basename(self.MicrographSelfile)

        dirname=self.ProjectDir+'/'+self.ImagesDir
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        self.allselfile = []
        self.allctfdatfile = []

        mysel=selfile.selfile()
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
                              
                    if (self.DoPhaseFlipping):
                        self.perform_phaseflip()

        if (self.DoSorting):
            self.perform_sort_junk()

    def process_all_pairs(self):
        import os
        print '*********************************************************************'
        print '*  Pre-processing all micrograph pairs in '+os.path.basename(self.MicrographSelfile)

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
        posname=self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.pos'
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
        rootname2=imgsubdir2+'/'+self.shortname2+'_'
        selname=self.downname+'.raw.sel' 
        selname2=self.downname2+'.raw.sel' 
        selnameb=self.shortname+'/'+self.downname+'.raw.sel' 
        selnameb2=self.shortname2+'/'+self.downname2+'.raw.sel' 
        posname=self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.pos'
        posname2=self.shortname2+'/'+self.downname2+'.raw.'+self.PosFile+'.pos'
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

        command='xmipp_micrograph_scissor -i ' + iname + ' -root ' + rootname + \
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
        self.remove_empty_images_pairs(selnameb,selnameb2)
        

    def remove_empty_images_pairs(self,selname,selname2):
        import os
        newsel=[]
        newsel2=[]

        count = 0
        # read old selfile
        fh  = open(selname,"r")
        sel = fh.readlines()
        fh.close()
        fh   = open(selname2,"r")
        sel2 = fh.readlines()
        fh.close()
        for i in range(len(sel)):
            args=sel[i].split()
            args2=sel2[i].split()
            if ((args[1].find('-1') > -1) or (args2[1].find('-1') > -1)):
                # remove empty images
                if (os.path.exists(args[0])):
                    os.remove(args[0])
                if (os.path.exists(args2[0])):
                    os.remove(args2[0])
                count += 1
            else:
                # or append to new selfile
                newsel.append(sel[i])
                newsel2.append(sel2[i])
                # For allselfiles, use relative paths wrt ProjectDir
                name= self.ImagesDir+'/'+self.shortname +'/'+os.path.basename(args[0]) +" 1\n"
                name2=self.ImagesDir+'/'+self.shortname2+'/'+os.path.basename(args2[0])+" 1\n"
                self.allselfile.append(name)
                self.allselfile2.append(name2)

        # write new selfiles
        fh=open(selname, 'w')
        fh.writelines(newsel)
        fh.close()
        fh=open(selname2, 'w')
        fh.writelines(newsel2)
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
        posname=self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.pos' 
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

        command='xmipp_micrograph_scissor -i ' + iname + ' -pos ' + posname + \
                 ' -root ' + rootname + ' -Xdim ' + str(size) + \
                 '|grep "corresponding image is set to blank"> ' + logname
        print '* ',command
        self.log.info(command)
        os.system(command)
        
        # Move selfile inside the subdirectory
        os.rename(selname,selnameb)

        # Remove particles near the border:
        self.remove_empty_images(selnameb)

    def remove_empty_images(self,selname):
        import os
        newsel=[]
        count = 0

        # read old selfile
        fh  = open(selname,"r")
        sel = fh.readlines()
        fh.close()
        for i in range(len(sel)):
            args=sel[i].split()
            if (args[1].find('-1') > -1):
                # remove empty image
                os.remove(args[0])
                count += 1
            else:
                # or append to new selfile
                newsel.append(sel[i])
                # Update all_images.sel     (relative paths wrt ProjectDir)
                name=self.ImagesDir+'/'+self.shortname+'/'+os.path.basename(args[0])
                self.allselfile.append(name+" 1\n")
                # Update all_images_ctf.dat (relative paths from ProjectDir)
                self.allctfdatfile.append(name + ' ' + \
                                          self.WorkingDir + '/' + \
                                          self.shortname + '/' + \
                                          self.downname + '_Periodogramavg.ctfparam\n')

        # write new selfile
        fh=open(selname, 'w')
        fh.writelines(newsel)
        fh.close()
        
        # Write updated all_images selfile
        outselfname=self.ProjectDir+'/'+self.OutSelFile
        fh=open(outselfname, 'w')
        fh.writelines(self.allselfile)
        fh.close()

        # Write updated all_images_ctf.dat
        ctfdatname=self.ProjectDir+'/'+self.OutSelFile.replace('.sel','.ctfdat')
        fh=open(ctfdatname,'w')
        fh.writelines(self.allctfdatfile)
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


    def perform_phaseflip(self):
        import os
        import spider_header, selfile

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
            
        # Make ctfdat file for this micrograph
        mysel=selfile.selfile()
        mysel.read(selname)
        ctfdatname=selname.replace('.sel','.ctfdat')
        mysel.write_ctfdat(ctfdatname,ctfparam)
        
        # Perform phase flipping operation
        print '*********************************************************************'
        print '*  Flipping phases for images in: '+ctfdatname
        command='xmipp_ctf_correct_phase -ctfdat '+ctfdatname
        print '* ',command
        self.log.info(command)
        os.system(command)

        # Re-window the images to their original size
        command='xmipp_window -i '+selname+' -size '+str(self.Size)
        print '* ',command
        self.log.info(command)
        os.system(command)
        
    def perform_sort_junk(self):
        import os
        print '*********************************************************************'
        print '*  Sorting images by statistics in: '+self.OutSelFile
        os.chdir(self.ProjectDir)
        command='xmipp_sort_by_statistics -i '+self.OutSelFile
        print '* ',command
        self.log.info(command)
        os.system(command)
        os.chdir(os.pardir)

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
                                                        IsPairList,
                                                        ProjectDir,
                                                        LogDir,
                                                        DoExtract,
                                                        PosFile,
                                                        Size,
                                                        ImagesDir,
                                                        DoNormalize,
                                                        BackGroundRadius,
                                                        DoUseRamp,
                                                        DoRemoveBlackDust,
                                                        DoRemoveWhiteDust,
                                                        DoPhaseFlipping,
                                                        DoSorting)

	# close 
	preprocess_particles.close()

