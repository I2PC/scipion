#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
#  - extraction of particles
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
""" Use the same directory where you executed xmipp_protocol_preprocess_micrographs.py
"""
WorkingDir='Preprocessing'
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
ProjectDir='/home/coss/temp/F22_cib'
# {expert} Directory name for logfiles:
LogDir='Logs'
#------------------------------------------------------------------------------------------------
# {section} Extract particles
#------------------------------------------------------------------------------------------------
# Extract particles from micrographs?
DoExtract=True
# Family name for the picked coordinates file
""" By default this is Common, and the posfiles are called <mic>.raw.Common.pos
    These files are supposed to be in the micrograph-related subdirectories
"""
PosFile='Common'
# Dimension of the particles to extract (in pix.)
Size=80
# {expert} Directory name for particle images:
""" This directory will be placed in the project directory
"""
ImagesDir='Images'
#------------------------------------------------------------------------------------------------
# {section} Normalization
#------------------------------------------------------------------------------------------------
# Perform particle normalization?
DoNormalize=False
# Pixels outside this circle are assumed to be noise and their stddev is set to 1.
# Radius for background circle definition (in pix.)
BackGroundRadius=30
# Perform ramping background correction?
""" Correct for inclined background densities by fitting a least-squares plane through the background pixels
"""
DoUseRamp=False
# Perform black dust particles removal?
""" Sets pixels with unusually low values to random values from a Gaussian with zero-mean and unity-standard deviation. For cryo, the default threshold value (see expert options) of 3.5 is a good value. For negative stain with high contrast signals, a higher value may be preferable.
"""
DoRemoveBlackDust=False
# Perform white dust particles removal?
""" Sets pixels with unusually high values to random values from a Gaussian with zero-mean and unity-standard deviation. For cryo, the default threshold value (see expert options) of 3.5 is a good value. For negative stain with high contrast signals, a higher value may be preferable.
"""
DoRemoveWhiteDust=False
# {expert} Threshold for dust removal:
""" Pixels with a signal higher or lower than this value times the standard devaition of the image will be affected. For cryo, 3.5 is a good value. For high-contrast negative stain, the signal itself may be affected so that a higher value may be preferable.
"""
DustRemovalThreshold=3.5
#------------------------------------------------------------------------------------------------
# {section} Particle sorting
#------------------------------------------------------------------------------------------------
# Perform particle sorting to identify junk particles?
DoSorting=False
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_preprocess_particles.py'
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
                 DustRemovalThreshold,
                 DoSorting,
                 ):
	     
        import os,sys
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,xmipp_protocol_preprocess_micrographs
        
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
        self.DustRemovalThreshold=DustRemovalThreshold
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
        self.Down=xmipp_protocol_preprocess_micrographs.Down
        
        # Check pairlist file
        if (not self.IsPairList):
            self.process_all_micrographs()
        else:
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
            self.shortname,self.allname=os.path.split(name)
            self.downname=self.allname.replace('.raw','')
            self.downname=self.downname.replace('.spi','')
            if (state.find('-1') < 0):
                if (self.check_have_marked()):
                    if (self.DoExtract):
                        self.perform_extract()

                    if (self.DoNormalize):
                        self.perform_normalize(self.shortname+'/'+self.allname+'.sel')
                             
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
        self.allselfileboth = []
        self.allctfdatfile = []

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
            self.shortname,self.allname=os.path.split(words[0])
            self.shortname2,self.allname2=os.path.split(words[1])
            self.downname=self.allname.replace('.raw','')
            self.downname2=self.allname2.replace('.raw','')
            self.downname=self.downname.replace('.spi','')
            self.downname2=self.downname2.replace('.spi','')
            state=words[2]
            if (state.find('-1') < 0):
                if (self.check_have_marked()):
                    if (self.DoExtract):
                        self.perform_extract_pairs()

                    if (self.DoNormalize):
                        self.perform_normalize(self.shortname+'/'+self.allname+'.sel')
                        self.perform_normalize(self.shortname2+'/'+self.allname2+'.sel')

        if (self.DoSorting):
            self.perform_sort_junk()

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
        posname=self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.auto.pos'
        if os.path.exists(posname):
            return True
        else:
            return False

    def perform_extract_pairs(self):
        import os,shutil
        import launch_job
        iname=self.shortname+'/'+self.allname
        iname2=self.shortname2+'/'+self.allname2
        imgsubdir=self.ProjectDir+'/'+self.ImagesDir+'/'+self.shortname
        imgsubdir2=self.ProjectDir+'/'+self.ImagesDir+'/'+self.shortname2
        rootname=imgsubdir+'/'+self.shortname+'_' 
        rootname2=imgsubdir2+'/'+self.shortname2+'_'
        selname=self.allname+'.sel' 
        selname2=self.allname2+'.sel' 
        selnameb=self.shortname+'/'+self.allname+'.sel' 
        selnameb2=self.shortname2+'/'+self.allname2+'.sel' 
        # micrograph_mark with tilt pairs expects <micrographname>.Common.pos files
        unpos   = self.shortname+'/'+self.allname+'.'+self.PosFile+'.pos'
        tilpos  = self.shortname2+'/'+self.allname2+'.'+self.PosFile+'.pos'
        posname = self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.pos'
        posname2= self.shortname2+'/'+self.downname2+'.raw.'+self.PosFile+'.pos'
        if (not os.path.exists(unpos)):
            shutil.copy(posname,unpos)
        if (not os.path.exists(tilpos)):
            shutil.copy(posname2,tilpos)
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

        command=' -i ' + iname + ' -root ' + rootname + \
                 ' -tilted ' + iname2 + ' -root_tilted ' + rootname2 + \
                 ' -Xdim ' + str(self.Size) + \
                 '|grep "corresponding image is set to blank"> ' + logname            
        launch_job.launch_job("xmipp_micrograph_scissor",
                              command,
                              self.log,
                              False,1,1,'')

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
                name= self.ImagesDir+'/'+self.shortname +'/'+os.path.basename(args[0]) 
                name2=self.ImagesDir+'/'+self.shortname2+'/'+os.path.basename(args2[0])
                self.allselfile.append(name + " 1\n")
                self.allselfile2.append(name2 + " 1\n")
                self.allselfileboth.append(name + " 1\n")
                self.allselfileboth.append(name2 + " 1\n")
                # Update all_images.ctfdat (relative paths from ProjectDir)
                self.allctfdatfile.append(name + ' ' + \
                                          self.WorkingDir + '/' + \
                                          self.shortname + '/' + \
                                          self.downname + '_Periodogramavg.ctfparam\n')
                self.allctfdatfile.append(name2 + ' ' + \
                                          self.WorkingDir + '/' + \
                                          self.shortname2 + '/' + \
                                          self.downname2 + '_Periodogramavg.ctfparam\n')



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

        # Also write out a self.OutSelFile for sorting
        fh=open(self.ProjectDir+'/'+self.OutSelFile, 'w')
        fh.writelines(self.allselfileboth)
        fh.close()

        # Write updated all_images.ctfdat
        ctfdatname=self.ProjectDir+'/'+self.OutSelFile.replace('.sel','.ctfdat')
        fh=open(ctfdatname,'w')
        fh.writelines(self.allctfdatfile)
        fh.close()
 
        # Output to screen
        message='Removed '+str(count)+' pairs from selfiles because at least one of the particles was too near to the border'
        print '* ',message
        self.log.info(message)


    def perform_extract(self):
        import os
        import launch_job
        iname=self.shortname+'/'+self.allname
        selname=self.allname+'.sel' 
        selnameb=self.shortname+'/'+self.allname+'.sel' 
        posname=self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.pos' 
        posnameauto=self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.auto.pos' 
        imgsubdir=self.ProjectDir+'/'+self.ImagesDir+'/'+self.shortname
        rootname=imgsubdir+'/'+self.shortname+'_'
        logname=self.shortname+'/scissor.log'
        size=self.Size

        # Make directory if necessary
        if not os.path.exists(imgsubdir):
            os.makedirs(imgsubdir)

        if os.path.exists(posname) and os.path.exists(posnameauto):
            launch_job.launch_job("cat",
                                  posname+" "+posnameauto+" > "+iname+"_all.pos",
                                  self.log,
                                  False,1,1,'')
        elif os.path.exists(posname):
            launch_job.launch_job("cp",
                                  posname+" "+iname+"_all.pos",
                                  self.log,
                                  False,1,1,'')
        elif os.path.exists(posnameauto):
            launch_job.launch_job("cp",
                                  posnameauto+" "+iname+"_all.pos",
                                  self.log,
                                  False,1,1,'')
        
        command= ' -i ' + iname + ' -pos ' + iname+"_all.pos" + \
                 ' -root ' + rootname + ' -Xdim ' + str(size) + \
                 '|grep "corresponding image is set to blank"> ' + logname
        launch_job.launch_job("xmipp_micrograph_scissor",
                              command,
                              self.log,
                              False,1,1,'')
        launch_job.launch_job("rm",
                              iname+"_all.pos",
                              self.log,
                              False,1,1,'')
        
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
                # Update all_images.ctfdat (relative paths from ProjectDir)
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

        # Write updated all_images.ctfdat
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
        import launch_job
        print '*********************************************************************'
        print '*  Normalize particles in: '+iname
        param=' -i ' +iname+' -background circle '+str(self.BackGroundRadius)
        if (self.DoUseRamp):
            param=param+' -method Ramp'
        if (self.DoRemoveBlackDust):
            param=param+'  -remove_black_dust -thr_black_dust -' + \
                str(self.DustRemovalThreshold)
        if (self.DoRemoveWhiteDust):
            param=param+'  -remove_white_dust -thr_white_dust ' + \
                str(self.DustRemovalThreshold)
        launch_job.launch_job("xmipp_normalize",
                              param,
                              self.log,
                              False,1,1,'')

    def perform_sort_junk(self):
        import os
        import launch_job
        print '*********************************************************************'
        print '*  Sorting images by statistics in: '+self.OutSelFile
        os.chdir(self.ProjectDir)
        command=' -i '+self.OutSelFile
        launch_job.launch_job("xmipp_sort_by_statistics",
                              command,
                              self.log,
                              False,1,1,'')
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
                                                        DustRemovalThreshold,
                                                        DoSorting)

	# close 
	preprocess_particles.close()

