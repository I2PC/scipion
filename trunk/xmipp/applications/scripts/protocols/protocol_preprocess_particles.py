#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
#  - phase flipping
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
#         Carlos Oscar, December 2010
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {dir} Working subdirectory:
""" Use the same directory where you executed xmipp_protocol_preprocess_micrographs.py
"""
WorkingDir='Images'

# {dir} Directory with the particle picking
PickingDir='ParticlePicking'

# {expert} Name for the output selfile:
""" This name should have extension .sel
"""
OutSelFile='all_images.sel'

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir='/home/coss/temp/F22_cib'

#------------------------------------------------------------------------------------------------
# {section} Processing parameters
#------------------------------------------------------------------------------------------------
# Box size of the particles to extract (in pix.)
Size=80

# Do phase flipping?
DoFlip=True

# {expert} Take Logarithm?
DoLog=False 

# Invert contrast?
DoInvert=False

# {expert} Background radius
"""Pixels outside this circle are assumed to be noise and their stddev is set to 1.
   Radius for background circle definition (in pix.).
   If this value is 0, then the same as the particle radius is used. """
BackGroundRadius=0

# Perform ramping background correction?
""" Correct for inclined background densities by fitting a least-squares plane through the background pixels
"""
DoUseRamp=True

# Perform dust particles removal?
""" Sets pixels with unusually large values to random values from a Gaussian with zero-mean and unity-standard deviation.
"""
DoRemoveDust=True

# {expert} Threshold for dust removal:
""" Pixels with a signal higher or lower than this value times the standard deviation of the image will be affected. For cryo, 3.5 is a good value. For high-contrast negative stain, the signal itself may be affected so that a higher value may be preferable.
"""
DustRemovalThreshold=3.5

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
AnalysisScript='visualize_preprocess_particles.py'
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#  {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#

class preprocess_particles_class:
    def getParameter(prm,filename):
        f = open(filename, 'r')
        lines=f.readlines()
        f.close()
        for line in lines:
            tokens=line.split('=')
            if tokens[0]==prm:
                return tokens[1]
        return ""

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

    #init variables
    def __init__(self,
                 WorkingDir,
                 PickingDir,
                 ProjectDir,
                 Size,
                 DoFlip,
                 DoLog, 
                 DoInvert,
                 BackGroundRadius,
                 DoUseRamp,
                 DoRemoveDust,
                 DustRemovalThreshold,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour
                 ):
	     
        import os,sys
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,xmipp_protocol_preprocess_micrographs
        
        self.WorkingDir=WorkingDir.strip()
        self.PickingDir=PickingDir.strip()
        self.ProjectDir=ProjectDir.strip()
        self.LogDir="Logs"
        self.PosFile="Common"
        self.Size=Size
        self.DoFlip=DoFlip
        self.DoLog=DoLog 
        self.DoInvert=DoInvert
        self.BackGroundRadius=BackGroundRadius
        self.DoUseRamp=DoUseRamp
        self.DoRemoveDust=DoRemoveBlackDust
        self.DustRemovalThreshold=DustRemovalThreshold
        self.OutSelFile=OutSelFile

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     self.LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Make working directory if it does not exist yet
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Save parameters and compare to possible previous runs
        self.saveAndCompareParameters([
                 "PickingDir",
                 "Size",
                 "BackGroundRadius",
                 "DoUseRamp",
                 "DoRemoveBlackDust",
                 "DoRemoveWhiteDust",
                 "DustRemovalThreshold"]);

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                       os.path.abspath(self.WorkingDir))
    
        # Preprocess paticles
        fnScript=self.WorkingDir+'/pickParticles.sh'
        self.fh_mpi=os.open(fnScript, os.O_WRONLY | os.O_TRUNC | os.O_CREAT, 0700)
        self.process_all_micrographs()
        os.close(self.fh_mpi)
        self.launchCommandFile(fnScript)

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

    def process_all_micrographs(self):
        import os
        print '*********************************************************************'
        print '*  Pre-processing micrographs in '+os.path.basename(self.PickingDir)

        fnPickingParameters=self.WorkingDir+"/protocolParameters.txt"
        isPairTilt=self.getParameter("IsPairTilt",fnPickingParameters)=="True"
        MicrographSelfile=self.getParameter("MicrographSelfile",fnPickingParameters)
        mD=xmipp.MetaData();
        xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, mD)
        for id in mD:
            micrograph=mD.getValue(xmipp.MDL_IMAGE)
            if isPairTilt:
                micrographTilted=mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1)
            
            # Phase flip
            command=""
            if self.DoFlip:
                command+="xmipp_"
            
                if (self.check_have_marked()):
                    if (self.DoExtract):
                        self.perform_extract()

                    if (self.DoNormalize):
                        self.perform_normalize(self.shortname+'/'+self.allname+'.sel')
                             
        self.perform_sort_junk()

    def process_all_pairs(self):
        import os
        print '*********************************************************************'
        print '*  Pre-processing all micrograph pairs in '+os.path.basename(self.PickingDir)

        dirname=self.ProjectDir+'/'+self.WorkingDir
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

        self.perform_sort_junk()

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
        imgsubdir=self.ProjectDir+'/'+self.WorkingDir+'/'+self.shortname
        imgsubdir2=self.ProjectDir+'/'+self.WorkingDir+'/'+self.shortname2
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
        
    def perform_extract(self):
        import os
        import launch_job
        iname=self.shortname+'/'+self.allname
        selname=self.allname+'.sel' 
        selnameb=self.shortname+'/'+self.allname+'.sel' 
        posname=self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.pos' 
        posnameauto=self.shortname+'/'+self.downname+'.raw.'+self.PosFile+'.auto.pos' 
        imgsubdir=self.ProjectDir+'/'+self.WorkingDir+'/'+self.shortname
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
    
    # Check that there is a valid list of micrographs
    if not os.path.exists(MicrographSelfile)>0:
        message="Cannot find "+MicrographSelfile
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    # Check that all micrographs exist
    import xmipp
    mD = xmipp.MetaData()
    xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, mD)
    message="Cannot find the following micrographs:\n"
    NnotFound=0
    for id in mD:
         micrograph = mD.getValue(xmipp.MDL_IMAGE)
         if not os.path.exists(micrograph):
            message+=micrograph+"\n"
            NnotFound=NnotFound+1
         if mD.containsLabel(xmipp.MDL_ASSOCIATED_IMAGE1):
             micrograph = mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1)
             if not os.path.exists(micrograph):
                 message+=micrograph+"\n"
                 NnotFound=NnotFound+1
    
    if not NnotFound>0:
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
   	# create preprocess_particles_class object
	preprocess_particles=preprocess_particles_class(
                 WorkingDir,
                 PickingDir,
                 ProjectDir,
                 Size,
                 DoFlip,
                 DoLog, 
                 DoInvert,
                 BackGroundRadius,
                 DoUseRamp,
                 DoRemoveDust,
                 DustRemovalThreshold,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour)
