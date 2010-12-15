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
def getParameter(prm,filename):
    f = open(filename, 'r')
    lines=f.readlines()
    f.close()
    for line in lines:
        tokens=line.split('=')
        if tokens[0]==prm:
            return tokens[1].strip()
    return ""

class preprocess_particles_class:

    def saveAndCompareParameters(self, listOfParameters):
        import os,shutil
        fnOut=self.WorkingDir + "/protocolParameters.txt"
        linesNew=[];
        for prm in listOfParameters:
            eval("linesNew.append('"+prm +"='+str("+prm+")+'\\n')")
        retval=False
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
                retval=True
        f = open(fnOut, 'w')
        f.writelines(linesNew)
        f.close()
        return retval

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
                 DoRemoveDust,
                 DustRemovalThreshold,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour
                 ):
	     
        import os,sys,time
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,xmipp,launch_job
        
        self.WorkingDir=WorkingDir.strip()
        self.PickingDir=PickingDir.strip()
        self.ProjectDir=ProjectDir.strip()
        self.LogDir="Logs"
        self.PosFile="Common"
        self.Size=Size
        self.DoFlip=DoFlip
        self.DoLog=DoLog 
        self.DoInvert=DoInvert
        if BackGroundRadius!=0:
            self.BackGroundRadius=BackGroundRadius
        else:
            self.BackGroundRadius=Size/2
        self.DoRemoveDust=DoRemoveDust
        self.DustRemovalThreshold=DustRemovalThreshold
        self.OutSelFile=OutSelFile
        self.DoParallel=DoParallel
        self.NumberOfMpiProcesses=NumberOfMpiProcesses
        self.SystemFlavour=SystemFlavour

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     self.LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Make working directory if it does not exist yet
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Save parameters and compare to possible previous runs
        deleted=self.saveAndCompareParameters([
                 "PickingDir",
                 "Size",
                 "DoFlip",
                 "DoInvert",
                 "DoLog",
                 "BackGroundRadius",
                 "DoRemoveDust",
                 "DustRemovalThreshold"]);

        # Update status
        fh=open(self.WorkingDir + "/status.txt", "a")
        fh.write("Step 0: Processed started at " + time.asctime() + "\n")
        fh.close()

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                       os.path.abspath(self.WorkingDir))
    
        # Preprocess paticles
        fnScript=self.WorkingDir+'/pickParticles.sh'
        self.fh_mpi=os.open(fnScript, os.O_WRONLY | os.O_TRUNC | os.O_CREAT, 0700)
        self.process_all_micrographs()
        os.close(self.fh_mpi)
        self.launchCommandFile(fnScript)

        # Join results
        generalMD=xmipp.MetaData()
        for selfile in self.outputSel:
            if not os.path.exists(selfile):
                fh=open(self.WorkingDir + "/status.txt", "a")
                fh.write("Step E: Cannot read "+selfile+". Finishing at " + time.asctime() + "\n")
                fh.close()
                sys.exit(1)
            MD=xmipp.MetaData(selfile)
            for id in MD:
                imageFrom=MD.getValue(xmipp.MDL_MICROGRAPH)
                MD.setValue(xmipp.MDL_MICROGRAPH,self.correspondingMicrograph[imageFrom])
                if (self.correspondingCTF[imageFrom]!=""):
                    MD.setValue(xmipp.MDL_CTFMODEL,self.correspondingCTF[imageFrom])
            generalMD.unionAll(MD)
        generalMD.write(self.OutSelFile)

        # Sort by statistics
        rootName,dummy=os.path.splitext(self.OutSelFile)
        launch_job.launch_job("xmipp_sort_by_statistics",
                              "-i "+self.OutSelFile+" -multivariate "+\
                              "-o "+rootName+"_sorted_by_score",
                              self.log,
                              False,1,1,'')

        # Remove intermediate selfiles
        for selfile in self.outputSel:
            if os.path.exists(selfile):
                os.remove(selfile)
    
        # Update status    
        if os.path.exists(rootName+"_sorted_by_score.sel"):
            fh=open(self.WorkingDir + "/status.txt", "a")
            fh.write("Step F: Processed finished at " + time.asctime() + "\n")
            fh.close()

    def launchCommandFile(self, commandFile):
        import launch_job, log, os
        log.cat(self.log, commandFile)
        if self.DoParallel:
            command=' -i ' + commandFile
            launch_job.launch_job("xmipp_run", command, self.log, True,
                  self.NumberOfMpiProcesses, 1, self.SystemFlavour)
        else:
            self.log.info(commandFile)     
            os.system(commandFile)     

    def process_all_micrographs(self):
        import os, xmipp
        print '*********************************************************************'
        print '*  Pre-processing micrographs in '+os.path.basename(self.PickingDir)

        self.outputSel=[]
        self.correspondingMicrograph={}
        self.correspondingCTF={}

        fnPickingParameters=self.PickingDir+"/protocolParameters.txt"
        isPairTilt=getParameter("IsPairTilt",fnPickingParameters)=="True"
        MicrographSelfile=getParameter("MicrographSelfile",fnPickingParameters)
        mD=xmipp.MetaData();
        xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, mD)
        preprocessingDir,dummy=os.path.split(MicrographSelfile)
        for id in mD:
            micrograph=mD.getValue(xmipp.MDL_IMAGE)
            dummy,micrographWithoutDirs=os.path.split(micrograph)
            if isPairTilt:
                micrographTilted=mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1)
            
            # Phase flip
            fnStack=self.WorkingDir+"/"+micrographWithoutDirs+".stk"
            command=''
            filesToDelete=[]
            fnToPick=preprocessingDir+"/"+micrograph
            if self.DoFlip and not isPairTilt:
                micrographName,micrographExt=os.path.splitext(micrographWithoutDirs)
                ctf=mD.getValue(xmipp.MDL_CTFMODEL)
                fnToPick=self.WorkingDir+"/"+micrographName+"_flipped.raw"
                command+="xmipp_micrograph_phase_flipping"+\
                         " -i "+preprocessingDir+"/"+micrograph+\
                         " -ctf "+preprocessingDir+"/"+ctf+\
                         " -o "+fnToPick + " ; "
                filesToDelete.append(fnToPick+"*")
            self.correspondingMicrograph[fnToPick]=preprocessingDir+"/"+micrograph
            if mD.containsLabel(xmipp.MDL_CTFMODEL):
                ctf=mD.getValue(xmipp.MDL_CTFMODEL)
                self.correspondingCTF[fnToPick]=preprocessingDir+"/"+ctf
            else:
                self.correspondingCTF[fnToPick]=""
            
            # Extract particles
            arguments=""
            if isPairTilt:
                pass
            else:
                posfile=""
                candidatePosFile=self.PickingDir+"/"+micrographWithoutDirs+"."+self.PosFile+".pos"
                if os.path.exists(candidatePosFile):
                    posfile=candidatePosFile
                candidatePosFile=self.PickingDir+"/"+micrographWithoutDirs+"."+self.PosFile+".auto.pos"
                if os.path.exists(candidatePosFile):
                    if posfile=="":
                        posfile=candidatePosFile
                    else:
                        MD1=xmipp.MetaData(posfile)
                        MD2=xmipp.MetaData(candidatePosFile)
                        MD1.unionAll(MD2)
                        candidatePosFile=self.PickingDir+"/"+micrographWithoutDirs+"."+self.PosFile+".both.pos"
                        MD1.write(candidatePosFile)
                        filesToDelete.append(candidatePosFile)
                        posfile=candidatePosFile
                if posfile!="":
                    fnStack=self.WorkingDir+"/"+micrographWithoutDirs+".stk"
                    self.outputSel.append(fnStack+".sel")
                    arguments="-i "+fnToPick+" --pos "+posfile+" -o "+fnStack
            if arguments=="":
                print "Cannot find positions for "+micrograph
                continue
            
            arguments+=" --Xdim "+str(self.Size)+" --rmStack"
            if self.DoInvert:
                arguments+=" --invert"
            if self.DoLog:
                arguments+=" --log"
            command+="xmipp_micrograph_scissor "+arguments+" ; "
            if isPairTilt:
                pass
            
            # Normalize particles
            normalizeArguments=\
                     ' -background circle '+str(self.BackGroundRadius)+\
                     ' -method Ramp'                
            if (self.DoRemoveDust):
                normalizeArguments+=' -thr_black_dust -' + str(self.DustRemovalThreshold)+\
                         ' -thr_white_dust ' + str(self.DustRemovalThreshold)
            command+='xmipp_normalize -i ' +fnStack+normalizeArguments
            if isPairTilt:
                pass

            # Remove temporary files
            for fileToDelete in filesToDelete:
                command+=" ; rm -f "+fileToDelete

            # Command done
            command += " ; if [ -e " + fnStack + ' ]; then ' + \
                        'echo "Step: '+micrograph+' processed " `date` >> ' + self.WorkingDir + "/status.txt; " + \
                       "fi"
            os.write(self.fh_mpi, command+"\n")
        
# Preconditions
def preconditions(gui):
    import os
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
    if not os.path.exists(PickingDir)>0:
        message="Cannot find "+PickingDir
        if gui:
            import tkMessageBox
            tkMessageBox.showerror("Error", message)
        else:
            print message
        retval=False
    
    # Check that all micrographs exist
    import xmipp
    fnPickingParameters=PickingDir+"/protocolParameters.txt"
    isPairTilt=getParameter("IsPairTilt",fnPickingParameters)=="True"
    MicrographSelfile=getParameter("MicrographSelfile",fnPickingParameters)
    print os.path.curdir
    mD=xmipp.MetaData();
    xmipp.readMetaDataWithTwoPossibleImages(MicrographSelfile, mD)
    preprocessingDir,dummy=os.path.split(MicrographSelfile)
    message="Cannot find the following micrographs:\n"
    NnotFound=0
    for id in mD:
        micrograph=mD.getValue(xmipp.MDL_IMAGE)
        if not os.path.exists(preprocessingDir+"/"+micrograph):
            message+=preprocessingDir+"/"+micrograph+"\n"
            NnotFound=NnotFound+1
        if isPairTilt:
            micrographTilted=mD.getValue(xmipp.MDL_ASSOCIATED_IMAGE1)
            if not os.path.exists(preprocessingDir+"/"+micrographTilted):
                message+=preprocessingDir+"/"+micrographTilted+"\n"
                NnotFound=NnotFound+1
    if NnotFound>0:
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
                 DoRemoveDust,
                 DustRemovalThreshold,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour)
