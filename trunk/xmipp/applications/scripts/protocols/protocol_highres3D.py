#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for High Resolution 3D reconstruction
#
#  - delete and create workdirectory
#  - alignment (MLalign2D)
#
# Example use:
# ./xmipp_protocol_highres3D.py
#
# Author: Carlos Oscar Sanchez Sorzano, April 2007
#
# required files: log.py
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------
# Selfile with the input images:
SelFileName='all_except_failures.sel'

# Reference file name (3D map)
ReferenceFileName="init_reference/LTA_rot_0.1_norm.vol"

# Working directory: 
WorkDirectory='Experiment1'

# Delete working directory if it already exists?
DoDeleteWorkingDir=False

# {expert} Root directory name for this project:
ProjectDir="/usr/scratch/cperez/CTD_LTA"

# {expert} Directory name for logfiles:
LogDir="Logs"

#-----------------------------------------------------------------------------
# {section} Particle description
#-----------------------------------------------------------------------------
# Particle radius (pixels)
ParticleRadius=18

# Particle mass (Daltons)
""" It is better to be conservative and provide a mass that is larger
    than the actual mass
"""
ParticleMass=800000

# Symmetry file
""" Use the syntax from xmipp_symmetry
"""
SymmetryFile="symmetryC6.txt"

# Sampling rate (Angstrom/pixel)
SamplingRate=5.6

#-----------------------------------------------------------------------------
# {section} Iteration parameters
#-----------------------------------------------------------------------------
# Pyramid levels
""" Scale=0 represents the plain images as they are. If Scale=1, then images
    are downsampled by 2 once. If Scale=2, they are downsampled twice, and so
    on. You must specify a pyramid level for each iteration. This can be done
    by a sequence of numbers (for instance, "2 2 1 1 1 0 0 0" specifies
    8 iterations, the first two with scale=2, then three with scale=1, and
    three with scale=0) or by a compact notation ("20x2 15x1 10x0", i.e.,
    10 iterations at value scale=2, 15 at scale=1, 10 at scale=0).
    
    Downsampled images should not be smaller than 64x64 pixels.
"""
PyramidLevels="35x1 15x0"

# Angular steps
""" Angular steps for each of the iterations. This parameter is used to build
    the discrete angular assignment library in the wavelet assignment. The
    smaller the angular step, the most accurate (but slow) it is. The
    syntax for this vector is the same as for the PyramidLevels
"""
AngularSteps="40x5 10x2"

# {expert} Use ART for reconstruction instead of WBP
UseART="0"

# {expert} Serial ART
""" Do serial ART even if parallel execution is available
"""
SerialART="1"

# {expert} ART lambda
""" ART relaxation factor for each of the iterations. There is a tradeoff
    between noisy reconstructions and speed of convergence. Too low relaxation
    factors result in too smooth reconstructions. On the other hand,
    too high relaxation factors result in too noisy reconstructions. As a
    a general rule, the noiser the projections, the smaller the relaxation
    factor.
"""
ARTLambda="20x0.001 10x0.001 15x0.0005"

# {expert} Discrete angular assignment
""" Especify for each iteration whether there is discrete assignment
    or not. Sometimes, it is desirable to do only continuous assignments
    in the last iterations. Set this variable to 0
    if no discrete assignment should be done, or to 1 if it should be done
"""
DiscreteAssignment="1"

# {expert} Continuous angular assignment
""" Especify for each iteration whether there is continuous assignment
    or not. It is not advisable to use continuous assignment if the
    reference volume is still too blurred. Set this variable to 0
    if no continuous assignment should be done, or to 1 if it should be done
"""
ContinuousAssignment="7x0 1"

# {expert} Compute resolution
""" Resolution is something slow to compute, do not abuse of it.
    The resolution is always computed in the last iteration
"""
ComputeResolution=""

# {expert} Resume at iteration
""" Provide the first iteration number for which a valid volume is not computed.
    If no iterations have been performed, set to 1.
"""
ResumeIteration=1

#-----------------------------------------------------------------------------
# {section} CTF Amplitude Correction (assuming previous phase correction)
#-----------------------------------------------------------------------------
# {expert} CTF dat file
""" This file specify the CTF parameters for each image
"""
CTFDat=""

# {expert} Amplitude correction
""" Specify whether amplitude correction is performed or not at each
    iteration.
"""
AmplitudeCorrection="0"

#-----------------------------------------------------------------------------
# {section} Post-processing
#-----------------------------------------------------------------------------
# {expert} Mask references
""" Masking the reference from one iteration to the next is important in
    order to reduce artefacts. The initial mask must be provided as a binary
    volume (0=background, 1=protein) of the same size as the reference volume.
 
    Masks for subsequent scales are generated by expansion, thresholding,
    opening, and closing of the initial mask.
"""
DoReferenceMask="1"

# Initial Reference Mask Volume
InitialReferenceMask="init_reference/initialMask.vol"

# {expert} Reference Lowpass filter
""" This vector specifies the frequency at which each reference volume
    will be filtered. The maximum frequency is 0.5, and 0.25 for all
    iterations is a reasonable value.
"""
FilterReference="0.25"

# {expert} Segment reference using particle mass
""" This vector specifies which iterations will use the particle mass
    to segment the reference for the next iteration. It is not
    adivsable to segment using the mass while the volume is still too
    blurred.
"""
SegmentUsingMass="0"

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel? (see Expert options)
DoParallel=True

# Number of processors to use:
MyNumberOfCPUs=8

# {expert} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
"""
MyMachineFile=""

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This variable serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_highres3D.py"

#-----------------------------------------------------------------------------
# {end-of-header} do not change anything bellow this line unless you know what you are doing
#-----------------------------------------------------------------------------
#===========================================================================
# Beginning of the protocol
#===========================================================================
import glob
import math
import os
import shutil
import string
import sys

scriptdir=os.path.expanduser('~')+'/scripts'
sys.path.append(scriptdir) # add default search path
import log
import SelFiles
import launch_parallel_job

class HighRes3DClass:
   #------------------------------------------------------------------------
   # Class constructor
   #------------------------------------------------------------------------
   def __init__(self,
      	        _SelFileName,
		_ReferenceFileName,
		_WorkDirectory,
		_DoDeleteWorkingDir,
		_ProjectDir,
		_LogDir,
		
		_ParticleRadius,
		_ParticleMass,
		_SymmetryFile,
		_SamplingRate,
		
		_PyramidLevels,
		_AngularSteps,
		_UseART,
		_SerialART,
		_ARTLambda,
		_DiscreteAssignment,
		_ContinuousAssignment,
		_ComputeResolution,
		_ResumeIteration,
		
		_CTFDat,
		_AmplitudeCorrection,
		
		_DoReferenceMask,
		_InitialReferenceMask,
		_FilterReference,
		_SegmentUsingMass,

		_DoParallel,
		_MyNumberOfCPUs,
		_MyMachineFile,
		
		_Verbose
                ):
       self.scriptFile=os.path.abspath(sys.argv[0])
       self.selFileName=os.path.abspath(_ProjectDir+"/"+_SelFileName)
       self.referenceFileName=os.path.abspath(_ProjectDir+"/"+_ReferenceFileName)
       self.workDirectory=os.getcwd()+'/'+_WorkDirectory
       self.doDeleteWorkingDir=_DoDeleteWorkingDir
       self.ProjectDir=_ProjectDir
       self.logDir=_LogDir
		
       self.particleRadius=_ParticleRadius
       self.particleMass=_ParticleMass
       if _SymmetryFile=="":
          self.symmetryFile=self.workDirectory+"/Src/symmetry.sym"
       else:
          self.symmetryFile=self.ProjectDir+"/"+_SymmetryFile
       self.samplingRate=_SamplingRate

       self.pyramidLevels="0 "+_PyramidLevels
       self.angularSteps="0 "+_AngularSteps
       self.useART="0 "+_UseART
       self.serialART="0 "+_SerialART
       self.ARTLambda="0 "+_ARTLambda
       self.discreteAssignment="0 "+_DiscreteAssignment
       self.continuousAssignment="0 "+_ContinuousAssignment
       self.computeResolution="0 "+_ComputeResolution
       self.resumeIteration=_ResumeIteration

       self.CTFDat=_CTFDat
       self.amplitudeCorrection="0 "+_AmplitudeCorrection
		
       self.doReferenceMask="0 "+_DoReferenceMask
       if not _InitialReferenceMask=="":
          self.initialReferenceMask=os.path.abspath(_ProjectDir+"/"+_InitialReferenceMask)
       else:
          self.initialReferenceMask=""
       self.filterReference="0 "+_FilterReference
       self.segmentUsingMass="0 "+_SegmentUsingMass

       self.doParallel=_DoParallel
       self.myNumberOfCPUs=_MyNumberOfCPUs
       self.myMachineFile=_MyMachineFile

       self.verbose=_Verbose

       if self.verbose:
	  self.mylog=log.init_log_system(_ProjectDir,
                                	 _LogDir,
                                	 sys.argv[0],
                                	 _WorkDirectory)
	  self.logLevel='info'
       
       # Produce side info
       self.SF=SelFiles.selfile()
       self.SF.read(self.selFileName)
                                      
   #------------------------------------------------------------------------
   # AdaptScaleFromPreviousIteration
   #------------------------------------------------------------------------
   def adaptScaleFromPreviousIteration(self,_iteration):
      self.log("# Adapting alignment and model from previous iteration ----------------")
      
      # Get the previous and current pyramid levels
      previousPyramidLevel=int(self.getPyramidLevel(_iteration-1))
      currentPyramidLevel=int(self.getPyramidLevel(_iteration))
      
      if previousPyramidLevel==currentPyramidLevel:
         # If they are the same, link the results from the previous
	 # _iteration in the current _iteration
	 self.linkFile("../"+self.getModelFilename(_iteration-1),
	               self.getModelFFilename(_iteration))
	 self.linkFile("../"+self.getAlignmentFilename(_iteration-1),
                       self.getAlignmentFFilename(_iteration))
      else:
         # If they are different
	 # Copy the last model and expand/reduce it
	 self.copyFile(self.getModelFilename(_iteration-1),
	               self.getModelFFilename(_iteration))
	 if previousPyramidLevel<currentPyramidLevel:
            self.execute("xmipp_pyramid -i "+\
	                 self.getModelFFilename(_iteration)+\
	        	 " -reduce -levels "+\
			 str(currentPyramidLevel-previousPyramidLevel))
         else:
            self.execute("xmipp_pyramid -i "+\
	                 self.getModelFFilename(_iteration)+\
	        	 " -expand -levels "+\
			 str(previousPyramidLevel-currentPyramidLevel))
	 
	 # Process the angle file and expand the shifts
	 if not previousPyramidLevel==currentPyramidLevel:
	    factor=math.pow(2.0,previousPyramidLevel-currentPyramidLevel)
	    self.execute(
               "export LC_ALL=C; export LANG=C; awk '$1==\";\" {print $0} "+\
      	       "$1!=\";\" {printf(\"%d %d %f %f %f %f %f %f\\n\",$1,6,$3,$4,$5,"+
	       str(factor)+"*$6,"+str(factor)+"*$7,$8);}' "+\
               " < "+self.getAlignmentFilename(_iteration-1)+
	       " | grep -v \"^0 0\" > "+\
	       self.getAlignmentFFilename(_iteration))
	 else:
	    self.linkFile(self.getAlignmentFilename(_iteration-1),
	                  self.getAlignmentFFilename(_iteration))
	 
	 # Check that the corresponding mask exists
	 if not self.initialReferenceMask=="":
	    if not os.path.exists(self.getMaskFilename(_iteration)):
	       if not currentPyramidLevel==0:
                  self.execute("xmipp_pyramid -i ../Src/referenceScaledMask.vol -o "+\
	                       self.getMaskFilename(_iteration)+\
	        	       " -reduce -levels "+\
			       str(currentPyramidLevel))
	       else:
	          self.linkFile("referenceScaledMask.vol",self.getMaskFilename(_iteration))
	       if not os.path.exists(self.getMaskFilename(_iteration)):
	          raise RuntimeError,"Cannot create "+self.getMaskFilename(_iteration)
	       self.execute("xmipp_threshold -i "+
	                    self.getMaskFilename(_iteration)+" -binarize -below 0.5")
	       self.execute("xmipp_morphology -i "+
	                    self.getMaskFilename(_iteration)+" -ope")
	       self.execute("xmipp_morphology -i "+
	                    self.getMaskFilename(_iteration)+" -clo")

   #------------------------------------------------------------------------
   # Angular assignment
   #------------------------------------------------------------------------
   def angularAssignment(self,_iteration):
       self.log("# Angular assignment --------------------------------------------------")
       
       # Take the last assignment to the image headers
       self.execute("xmipp_headerinfo -assign -i "+\
                    self.getAlignmentFFilename(_iteration)+\
		    " -o preproc_assign.sel -force")
       
       # Perform a discrete angular assignment
       if self.getDiscreteAssignment(_iteration)=="1":
	  params="-i preproc_assign.sel "+\
        	 "-ref "+self.getModelFFilename(_iteration)+" "+\
		 "-oang "+self.getDiscreteAnglesFilename(_iteration)+" "+\
		 "-proj_step "+self.getAngularSteps(_iteration)+" "+\
		 "-psi_step "+self.getAngularSteps(_iteration)+" "+\
                 "-sym "+self.symmetryFile
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_angular_predict",
        			       "xmipp_mpi_angular_predict",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
       else:
          self.linkFile(removeDirectories(self.getAlignmentFFilename(_iteration)),
	                self.getDiscreteAnglesFilename(_iteration))
       if not os.path.exists(self.getDiscreteAnglesFilename(_iteration)):
          raise RuntimeError,"There is a problem with the discrete assignment"

       # Perform a continuous angular assignment
       if self.getContinuousAssignment(_iteration)=="1":
	  params="-i preproc_assign.sel "+\
        	 "-ref "+self.getModelFFilename(_iteration)+" "+\
		 "-ang "+self.getDiscreteAnglesFilename(_iteration)+" "+\
		 "-oang "+self.getContinuousAnglesFilename(_iteration)
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_angular_predict_continuous",
        			       "xmipp_mpi_angular_predict_continuous",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
       else:
          self.linkFile(removeDirectories(self.getDiscreteAnglesFilename(_iteration)),
	                self.getContinuousAnglesFilename(_iteration))
       if not os.path.exists(self.getContinuousAnglesFilename(_iteration)):
          raise RuntimeError,"There is a problem with the continuous assignment"

       # Measure distance between the two angular assignments
       self.execute("xmipp_angular_distance -ang1 "+
                    self.getAlignmentFFilename(_iteration)+" -ang2 "+\
      	            self.getAlignmentFilename(_iteration)+\
		    " -o Iteration"+itoa(_iteration,2)+"/diff_angles"+\
		    itoa(_iteration,2)+" -sym "+self.symmetryFile+\
                    " -check_mirrors | "+\
		    " sed -n 's/Global distance = /"+str(_iteration)+" /p' "+\
		    ">> angle_convergence.txt")

   #------------------------------------------------------------------------
   # Backup protocol
   #------------------------------------------------------------------------
   def backupProtocol(self):
       self.log("# Creating backup")
       self.copyFile(self.scriptFile,self.workDirectory+"/Src/"+\
          removeDirectories(sys.argv[0]))

   #------------------------------------------------------------------------
   # Change directory
   #------------------------------------------------------------------------
   def changeDirectory(self,_directory):
       self.log("cd "+_directory)
       os.chdir(_directory)

   #------------------------------------------------------------------------
   # Copy file
   #------------------------------------------------------------------------
   def copyFile(self,_source,_target):
       if os.path.exists(_target):
          self.execute("rm -f "+_target)
          if os.path.exists(_target):
	     raise RuntimeError,"There is an error deleting "+_target
       self.execute("cp "+_source+" "+_target)
       if not os.path.exists(_target):
	  raise RuntimeError,"There is an error copying "+_target

   #------------------------------------------------------------------------
   # Create directory
   #------------------------------------------------------------------------
   def createDirectory(self,_directory):
       if not os.path.exists(_directory):
          self.log("mkdir "+_directory)
          os.makedirs(_directory)
          if not os.path.exists(_directory):
	     raise RuntimeError,"There is an error creating "+_directory
	  self.execute("chmod -R u+rw "+_directory)

   #------------------------------------------------------------------------
   # Delete directory
   #------------------------------------------------------------------------
   def deleteDirectory(self,_directory):
       if os.path.exists(_directory):
          self.log("rm -rf "+_directory)
          shutil.rmtree(_directory)
          if os.path.exists(_directory):
	     raise RuntimeError,"There is an error deleting "+_directory

   #------------------------------------------------------------------------
   # Delete file
   #------------------------------------------------------------------------
   def deleteFile(self,_filename):
       if os.path.exists(_filename):
          self.execute("rm -f "+_filename)
          if os.path.exists(_filename):
	     raise RuntimeError,"There is an error deleting "+_filename

   #------------------------------------------------------------------------
   # Execute command
   #------------------------------------------------------------------------
   def execute(self,_command):
      self.log(_command)
      os.system(_command)

   #------------------------------------------------------------------------
   # Generate images whose size is a power of 2
   #------------------------------------------------------------------------
   def generateImagesPower2(self):
       self.log("Generating a set of images whose size is a power of 2",
                'info',True)

       # Compute the apropriate power of 2 size
       self.getImageSize()
       
       # Compute the apropriate particle radius
       self.particleWorkingRadius=math.ceil(
          self.particleRadius*self.workXDim/self.xDim)

       # Generate the images if not already generated
       if not os.path.exists(self.workDirectory+"/ScaledImgs"):
          self.createDirectory(self.workDirectory+"/ScaledImgs")

          # Copy the images as they are and change name
      	  self.changeDirectory(self.ProjectDir)
      	  self.execute("xmipp_adapt_for_spider rename -i "+self.selFileName+ \
                       " -oroot "+self.workDirectory+"/ScaledImgs/img"+\
	               " -o "+self.workDirectory+"/imgs.sel")

          # Rescale the images
	  self.changeDirectory(self.workDirectory)
	  if not self.xDim==self.workXDim:
	     self.execute("xmipp_scale -i imgs.sel -xdim "+str(self.workXDim))
	  
	  # Rescale the reference volume and mask
	  self.copyFile(self.referenceFileName,"Src/referenceScaledVolume.vol")
          if not self.initialReferenceMask=="":
	     self.copyFile(self.initialReferenceMask,"Src/referenceScaledMask.vol")
	  if not self.xDim==self.workXDim:
	     self.execute("xmipp_scale -i Src/referenceScaledVolume.vol -xdim "+str(self.workXDim))
             if not self.initialReferenceMask=="":
	        self.execute("xmipp_scale -i Src/referenceScaledMask.vol -xdim "+str(self.workXDim))

   #------------------------------------------------------------------------
   # Generate images for this iteration
   #------------------------------------------------------------------------
   def generateImagesForThisIteration(self,_iteration):
       self.log("# Generating images for this iteration --------------------------------")

       if not os.path.exists("preproc_recons.sel") \
          or _iteration==self.resumeIteration \
          or not self.getPyramidLevel(_iteration)==self.getPyramidLevel(_iteration-1) \
          or self.getAmplitudeCorrection(_iteration)=="1" \
	  or True:
          
          # Generate plain set of images
          self.deleteDirectory("preproc")
          self.createDirectory("preproc")
          self.execute("xmipp_adapt_for_spider rename -i ../imgs.sel -oroot preproc/preproc -o preproc.sel")
	  self.execute("chmod -R u+w preproc*")

          # Correct for the amplitude
          if self.getAmplitudeCorrection(_iteration)=="1" \
	     and self.getPyramidLevel(_iteration)=="0":
	     self.execute("xmipp_headerinfo -assign -i "+\
	                  self.getAlignmentFFilename(_iteration)+" -o preproc.sel"+\
		          " -force")
	     self.execute("xmipp_idr_art -exp preproc.sel -vol "+\
	                  self.getModelFilename(_iteration)+" -ctf "+\
		          self.CTFDat+\
		          " -oroot preproc/preproc -adjust_gray_levels")
	     # COSS: *** IDR no funciona con el tipo de ficheros que estoy pasando

          # Generate images for assignment
	  factor=pow(2,-int(self.getPyramidLevel(_iteration)))
          self.execute("xmipp_adapt_for_spider rename -i preproc.sel -oroot preproc/preproc_assign -o preproc_assign.sel")
	  self.execute("chmod -R u+w preproc*")
          if not self.getPyramidLevel(_iteration)=="0":
	     self.execute("xmipp_pyramid -i preproc_assign.sel -reduce -levels "+\
	                  str(self.getPyramidLevel(_iteration)))
          self.execute("xmipp_normalize -i preproc_assign.sel -method NewXmipp -background circle "+\
	               str(math.ceil(self.particleWorkingRadius*1.1*factor)))

          # Remove useless images
          self.execute("xmipp_rmsel preproc.sel")

          # Generate images for reconstruction
          self.execute("xmipp_adapt_for_spider rename -i preproc_assign.sel -oroot preproc/preproc_recons -o preproc_recons.sel")
	  self.execute("chmod -R u+w preproc*")
          self.execute("xmipp_mask -i preproc_recons.sel -mask raised_cosine "+\
	               str(-(math.ceil(self.particleWorkingRadius*factor)))+" "+\
		       str(-(math.ceil(self.particleWorkingRadius*1.1*factor))))
	  
   #------------------------------------------------------------------------
   # Generate Next Reference
   #------------------------------------------------------------------------
   def generateNextReference(self,_iteration):
      self.log("# Generating next reference -------------------------------------------")
      self.copyFile(self.getReconstructionRootname(_iteration)+".vol",
                    self.getModelFilename(_iteration))
		  
      # Filter
      if not self.getFilterReference(_iteration)=="0":
         self.execute("xmipp_fourierfilter -i "+\
	              self.getModelFilename(_iteration)+" "+\
		      "-low_pass "+self.getFilterReference(_iteration)+" "+\
			 "-fourier_mask raised_cosine 0.05")
	 if self.getDoReferenceMask(_iteration)=="1":
            self.execute("xmipp_mask -i "+self.getModelFilename(_iteration)+" "+\
	        	 "-mask "+self.getMaskFilename(_iteration))

      # Segment
      if self.getSegmentUsingMass(_iteration)=="1":
         self.execute("xmipp_segment -i "+self.getModelFilename(_iteration)+" "+\
	              "-dalton_mass "+str(self.particleMass)+" "+\
		      "-sampling_rate "+
		      str(self.workingSamplingRate*pow(2,int(self.getPyramidLevel(_iteration))))+" "+\
		      "-o temp_mask.vol")
	 self.execute("xmipp_mask -i "+self.getModelFilename(_iteration)+" "+\
	              "-mask temp_mask.vol")
         self.deleteFile("temp_mask.vol")

      # Move the center of mass to 0
      self.execute("xmipp_findcenter3D -i "+self.getModelFilename(_iteration)+" "+\
                   "-center_volume")

   #------------------------------------------------------------------------
   # Get
   #------------------------------------------------------------------------
   def getAlignmentFilename(self,_iteration):
      return "Iteration"+itoa(_iteration,2)+"/angles"+itoa(_iteration,2)+"b.txt"
   def getAlignmentFFilename(self,_iteration):
      return "Iteration"+itoa(_iteration,2)+"/angles"+itoa(_iteration-1,2)+"F.txt"
   def getAmplitudeCorrection(self,_iteration):
      return getComponentFromVector(self.amplitudeCorrection,_iteration)
   def getAngularSteps(self,_iteration):
      return getComponentFromVector(self.angularSteps,_iteration)
   def getARTLambda(self,_iteration):
      return getComponentFromVector(self.ARTLambda,_iteration)
   def getComputeResolution(self,_iteration):
      return getComponentFromVector(self.computeResolution,_iteration)
   def getContinuousAnglesFilename(self,_iteration):
      return self.getAlignmentFilename(_iteration)
   def getContinuousAssignment(self,_iteration):
      return getComponentFromVector(self.continuousAssignment,_iteration)
   def getDiscreteAnglesFilename(self,_iteration):
      return "Iteration"+itoa(_iteration,2)+"/angles"+itoa(_iteration,2)+".txt"
   def getDiscreteAssignment(self,_iteration):
      return getComponentFromVector(self.discreteAssignment,_iteration)
   def getDoReferenceMask(self,_iteration):
      return getComponentFromVector(self.doReferenceMask,_iteration)
   def getFilterReference(self,_iteration):
      return getComponentFromVector(self.filterReference,_iteration)
   def getMaskFilename(self,_iteration):
      return "../Src/maskPyramidLevel"+str(self.getPyramidLevel(_iteration))+".vol"
   def getModelFilename(self,_iteration):
      return "Iteration"+itoa(_iteration,2)+"/model"+itoa(_iteration,2)+".vol"
   def getModelFFilename(self,_iteration):
      return "Iteration"+itoa(_iteration,2)+"/model"+itoa(_iteration-1,2)+"F.vol"
   def getReconstructionRootname(self,_iteration):
      return "Iteration"+itoa(_iteration,2)+"/volume"+itoa(_iteration,2)
   def getPyramidLevel(self,_iteration):
      return getComponentFromVector(self.pyramidLevels,_iteration)
   def getSegmentUsingMass(self,_iteration):
      return getComponentFromVector(self.segmentUsingMass,_iteration)
   def getSerialART(self,_iteration):
      return getComponentFromVector(self.serialART,_iteration)
   def getUseART(self,_iteration):
      return getComponentFromVector(self.useART,_iteration)

   #------------------------------------------------------------------------
   # Get image size
   #------------------------------------------------------------------------
   def getImageSize(self):
      SFaux=self.SF.add_1directory_begin(self.ProjectDir)
      self.xDim,self.yDim=SFaux.imgSize()
      self.workXDim=pow(2.0,math.ceil(math.log10(self.xDim)/math.log10(2.0)));
      self.workingSamplingRate=self.samplingRate*self.xDim/self.workXDim
      self.log("# Current image size="+str(self.xDim))
      self.log("# Working image size="+str(self.workXDim))
      self.log("# Working sampling rate="+\
               str(self.workingSamplingRate))

   #------------------------------------------------------------------------
   # Initialize Directories
   #------------------------------------------------------------------------
   def initDirectories(self):
       self.log("Initializing directory structure","info",True)

       # Delete working dir if necessary
       if self.doDeleteWorkingDir or not os.path.exists(self.workDirectory):
          self.deleteDirectory(self.workDirectory)
	  self.createDirectory(self.workDirectory)
	  self.changeDirectory(self.workDirectory)
	  self.createDirectory('Src')
	  self.createDirectory('Results')
	  if not os.path.exists(self.symmetryFile):
             self.touchFile(self.symmetryFile)

       # Backup the protocol
       self.backupProtocol()

   #------------------------------------------------------------------------
   # Link file
   #------------------------------------------------------------------------
   def linkFile(self,_source,_target):
       if os.path.exists(_target):
          self.execute("rm -f "+_target)
          if os.path.exists(_target):
	     raise RuntimeError,"There is an error deleting "+_target
       self.execute("ln -s "+_source+" "+_target)
       if not os.path.exists(_target):
	  raise RuntimeError,"There is an error linking "+_target

   #------------------------------------------------------------------------
   # Log messages
   #------------------------------------------------------------------------
   def log(self, _message, _level='info', _frame=False):
       if not self.verbose: return
       if _level=='debug' and self.logLevel=='info':
          return
       if _frame:
          print '# '+'*'*70
	  self.mylog.info('# '+'*'*148)
       if len(_message)<150:
          _message+=' '*(150-len(_message))
       print _message
       if _level=='info':
          self.mylog.info(_message)
       elif _level=='debug':
          self.mylog.debug(_message)
       if _frame:
          print '# '+'*'*70
	  self.mylog.info('# '+'*'*148)

   #------------------------------------------------------------------------
   # Prealignment
   #------------------------------------------------------------------------
   def prealignment(self):
       self.log("Prealigning","info",True);
       if not os.path.exists(self.workDirectory+"/Src/prealignment.txt"):
          self.changeDirectory(self.workDirectory+"/Results")
	  self.createDirectory("preproc")
	  self.execute("xmipp_adapt_for_spider rename -i ../imgs.sel "+
                       "-oroot preproc/preproc -o preproc.sel");
	  self.execute("xmipp_fourierfilter -i preproc.sel "+
                       "-low_pass 0.25 -fourier_mask raised_cosine 0.02");
	  params="-i preproc.sel -nref 1 -output_docfile -fast"
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_MLalign2D",
        			       "xmipp_mpi_MLalign2D",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
	  self.copyFile("ml2d.doc","../Src/prealignment.txt")
	  self.execute("rm -rf MLalign2D_* ml2d* preproc*")

   #------------------------------------------------------------------------
   # Prepare for first iteration
   #------------------------------------------------------------------------
   def prepareForFirstIteration(self):
       self.log("Preparing for first iteration","info",True);
       self.changeDirectory(self.workDirectory+"/Results")
       if self.resumeIteration==1:
          self.deleteDirectory("Iteration00")
          self.createDirectory("Iteration00")
	  self.linkFile("../../Src/prealignment.txt",self.getAlignmentFilename(0))
	  self.linkFile("../../Src/referenceScaledVolume.vol",self.getModelFilename(0))
	  self.deleteFile("angle_convergence.txt")
	  self.touchFile("angle_convergence.txt")
       else:
          if not os.path.exists("angle_convergence.txt"):
	     raise RuntimeError,"File angle_convergence.txt does not exist"
          self.execute("head -"+str(self.resumeIteration-1)+" angle_convergence.txt > inter")
	  self.execute("mv inter angle_convergence.txt")
       for it in range(self.resumeIteration,99):
	  self.deleteDirectory("Iteration"+itoa(it,2))
       
   #------------------------------------------------------------------------
   # Reconstruct
   #------------------------------------------------------------------------
   def reconstruct(self,_iteration):
       self.log("# 3D reconstruction ---------------------------------------------------")
       
       # Take the last assignment to the image headers
       self.execute("xmipp_headerinfo -assign -i "+\
                    self.getAlignmentFilename(_iteration)+\
		    " -o preproc_recons.sel -force")
       
       # Reconstruct
       if self.getUseART(_iteration)=="1":
          params="-i preproc_recons.sel "+\
	         "-o "+self.getReconstructionRootname(_iteration)+" "+\
		 "-l "+self.getARTLambda(_iteration)+" "+\
		 "-sym "+self.symmetryFile+" "+\
		 "-R "+str(math.ceil(self.particleWorkingRadius*1.1))
	  doParallel=self.doParallel
	  if self.getSerialART(_iteration)=="1":
	     doParallel=False
	  launch_parallel_job.launch_job(doParallel,
        			       "xmipp_art",
        			       "xmipp_mpi_art",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
	  self.deleteFile(self.getReconstructionRootname(_iteration)+".hist")
       else:
          params="-i preproc_recons.sel "+\
	         "-o "+self.getReconstructionRootname(_iteration)+".vol "+\
		 "-sym "+self.symmetryFile+" "+\
		 "-use_each_image"
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_wbp",
        			       "xmipp_mpi_wbp",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
       if not os.path.exists(self.getReconstructionRootname(_iteration)+".vol"):
	  raise RuntimeError,"There is a problem when reconstructing"
       
       # Apply raised cosine mask
       self.execute("xmipp_mask -i "+self.getReconstructionRootname(_iteration)+\
                    ".vol -mask raised_cosine "+str(-self.particleWorkingRadius)+\
		    " "+str(-math.ceil(self.particleWorkingRadius*1.1)))
       
       # Symmetrize volume
       self.execute("xmipp_symmetrize -i "+\
                    self.getReconstructionRootname(_iteration)+".vol "+\
		    "-sym "+self.symmetryFile)
      
       # Mask volume
       if not self.initialReferenceMask=="":
          self.execute("xmipp_mask -i "+\
                       self.getReconstructionRootname(_iteration)+".vol "+\
	               "-mask "+self.getMaskFilename(_iteration))

   #------------------------------------------------------------------------
   # Reconstruction iteration
   #------------------------------------------------------------------------
   def reconstructionIteration(self,_iteration):
      self.changeDirectory(self.workDirectory+"/Results")
      self.createDirectory("Iteration"+itoa(_iteration,2))
      self.adaptScaleFromPreviousIteration(_iteration)
      self.generateImagesForThisIteration(_iteration)
      self.angularAssignment(_iteration)
      self.reconstruct(_iteration)
      
   #------------------------------------------------------------------------
   # Run
   #------------------------------------------------------------------------
   def run(self):
       self.initDirectories()
       self.generateImagesPower2()
       self.prealignment()
       self.prepareForFirstIteration()
       for it in range(len(getListFromVector(self.pyramidLevels))):
           if it>=self.resumeIteration:
              self.log("Iteration "+str(it),'info',True)
	      self.reconstructionIteration(it)
	      self.generateNextReference(it)
       
   #------------------------------------------------------------------------
   # Touch file
   #------------------------------------------------------------------------
   def touchFile(self,_filename):
       self.execute("touch "+_filename)
       if not os.path.exists(_filename):
	  raise RuntimeError,"There is an error touching "+_filename

#===========================================================================
# Useful functions
#===========================================================================
#---------------------------------------------------------------------------
# getComponentFromVector
#---------------------------------------------------------------------------
def getComponentFromVector(_vector,_iteration):
   listValues=getListFromVector(_vector)
   if _iteration<0: _iteration=0
   if _iteration<len(listValues): return listValues[_iteration]
   else:                          return listValues[len(listValues)-1]

#---------------------------------------------------------------------------
# getListFromVector
#---------------------------------------------------------------------------
def getListFromVector(_vector):
   intervalos=string.split(_vector)
   if len(intervalos)==0:
      raise RuntimeError,"Empty vector"
   listValues=[]
   for i in range(len(intervalos)):
      intervalo=intervalos[i]
      listaIntervalo=string.split(intervalo,'x')
      if len(listaIntervalo)==1:
	 listValues+=listaIntervalo
      elif len(listaIntervalo)==2:
         listValues+=[ listaIntervalo[1] ] * string.atoi(listaIntervalo[0])
      else:
         raise RuntimeError,"Unknown syntax: "+intervalos
   return listValues

#---------------------------------------------------------------------------
# itoa
#---------------------------------------------------------------------------
def itoa(_number,_length):
   return ("%0"+str(_length)+"d") % _number

#---------------------------------------------------------------------------
# Remove directories from filename
#---------------------------------------------------------------------------
def removeDirectories(_filename):
    lista=string.split(_filename,'/')
    return lista[len(lista)-1]

#===========================================================================
# Main
#===========================================================================
if __name__ == '__main__':
   myHighRes3D=HighRes3DClass(
      	        SelFileName,
		ReferenceFileName,
		WorkDirectory,
		DoDeleteWorkingDir,
		ProjectDir,
		LogDir,
		
		ParticleRadius,
		ParticleMass,
		SymmetryFile,
		SamplingRate,
		
		PyramidLevels,
		AngularSteps,
		UseART,
		SerialART,
		ARTLambda,
		DiscreteAssignment,
		ContinuousAssignment,
		ComputeResolution,
		ResumeIteration,
		
		CTFDat,
		AmplitudeCorrection,
		
		DoReferenceMask,
		InitialReferenceMask,
		FilterReference,
		SegmentUsingMass,

		DoParallel,
		MyNumberOfCPUs,
		MyMachineFile,
		
		True
              )
   myHighRes3D.run()
