#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Xmipp protocol for multi-resolution refinement
#
# {please cite} C.O.S. Sorzano et al. Journal of Structural Biology 148: 194-204 (2004)
# {please cite} S. Jonic et al. Ultramicroscopy,  103:  303-317 (2005)
#
# Example use:
# ./xmipp_protocol_multires.py
#
# Author: Carlos Oscar Sanchez Sorzano, June 2007
#
#-----------------------------------------------------------------------------
# {section} Global parameters
#-----------------------------------------------------------------------------

# {file} Selfile with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
SelFileName='img.sel'

# {file} Initial 3D reference map:
ReferenceFileName='Src/initial_reference.vol'

# {dir} Working subdirectory: 
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkDirectory='MultiRes/Test1'

# Delete working directory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir=False

# Number of iterations to perform
NumberofIterations=10

# Resume at iteration
""" This option may be used to finish a previously performed run.
    Set to 1 to start a new run 
    Note: Do NOT delete working directory if this option is not set to 1
"""
ResumeIteration=2

# {expert} {dir} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/media/usbdisk/Experiments/TestMuySencillo'

# {expert} {dir} Directory name for logfiles:
LogDir='Logs'

# {expert} Skip prealignment:
SkipPrealignment=True

#-----------------------------------------------------------------------------
# {section} Particle description
#-----------------------------------------------------------------------------
# Particle radius (pixels)
ParticleRadius=32

# Particle mass (Daltons)
ParticleMass=''

# Symmetry group
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the symmetry groups format
    If no symmetry is present, give c1
"""
SymmetryGroup='c1'

# Sampling rate (Angstrom/pixel)
SamplingRate=1

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
    
    More than scale=2 should not be used.
    
    Scaling is done via spline pyramids, please visit:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Pyramid
"""
PyramidLevels='0'

# Angular steps
""" Angular steps for each of the iterations. This parameter is used to build
    the discrete angular assignment library in the wavelet assignment. The
    smaller the angular step, the most accurate (but slow) it is. The
    syntax for this vector is the same as for the PyramidLevels.
    
    The discrete angular assignment is done with xmipp_angular_predict:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Angular_predict
"""
AngularSteps='5'

# {expert} Quality percentil
""" The quality percentil indicates what is the percentil of images to be
    discarded according to several criteria. You can use different quality
    criteria along iterations. For instance, 2 5 8x10 means that in the
    first iteration 2% of the images will be discarded, then 5%, and 10%
    in the following 8 iterations. Use 0% for not removing any image.
"""
QualityPercentil='10'

# {expert} Quality angle movement
""" If an image moves from one iteration to the next more than this
    threshold (expressed in degrees), it will not be considered for this
    reconstruction. However, the image is not removed from the dataset
    and it might be reused in a latter iteration.
    You may use different thresholds for the different iterations.
    If you don't want to use this feature, set it to 360.
"""
QualityAngleMovement='360'

# {expert} Quality shift movement
""" If an image moves from one iteration to the next more than this
    threshold (expressed as a percentage of the image size),
    it will not be considered for this reconstruction. However, the image
    is not removed from the dataset and it might be reused in a latter
    iteration.
    You may use different thresholds for the different iterations.
    If you don't want to use this feature, set it to 100.
"""
QualityShiftMovement='100'

# {expert} Reconstruction method
""" Choose between fourier, wbp or art
    You must specify this option for each iteration. 
    This can be done by a sequence of numbers (for instance, "wbp wbp wbp art " 
    specifies 4 iterations, the first three set the value to wbp (no restriction)
    and the last  to art. An alternative compact notation 
    is ("3xwbp 1xart", i.e.,
    3 iterations with wbp, and 1 with art).
    Note: if there are less values than iterations the last value is reused
    Note: if there are more values than iterations the extra value are ignored
"""
ReconstructionMethod='fourier'

# {expert} Serial ART
""" Do serial ART even if parallel execution is available. This parameter
    is a vector specifying whether serial ART is used or not.
"""
SerialART=''

# {expert} ART lambda
""" ART relaxation factor for each of the iterations. There is a tradeoff
    between noisy reconstructions and speed of convergence. Too low relaxation
    factors result in too smooth reconstructions. On the other hand,
    too high relaxation factors result in too noisy reconstructions. As a
    a general rule, the noiser the projections, the smaller the relaxation
    factor.
    
    For more information about art, visit:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Art
"""
ARTLambda=''

# {expert} Discrete angular assignment
""" Especify for each iteration whether there is discrete assignment
    or not. Sometimes, it is desirable to do only continuous assignments
    in the last iterations. Set this variable to 0
    if no discrete assignment should be done, or to 1 if it should be done

    The discrete angular assignment is done with xmipp_angular_predict:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Angular_predict
"""
DiscreteAssignment='1'

# {expert} Continuous angular assignment
""" Especify for each iteration whether there is continuous assignment
    or not. It is not advisable to use continuous assignment if the
    reference volume is still too blurred. Set this variable to 0
    if no continuous assignment should be done, or to 1 if it should be done

    The discrete angular assignment is done with xmipp_angular_predict:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Angular_predict_continuous
"""
ContinuousAssignment='1'

# {expert} Compute resolution using FSC
""" The resolution is always computed in the last iteration
"""
DoComputeFSC='1'

# {expert} Compute resolution using SSNR
""" Computation of the spectral signal-to-noise ratio is slow, do not abuse it.
    The resolution is always computed in the last iteration
"""
DoComputeSSNR='0'

#-----------------------------------------------------------------------------
# {section} CTF Amplitude Correction (assuming previous phase correction)
#-----------------------------------------------------------------------------
# {expert} {file} CTF dat file
""" This file specify the CTF parameters for each image.
    The file is a two-column text file: the first column corresponds to
    the filename of a projection, and the second column is the
    corresponding CTF file.
    
    This is the format output by the preprocessing protocol so that
    if you have run this other protocol, you simply have to provide
    the name generated at that stage.
"""
CTFDat=''

# {expert} Amplitude correction
""" Specify whether amplitude correction is performed or not at each
    iteration.
    E.g. 45x0 5x1
"""
AmplitudeCorrection=''

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
DoReferenceMask='1'

# {file} Initial Reference Mask Volume
InitialReferenceMask='Src/initial_mask.vol'

# {expert} Reference Lowpass filter (Normalized digital freq.)
""" This vector specifies the frequency at which each reference volume
    will be filtered. The maximum frequency is 0.5, and 0.25 for all
    iterations is a reasonable value.

    For more information about Fourier filtering, please visit:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/FourierFilter
"""
FilterLowPassReference='50x0.25'

# {expert} Reference Highpass filter (Normalized digital freq.)
""" This vector specifies the frequency at which each reference volume
    will be filtered. The maximum frequency is 0.5, and 0.05 for all
    iterations is a reasonable value.

    For more information about Fourier filtering, please visit:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/FourierFilter
"""
FilterHighPassReference='0'

# {expert} Segment reference using particle mass
""" This vector specifies which iterations will use the particle mass
    to segment the reference for the next iteration. It is not
    advisable to segment using the mass while the volume is still too
    blurred.
    
    Segmentation is done with xmipp_segment:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Segment
"""
SegmentUsingMass='50x0'

# {expert} Recenter reference using center of mass
""" Specify whether the reference for the next iteration must be recentered
    or not after each iteration.
"""
Recenter=False

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel?
DoParallel=True

# Number of processors to use:
MyNumberOfCPUs=3

# Number of processors to use by large memory demanding algorithms:
""" In fact only the Fourier reconstruction method is so large memory demanding
"""
MyNumberOfCPUsReduced=2

# Number of threads available on a single node
""" Maximum number of threads that can be launched in a single node
"""
MyNumberOfThreads=2

# {file} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
    if your queueing system using an environment variable, give it here (with the leading $, e.g. $PBS_NODEFILE
"""
MyMachineFile=''

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_multires.py'
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

scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
sys.path.append(scriptdir) # add default search path
import ctfdat
import launch_parallel_job
import log
import selfile

debugSteps = False

class MultiResClass:
   #------------------------------------------------------------------------
   # Class constructor
   #------------------------------------------------------------------------
   def __init__(self,
      	        _SelFileName,
		_ReferenceFileName,
		_WorkDirectory,
		_DoDeleteWorkingDir,
		_NumberofIterations,
		_ProjectDir,
		_LogDir,
                _SkipPrealignment,
		
		_ParticleRadius,
		_ParticleMass,
		_SymmetryGroup,
		_SamplingRate,
		
		_PyramidLevels,
		_AngularSteps,
                _QualityPercentil,
                _QualityAngleMovement,
                _QualityShiftMovement,
		_ReconstructionMethod,
		_SerialART,
		_ARTLambda,
		_DiscreteAssignment,
		_ContinuousAssignment,
		_DoComputeFSC,
		_DoComputeSSNR,
		_ResumeIteration,
		
		_CTFDat,
		_AmplitudeCorrection,
		
		_DoReferenceMask,
		_InitialReferenceMask,
		_FilterLowPassReference,
		_FilterHighPassReference,
		_SegmentUsingMass,
		_Recenter,

		_DoParallel,
		_MyNumberOfCPUs,
		_MyNumberOfCPUsReduced,
                _MyNumberOfThreads,
		_MyMachineFile,
		
		_Verbose
                ):
       self.projectDir=_ProjectDir
       self.scriptFile=os.path.abspath(sys.argv[0])
       self.selFileName=os.path.abspath(_SelFileName)
       self.referenceFileName=os.path.abspath(_ReferenceFileName)
       self.workDirectory=os.path.abspath(_WorkDirectory)
       self.doDeleteWorkingDir=_DoDeleteWorkingDir
       self.numberOfIterations=_NumberofIterations
       self.logDir=os.path.abspath(_LogDir)
       self.skipPrealignment=_SkipPrealignment
		
       self.particleRadius=_ParticleRadius
       self.particleMass=_ParticleMass
       if _SymmetryGroup=="":
          self.symmetryGroup="c1"
       else:
          self.symmetryGroup=_SymmetryGroup
       self.samplingRate=_SamplingRate

       self.pyramidLevels="0 "+_PyramidLevels
       self.angularSteps="0 "+_AngularSteps
       self.qualityPercentil="0 "+_QualityPercentil
       self.qualityAngleMovement="0 "+_QualityAngleMovement
       self.qualityShiftMovement="0 "+_QualityShiftMovement
       self.reconstructionMethod="null "+_ReconstructionMethod
       self.serialART="0 "+_SerialART
       self.ARTLambda="0 "+_ARTLambda
       self.discreteAssignment="0 "+_DiscreteAssignment
       self.continuousAssignment="0 "+_ContinuousAssignment
       self.doComputeFSC="0 "+_DoComputeFSC
       self.doComputeSSNR="0 "+_DoComputeSSNR
       self.resumeIteration=_ResumeIteration

       if _CTFDat!="":
          self.CTFDat=os.path.abspath(_CTFDat)
       else:
          self.CTFDat=""
       self.amplitudeCorrection="0 "+_AmplitudeCorrection
		
       self.doReferenceMask="0 "+_DoReferenceMask
       if not _InitialReferenceMask=="":
          self.initialReferenceMask=os.path.abspath(_InitialReferenceMask)
       else:
          self.initialReferenceMask=""
	  self.doReferenceMask="0"
       self.filterLowPassReference="0 "+_FilterLowPassReference
       self.filterHighPassReference="0 "+_FilterHighPassReference
       self.segmentUsingMass="0 "+_SegmentUsingMass
       self.recenter=_Recenter

       self.doParallel=_DoParallel
       self.myNumberOfCPUs=_MyNumberOfCPUs
       self.myNumberOfCPUsReduced=_MyNumberOfCPUsReduced
       self.myNumberOfThreads=_MyNumberOfThreads
       if _MyMachineFile=='':
          self.myMachineFile=''
       elif _MyMachineFile[0]=='/':
          self.myMachineFile=os.path.abspath(_MyMachineFile)
       elif _MyMachineFile[0]=='$':
          self.myMachineFile=_MyMachineFile
       else:
          self.myMachineFile=os.path.abspath(self.projectDir+"/"+_MyMachineFile)

       self.verbose=_Verbose
       if self.verbose:
	  self.mylog=log.init_log_system(_ProjectDir,
                                	 _LogDir,
                                	 sys.argv[0],
                                	 _WorkDirectory)
	  self.logLevel='info'
       
       # Produce side info
       self.SF=selfile.selfile()
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
            self.execute("xmipp_scale_pyramid -i "+\
	                 self.getModelFFilename(_iteration)+\
	        	 " -reduce -levels "+\
			 str(currentPyramidLevel-previousPyramidLevel))
         else:
            self.execute("xmipp_scale_pyramid -i "+\
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
               self.execute("xmipp_scale_pyramid -i ../Src/referenceScaledMask.vol -o "+\
        		    self.getMaskFilename(_iteration)+\
        		    " -reduce -levels "+\
	     		    str(currentPyramidLevel))
	       self.execute("xmipp_threshold -i "+
        		    self.getMaskFilename(_iteration)+" -binarize -below 0.5")
	       self.execute("xmipp_morphology -i "+
        		    self.getMaskFilename(_iteration)+" -ope")
	       self.execute("xmipp_morphology -i "+
        		    self.getMaskFilename(_iteration)+" -clo")
            else:
               self.linkFile("referenceScaledMask.vol",self.getMaskFilename(_iteration))
            if not os.path.exists(self.getMaskFilename(_iteration)):
               raise RuntimeError,"Cannot create "+self.getMaskFilename(_iteration)

   #------------------------------------------------------------------------
   # Angular assignment
   #------------------------------------------------------------------------
   def angularAssignment(self,_iteration):
       self.log("# Angular assignment --------------------------------------------------")
       
       # Perform a discrete angular assignment
       if self.getDiscreteAssignment(_iteration)=="1":
          params0="-i "+self.getModelFFilename(_iteration)+" "+\
                   " -sampling_rate "+self.getAngularSteps(_iteration)+" "+\
                   "-o ref -sym "+self.symmetryGroup
             
          params="-i "+self.getAlignmentFFilename(_iteration)+" "+\
                  "-ref ref.sel "+\
                  "-oang "+self.getDiscreteAnglesFilename(_iteration)+" "+\
                  "-psi_step "+self.getAngularSteps(_iteration)+" "+\
		  "-max_shift_change "+str(self.particleWorkingRadius/5)
          if (_iteration==1):
             params+=" -5D -shift_step 2"
          if not self.symmetryGroup=="c1":
             params+=" -sym "+self.symmetryGroup

          launch_parallel_job.launch_job(self.doParallel,
                                         "xmipp_angular_project_library",
                                         "xmipp_mpi_angular_project_library",
                                         params0,
                                         self.mylog,
                                         self.myNumberOfCPUs,
                                         self.myMachineFile,
                                         False)
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_angular_discrete_assign",
        			       "xmipp_mpi_angular_discrete_assign",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
      	  self.execute("find . -name \"ref*\" -exec rm -f {} \; &")
       else:
          self.linkFile(removeDirectories(self.getAlignmentFFilename(_iteration)),
	                self.getDiscreteAnglesFilename(_iteration))
       if not os.path.exists(self.getDiscreteAnglesFilename(_iteration)):
          raise RuntimeError,"There is a problem with the discrete assignment"

       # Perform a continuous angular assignment
       if self.getContinuousAssignment(_iteration)=="1":
	  params="-ref "+self.getModelFFilename(_iteration)+" "+\
		 "-ang "+self.getDiscreteAnglesFilename(_iteration)+" "+\
		 "-oang "+self.getContinuousAnglesFilename(_iteration)+" "+\
		 "-max_shift "+str(self.particleWorkingRadius/5)
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_angular_continuous_assign",
        			       "xmipp_mpi_angular_continuous_assign",
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
       angular_distance_command="xmipp_angular_distance -ang1 "+\
                    self.getAlignmentFFilename(_iteration)+" -ang2 "+\
      	            self.getAlignmentFilename(_iteration)+\
		    " -o Iteration"+itoa(_iteration,2)+"/diff_angles"+\
		    itoa(_iteration,2)+" -check_mirrors"
       if not self.symmetryGroup=="":
            angular_distance_command+=" -sym "+self.symmetryGroup
       self.execute(angular_distance_command+" > inter ");
       self.execute(" sed -n 's/Global angular distance = /"+str(_iteration)+" /p' "+\
		    "< inter >> angle_convergence.txt")
       self.execute(" sed -n 's/Global shift   distance = /"+str(_iteration)+" /p' "+\
		    "< inter >> shift_convergence.txt")
       self.execute("rm -f inter")
       self.execute("mv Iteration"+itoa(_iteration,2)+"/diff_angles"+\
		    itoa(_iteration,2)+"_vec_diff_hist.txt Iteration"+\
                    itoa(_iteration,2)+"/angular_change_histogram.txt")
       self.execute("rm Iteration"+itoa(_iteration,2)+"/diff_angles*")

   #------------------------------------------------------------------------
   # Change directory
   #------------------------------------------------------------------------
   def changeDirectory(self,_directory):
       self.log("cd "+_directory)
       os.chdir(_directory)

   #------------------------------------------------------------------------
   # Compute resolution
   #------------------------------------------------------------------------
   def computeResolution(self,_iteration):
       self.log("# Computing resolution ------------------------------------------------")
       
       # FSC ...............................................................
       if self.getComputeFSC(_iteration)=="1" and \
          not self.getReconstructionMethod(_iteration)=="fourier":
          # Split the data in two halves
          self.execute("xmipp_selfile_split -i preproc_recons.sel -n 2 -o preproc_recons")

          # Make the two reconstructions
          self.runReconstructionAlgorithm(_iteration,"preproc_recons_1.sel",
	     self.getReconstructionRootname(_iteration)+"_1",False)
          self.runReconstructionAlgorithm(_iteration,"preproc_recons_2.sel",
	     self.getReconstructionRootname(_iteration)+"_2",False)

          # Compute the FSC
          self.execute("xmipp_resolution_fsc "+\
             "-ref "+self.getReconstructionRootname(_iteration)+"_1.vol "+\
             "-i "+self.getReconstructionRootname(_iteration)+"_2.vol "+\
	     "-sam "+str(self.workingSamplingRate*
	                 pow(2,int(self.getPyramidLevel(_iteration)))))
          self.execute("mv "+self.getReconstructionRootname(_iteration)+\
            "_2.vol.fsc "+self.getReconstructionRootname(_iteration)+".fsc")

          # Remove unnecessary files
          self.execute("rm -f preproc_recons_?.sel "+\
             self.getReconstructionRootname(_iteration)+"_?.vol")

       # SSNR ..............................................................
       if self.getComputeSSNR(_iteration)=="1":
          # Reconstruct the noise volume
          self.copySelFile("preproc_recons.sel","preproc_noise")
          self.execute("xmipp_mask -i preproc_noise.sel -mask circular "+\
             str(self.workXDim))
          self.execute("xmipp_add_noise -i preproc_noise.sel -gaussian 1")
          self.runReconstructionAlgorithm(_iteration,"preproc_noise",
	     self.getReconstructionRootname(_iteration)+"_noise",False)

          # Compute the SSNR
          self.execute("xmipp_resolution_ssnr "+\
             "-S "+self.getReconstructionRootname(_iteration)+".vol "+\
             "-N "+self.getReconstructionRootname(_iteration)+"_noise.vol "+\
	     "-selS preproc_recons.sel -selN preproc_noise.sel "+\
	     "-sampling_rate "+str(self.workingSamplingRate*
	        pow(2,int(self.getPyramidLevel(_iteration))))+" "+\
	     "-o "+self.getReconstructionRootname(_iteration)+".ssnr")

          # Remove unnecessary files
          self.execute("rm -rf preproc_noise*")

   #------------------------------------------------------------------------
   # Copy CTFs
   #------------------------------------------------------------------------
   def copyCTFs(self):
       if self.CTFDat=="":
          return
       self.log("Adapting the CTFs",'info',True)
       self.changeDirectory(self.workDirectory)
       self.createDirectory("ScaledCTFs")
       CTFs=ctfdat.ctfdat()
       CTFs.read(self.CTFDat)
       self.scaledCTFs=ctfdat.ctfdat()
       for line in CTFs.lines:
           # Read a line from the ctfdat
           aux=line.split()
	   fnProjection=aux[0]
	   fnCTF=aux[1]
           splitFnProjection=fnProjection.split('/')
           splitFnCTF=fnCTF.split('/')
	   
	   # Copy the CTF file to ScaledCTFs
	   if not os.path.exists("ScaledCTFs/"+splitFnCTF[-1]):
	       if fnCTF[0]=='/':
        	   self.copyFile(fnCTF,"ScaledCTFs/"+splitFnCTF[-1])
	       else:
        	   self.copyFile(self.projectDir+'/'+fnCTF,
	              "ScaledCTFs/"+splitFnCTF[-1])

    	       # Change the sampling rate
	       if not self.xDim==self.workXDim:
		   self.execute("grep -v sampling_rate ScaledCTFs/"+\
	              splitFnCTF[-1]+" > inter.txt")
		   self.execute("echo sampling_rate="+\
	              str(self.samplingRate*self.xDim/self.workXDim)+\
		      " > ScaledCTFs/"+splitFnCTF[-1])
		   self.execute("cat inter.txt >> ScaledCTFs/"+splitFnCTF[-1])
		   self.execute("rm -f inter.txt")

	   # Create corresponding line in the scaled ctfdat file
	   self.scaledCTFs.append("ScaledImgs/"+splitFnProjection[-1],
	      "ScaledCTFs/"+splitFnCTF[-1])

       # Write the scaled ctfdat file
       self.log("Creating file ctfdat.txt")
       self.scaledCTFs.write("ctfdat.txt")
       
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
   # Copy selfile
   #------------------------------------------------------------------------
   def copySelFile(self,_sourceSelFile,_targetDirectory,_targetSelFile="",
       _baseDirectory="",_pathToSelFile=""):
       self.log("# Copying "+_sourceSelFile+" into "+_targetDirectory);
       SFin=selfile.selfile()
       SFin.read(_sourceSelFile)
       if not _pathToSelFile=="":
          SFin=SFin.add_1directory_begin(_pathToSelFile)
       SFout=SFin.copy_sel(_targetDirectory)
       if not _baseDirectory=="":
          SFout=SFout.replace_string(_baseDirectory+"/","")
       if not os.path.exists(_targetDirectory):
	  raise RuntimeError,"There is an error copying "+_sourceSelFile
       if _targetSelFile=="":
          _targetSelFile=_targetDirectory+".sel"
       self.log("# Creating "+_targetSelFile)
       SFout.write(_targetSelFile)
       if not os.path.exists(_targetSelFile):
	  raise RuntimeError,"There is an error creating "+_targetSelFile
       self.execute("chmod -R u+w "+_targetDirectory)

   #------------------------------------------------------------------------
   # Create directory
   #------------------------------------------------------------------------
   def createDirectory(self,_directory):
       if not os.path.exists(_directory):
          self.log("mkdir -p "+_directory)
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
      	  self.changeDirectory(self.projectDir)
	  self.copySelFile(self.selFileName,self.workDirectory+"/ScaledImgs",
	                   self.workDirectory+"/imgs.sel", self.workDirectory)

          # Rescale the images
	  self.changeDirectory(self.workDirectory)
	  if not self.xDim==self.workXDim:
	     self.execute("xmipp_scale -i imgs.sel -xdim "+str(self.workXDim))
	  
          # Normalize images
          self.execute("xmipp_normalize -i imgs.sel -method Ramp -background circle "+\
        	       str(math.ceil(self.particleWorkingRadius*1.1)))

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
       factor=pow(2,-int(self.getPyramidLevel(_iteration)))

       # Generate images for assignment
       self.copySelFile("../imgs.sel","preproc_assign","preproc_assign.sel",
          self.workDirectory+"/Results","..")

       # Correct for the amplitude
       if self.CTFDat!="" and self.getAmplitudeCorrection(_iteration)=="1" \
              and self.getPyramidLevel(_iteration)=="0":
	  CTFs=ctfdat.ctfdat()
	  CTFs.read("../ctfdat.txt")
	  preprocCTFs=CTFs.changeDirectory("preproc_assign","../ScaledCTFs")
	  self.log("# Creating a ctfdat file for preproc_assign")
          preprocCTFs.write("preproc_assign_ctfdat.txt")
          self.execute("xmipp_header_assign -i "+\
        	       self.getAlignmentFFilename(_iteration)+" -o preproc_assign.sel"+\
        	       " -force")
          self.createDirectory("preproc_assign_IDR");
          params="-vol "+\
        	 self.getModelFFilename(_iteration)+\
		 " -ctfdat preproc_assign_ctfdat.txt"+\
		 " -oroot preproc_assign_IDR/preproc_assign_IDR_";
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_ctf_correct_idr",
        			       "xmipp_mpi_ctf_correct_idr",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
          self.execute('xmipp_selfile_create "preproc_assign_IDR/*" > preproc_assign.sel')

       # Scale if necessary
       if (not self.getPyramidLevel(_iteration)=="0" and
          (_iteration==1 or _iteration>1 and
             not self.getPyramidLevel(_iteration)==
                self.getPyramidLevel(_iteration))):
          self.execute("xmipp_scale_pyramid -i preproc_assign.sel -reduce -levels "+\
        	       str(self.getPyramidLevel(_iteration)))
          self.execute("xmipp_normalize -i preproc_assign.sel -method NewXmipp -background circle "+\
        	       str(math.ceil(self.particleWorkingRadius*1.1*factor)))

   #------------------------------------------------------------------------
   # Generate Next Reference
   #------------------------------------------------------------------------
   def generateNextReference(self,_iteration):
      self.log("# Generating next reference -------------------------------------------")
      self.copyFile(self.getReconstructionRootname(_iteration)+".vol",
                    self.getModelFilename(_iteration))
		  
      # Low Pass Filter
      if not self.getFilterLowPassReference(_iteration)=="0":
         self.execute("xmipp_fourier_filter -i "+\
	              self.getModelFilename(_iteration)+" "+\
		      "-low_pass "+self.getFilterLowPassReference(_iteration)+" "+\
			 "-fourier_mask raised_cosine 0.05")
	 if self.getDoReferenceMask(_iteration)=="1":
            self.execute("xmipp_mask -i "+self.getModelFilename(_iteration)+" "+\
	        	 "-mask "+self.getMaskFilename(_iteration))

      # Segment
      if self.getSegmentUsingMass(_iteration)=="1":
         self.execute("xmipp_volume_segment -i "+self.getModelFilename(_iteration)+" "+\
	              "-dalton_mass "+str(self.particleMass*2)+" "+\
		      "-sampling_rate "+
		      str(self.workingSamplingRate*pow(2,int(self.getPyramidLevel(_iteration))))+" "+\
		      "-o temp_mask.vol")
	 self.execute("xmipp_mask -i "+self.getModelFilename(_iteration)+" "+\
	              "-mask temp_mask.vol")
         self.deleteFile("temp_mask.vol")

      # Move the center of mass to 0
      if self.recenter:
          self.execute("xmipp_find_center3d -i "+self.getModelFilename(_iteration)+" "+\
                       "-center_volume")

      # High Pass Filter
      if not self.getFilterHighPassReference(_iteration)=="0":
         self.execute("xmipp_fourier_filter -i "+\
	              self.getModelFilename(_iteration)+" "+\
		      "-high_pass "+self.getFilterHighPassReference(_iteration)+" "+\
			 "-fourier_mask raised_cosine 0.02")
	 if self.getDoReferenceMask(_iteration)=="1":
            self.execute("xmipp_mask -i "+self.getModelFilename(_iteration)+" "+\
	        	 "-mask "+self.getMaskFilename(_iteration))
      if (debugSteps):
          a=input("Volume for next iteration prepared. Press a key")

   #------------------------------------------------------------------------
   # Get
   #------------------------------------------------------------------------
   def getAlignmentFilename(self,_iteration):
      return self.getIterationDirectory(_iteration)+"/angles"+itoa(_iteration,2)+"b.txt"
   def getAlignmentFFilename(self,_iteration):
      return self.getIterationDirectory(_iteration)+"/angles"+itoa(_iteration-1,2)+"F.txt"
   def getAmplitudeCorrection(self,_iteration):
      return getComponentFromVector(self.amplitudeCorrection,_iteration)
   def getAngularSteps(self,_iteration):
      return getComponentFromVector(self.angularSteps,_iteration)
   def getARTLambda(self,_iteration):
      return getComponentFromVector(self.ARTLambda,_iteration)
   def getComputeFSC(self,_iteration):
      return getComponentFromVector(self.doComputeFSC,_iteration)
   def getComputeResolution(self,_iteration):
      if self.getComputeFSC(_iteration)=="1": return "1"
      if self.getComputeSSNR(_iteration)=="1": return "1"
      return "0"
   def getComputeSSNR(self,_iteration):
      return getComponentFromVector(self.doComputeSSNR,_iteration)
   def getContinuousAnglesFilename(self,_iteration):
      return self.getAlignmentFilename(_iteration)
   def getContinuousAssignment(self,_iteration):
      return getComponentFromVector(self.continuousAssignment,_iteration)
   def getDiscreteAnglesFilename(self,_iteration):
      return self.getIterationDirectory(_iteration)+"/angles"+itoa(_iteration,2)+".txt"
   def getDiscreteAssignment(self,_iteration):
      return getComponentFromVector(self.discreteAssignment,_iteration)
   def getDoReferenceMask(self,_iteration):
      return getComponentFromVector(self.doReferenceMask,_iteration)
   def getFilterLowPassReference(self,_iteration):
      return getComponentFromVector(self.filterLowPassReference,_iteration)
   def getFilterHighPassReference(self,_iteration):
      return getComponentFromVector(self.filterHighPassReference,_iteration)
   def getIterationDirectory(self,_iteration):
      return "Iteration"+itoa(_iteration,2)
   def getMaskFilename(self,_iteration):
      return "../Src/maskPyramidLevel"+str(self.getPyramidLevel(_iteration))+".vol"
   def getModelFilename(self,_iteration):
      return self.getIterationDirectory(_iteration)+"/model"+itoa(_iteration,2)+".vol"
   def getModelFFilename(self,_iteration):
      return self.getIterationDirectory(_iteration)+"/model"+itoa(_iteration-1,2)+"F.vol"
   def getReconstructionRootname(self,_iteration):
      return self.getIterationDirectory(_iteration)+"/volume"+itoa(_iteration,2)
   def getPyramidLevel(self,_iteration):
      return getComponentFromVector(self.pyramidLevels,_iteration)
   def getQualityAngleMovement(self,_iteration):
      return getComponentFromVector(self.qualityAngleMovement,_iteration)
   def getQualityPercentil(self,_iteration):
      return getComponentFromVector(self.qualityPercentil,_iteration)
   def getQualityShiftMovement(self,_iteration):
      return getComponentFromVector(self.qualityShiftMovement,_iteration)
   def getReconstructionAnglesFilename(self,_iteration):
      return self.getIterationDirectory(_iteration)+"/angles_reconstruction"+itoa(_iteration,2)+".txt"
   def getSegmentUsingMass(self,_iteration):
      return getComponentFromVector(self.segmentUsingMass,_iteration)
   def getSerialART(self,_iteration):
      return getComponentFromVector(self.serialART,_iteration)
   def getReconstructionMethod(self,_iteration):
      return getComponentFromVector(self.reconstructionMethod,_iteration)

   #------------------------------------------------------------------------
   # Get image size
   #------------------------------------------------------------------------
   def getImageSize(self):
      [name,state]=self.SF.find_first_active_image()
      if name[0]=='/':
         SFaux=self.SF
      else:
         SFaux=self.SF.add_1directory_begin(self.projectDir)
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

       # Backup the protocol
       self.changeDirectory(self.projectDir)
       log.make_backup_of_script_file(sys.argv[0],
                                      os.path.abspath(self.workDirectory))

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
	  self.mylog.info('# '+'*'*70)
       print _message
       if _level=='info':
          self.mylog.info(_message)
       elif _level=='debug':
          self.mylog.debug(_message)
       if _frame:
          print '# '+'*'*70
	  self.mylog.info('# '+'*'*70)

   #------------------------------------------------------------------------
   # Prealignment
   #------------------------------------------------------------------------
   def prealignment(self):
       self.log("Prealigning","info",True);
       if not os.path.exists(self.workDirectory+"/Src/prealignment.txt"):
          if not self.skipPrealignment:
             self.changeDirectory(self.workDirectory+"/Results")
	     self.copySelFile("../imgs.sel","preproc","preproc.sel",
	        self.workDirectory+"/Results","..")
	     self.execute("xmipp_average -i preproc.sel")
	     self.execute("xmipp_align2d -i preproc.sel -ref preproc.med.xmp "+\
	                  "-iter 4 -only_trans")
	     self.execute("xmipp_header_extract -i preproc.sel -o ../Src/prealignment.txt")
	     self.execute("rm -rf preproc*")
          else:
             self.changeDirectory(self.workDirectory)
#             self.execute("xmipp_header_reset -i imgs.sel")
	     self.execute("xmipp_header_extract -i imgs.sel -o Src/prealignment.txt")

   #------------------------------------------------------------------------
   # Prepare for first iteration
   #------------------------------------------------------------------------
   def prepareForFirstIteration(self):
       self.log("Preparing for first iteration","info",True);
       self.changeDirectory(self.workDirectory+"/Results")
       if self.resumeIteration==1:
          self.deleteDirectory("Iteration00")
          self.createDirectory("Iteration00")
	  self.execute("sed 's/preproc/preproc_assign/' < ../Src/prealignment.txt > "+self.getAlignmentFilename(0))
	  self.linkFile("../../Src/referenceScaledVolume.vol",self.getModelFilename(0))
	  self.deleteFile("angle_convergence.txt")
	  self.touchFile("angle_convergence.txt")
	  self.deleteFile("shift_convergence.txt")
	  self.touchFile("shift_convergence.txt")
       else:
          if not os.path.exists("angle_convergence.txt"):
	     raise RuntimeError,"File angle_convergence.txt does not exist"
          self.execute("head -"+str(self.resumeIteration-1)+" angle_convergence.txt > inter")
	  self.execute("mv inter angle_convergence.txt")
          if not os.path.exists("shift_convergence.txt"):
	     raise RuntimeError,"File shift_convergence.txt does not exist"
          self.execute("head -"+str(self.resumeIteration-1)+" shift_convergence.txt > inter")
	  self.execute("mv inter shift_convergence.txt")
       for it in range(self.resumeIteration,99):
	  self.deleteDirectory("Iteration"+itoa(it,2))
       
   #------------------------------------------------------------------------
   # Reconstruct
   #------------------------------------------------------------------------
   def reconstruct(self,_iteration):
       self.log("# 3D reconstruction ---------------------------------------------------")
       factor=pow(2,-int(self.getPyramidLevel(_iteration)))
       
       # Apply a quality criteria
       cmd="xmipp_filter_projections"+\
          " -i "+self.getAlignmentFilename(_iteration)+\
          " -o preproc_recons -thr "+str(self.myNumberOfThreads)

       if self.getQualityPercentil(_iteration)>0 and \
          self.getDiscreteAssignment(_iteration):
          cmd+=" -filter_score "+self.getDiscreteAnglesFilename(_iteration)+" "+\
             self.getQualityPercentil(_iteration)

       if self.getQualityPercentil(_iteration)>0 and \
          self.getContinuousAssignment(_iteration):
          cmd+=" -filter_cost "+self.getContinuousAnglesFilename(_iteration)+" "+\
             self.getQualityPercentil(_iteration)

       if self.getQualityPercentil(_iteration)>0:
          cmd+=" -filter_normalization "+\
               self.getModelFFilename(_iteration)+" "+\
               str(math.ceil(self.particleWorkingRadius*factor))+" "+\
               str(math.ceil(self.particleWorkingRadius*factor*1.1))+" "+\
               self.getQualityPercentil(_iteration)

       if _iteration>1:
          cmd+=" -filter_movement "+\
               self.getAlignmentFFilename(_iteration)+" "+\
               self.getQualityAngleMovement(_iteration)+" "+\
               str(math.ceil(self.particleWorkingRadius*factor*\
                  float(self.getQualityShiftMovement(_iteration))/100.0))

       self.execute(cmd)
       self.execute("cp preproc_recons.doc "+\
          self.getReconstructionAnglesFilename(_iteration))
       
       self.runReconstructionAlgorithm(_iteration,"preproc_recons",
	  self.getReconstructionRootname(_iteration),True)

   #------------------------------------------------------------------------
   # Reconstruction iteration
   #------------------------------------------------------------------------
   def reconstructionIteration(self,_iteration):
      self.changeDirectory(self.workDirectory+"/Results")
      self.createDirectory("Iteration"+itoa(_iteration,2))
      self.adaptScaleFromPreviousIteration(_iteration)
      self.generateImagesForThisIteration(_iteration)
      if (debugSteps):
          a=input("Images prepared. Press a key")
      self.angularAssignment(_iteration)
      if (debugSteps):
          a=input("Angular assignment done. Press a key")
      self.reconstruct(_iteration)
      if (debugSteps):
          a=input("Reconstruction done. Press a key")
      if self.getComputeResolution(_iteration)=="1":
         self.computeResolution(_iteration)
      
   #------------------------------------------------------------------------
   # Run
   #------------------------------------------------------------------------
   def run(self):
       self.initDirectories()
       self.generateImagesPower2()
       self.copyCTFs()
       self.prealignment()
       self.prepareForFirstIteration()
       if (debugSteps):
           a=input("Prepared for first iteration. Press a key")
       for it in range(self.numberOfIterations):
           if it>=self.resumeIteration:
              self.log("Iteration "+str(it),'info',True)
	      self.reconstructionIteration(it)
	      self.generateNextReference(it)
       if not self.getComputeResolution(self.numberOfIterations-1)=="1":
          self.computeResolution(self.numberOfIterations-1)
       
   #------------------------------------------------------------------------
   # Run Reconstruction Algorithm
   #------------------------------------------------------------------------
   def runReconstructionAlgorithm(self,_iteration,_rootname,_outputRootName,
       _applyMasks):

       docfile = _rootname+".doc"
       selfile = _rootname+".sel"

       # Reconstruct
       if self.getReconstructionMethod(_iteration)=="art":
          params="-i "+selfile+" "+\
	         "-o "+_outputRootName+" "+\
		 "-l "+self.getARTLambda(_iteration)
          if not self.symmetryGroup=="":
             params+=" -sym "+self.symmetryGroup
	  if _applyMasks:
	     params+="-R "+str(math.ceil(self.particleWorkingRadius*1.1))
	  doParallel=self.doParallel
	  if self.getSerialART(_iteration)=="1":
	     doParallel=False
	  launch_parallel_job.launch_job(doParallel,
        			       "xmipp_reconstruct_art",
        			       "xmipp_mpi_reconstruct_art",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
	  self.deleteFile(_outputRootName+".hist")
       elif self.getReconstructionMethod(_iteration)=="wbp":
          params="-i "+selfile+" "+\
	         "-o "+_outputRootName+".vol "+\
		 "-use_each_image"
          if not self.symmetryGroup=="":
	     params+="-sym "+self.symmetryGroup+" "
	  if not _applyMasks:
	     params+=" -radius "+str(0.5*math.floor(math.sqrt(2.0)*
	        self.workXDim/
	        math.pow(2.0,int(self.getPyramidLevel(_iteration)))))
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_reconstruct_wbp",
        			       "xmipp_mpi_reconstruct_wbp",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUs,
        			       self.myMachineFile,
        			       False)
       elif self.getReconstructionMethod(_iteration)=="fourier":
          params="-i "+selfile+" "+\
	         "-o "+_outputRootName+".vol "
          if not self.symmetryGroup=="":
	     params+="-sym "+self.symmetryGroup+" "
          if self.getComputeFSC(_iteration):
            params+="-prepare_fsc "+self.getReconstructionRootname(_iteration)+\
               ".fsc "
	  launch_parallel_job.launch_job(self.doParallel,
        			       "xmipp_reconstruct_fourier",
        			       "xmipp_mpi_reconstruct_fourier",
        			       params,
        			       self.mylog,
        			       self.myNumberOfCPUsReduced,
        			       self.myMachineFile,
        			       False)
       else:
          raise RuntimeError,"Unknown reconstruction method"+\
	        self.getReconstructionMethod(_iteration)
       if not os.path.exists(_outputRootName+".vol"):
	  raise RuntimeError,"There is a problem when reconstructing"

       # Apply raised cosine mask
       if _applyMasks:
          self.execute("xmipp_mask -i "+_outputRootName+".vol "\
                       "-mask raised_cosine "+str(-self.particleWorkingRadius)+\
	               " "+str(-math.ceil(self.particleWorkingRadius*1.1)))
       
       # Symmetrize volume
       if not self.symmetryGroup=="":
          self.execute("xmipp_symmetrize -i "+_outputRootName+".vol "+\
	               "-sym "+self.symmetryGroup)
      
       # Mask volume
       if (not self.initialReferenceMask=="") and _applyMasks:
          self.execute("xmipp_mask -i "+_outputRootName+".vol "+\
	               "-mask "+self.getMaskFilename(_iteration))

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
   myMultiRes=MultiResClass(
      	        SelFileName,
		ReferenceFileName,
		WorkDirectory,
		DoDeleteWorkingDir,
		NumberofIterations,
		ProjectDir,
		LogDir,
                SkipPrealignment,
		
		ParticleRadius,
		ParticleMass,
		SymmetryGroup,
		SamplingRate,
		
		PyramidLevels,
		AngularSteps,
                QualityPercentil,
                QualityAngleMovement,
                QualityShiftMovement,
		ReconstructionMethod,
		SerialART,
		ARTLambda,
		DiscreteAssignment,
		ContinuousAssignment,
		DoComputeFSC,
		DoComputeSSNR,
		ResumeIteration,
		
		CTFDat,
		AmplitudeCorrection,
		
		DoReferenceMask,
		InitialReferenceMask,
		FilterLowPassReference,
		FilterHighPassReference,
		SegmentUsingMass,
		Recenter,

		DoParallel,
		MyNumberOfCPUs,
		MyNumberOfCPUsReduced,
                MyNumberOfThreads,
		MyMachineFile,
		
		True
              )
   myMultiRes.run()
