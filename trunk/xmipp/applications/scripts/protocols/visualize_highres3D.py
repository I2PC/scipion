#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_highres.py
#
# Example use:
# python visualize_highres3D.py
# python protocol_gui.py visualize_highres3D.py
#
# Author: Carlos Oscar Sanchez Sorzano, April 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Iterations
Iterations="l15-19"
# Show reconstructed volume
DoShowVolume=1
# Show model
DoShowModel=1
# Show angle convergence
DoShowAngleConvergence=1
# Show vector differences
DoShowVectorDifferences=1
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
import os
import shutil
import string
import sys

scriptdir=os.path.expanduser('~')+'/trunk/Xmipp/Applications/Batches/Protocols_lib'
sys.path.append(scriptdir) # add default search path
import protocol_highres3D

#===========================================================================
# Beginning of the protocol
#===========================================================================
class VisualizeHighres3DClass:
   #------------------------------------------------------------------------
   # Class constructor
   #------------------------------------------------------------------------
   def __init__(self,
       _Iterations,
       _DoShowVolume,
       _DoShowModel,
       _DoShowAngleConvergence,
       _DoShowVectorDifferences):
       self.iterations=_Iterations
       self.doShowVolume=_DoShowVolume
       self.doShowModel=_DoShowModel
       self.doShowAngleConvergence=_DoShowAngleConvergence
       self.doShowVectorDifferences=_DoShowVectorDifferences

       # Produce side info
       self.myHighRes3D=protocol_highres3D.HighRes3DClass(
      	        protocol_highres3D.SelFileName,
		protocol_highres3D.ReferenceFileName,
		protocol_highres3D.WorkDirectory,
		protocol_highres3D.DoDeleteWorkingDir,
		protocol_highres3D.projectDirectory,
		protocol_highres3D.LogDir,
		
		protocol_highres3D.ParticleRadius,
		protocol_highres3D.ParticleMass,
		protocol_highres3D.SymmetryFile,
		protocol_highres3D.SamplingRate,
		
		protocol_highres3D.PyramidLevels,
		protocol_highres3D.AngularSteps,
		protocol_highres3D.UseART,
		protocol_highres3D.SerialART,
		protocol_highres3D.ARTLambda,
		protocol_highres3D.DiscreteAssignment,
		protocol_highres3D.ContinuousAssignment,
		protocol_highres3D.GaussianLowPass,
		protocol_highres3D.ComputeResolution,
		protocol_highres3D.ResumeIteration,
		
		protocol_highres3D.CTFDat,
		protocol_highres3D.AmplitudeCorrection,
		
		protocol_highres3D.DoReferenceMask,
		protocol_highres3D.InitialReferenceMask,
		protocol_highres3D.FilterReference,
		protocol_highres3D.SegmentUsingMass,

		protocol_highres3D.DoParallel,
		protocol_highres3D.MyNumberOfCPUs,
		protocol_highres3D.MyMachineFile,
		
		False
              )
       os.chdir(self.myHighRes3D.workDirectory+"/Results")
       self.expandIterations()

   #------------------------------------------------------------------------
   # Execute command
   #------------------------------------------------------------------------
   def execute(self,_command):
      os.system(_command)

   #------------------------------------------------------------------------
   # Expand iterations
   #------------------------------------------------------------------------
   def expandIterations(self):
      lista=string.split(self.iterations," ")
      if len(lista)==0: lista="last"
      self.iterationList=[]
      for i in range(len(lista)):
         listaIntervalo=string.split(lista[i],'-')
	 if len(listaIntervalo)==1:
	    if listaIntervalo[0]=="last":
	       lastIter=-1
	       for iteration in range(99):
	          if os.path.exists("Iteration"+itoa(iteration,2)):
		     lastIter=iteration
	       if lastIter>=0: self.iterationList+=[lastIter]
	    else:
	       self.iterationList+=[int(listaIntervalo[0])]
	 else:
	    self.iterationList+=range(int(listaIntervalo[0]),
	                              int(listaIntervalo[1]))

   #------------------------------------------------------------------------
   # Run
   #------------------------------------------------------------------------
   def run(self):
       if self.doShowVolume: self.showVolumes()
       if self.doShowModel:  self.showModels()
       if self.doShowAngleConvergence: self.showAngleConvergence()
       if self.doShowVectorDifferences: self.showVectorDifferences()

   #------------------------------------------------------------------------
   # Show Angle Convergence
   #------------------------------------------------------------------------
   def showAngleConvergence(self):
       command="echo plot \\\"angle_convergence.txt\\\" u 1:2 w l"
       command+=" \; pause 300 | gnuplot &"
       self.execute(command)
	     
   #------------------------------------------------------------------------
   # Show Models
   #------------------------------------------------------------------------
   def showModels(self):
       if not len(self.iterationList)==0:
          command="xmipp_show -vol "
	  for i in range(len(self.iterationList)):
	      command+=self.myHighRes3D.getReconstructionRootname(
	         self.iterationList[i])+".vol "
	  command+="&"
	  self.execute(command)
	     
   #------------------------------------------------------------------------
   # Show Vector Differences
   #------------------------------------------------------------------------
   def showVectorDifferences(self):
       if not len(self.iterationList)==0:
          command="echo plot "
	  for i in range(len(self.iterationList)):
	      command+="\\\"Iteration"+itoa(self.iterationList[i],2)+\
	               "/diff_angles"+itoa(self.iterationList[i],2)+\
		       "_vec_diff_hist.txt\\\" u 1:2 w st"
	      if not i==len(self.iterationList)-1:
	         command+=", "
          command+=" \; pause 300 | gnuplot &"
	  print command
	  self.execute(command)

   #------------------------------------------------------------------------
   # Show Volumes
   #------------------------------------------------------------------------
   def showVolumes(self):
       if not len(self.iterationList)==0:
          command="xmipp_show -vol "
	  for i in range(len(self.iterationList)):
	      command+=self.myHighRes3D.getReconstructionRootname(
	         self.iterationList[i])+".vol "
	  command+="&"
	  self.execute(command)
	     
#---------------------------------------------------------------------------
# itoa
#---------------------------------------------------------------------------
def itoa(_number,_length):
   return ("%0"+str(_length)+"d") % _number

#===========================================================================
# Main
#===========================================================================
if __name__ == '__main__':
    visualizeHighres3D=VisualizeHighres3DClass(
       Iterations,
       DoShowVolume,
       DoShowModel,
       DoShowAngleConvergence,
       DoShowVectorDifferences)
    visualizeHighres3D.run()

