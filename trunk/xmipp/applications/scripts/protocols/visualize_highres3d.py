#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_highres.py
#
# Example use:
# ./visualize_highres3d.py
#
# Author: Carlos Oscar Sanchez Sorzano, April 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Iterations
Iterations='1-last'
# Show reconstructed volume
DoShowVolume=True
# Show model
DoShowModel=True
# Show angle convergence
DoShowAngleConvergence=True
# Show vector difference
DoShowVectorDifferences=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
import os
import shutil
import string
import sys

scriptdir=os.path.expanduser('~')+'/scripts'
sys.path.append(scriptdir) # add default search path
import protocol_highres3d

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
       self.myHighRes3D=protocol_highres3d.HighRes3DClass(
      	        protocol_highres3d.SelFileName,
		protocol_highres3d.ReferenceFileName,
		protocol_highres3d.WorkDirectory,
		protocol_highres3d.DoDeleteWorkingDir,
		protocol_highres3d.ProjectDir,
		protocol_highres3d.LogDir,
		
		protocol_highres3d.ParticleRadius,
		protocol_highres3d.ParticleMass,
		protocol_highres3d.SymmetryFile,
		protocol_highres3d.SamplingRate,
		
		protocol_highres3d.PyramidLevels,
		protocol_highres3d.AngularSteps,
		protocol_highres3d.UseART,
		protocol_highres3d.SerialART,
		protocol_highres3d.ARTLambda,
		protocol_highres3d.DiscreteAssignment,
		protocol_highres3d.ContinuousAssignment,
		protocol_highres3d.ComputeResolution,
		protocol_highres3d.ResumeIteration,
		
		protocol_highres3d.CTFDat,
		protocol_highres3d.AmplitudeCorrection,
		
		protocol_highres3d.DoReferenceMask,
		protocol_highres3d.InitialReferenceMask,
		protocol_highres3d.FilterReference,
		protocol_highres3d.SegmentUsingMass,

		protocol_highres3d.DoParallel,
		protocol_highres3d.MyNumberOfCPUs,
		protocol_highres3d.MyMachineFile,
		
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
	       lastIter=self.lastIteration()
	       if lastIter>=0: self.iterationList+=[lastIter]
	    else:
	       self.iterationList+=[int(listaIntervalo[0])]
	 else:
	    if listaIntervalo[1]=="last":
	       lastIter=self.lastIteration()
	       if lastIter>=0: self.iterationList+=[lastIter]
	       else: lastIter=listaIntervalo[0]
	       listaIntervalo[1]=lastIter
	    self.iterationList+=range(int(listaIntervalo[0]),
	                              int(listaIntervalo[1]))

   #------------------------------------------------------------------------
   # Last iteration
   #------------------------------------------------------------------------
   def lastIteration(self):
      lastIter=-1
      for iteration in range(99):
	 if os.path.exists(self.myHighRes3D.getModelFilename(iteration)):
	    lastIter=iteration
      return lastIter

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
	      command+=self.myHighRes3D.getModelFilename(
	         self.iterationList[i])+" "
	  command+="&"
	  print command
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

