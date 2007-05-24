#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_highres.py
#
# Example use:
# ./visualize_highres3d.py
#
# Author: Carlos Oscar Sanchez Sorzano, April 2007
# QUEDA POR HACER: Mostrar la distribucion angular del summary
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Iterations
Iterations='1-2'
# Show reconstructed volume
DoShowVolume=False
# Show model
DoShowModel=False
# Show angle convergence
DoShowAngleConvergence=False
# Show vector difference
DoShowVectorDifferences=False
# Show discrete angular assignment summary
DoShowDiscreteSummary=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
import os
import shutil
import string
import sys

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
                _DoShowVectorDifferences,
                _DoShowDiscreteSummary,
                _ProtocolName):

       self.iterations=_Iterations
       self.doShowVolume=_DoShowVolume
       self.doShowModel=_DoShowModel
       self.doShowAngleConvergence=_DoShowAngleConvergence
       self.doShowVectorDifferences=_DoShowVectorDifferences
       self.doShowDiscreteSummary=DoShowDiscreteSummary

       # Import the corresponding protocol
       pardir=os.path.abspath(os.getcwd())
       shutil.copy(ProtocolName,'protocol.py')
       import protocol

       # Produce side info
       self.myHighRes3D=protocol.HighRes3DClass(
      	        protocol.SelFileName,
		protocol.ReferenceFileName,
		protocol.WorkDirectory,
		protocol.DoDeleteWorkingDir,
		protocol.NumberofIterations,
		protocol.ProjectDir,
		protocol.LogDir,
		protocol.ParticleRadius,
		protocol.ParticleMass,
		protocol.SymmetryFile,
		protocol.SamplingRate,
		protocol.PyramidLevels,
		protocol.AngularSteps,
		protocol.UseART,
		protocol.SerialART,
		protocol.ARTLambda,
		protocol.DiscreteAssignment,
		protocol.ContinuousAssignment,
		protocol.ComputeResolution,
		protocol.ResumeIteration,
		protocol.CTFDat,
		protocol.AmplitudeCorrection,
		protocol.DoReferenceMask,
		protocol.InitialReferenceMask,
		protocol.FilterReference,
		protocol.SegmentUsingMass,
		protocol.DoParallel,
		protocol.MyNumberOfCPUs,
		protocol.MyMachineFile,
		
		False
              )
       # Remove protocol.py(c)
       if (os.path.exists('protocol.py')):
           os.remove('protocol.py')
       if (os.path.exists('protocol.pyc')):
           os.remove('protocol.pyc')

       os.chdir(self.myHighRes3D.workDirectory+"/Results")
       self.expandIterations()

   #------------------------------------------------------------------------
   # Execute command
   #------------------------------------------------------------------------
   def execute(self,_command):
      print _command
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
	                              int(listaIntervalo[1])+1)

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
       if self.doShowDiscreteSummary: self.showDiscreteSummary()

   #------------------------------------------------------------------------
   # Show Angle Convergence
   #------------------------------------------------------------------------
   def showAngleConvergence(self):
       command="echo plot \\\"angle_convergence.txt\\\" u 1:2 w l"
       command+=" \; pause 300 | gnuplot &"
       self.execute(command)
	     
   #------------------------------------------------------------------------
   # Show Angular Assignment Discrete Summary
   #------------------------------------------------------------------------
   def showDiscreteSummary(self):
       if not len(self.iterationList)==0:
 	  import visualize_projmatch
	  for i in range(len(self.iterationList)):
              ShowPlots=[]
	      ShowPlots.append(self.myHighRes3D.getDiscreteAnglesSummaryDir(
 	         self.iterationList[i])+"/summary_summary.doc")
	      title='Angular distribution for iteration '+\
        	     str(self.iterationList[i])
      	      visualize_projmatch.show_ang_distribution(
	         ShowPlots,self.iterationList[i],title)
	      print "Calling",i

	  for i in range(len(self.iterationList)):
              command="( cd "+self.myHighRes3D.workDirectory+"/Results/"+\
	         self.myHighRes3D.getDiscreteAnglesSummaryDir(
	            self.iterationList[i])+\
		    " ; xmipp_show -sel summary_comparison.sel ) &"
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

    import sys
    ProtocolName=sys.argv[1]

    visualizeHighres3D=VisualizeHighres3DClass(Iterations,
                                               DoShowVolume,
                                               DoShowModel,
                                               DoShowAngleConvergence,
                                               DoShowVectorDifferences,
                                               DoShowDiscreteSummary,
                                               ProtocolName)

    visualizeHighres3D.run()

