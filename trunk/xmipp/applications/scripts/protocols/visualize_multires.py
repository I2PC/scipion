#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_multires.py
#
# Example use:
# ./visualize_multires.py
#
# Author: Carlos Oscar Sanchez Sorzano, April 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Iterations
Iterations='1'
# Show model for angular assignment
DoShowModelF=False
# Show reconstructed volume
DoShowVolume=False
# Show model after postprocessing
DoShowModel=False
# Show angle convergence
DoShowAngleConvergence=False
# Show vector difference
DoShowVectorDifferences=False
# Show discrete angular assignment summary
DoShowDiscreteSummary=False
# Show resolution
DoShowResolution=True
#------------------------------------------------------------------------------------------------
# {section} Volume visualization
#------------------------------------------------------------------------------------------------
# Visualize volumes in slices along Z?
VisualizeVolZ=True
# Visualize volumes in slices along X?
VisualizeVolX=False
# Visualize volumes in slices along Y?
VisualizeVolY=False
# Visualize volumes in UCSF Chimera?
""" For this to work, you need to have chimera installed!
"""
VisualizeVolChimera=False
# {expert} Width of xmipp_show (multiple of 3!):
MatrixWidth=9
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
import os
import shutil
import string
import sys

scriptdir=os.path.expanduser('~')+'/scripts/'
sys.path.append(scriptdir) # add default search path
import visualization

#===========================================================================
# Beginning of the protocol
#===========================================================================
class VisualizeMultires3DClass:
   #------------------------------------------------------------------------
   # Class constructor
   #------------------------------------------------------------------------
   def __init__(self,
                _Iterations,
                _DoShowModelF,
                _DoShowVolume,
                _DoShowModel,
                _DoShowAngleConvergence,
                _DoShowVectorDifferences,
                _DoShowDiscreteSummary,
		_DoShowResolution,
		_VisualizeVolZ,
		_VisualizeVolX,
		_VisualizeVolY,
		_VisualizeVolChimera,
		_MatrixWidth,
                _ProtocolName):

       self.iterations=_Iterations
       self.doShowModelF=_DoShowModelF
       self.doShowVolume=_DoShowVolume
       self.doShowModel=_DoShowModel
       self.doShowAngleConvergence=_DoShowAngleConvergence
       self.doShowVectorDifferences=_DoShowVectorDifferences
       self.doShowDiscreteSummary=DoShowDiscreteSummary
       self.doShowResolution=_DoShowResolution
       self.visualizeVolZ=_VisualizeVolZ
       self.visualizeVolX=_VisualizeVolX
       self.visualizeVolY=_VisualizeVolY
       self.visualizeVolChimera=_VisualizeVolChimera
       self.matrixWidth=_MatrixWidth

       # Import the corresponding protocol
       import protocol_highres3d

       # Produce side info
       self.myMultiRes=protocol_highres3d.MultiResClass(
      	        protocol_highres3d.SelFileName,
		protocol_highres3d.ReferenceFileName,
		protocol_highres3d.WorkDirectory,
		protocol_highres3d.DoDeleteWorkingDir,
		protocol_highres3d.NumberofIterations,
		protocol_highres3d.ProjectDir,
		protocol_highres3d.LogDir,
		protocol_highres3d.ParticleRadius,
		protocol_highres3d.ParticleMass,
		protocol_highres3d.SymmetryFile,
		protocol_highres3d.SamplingRate,
		protocol_highres3d.PyramidLevels,
		protocol_highres3d.AngularSteps,
		protocol_highres3d.ReconstructionMethod,
		protocol_highres3d.SerialART,
		protocol_highres3d.ARTLambda,
		protocol_highres3d.DiscreteAssignment,
		protocol_highres3d.ContinuousAssignment,
		protocol_highres3d.DoComputeResolution,
		protocol_highres3d.ResumeIteration,
		protocol_highres3d.CTFDat,
		protocol_highres3d.PhaseCorrection,
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

       os.chdir(self.myMultiRes.workDirectory+"/Results")
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
	 if os.path.exists(self.myMultiRes.getModelFilename(iteration)):
	    lastIter=iteration
      return lastIter

   #------------------------------------------------------------------------
   # Run
   #------------------------------------------------------------------------
   def run(self):
       if self.doShowModelF: self.showModelsF()
       if self.doShowVolume: self.showVolumes()
       if self.doShowModel:  self.showModels()
       if self.doShowAngleConvergence: self.showAngleConvergence()
       if self.doShowVectorDifferences: self.showVectorDifferences()
       if self.doShowDiscreteSummary: self.showDiscreteSummary()
       if self.doShowResolution: self.showResolution()

   #------------------------------------------------------------------------
   # Show Angle Convergence
   #------------------------------------------------------------------------
   def showAngleConvergence(self):
       plot=visualization.gnuplot()
       plot.plot_xy_file("angle_convergence.txt",
                         Title="Angular assignment convergence",
                         X_Label="Iteration",
                         Y_Label="Degrees",
                         X_col=1,
                         Y_col=2)
	     
   #------------------------------------------------------------------------
   # Show Angular Assignment Discrete Summary
   #------------------------------------------------------------------------
   def showDiscreteSummary(self):
       if not len(self.iterationList)==0:
 	  import visualize_projmatch
	  for i in range(len(self.iterationList)):
              ShowPlots=[]
	      ShowPlots.append(self.myMultiRes.getDiscreteAnglesSummaryDir(
 	         self.iterationList[i])+"/summary_summary.doc")
	      title=[]
	      title.append('Angular distribution for iteration '+\
        	 str(self.iterationList[i]))
      	      visualize_projmatch.show_ang_distribution(
	         ShowPlots,self.iterationList[i],title)

	  for i in range(len(self.iterationList)):
              command="( cd "+self.myMultiRes.workDirectory+"/Results/"+\
	         self.myMultiRes.getDiscreteAnglesSummaryDir(
	            self.iterationList[i])+\
		    " ; xmipp_show -sel summary_comparison.sel ) &"
 	      self.execute(command)
	     
   #------------------------------------------------------------------------
   # Show Models
   #------------------------------------------------------------------------
   def showModels(self):
       if not len(self.iterationList)==0:
          volumeList=[]
	  for i in range(len(self.iterationList)):
	      volumeList.append(self.myMultiRes.getModelFilename(
	         self.iterationList[i]))
          visualization.visualize_volumes(volumeList,
                                          self.visualizeVolZ,
                                          self.visualizeVolX,
                                          self.visualizeVolY,
                                          self.visualizeVolChimera)
	     
   #------------------------------------------------------------------------
   # Show ModelsF
   #------------------------------------------------------------------------
   def showModelsF(self):
       if not len(self.iterationList)==0:
          volumeList=[]
	  for i in range(len(self.iterationList)):
	      volumeList.append(self.myMultiRes.getModelFFilename(
	         self.iterationList[i]))
          visualization.visualize_volumes(volumeList,
                                          self.visualizeVolZ,
                                          self.visualizeVolX,
                                          self.visualizeVolY,
                                          self.visualizeVolChimera)
	     
   #------------------------------------------------------------------------
   # Show Resolution
   #------------------------------------------------------------------------
   def showResolution(self):
       if not len(self.iterationList)==0:
          launchGnuplot=False
          command="echo plot "
	  for i in range(len(self.iterationList)):
	      if os.path.exists(self.myMultiRes.getReconstructionRootname(
	         self.iterationList[i])+"_2.vol.frc"):
		 if launchGnuplot: command+=", "
         	 command+="\\\""+self.myMultiRes.getReconstructionRootname(
	                  self.iterationList[i])+\
			  '_2.vol.frc\\\" u 1:2 title \\"FSC iteration '+\
			  str(self.iterationList[i])+'\\" w l'
		 if not launchGnuplot:
		    command+=", \\\""+self.myMultiRes.getReconstructionRootname(
                             self.iterationList[i])+\
		             '_2.vol.frc\\\" u 1:3 title \\"Noise threshold\\" w l'
		 launchGnuplot=True
    	  if launchGnuplot:
             command+=" \; set xlabel \\\"1/Angstrom\\\" \; replot \; "+\
	              "pause 300 | gnuplot &"
	     self.execute(command)

       if not len(self.iterationList)==0:
          launchGnuplot=False
          command="echo plot "
	  for i in range(len(self.iterationList)):
	      if os.path.exists(self.myMultiRes.getReconstructionRootname(
	         self.iterationList[i])+".vol.ssnr"):
		 if launchGnuplot: command+=", "
         	 command+="\\\""+self.myMultiRes.getReconstructionRootname(
	                  self.iterationList[i])+\
			  '.vol.ssnr\\\" u 2:3 title \\"SSNR iteration '+\
			  str(self.iterationList[i])+'\\" w l'
		 if not launchGnuplot:
		    command+=", 0"
		 launchGnuplot=True
    	  if launchGnuplot:
             command+=" \; set xlabel \\\"1/Angstrom\\\" \; "+\
	              "set yrange [-5:30] \; replot \; pause 300 | gnuplot &"
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
          volumeList=[]
	  for i in range(len(self.iterationList)):
	      volumeList.append(self.myMultiRes.getReconstructionRootname(
	         self.iterationList[i])+".vol")
          visualization.visualize_volumes(volumeList,
                                          self.visualizeVolZ,
                                          self.visualizeVolX,
                                          self.visualizeVolY,
                                          self.visualizeVolChimera)
	     
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

    visualizeMultires3D=VisualizeMultires3DClass(Iterations,
                			       DoShowModelF,
                			       DoShowVolume,
                			       DoShowModel,
                			       DoShowAngleConvergence,
                			       DoShowVectorDifferences,
                			       DoShowDiscreteSummary,
					       DoShowResolution,
					       VisualizeVolZ,
					       VisualizeVolX,
					       VisualizeVolY,
					       VisualizeVolChimera,
					       MatrixWidth,
                			       ProtocolName)

    visualizeMultires3D.run()
