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
Iterations='1 2 3'
# Show model for angular assignment
DoShowModelF=False
# Show reconstructed volume
DoShowVolume=True
# Show model after postprocessing
DoShowModel=False
# Show angle convergence
DoShowAngleConvergence=False
# Show vector difference
DoShowVectorDifferences=False
# Show discrete angular assignment summary
DoShowDiscreteSummary=False
# Show resolution
DoShowResolution=False
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

scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
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
       import protocol_multires

       # Produce side info
       self.myMultiRes=protocol_multires.MultiResClass(
      	        protocol_multires.SelFileName,
		protocol_multires.ReferenceFileName,
		protocol_multires.WorkDirectory,
		protocol_multires.DoDeleteWorkingDir,
		protocol_multires.NumberofIterations,
		protocol_multires.ProjectDir,
		protocol_multires.LogDir,
		protocol_multires.ParticleRadius,
		protocol_multires.ParticleMass,
		protocol_multires.SymmetryFile,
		protocol_multires.SamplingRate,
		protocol_multires.PyramidLevels,
		protocol_multires.AngularSteps,
		protocol_multires.ReconstructionMethod,
		protocol_multires.SerialART,
		protocol_multires.ARTLambda,
		protocol_multires.DiscreteAssignment,
		protocol_multires.ContinuousAssignment,
		protocol_multires.DoComputeResolution,
		protocol_multires.ResumeIteration,
		protocol_multires.CTFDat,
		protocol_multires.PhaseCorrection,
		protocol_multires.AmplitudeCorrection,
		protocol_multires.DoReferenceMask,
		protocol_multires.InitialReferenceMask,
		protocol_multires.FilterReference,
		protocol_multires.SegmentUsingMass,
		protocol_multires.DoParallel,
		protocol_multires.MyNumberOfCPUs,
		protocol_multires.MyMachineFile,
		
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
          plotFSC=visualization.gnuplot()
          plotFSC.prepare_empty_plot("Resolution measure by FSC",
             "Armstrong^-1","Fourier Shell Correlation")
	  for i in range(len(self.iterationList)):
              fscFile=self.myMultiRes.getReconstructionRootname(
	         self.iterationList[i])+"_2.vol.frc"
	      if os.path.exists(fscFile):
                 if i==0:
                    plotFSC.send(" plot '" + fscFile +
                       "' using 1:2 title 'Iteration "+
                       str(self.iterationList[i])+"' with lines")
                    plotFSC.send(" replot '" + fscFile +
                       "' using 1:3 title 'Noise level' with lines")
                 else:
                    plotFSC.send(" replot '" + fscFile +
                       "' using 1:2 title 'Iteration "+
                       str(self.iterationList[i])+"' with lines")

       if not len(self.iterationList)==0:
          plotSSNR=visualization.gnuplot()
          plotSSNR.prepare_empty_plot("Resolution measured by SSNR",
             "Armstrong^-1","Spectral Signal-to-Noise Ratio")
	  for i in range(len(self.iterationList)):
              ssnrFile=self.myMultiRes.getReconstructionRootname(
	         self.iterationList[i])+".vol.ssnr"
	      if os.path.exists(ssnrFile):
                 if i==0:
                    plotSSNR.send(" plot '" + ssnrFile +
                       "' using 2:3 title 'Iteration "+
                       str(self.iterationList[i])+"' with lines")
                    plotSSNR.send(" replot 0 title 'Noise level' with lines")
                 else:
                    plotSSNR.send(" replot '" + fscFile +
                       "' using 1:2 title 'Iteration "+
                       str(self.iterationList[i])+"' with lines")
          plotSSNR.send("set yrange [-5:30]")
          plotSSNR.send("replot")

   #------------------------------------------------------------------------
   # Show Vector Differences
   #------------------------------------------------------------------------
   def showVectorDifferences(self):
       if not len(self.iterationList)==0:
          plot=visualization.gnuplot()
          plot.prepare_empty_plot("Vector differences","Degrees","Histogram")
	  for i in range(len(self.iterationList)):
              diffFile="Iteration"+itoa(self.iterationList[i],2)+\
	               "/diff_angles"+itoa(self.iterationList[i],2)+\
		       "_vec_diff_hist.txt"
              if i==0:
                 plot.send(" plot '" + diffFile +
                    "' using 1:2 title 'Iteration "+
                    str(self.iterationList[i])+"' with steps")
              else:
                 plot.send(" replot '" + diffFile +
                    "' using 1:2 title 'Iteration "+
                    str(self.iterationList[i])+"' with steps")

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
