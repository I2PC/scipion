#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based ML-refinement of subtomograms or RCT reconstructions, according to:
# {please cite} Scheres et al. (2009) submitted
#
# Example use:
# ./xmipp_protocol_mltomo.py
#
# Author: Sjors Scheres, August 2009
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------

# {file} Selfile with the input images:
""" This selfile points to the spider single-file format 3D-images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile in the PROJECTDIR. 
"""
InSelFile='images.sel'

# {file} Docfile with parameters for all input images
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Ml_tomo for a detailed explanation of this format. Note that all images in the selfile should have an entry in this docfile (with exactly the same image name), but the docfile may contain more images than the selfile.
"""
InDocFile='images.doc'

# {file} Docfile with the definition of the missing data regions
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Ml_tomo for a detailed explanation of this format. If your data have no missing regions in Fourier space, leave this field empty 
"""
MissingDocFile='wedges.doc'

# Working subdirectory:
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='MLtomo/run1'

# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir=False

# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/home/scheres/tomo'

# {expert} Directory name for logfiles:
LogDir='Logs'

#------------------------------------------------------------------------------------------------
# {section} General refinement parameters
#------------------------------------------------------------------------------------------------
# Output rootname for this run
""" Use different names if you want to run multiple runs in the same WorkingDir
"""
OutputRootName='mltomo'

# Number of ML iterations to perform:
NumberOfIterations=25

# Keep the angles from the input docfile?
""" If set to false, random assignments will be used for all three Euler angles. This option is usually only set to True after a previous run of this program.
""" 
DoKeepAngles=False

# {expert} Restart after iteration:
""" For previous runs that stopped before convergence,
    resume the calculations after the completely finished iteration. 
    Set to zero to start a new run.
"""
RestartIter=0

# Point group symmetry:
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the point group symmetries
    Give c1 if no symmetry is present
"""
Symmetry='c1'

# Dimension of the image to be used
""" Internally, the images will be resized to this value (in pixels). Use a negative value to use the original image size. 
"""
Dimension=-1

# {expert} Maximum resolution (in pixel^-1) to be used
""" The maximum (Nyquist) resolution is 0.5. Use smaller values, e.g. 0.45, to prevent high-resolution artifacts.
"""
MaximumResolution=0.45

# {expert} Additional xmipp_ml_tomo parameters:
""" For a complete description see the manual pages:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Ml_tomo
"""
ExtraParamsMLtomo=''

#------------------------------------------------------------------------------------------------
# {section} Alignment parameters
#------------------------------------------------------------------------------------------------
# Allow images to change orientation?
""" Only set this option to False when you want to classify without changing the orientations of the images.
"""
DoAlign=True

# Angular sampling for ML refinement
""" Fine samplings take huge amounts of CPU and memory. Therefore, one often performs a first run with coarse samplings and wide angular search ranges (e.g. exhaustive searches), and in subsequent runs one then uses finer samplings with reduced search ranges. 
"""
AngularSampling=15

# Angular search range for ML refinement
""" All angles within the value given here around the angles in the input docfile will be searched for every input image. Use a large value (e.g. 1000) for exhaustive angular searches. 
"""
AngularSearchRange=1000

# {expert} Apply small random perturbations on the angular sampling?
""" This option is recommended, as it will prevent the algorithm from getting stuck in local minima.
"""
DoPerturb=True

#------------------------------------------------------------------------------------------------
# {section} Classification parameters
#------------------------------------------------------------------------------------------------

# Number of classes to use in the refinement:
""" Use 1 for a run with only aligment and no classification
"""
NumberOfReferences=3

# {expert} {file} Provide a selfile with user-defined references:
""" By default, the references will be created from random assignments of the angles and classes (or only of the classes, when the angles from a previous run are kept). By providing a selfile with references here, one may supervise the classification. Note that these references should be on the correct intensity greyscale
"""
SeedsSelfile=''

# Initial regularization parameter
""" We have found that imposing similarity on the difference references during the initial iterations of ML refinement may prevent the classification from getting stuck into local minima.
"""
InitialRegularization=5

# Number of iterations to use regularization
""" Useful values are 3-5. Note that during these iterations, the references will be enforced to be similar, so do not expect any useful classification until after these iterations. 
"""
NumberRegularizationSteps=5

# {expert} Mask for focussed classification 
""" A mask can only be used when the images are NOT allowed to change their orientation (see Alignment parameters). Leave this field empty to use the entire images. 
"""
MaskName=''

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------

# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads=8

# Use distributed-memory parallelization (MPI)?
""" This option provides parallelization on clusters with distributed memory architecture.
    It requires mpi to be installed.
"""
DoParallel=True

# Number of MPI processes to use:
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
NumberOfMpiProcesses=2

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, 
    different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. Or read: xxx
"""
SystemFlavour='TORQUE-OPENMPI'

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_mltomo.py'

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class MLtomo_class:

    #init variables
    def __init__(self,
                 InSelFile,
                 InDocFile,
                 MissingDocFile,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 ProjectDir,
                 LogDir,
                 OutputRootName,
                 NumberOfIterations,
                 DoKeepAngles,
                 RestartIter,
                 Symmetry,
                 Dimension,
                 MaximumResolution,
                 ExtraParamsMLtomo,
                 DoAlign,
                 AngularSampling,
                 AngularSearchRange,
                 DoPerturb,
                 NumberOfReferences,
                 SeedsSelfile,
                 InitialRegularization,
                 NumberRegularizationSteps,
                 MaskName,
                 NumberOfThreads,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,selfile

        self.InSelFile=InSelFile
        self.InDocFile=InDocFile
        self.MissingDocFile=MissingDocFile
        self.WorkingDir=WorkingDir
        self.DoDeleteWorkingDir=DoDeleteWorkingDir
        self.ProjectDir=ProjectDir
        self.LogDir=LogDir
        self.OutputRootName=OutputRootName
        self.NumberOfIterations=NumberOfIterations
        self.DoKeepAngles=DoKeepAngles
        self.RestartIter=RestartIter
        self.Symmetry=Symmetry
        self.Dimension=Dimension
        self.MaximumResolution=MaximumResolution
        self.ExtraParamsMLtomo=ExtraParamsMLtomo
        self.DoAlign=DoAlign
        self.AngularSampling=AngularSampling
        self.AngularSearchRange=AngularSearchRange
        self.DoPerturb=DoPerturb
        self.NumberOfReferences=NumberOfReferences
        self.SeedsSelfile=SeedsSelfile
        self.InitialRegularization=InitialRegularization
        self.NumberRegularizationSteps=NumberRegularizationSteps
        self.MaskName=MaskName
        self.NumberOfThreads=NumberOfThreads
        self.DoParallel=DoParallel
        self.NumberOfMpiProcesses=NumberOfMpiProcesses
        self.SystemFlavour=SystemFlavour

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     self.LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)

        # This is not a restart
        if (self.RestartIter < 1):
            # Delete working directory if it exists, make a new one
            if (DoDeleteWorkingDir): 
                if os.path.exists(self.WorkingDir):
                    shutil.rmtree(self.WorkingDir)
            if not os.path.exists(self.WorkingDir):
                os.makedirs(self.WorkingDir)


        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                       os.path.abspath(self.WorkingDir))
    
        # This protocol is executed from the ProjectDir
        self.execute_MLtomo()


    # Execute the actual program
    def execute_MLtomo(self):
        import os
        import launch_job, utils_xmipp

        if (self.RestartIter>0):
            print '*********************************************************************'
            print '*  Restarting ml_tomo program :' 
            restartname = utils_xmipp.composeFileName('ml3d_it',iter,'log')
            params= ' -restart ' + restartname
        else:
        
            print '*********************************************************************'
            print '*  Executing ml_tomo program :' 

            # General parameters
            params= ' -i '    + str(self.InSelFile) + \
                    ' -doc '  + str(self.InDocFile) + \
                    ' -o '    + str(self.WorkingDir) + "/" + str(self.OutputRootName) + \
                    ' -iter ' + str(self.NumberOfIterations) + \
                    ' -sym '  + str(self.Symmetry) +\
                    ' -maxres ' + str(self.MaximumResolution)
            if (len(self.MissingDocFile) > 0):
                params+= ' -missing ' + str(self.MissingDocFile)
            if (self.Dimension > 0):
                params+= ' -dim ' + str(self.Dimension)
            if (self.DoKeepAngles):
                params+= ' -keep_angles '

            # Alignment parameters
            if (self.DoAlign):
                params+= ' -ang ' + str(self.AngularSampling)
                if (self.AngularSearchRange<360):
                    params+= ' -ang_search ' + str(self.AngularSearchRange)
                if (self.DoPerturb):
                    params+= ' -perturb '
            else:
                params+= ' -dont_align '

            # Classification parameters
            if (len(self.SeedsSelfile)>0):
                params+= ' -ref ' + str(self.SeedsSelfile)
            else:
                params+= ' -nref ' + str(self.NumberOfReferences)
            if (self.InitialRegularization>0):
                params+= ' -reg0 ' + str(self.InitialRegularization)
                params+= ' -reg_steps ' + str(self.NumberRegularizationSteps)
            if (len(self.MaskName) > 0):
                if (self.DoAlign):
                    print '********************************************************************************'
                    print '*  ERROR: you can only provide a mask if you do NOT allow orientations to change'
                    self.close()
                else:
                    params+= ' -mask ' + str(self.MaskName)

            # Extra parameters
            if (len(self.ExtraParamsMLtomo) > 0):
                params+=' ' + str (self.ExtraParamsMLtomo)

            # Thread parallelization
            if (self.NumberOfThreads > 1):
                params+=' -thr ' + str(self.NumberOfThreads)

                
        launch_job.launch_job("xmipp_ml_tomo",
                              params,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpiProcesses,
                              self.NumberOfThreads,
                              self.SystemFlavour)


    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    MLtomo=MLtomo_class(InSelFile,
                        InDocFile,
                        MissingDocFile,
                        WorkingDir,
                        DoDeleteWorkingDir,
                        ProjectDir,
                        LogDir,
                        OutputRootName,
                        NumberOfIterations,
                        DoKeepAngles,
                        RestartIter,
                        Symmetry,
                        Dimension,
                        MaximumResolution,
                        ExtraParamsMLtomo,
                        DoAlign,
                        AngularSampling,
                        AngularSearchRange,
                        DoPerturb,
                        NumberOfReferences,
                        SeedsSelfile,
                        InitialRegularization,
                        NumberRegularizationSteps,
                        MaskName,
                        NumberOfThreads,
                        DoParallel,
                        NumberOfMpiProcesses,
                        SystemFlavour)

    
    # close 
    MLtomo.close()

