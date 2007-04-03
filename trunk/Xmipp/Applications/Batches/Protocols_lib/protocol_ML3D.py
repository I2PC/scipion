#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based ML3D classification, according to:
# - Scheres et al. (2007), Nature Methods, 4, 27-29
#
# Example use:
# ./protocol_ML3D.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Selfile with the input images (relative path from ProjectDir):
InSelFile="all_images.sel"
# Working subdirectory:
WorkingDir="ML3ref"
# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
""" All logfiles will be stored in $ProjectDir/$LogDir
"""
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Correct absolute grey scale of initial reference
#------------------------------------------------------------------------------------------------
# Correct the absolute grey-scale of the initial reference map?
""" The probabilities are based on squared differences, so that the absolute grey scale is important.
"""
DoCorrectGreyScale=True
# {expert} Angular sampling for a quick projection matching to obtain right grey scale
""" As the resolution of the intial reference should be low, this sampling can be relatively crude, e.g. 15
"""
ProjMatchSampling=15
# {expert} Threshold for a quick weighted back-projection to obtain right grey scale
""" As the resolution of the intial reference should be low, this value can be relatively high, e.g. 0.02
"""
WbpThreshold=0.02
#------------------------------------------------------------------------------------------------
# {section} Seed generation
#------------------------------------------------------------------------------------------------
# Generate seeds from a single initial reference map?
DoGenerateSeeds=True
# Number of seeds to be generated (and later on used in ML3D-classification):
NumberOfReferences=3
# Initial reference map to create seeds:
InitialReference="ref.vol"
#------------------------------------------------------------------------------------------------
# {section} ML3D-classification parameters
#------------------------------------------------------------------------------------------------
# Perform ML3D classification run?
DoML3DClassification=False
# Angular sampling for ML3D classification:
""" Fine samplings take huge amounts of CPU and memory.
    Therefore, in general, dont use samplings finer than 10 degrees.
"""
AngularSampling=10
# Number of ML3D iterations to perform:
NumberOfIterations=25
# Symmetry description file:
""" See WIKI link for a description of the symmetry file format
    Dont give anything, if no symmetry is present
"""
SymmetryFile=""
# {expert} Additional xmipp_MLalign2D parameters:
ExtraParamsMLrefine3D=""
# {expert} Selfile with user-provided seeds (from ProjectDir):
""" This option is NOT recommended!
    The seeds should already be on the correct absolute greyscale!
"""
SeedsSelfile=""
#------------------------------------------------------------------------------------------------
# {section} Job submission issues
#------------------------------------------------------------------------------------------------
# Batch submission command (use "" to launch without batch submission):
""" This will depend on your queueing system., ask your system administrator...

    Examples: LaunchJobCommand=\"bsub -q 1day\"
      or, if you do not use a queueing system: LaunchJobCommand=\"\"
"""
LaunchJobCommand="" 
# Use multiple processors in parallel? (see Expert options)
DoParallel=False
# Number of processors to use:
MyNumberOfCPUs=10
# {expert} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
"""
MyMachineFile="/home2/bioinfo/scheres/machines.dat"
# {expert} Standard script to launch MPI-jobs:
""" This will also depend on your system...
    The simplest script consists of the following two lines:

    #!/usr/bin/env sh
    `which mpirun` -np MyNumberOfCPUs -machinefile MyMachineFile ~/machines.dat \

    Note that the next line with the xmipp_mpi_MLalign2D command will be
    generated automatically, and the variables MyNumberOfCPUs and MyMachineFile
    will be replaced by the corresponding values given here above.
    More scripts for different batch systems can be found at:
    [Wiki link]
"""
ParallelScript="/home2/bioinfo/scheres/submit_mpi_job.sh"
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_rct.py"
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class ML2D_class:

    #init variables
    def __init__(self,
                 InSelFile,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 ProjectDir,
                 LogDir,
                 DoCorrectGreyScale,
                 ProjMatchSampling,
                 WbpThreshold,
                 DoGenerateSeeds,
                 NumberOfReferences,
                 InitialReference,
                 DoML3DClassification,
                 AngularSampling,
                 NumberOfIterations,
                 SymmetryFile,
                 ExtraParamsMLrefine3D,
                 SeedsSelfile,
                 LaunchJobCommand,
                 DoParallel,
                 MyNumberOfCPUs,
                 MyMachineFile,
                 ParallelScript):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self.InSelFile=self.ProjectDir+'/'+str(InSelFile)
        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.NumberOfReferences=NumberOfReferences
        self.DoGenerateSeeds=DoGenerateSeeds
        # STILL SORT OUT ABS PATHS FOR THIS ONE!!
        self.InitialReference=InitialReference
        self.SeedsSelfile=SeedsSelfile
        self.AngularSampling=AngularSampling
        self.NumberOfIterations=NumberOfIterations
        self.SymmetryFile=SymmetryFile
        self.ExtraParamsMLrefine3D=ExtraParamsMLrefine3D
        self.DoCorrectGreyScale=DoCorrectGreyScale
        self.ProjMatchSampling=ProjMatchSampling
        self.WbpThreshold=WbpThreshold
        self.DoParallel=DoParallel
        self.MyNumberOfCPUs=MyNumberOfCPUs
        self.MyMachineFile=MyMachineFile
        self.LaunchJobCommand=LaunchJobCommand
        self.ParallelScript=ParallelScript

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Delete working directory if it exists, make a new one, and go there
        if (DoDeleteWorkingDir and DoML2D): 
            if os.path.exists(self.WorkingDir):
                shutil.rmtree(self.WorkingDir)
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Execute MLalign2D in the working directory
        os.chdir(self.WorkingDir)

        if self.DoCorrectGreyScale:
            self.correct_greyscale()

        if self.DoGenerateSeeds:
            self.generate_seeds()

        if self.DoML3DClassification:
            self.execute_ML3D_classification()
        
        # Return to parent dir
        os.chdir(os.pardir)


    # Crude correction of grey-scale, by performing a single iteration of projection matching and wbp
    def correct_greyscale(self):
        import os
        import launch_parallel_job
        print '*********************************************************************'
        print '*  Correcting absolute grey scale of initial reference:'

        # A single cycle of projection matching
        basename='corrected_'+os.path.namebase(self.InitialReference)
        params= ' -i '    + str(self.InSelFile) + \
                ' -o '    + str(basename) + \
                ' -vol '  + str(self.InitialReference) + \
                ' -sam ' + str(ProjMatchSampling)
        if not self.SymmetryFile=="":
            param+= ' -sym '+str(self.SymmetryFile)
        params+=' -dont_modify_header -output_classes'
                
        launch_parallel_job.launch_job(self.DoParallel,
                                       "xmipp_projection_matching",
                                       "xmipp_mpi_projection_matching",
                                       params,
                                       self.ParallelScript,
                                       self.LaunchJobCommand,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       self.WorkingDir,
                                       False)

        # Followed by a weighted back-projection reconstruction
        iname=basename+'_classes.sel'
        outname=basename+'.vol'
        params= ' -i '    + str(iname) + \
                ' -o '    + str(outname) + \
                ' -threshold '+str(self.WbpThreshold) + \
                '  -use_each_image -weight '
        if not self.SymmetryFile=="":
            param+= ' -sym '+str(self.SymmetryFile)
           
        launch_parallel_job.launch_job(self.DoParallel,
                                       "xmipp_wbp",
                                       "xmipp_mpi_wbp",
                                       params,
                                       self.ParallelScript,
                                       self.LaunchJobCommand,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       self.WorkingDir,
                                       False)

        return outname


    # Splits selfile and performs a sinmgle cycle of ML3D-classification for each subset
    def generate_seeds(self):
        import os
        import SelFiles
        print '*********************************************************************'
        print '*  Generating seeds:' 
        newsel=SelFiles.selfile()

        # Split selfiles
        command='xmipp_split_selfile -o seeds_split'+ \
                 ' -i ' + str(self.InSelFile) + \
                 ' -n ' + str(self.NumberOfReferences) +' \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

        corrvol='corrected_'+os.path.namebase(self.InitialReference)+'.vol'
        if os.path.exists(corrvol):
            reference=corrvol
        else:
            reference=self.InitialReference

        # Launch MLrefine3D with output to subdirectories
        for i in range(self.NumberOfReferences):
            inselfile='seeds_split_'+str(i+1)+'.sel'
            dirname='MLseeds_split_'+str(i+1)'/'
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            outname=dirname+'MLseeds_split_'+str(i+1)
            self.execute_MLrefine3D(inselfile,
                                    outname,
                                    reference,
                                    self.AngularSampling,
                                    1,
                                    self.SymmetryFile,
                                    self.ExtraParamsMLrefine3D)
            newsel.insert(dirname+'MLseeds_split_'+str(i+1)+'_it00001.vol')
        newsel.write('ml3d_seeds.sel')


    # Perform the actual classification
    def execute_ML3D_classification(self):
        import os
        import launch_parallel_job
 
        outname='ml3d'
        if self.SeedsSelfile=="":
            volname='ml3d_seeds.sel'
        else:
            volname=self.SeedsSelfile
        self.execute_MLrefine3D(self.InSelFile,
                                outname,
                                volname,
                                self.AngularSampling,
                                self.NumberOfIterations,
                                self.SymmetryFile,
                                self.ExtraParamsMLrefine3D)


    # Either for seeds generation or for ML3D-classification
    def execute_MLrefine3D(self,inselfile,outname,volname,sampling,iter,symfile,extraparam):
        import os
        import launch_parallel_job

        print '*********************************************************************'
        print '*  Executing MLrefine3D program :' 
        params= ' -i '    + str(inselfile) + \
                ' -o '    + str(outname) + \
                ' -vol '  + str(volname) + \
                ' -iter ' + str(iter) + \
                ' -sam '  + str(sampling) + \
                
        if not symfile=="":
            params+= ' -sym '+str(symfile)
        params+=' '+extraparam

        launch_parallel_job.launch_job(self.DoParallel,
                                       "xmipp_MLrefine3D",
                                       "xmipp_mpi_MLrefine3D",
                                       params,
                                       self.ParallelScript,
                                       self.LaunchJobCommand,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       self.WorkingDir,
                                       False)

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    ML3D=ML3D_class(InSelFile,
                    WorkingDir,
                    DoDeleteWorkingDir,
                    ProjectDir,
                    LogDir,
                    DoCorrectGreyScale,
                    ProjMatchSampling,
                    WbpThreshold,
                    DoGenerateSeeds,
                    NumberOfReferences,
                    InitialReference,
                    DoML3DClassification,
                    AngularSampling,
                    NumberOfIterations,
                    SymmetryFile,
                    ExtraParamsMLrefine3D,
                    SeedsSelfile,
                    LaunchJobCommand,
                    DoParallel,
                    MyNumberOfCPUs,
                    MyMachineFile,
                    ParallelScript)

    # close 
    ML3D.close()

