#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based ML3D/MLF3D classification, according to:
# {please cite} for ML3D:  Scheres et al. (2007) Nature Methods, 4, 27-29
# {please cite} for MLF3D: Scheres et al. (2007) Structure, 15, 1167-1177
#
# Example use:
# ./xmipp_protocol_ml3d.py
#
# Author: Sjors Scheres, December 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {file} Selfile with the input images:
InSelFile="all_images.sel"
# {file} Initial 3D reference map:
InitialReference="my_ref.vol"
# Working subdirectory:
WorkingDir="ML3D/test1"
# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir="/home/scheres"
# {expert} Directory name for logfiles:
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} MLF-specific parameters
#------------------------------------------------------------------------------------------------
# Perform MLF3D instead of ML3D classification?
DoMlf=False
# {file} CTFdat file with the input images:
""" The names of both the images and the ctf-parameter files should be with absolute paths.
"""
InCtfDatFile="all_images.ctfdat"
# Is the initial 3D reference map CTF-amplitude corrected?
""" If coming from programs other than xmipp_mlf_refine3d this is usually not the case. If you will perform a grey-scale correction, this parameter becomes irrelevant as the output maps never have the CTF-amplitudes corrected.
"""
InitialMapIsAmplitudeCorrected=True
# {expert} Are the seeds  CTF-amplitude corrected?
""" This option is only relevant if you provide your own seeds! If the seeds are generated automatically, this parameter becomes irrelevant as they will always be amplitude-corrected
"""
SeedsAreAmplitudeCorrected=False
# High-resolution limit (in Angstroms)
""" No frequencies higher than this limit will be taken into account. If zero is given, no limit is imposed
"""
HighResLimit=20
#------------------------------------------------------------------------------------------------
# {section} Correct absolute grey scale of initial reference
#------------------------------------------------------------------------------------------------
# Correct the absolute grey-scale of the initial reference map?
""" The probabilities are based on squared differences, so that the absolute grey scale is important.
"""
DoCorrectGreyScale=False
# {expert} Angular sampling for a quick projection matching to obtain right grey scale
""" As the resolution of the intial reference should be low, this sampling can be relatively crude, e.g. 15
"""
ProjMatchSampling=15
# {expert} Threshold for a quick weighted back-projection to obtain right grey scale
""" As the resolution of the intial reference should be low, this value can be relatively high, e.g. 0.02
"""
WbpThreshold=0.02
#------------------------------------------------------------------------------------------------
# {section} Low-pass filter initial reference?
#------------------------------------------------------------------------------------------------
# Low-pass filter the initial reference?
""" It is highly recommended to low-pass filter your initial reference volume as much as you can.
"""
DoLowPassFilterReference=True
# Resolution of the low-pass filter (in Angstroms):
LowPassFilter=50
# Pixel size (in Angstroms):
PixelSize=5.6
#------------------------------------------------------------------------------------------------
# {section} Seed generation
#------------------------------------------------------------------------------------------------
# Generate unbiased seeds from a single initial reference map?
DoGenerateSeeds=True
# Number of seeds to be generated (and later on used in ml3d-classification):
NumberOfReferences=3
# {expert} Alternative 1: don't classify anything, just refine my initial map
""" No multiple seeds will be generated, only the initial 3D reference
    map will be refined. One may still perform grey-scale correction
    and/or low-pass filtering of this map.
    Note that the option to generate unbiased seeds should be set to False
    
"""
DoJustRefine=False
# {expert} {file} Alternative 2: provide a selfile with user-defined seeds:
""" Automated (unbiased!) seed generation is highly recommended...
    But you may use this option to provide your own seeds.
    The seeds should already be on the correct absolute greyscale!
    Note that the options to generate unbiased seeds and to just
    refine the initial map should be set to False
"""
SeedsSelfile=""
#------------------------------------------------------------------------------------------------
# {section} ML3D classification
#------------------------------------------------------------------------------------------------
# Perform ML(F)3D classification run?
DoML3DClassification=True
# Angular sampling for ML(F)3D classification:
""" Fine samplings take huge amounts of CPU and memory.
    Therefore, in general, dont use samplings finer than 10 degrees.
"""
AngularSampling=10
# Number of ML(F)3D iterations to perform:
NumberOfIterations=25
# {file} In the case of symmetry supply the symmetry description file:
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetrize
    for a description of the symmetry file format
    dont give anything, if no symmetry is present
"""
SymmetryFile=""
# {expert} Restart after iteration:
""" For previous runs that stopped before convergence,
    resume the calculations after the completely finished iteration,
    i.e. including all 3D reconstructions.
    Note that all flags about grey-scale correction, filtering and
    seed generation will be ignored if a value larger than 0 is given,
    since this option only concerns the ML3D classification part
"""
RestartIter=0
# {expert} Additional xmipp_ml(f)_refine3d parameters:
""" For a complete description see the manual pages:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/MLrefine3D
"""
ExtraParamsMLrefine3D=""
#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Use multiple processors in parallel?
DoParallel=True
# Number of processors to use:
MyNumberOfCPUs=5
# {file} A list of all available CPUs (the MPI-machinefile):
""" Depending on your system, your standard script to launch MPI-jobs may require this
    if your queueing system using an environment variable, give it here (with the leading $, e.g. $PBS_NODEFILE
"""
MyMachineFile="/home2/bioinfo/scheres/machines.dat"
# {expert} Control file
""" This is an ugly solution to have improved killing control over the mpi jobs.
    The script will create this file, and any mpi-job will be killed as soon as this file doesn't exist anymore. This is required with certain queueing systems.
"""
MyControlFile=""
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_ml3d.py"
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class ML3D_class:

    #init variables
    def __init__(self,
                 InSelFile,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 ProjectDir,
                 LogDir,
                 DoMlf,
                 InCtfDatFile,
                 InitialMapIsAmplitudeCorrected,
                 SeedsAreAmplitudeCorrected,
                 HighResLimit,
                 DoCorrectGreyScale,
                 ProjMatchSampling,
                 WbpThreshold,
                 DoLowPassFilterReference,
                 LowPassFilter,
                 PixelSize,
                 DoGenerateSeeds,
                 NumberOfReferences,
                 DoJustRefine,
                 InitialReference,
                 DoML3DClassification,
                 AngularSampling,
                 NumberOfIterations,
                 Symmetry,
                 RestartIter,
                 ExtraParamsMLrefine3D,
                 SeedsSelfile,
                 DoParallel,
                 MyNumberOfCPUs,
                 MyMachineFile,
                 MyControlFile):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,selfile

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.DoMlf=DoMlf
        self.RefForSeedsIsAmplitudeCorrected=InitialMapIsAmplitudeCorrected
        self.SeedsAreAmplitudeCorrected=SeedsAreAmplitudeCorrected
        self.HighResLimit=HighResLimit
        self.NumberOfReferences=NumberOfReferences
        self.DoJustRefine=DoJustRefine
        self.DoGenerateSeeds=DoGenerateSeeds
        self.InitialReference=os.path.abspath(InitialReference)
        if (len(SeedsSelfile)>0):
            self.SeedsSelfile=os.path.abspath(SeedsSelfile)
        else:
            self.SeedsSelfile=SeedsSelfile
        self.AngularSampling=AngularSampling
        self.LowPassFilter=LowPassFilter
        self.PixelSize=PixelSize
        self.NumberOfIterations=NumberOfIterations
        self.Symmetry=Symmetry
        self.ExtraParamsMLrefine3D=ExtraParamsMLrefine3D
        self.ProjMatchSampling=ProjMatchSampling
        self.WbpThreshold=WbpThreshold
        self.DoParallel=DoParallel
        self.MyNumberOfCPUs=MyNumberOfCPUs
        if (MyMachineFile[0]=="$"):
            self.MyMachineFile=MyMachineFile
        else:
            self.MyMachineFile=os.path.abspath(MyMachineFile)
        if (MyControlFile==""):
            self.DoControl=False
        else:
            self.DoControl=True

        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)

        # This is not a restart
        if (RestartIter < 1):
            # Delete working directory if it exists, make a new one
            if (DoDeleteWorkingDir): 
                if os.path.exists(self.WorkingDir):
                    shutil.rmtree(self.WorkingDir)
            if not os.path.exists(self.WorkingDir):
                os.makedirs(self.WorkingDir)

            # Create a CONTROL file for improved killing control
            if (self.DoControl):
                self.MyControlFile=os.path.abspath(self.WorkingDir+
                                                   '/'+MyControlFile)
                FILE = open(self.MyControlFile,"w")
                FILE.write("Delete this file to kill all current processes\n")
                FILE.close()

            # Create a selfile with absolute pathname in the WorkingDir
            mysel=selfile.selfile()
            mysel.read(InSelFile)
            newsel=mysel.make_abspath()
            self.InSelFile=os.path.abspath(self.WorkingDir+'/'+InSelFile)
            newsel.write(self.InSelFile)

            if (self.DoMlf):
                # Copy CTFdat to the workingdir as well
                shutil.copy(InCtfDatFile,self.WorkingDir+'/my.ctfdat')

            # Backup script
            log.make_backup_of_script_file(sys.argv[0],
                                           os.path.abspath(self.WorkingDir))
    
            # Execute in the working directory
            os.chdir(self.WorkingDir)

            self.copy_initial_refs()
        
            if DoCorrectGreyScale:
                self.correct_greyscale()

            if DoLowPassFilterReference:
                self.filter_reference()

            if DoGenerateSeeds:
                self.generate_seeds()

            if DoML3DClassification:
                self.execute_ML3D_classification()

        # Restarting a previous run...
        else:
            # Create a CONTROL file for improved killing control
            if (self.DoControl):
                self.MyControlFile=os.path.abspath(self.WorkingDir+
                                                   '/'+MyControlFile)
                FILE = open(self.MyControlFile,"w")
                FILE.write("Delete this file to kill current processes\n")
                FILE.close()

            # Execute protocol in the working directory
            os.chdir(self.WorkingDir)
            self.restart_MLrefine3D(RestartIter)
        
        # Return to parent dir
        os.chdir(os.pardir)

        # Delete the CONTROL file for improved killing control
        if (self.DoControl):
            os.remove(self.MyControlFile)

    # Make a copy of the initial references:
    def copy_initial_refs(self):
        import os,shutil
        if os.path.exists(self.InitialReference):
            shutil.copy(self.InitialReference,'initial_reference.vol')
        if (self.DoGenerateSeeds==False and DoJustRefine==False):
            fh=open(self.SeedsSelfile,'r')
            lines=fh.readlines()
            fh.close()
            i = 0
            newlines=[]
            for line in lines:
                i = i + 1
                oname='initial_seed'+str(i).zfill(5)+'.vol'
                words=line.split()
                if words[0][0]=="/":
                    shutil.copy(words[0],oname)
                else:
                    dirname=os.path.dirname(self.SeedsSelfile)
                    iname=dirname+'/'+words[0]
                    shutil.copy(iname,oname)
                newlines.append(oname+' 1\n')
            fh=open('ml3d_seeds.sel','w')
            fh.writelines(newlines)
            fh.close()

    # Crude correction of grey-scale, by performing a single iteration of projection matching and wbp
    def correct_greyscale(self):
        import os
        import launch_parallel_job
        print '*********************************************************************'
        print '*  Correcting absolute grey scale of initial reference:'

        dirname='CorrectGreyscale/'
        basename='corrected_reference'
        docfile='original_angles.doc'
        refname='ref'
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        # Grey-scale correction always leads to an amplitude uncorrected map
        self.RefForSeedsIsAmplitudeCorrected=False
            
        print '*********************************************************************'
        print '* Create initial docfile'
        command='xmipp_header_extract -i ' + str(self.InSelFile) + \
            ' -o ' + dirname + docfile
        self.log.info(command)
        os.system(command)

        print '*********************************************************************'
        print '* Create projection library'
        parameters= ' -i '                   + str(self.InitialReference) + \
                    ' -experimental_images ' + dirname + docfile + \
                    ' -o '                   + dirname + refname + \
                    ' -sampling_rate '       + str(self.ProjMatchSampling)  + \
                    ' -sym '                 + self.Symmetry + 'h' + \
                    ' -compute_neighbors -angular_distance -1 '

        launch_parallel_job.launch_job(self.DoParallel,
                                       "xmipp_angular_project_library",
                                       "xmipp_mpi_angular_project_library",
                                       parameters,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       False)


        print '*********************************************************************'
        print '* Perform projection matching'
        parameters= ' -i '              + dirname + docfile + \
                    ' -o '              + dirname + basename + \
                    ' -ref '            + dirname + refname

        launch_parallel_job.launch_job(self.DoParallel,
                                       'xmipp_angular_projection_matching',
                                       'xmipp_mpi_angular_projection_matching',
                                       parameters,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       False)

        print '*********************************************************************'
        print '* Make the class averages '
        parameters =  ' -i '      + dirname + basename + '.doc'  + \
                      ' -lib '    + dirname + refname  + '_angles.doc' + \
                      ' -o '      + dirname + basename 

        launch_parallel_job.launch_job(self.DoParallel,
                                       'xmipp_angular_class_average',
                                       'xmipp_mpi_angular_class_average',
                                       parameters,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       False)


        print '*********************************************************************'
        print '* Perform WBP reconstruction '
        iname=dirname+basename+'_classes.sel'
        outname=basename+'.vol'
        params= ' -i '    + str(iname) + \
                ' -o '    + str(outname) + \
                ' -threshold '+str(self.WbpThreshold) + \
                ' -sym '+ self.Symmetry + \
                '  -use_each_image -weight '
        if (self.DoControl):
            params+=' -control ' + self.MyControlFile
           
        launch_parallel_job.launch_job(self.DoParallel,
                                       "xmipp_reconstruct_wbp",
                                       "xmipp_mpi_reconstruct_wbp",
                                       params,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       False)

        return outname

    # Low-pass filter
    def filter_reference(self):
        import os
        print '*********************************************************************'
        print '*  Low-pass filtering of the initial reference:'


        if os.path.exists('corrected_reference.vol'):
            reference='corrected_reference.vol'
        else:
            reference=self.InitialReference
        command='xmipp_fourier_filter -o filtered_reference.vol' + \
                 ' -i ' + reference  + \
                 ' -sampling ' + str(self.PixelSize) + \
                 ' -low_pass ' + str(self.LowPassFilter)

        print '* ',command
        self.log.info(command)
        os.system(command)

    # Splits selfile and performs a single cycle of ML3D-classification for each subset
    def generate_seeds(self):
        import os
        import selfile
        print '*********************************************************************'
        print '*  Generating seeds:' 
        newsel=selfile.selfile()

        # Split selfiles
        command='xmipp_selfile_split -o seeds_split'+ \
                 ' -i ' + str(self.InSelFile) + \
                 ' -n ' + str(self.NumberOfReferences) +' \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

        if os.path.exists('filtered_reference.vol'):
            reference='filtered_reference.vol'
        elif os.path.exists('corrected_reference.vol'):
            reference='corrected_reference.vol'
        else:
            reference='initial_reference.vol'

        # Launch MLrefine3D with output to subdirectories
        for i in range(self.NumberOfReferences):
            inselfile='seeds_split_'+str(i+1)+'.sel'
            dirname='GenerateSeed_'+str(i+1)+'/'
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            outname=dirname+'seeds_split_'+str(i+1)
            self.execute_MLrefine3D(inselfile,
                                    outname,
                                    reference,
                                    self.AngularSampling,
                                    1,
                                    self.Symmetry,
                                    self.RefForSeedsIsAmplitudeCorrected,
                                    self.ExtraParamsMLrefine3D)
            newsel.insert(outname+'_it00001.vol','1')
        newsel.write('ml3d_seeds.sel')

        # Seed generation with MLF always does amplitude correction
        self.SeedsAreAmplitudeCorrected=True

    # Perform the actual classification
    def execute_ML3D_classification(self):
        import os
        import launch_parallel_job
 
        dirname='RunML3D/'
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        if (self.DoJustRefine):
            if os.path.exists('filtered_reference.vol'):
                reference='filtered_reference.vol'
            elif os.path.exists('corrected_reference.vol'):
                reference='corrected_reference.vol'
            else:
                reference='initial_reference.vol'
        else:
            reference='ml3d_seeds.sel'

        self.execute_MLrefine3D(self.InSelFile,
                                dirname+'ml3d',
                                reference,
                                self.AngularSampling,
                                self.NumberOfIterations,
                                self.Symmetry,
                                self.SeedsAreAmplitudeCorrected,
                                self.ExtraParamsMLrefine3D)


    # Either for seeds generation or for ML3D-classification
    def execute_MLrefine3D(self,inselfile,outname,
                           volname,sampling,iter,symmetry,amplitude_corrected,extraparam):
        import os
        import launch_parallel_job

        print '*********************************************************************'
        print '*  Executing ml_refine3d program :' 
        params= ' -i '    + str(inselfile) + \
                ' -o '    + str(outname) + \
                ' -vol '  + str(volname) + \
                ' -iter ' + str(iter) + \
                ' -sym '  + symmetry +\
                ' -ang '  + str(sampling)
        params+=' '+extraparam
        if (self.DoControl):
            params+=' -control ' + self.MyControlFile
        if (self.DoMlf):
            params+= ' -ctfdat my.ctfdat'
            if (not amplitude_corrected):
                params+= ' -ctf_affected_refs'
            if (self.HighResLimit > 0):
                params += ' -high ' + str(self.HighResLimit)

        if (self.DoMlf):
            program="xmipp_mlf_refine3d"
            mpiprogram="xmipp_mpi_mlf_refine3d"
        else:
            program="xmipp_ml_refine3d"
            mpiprogram="xmipp_mpi_ml_refine3d"
           
        launch_parallel_job.launch_job(self.DoParallel,
                                       program,
                                       mpiprogram,
                                       params,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
                                       False)

    # Either for seeds generation or for ML3D-classification
    def restart_MLrefine3D(self,iter):
        import os
        import launch_parallel_job

        print '*********************************************************************'
        print '*  Restarting ml(f)_refine3d program :' 
        params= ' -restart RunML3D/ml3d_it' + str(iter).zfill(5) + '.log'

        if (self.DoMlf):
            program="xmipp_mlf_refine3d"
            mpiprogram="xmipp_mpi_mlf_refine3d"
        else:
            program="xmipp_ml_refine3d"
            mpiprogram="xmipp_mpi_ml_refine3d"

        launch_parallel_job.launch_job(self.DoParallel,
                                       program,
                                       mpiprogram,
                                       params,
                                       self.log,
                                       self.MyNumberOfCPUs,
                                       self.MyMachineFile,
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
                    DoMlf,
                    InCtfDatFile,
                    InitialMapIsAmplitudeCorrected,
                    SeedsAreAmplitudeCorrected,
                    HighResLimit,
                    DoCorrectGreyScale,
                    ProjMatchSampling,
                    WbpThreshold,
                    DoLowPassFilterReference,
                    LowPassFilter,
                    PixelSize,
                    DoGenerateSeeds,
                    NumberOfReferences,
                    DoJustRefine,
                    InitialReference,
                    DoML3DClassification,
                    AngularSampling,
                    NumberOfIterations,
                    Symmetry,
                    RestartIter,
                    ExtraParamsMLrefine3D,
                    SeedsSelfile,
                    DoParallel,
                    MyNumberOfCPUs,
                    MyMachineFile,
                    MyControlFile)
    
    # close 
    ML3D.close()

