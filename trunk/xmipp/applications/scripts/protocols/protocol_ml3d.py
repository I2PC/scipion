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
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
InSelFile='all_images.sel'
# {file} Initial 3D reference map:
""" Spider-format 3D density map with the same dimensions as your particles. This file may be in any directory.
"""
InitialReference='my_ref.vol'
# Working subdirectory:
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='ML3D/test1'
# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/home/scheres/proteinA/dataset1'
# {expert} Directory name for logfiles:
LogDir='Logs'
#------------------------------------------------------------------------------------------------
# {section} MLF-specific parameters
#------------------------------------------------------------------------------------------------
# Perform MLF3D instead of ML3D classification?
DoMlf=False
# Use CTF-amplitude correction inside MLF?
""" If set to true, provide the ctfdat file in the field below. If set to false, the ctfdat field can be ignored"""
DoCorrectAmplitudes=True
# {file} CTFdat file with the input images:
""" The names of both the images and the ctf-parameter files should be with absolute paths. If you want to use this, make sure also the images in the input selfile (see above) are with absolute paths.
"""
InCtfDatFile='all_images.ctfdat'
# High-resolution limit (in Angstroms)
""" No frequencies higher than this limit will be taken into account. If zero is given, no limit is imposed
"""
HighResLimit=20
# Are the images CTF phase flipped?
""" You can run MLF with or without having phase flipped the images.
"""
ImagesArePhaseFlipped=True
# Is the initial 3D reference map CTF-amplitude corrected?
""" If coming from programs other than xmipp_mlf_refine3d this is usually not the case. If you will perform a grey-scale correction, this parameter becomes irrelevant as the output maps never have the CTF-amplitudes corrected.
"""
InitialMapIsAmplitudeCorrected=False
# {expert} Are the seeds  CTF-amplitude corrected?
""" This option is only relevant if you provide your own seeds! If the seeds are generated automatically, this parameter becomes irrelevant as they will always be amplitude-corrected
"""
SeedsAreAmplitudeCorrected=False
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
SeedsSelfile=''
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
# Point group symmetry:
""" See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for a description of the point group symmetries
    Give c1 if no symmetry is present
"""
Symmetry='c1'
# Refine the normalization parameters for each image?
""" This variant of the algorithm deals with normalization errors.
    For more info see (and please cite) Scheres et. al. (2009) J. Struc. Biol., in press
"""
DoNorm=False
# {expert} Use Fourier-interpolation instead of wlsART?
""" The Fourier-interpolation reconstruction method is much faster than wlsART 
    and may give similar results. It however is not guaranteed to optimize the 
    likelihood function. This is an experimental feature. One may limit the 
    maximum resolution of the fourier-interpolation using "-max_resolution 0.3"
    (to 0.3 digital frequency). Use the extra parameter entry below for that. 
"""
DoFourier=False
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
ExtraParamsMLrefine3D=''
#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads=1
# Use distributed-memory parallelization (MPI)?
""" This option provides parallelization on clusters with distributed memory architecture.
    It requires mpi to be installed.
"""
DoParallel=True
# Number of MPI processes to use:
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
NumberOfMpiProcesses=5
# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. 
"""
SystemFlavour=''

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_ml3d.py'
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
                 DoCorrectAmplitudes,
                 InCtfDatFile,
                 ImagesArePhaseFlipped,
                 InitialMapIsAmplitudeCorrected,
                 SeedsAreAmplitudeCorrected,
                 HighResLimit,
                 DoCorrectGreyScale,
                 ProjMatchSampling,
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
                 DoNorm,
                 DoFourier,
                 RestartIter,
                 ExtraParamsMLrefine3D,
                 SeedsSelfile,
                 NumberOfThreads,
                 DoParallel,
                 NumberOfMpiProcesses,
                 SystemFlavour):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,selfile

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.DoMlf=DoMlf
        self.DoCorrectAmplitudes=DoCorrectAmplitudes
        self.RefForSeedsIsAmplitudeCorrected=InitialMapIsAmplitudeCorrected
        self.ImagesArePhaseFlipped=ImagesArePhaseFlipped
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
        self.DoNorm=DoNorm
        self.DoFourier=DoFourier
        self.ExtraParamsMLrefine3D=ExtraParamsMLrefine3D
        self.ProjMatchSampling=ProjMatchSampling
        self.NumberOfThreads=NumberOfThreads
        self.DoParallel=DoParallel
        self.NumberOfMpiProcesses=NumberOfMpiProcesses
        self.SystemFlavour=SystemFlavour

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

            # Create a selfile with absolute pathname in the WorkingDir
            mysel=selfile.selfile()
            mysel.read(InSelFile)
            newsel=mysel.make_abspath()
            self.InSelFile=os.path.abspath(self.WorkingDir+'/'+InSelFile)
            newsel.write(self.InSelFile)

            if (self.DoMlf and self.DoCorrectAmplitudes):
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
            # Execute protocol in the working directory
            os.chdir(self.WorkingDir)
            self.restart_MLrefine3D(RestartIter)
        
        # Return to parent dir
        os.chdir(os.pardir)


    # Make a copy of the initial references:
    def copy_initial_refs(self):
        import os,shutil
        import utils_xmipp
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
                oname=utils_xmipp.composeFileName('initial_seed',i,'vol')
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

    # Crude correction of grey-scale, by performing a single iteration of 
    # projection matching and fourier reconstruction
    def correct_greyscale(self):
        import os
        import launch_job
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
        params= ' -i ' + str(self.InSelFile) + \
                ' -o ' + dirname + docfile
        launch_job.launch_job("xmipp_header_extract",
                              params,
                              self.log,
                              False,1,1,'')

        print '*********************************************************************'
        print '* Create projection library'
        parameters= ' -i '                   + str(self.InitialReference) + \
                    ' -experimental_images ' + dirname + docfile + \
                    ' -o '                   + dirname + refname + \
                    ' -sampling_rate '       + str(self.ProjMatchSampling)  + \
                    ' -sym '                 + self.Symmetry + 'h' + \
                    ' -compute_neighbors -angular_distance -1 -shears'

        launch_job.launch_job("xmipp_angular_project_library",
                              parameters,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpiProcesses,
                              self.NumberOfThreads,
                              self.SystemFlavour)


        print '*********************************************************************'
        print '* Perform projection matching'
        parameters= ' -i '              + dirname + docfile + \
                    ' -o '              + dirname + basename + \
                    ' -ref '            + dirname + refname

        launch_job.launch_job('xmipp_angular_projection_matching',
                              parameters,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpiProcesses,
                              self.NumberOfThreads,
                              self.SystemFlavour)

        print '*********************************************************************'
        print '* Make the class averages '
        parameters =  ' -i '      + dirname + basename + '.doc'  + \
                      ' -lib '    + dirname + refname  + '_angles.doc' + \
                      ' -o '      + dirname + basename 

        launch_job.launch_job('xmipp_angular_class_average',
                              parameters,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpiProcesses,
                              self.NumberOfThreads,
                              self.SystemFlavour)


        print '*********************************************************************'
        print '* Perform Fourier-interpolation reconstruction '
        iname=dirname+basename+'_classes.sel'
        outname=basename+'.vol'
        parameters= ' -i '    + str(iname) + \
                    ' -o '    + str(outname) + \
                    ' -sym '+ self.Symmetry + \
                    '  -weight '
        if (self.NumberOfThreads>1):
            parameters += ' -thr ' + str(self.NumberOfThreads)
           
        launch_job.launch_job("xmipp_reconstruct_fourier",
                              parameters,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpiProcesses,
                              self.NumberOfThreads,
                              self.SystemFlavour)
 
        return outname

    # Low-pass filter
    def filter_reference(self):
        import os
        import launch_job
        print '*********************************************************************'
        print '*  Low-pass filtering of the initial reference:'


        if os.path.exists('corrected_reference.vol'):
            reference='corrected_reference.vol'
        else:
            reference=self.InitialReference
        params=' -o filtered_reference.vol' + \
               ' -i ' + reference  + \
               ' -sampling ' + str(self.PixelSize) + \
               ' -low_pass ' + str(self.LowPassFilter)
        launch_job.launch_job("xmipp_fourier_filter",
                              params,
                              self.log,
                              False,1,1,'')

    # Splits selfile and performs a single cycle of ML3D-classification for each subset
    def generate_seeds(self):
        import os
        import launch_job
        import selfile, utils_xmipp
        print '*********************************************************************'
        print '*  Generating seeds:' 
        newsel=selfile.selfile()

        # Split selfiles
        params=' -o seeds_split'+ \
               ' -i ' + str(self.InSelFile) + \
               ' -n ' + str(self.NumberOfReferences) +' \n'
        launch_job.launch_job("xmipp_selfile_split",
                              params,
                              self.log,
                              False,1,1,'')

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
                                    self.ImagesArePhaseFlipped,
                                    self.RefForSeedsIsAmplitudeCorrected,
                                    self.ExtraParamsMLrefine3D)
            seedname=utils_xmipp.composeFileName(outname+'_it',1,'vol')
            newsel.insert(seedname,'1')
        newsel.write('ml3d_seeds.sel')

        # Seed generation with MLF always does amplitude correction
        self.SeedsAreAmplitudeCorrected=True

    # Perform the actual classification
    def execute_ML3D_classification(self):
        import os
        import launch_job
 
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
                                self.ImagesArePhaseFlipped,
                                self.SeedsAreAmplitudeCorrected,
                                self.ExtraParamsMLrefine3D)


    # Either for seeds generation or for ML3D-classification
    def execute_MLrefine3D(self,inselfile,outname,
                           volname,sampling,iter,symmetry,phase_flipped,amplitude_corrected,extraparam):
        import os
        import launch_job

        print '*********************************************************************'
        print '*  Executing ml_refine3d program :' 
        params= ' -i '    + str(inselfile) + \
                ' -o '    + str(outname) + \
                ' -vol '  + str(volname) + \
                ' -iter ' + str(iter) + \
                ' -sym '  + symmetry +\
                ' -ang '  + str(sampling)
        params+=' '+extraparam
        if (self.NumberOfThreads > 1):
            params+=' -thr ' + str(self.NumberOfThreads)
        if (self.DoNorm):
            params+=' -norm '
        if (self.DoFourier):
            params+=' -fourier '
        if (self.DoMlf):
            if (self.DoCorrectAmplitudes):
                params+= ' -ctfdat my.ctfdat'
            else:
                params+= ' -no_ctf -pixel_size ' + str(self.PixelSize)
            if (not phase_flipped):
                params+= ' -not_phase_flipped'
            if (not amplitude_corrected):
                params+= ' -ctf_affected_refs'
            if (self.HighResLimit > 0):
                params += ' -high ' + str(self.HighResLimit)

        if (self.DoMlf):
            program="xmipp_mlf_refine3d"
        else:
            program="xmipp_ml_refine3d"
           
        launch_job.launch_job(program,
                              params,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpiProcesses,
                              self.NumberOfThreads,
                              self.SystemFlavour)


    # Either for seeds generation or for ML3D-classification
    def restart_MLrefine3D(self,iter):
        import os
        import launch_job, utils_xmipp

        print '*********************************************************************'
        print '*  Restarting ml(f)_refine3d program :' 
        restartname = utils_xmipp.composeFileName('RunML3D/ml3d_it',iter,'log')
        params= ' -restart ' + restartname

        if (self.DoMlf):
            program="xmipp_mlf_refine3d"
        else:
            program="xmipp_ml_refine3d"

        launch_job.launch_job(program,
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

    ML3D=ML3D_class(InSelFile,
                    WorkingDir,
                    DoDeleteWorkingDir,
                    ProjectDir,
                    LogDir,
                    DoMlf,
                    DoCorrectAmplitudes,
                    InCtfDatFile,
                    ImagesArePhaseFlipped,
                    InitialMapIsAmplitudeCorrected,
                    SeedsAreAmplitudeCorrected,
                    HighResLimit,
                    DoCorrectGreyScale,
                    ProjMatchSampling,
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
                    DoNorm,
                    DoFourier,
                    RestartIter,
                    ExtraParamsMLrefine3D,
                    SeedsSelfile,
                    NumberOfThreads,
                    DoParallel,
                    NumberOfMpiProcesses,
                    SystemFlavour)
    
    # close 
    ML3D.close()

