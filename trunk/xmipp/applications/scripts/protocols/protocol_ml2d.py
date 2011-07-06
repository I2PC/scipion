#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles, according to:
# {please cite} for ML2D:  Scheres et al. (2005) J.Mol.Biol 348, 139-149 
# {please cite} for MLF2D: Scheres et al. (2007) Structure 15, 1167-1177
#
# Example use:
# ./xmipp_protocol_ml2d.py
#
#  Author:  Sjors Scheres, January 2008
# Updated:  J. M. de la Rosa Trevin July 2011
#
# {begin_of_header} 
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------

# Working subdirectory: 
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
RunName ='run_001'

# {file} Selfile with the input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
InSelFile="all_images.sel"

# Perform MLF2D instead of ML2D classification?
DoMlf=False

#------------------------------------------------------------------------------------------------
# {section} {condition}(DoMlf=True) MLF-specific parameters
#------------------------------------------------------------------------------------------------

# Use CTF-amplitude correction inside MLF?
""" If set to true, provide the ctfdat file in the field below. If set to false, one can ignore the ctfdat field, but has to provide the image pixel size in Angstrom"""
DoCorrectAmplitudes=True

# {file} CTFdat file with the input images:
""" The names of both the images and the ctf-parameter files should be with absolute paths.
"""
InCtfDatFile="all_images.ctfdat"

# Image pixel size (in Angstroms)
PixelSize=4.2

# Are the images CTF phase flipped?
""" You can run MLF with or without having phase flipped the images.
"""
ImagesArePhaseFlipped=True

# High-resolution limit (in Angstroms)
""" No frequencies higher than this limit will be taken into account. If zero is given, no limit is imposed
"""
HighResLimit=20

#------------------------------------------------------------------------------------------------
# {section}{expert} ml(f)_align2d parameters
#------------------------------------------------------------------------------------------------

# Perform 2D maximum-likelihood refinement?
DoML2D=True

# Number of references (or classes) to be used:
NumberOfReferences=3

# Also include mirror transformation in the alignment?
"""  Including the mirror transformation is useful if your particles have a handedness and may fall either face-up or face-down on the grid.
 Note that when you want to use this ML2D run for later RCT reconstruction, you can NOT include the mirror transformation here.
"""
DoMirror=False

# Use the fast version of this algorithm?
""" See Scheres et al., Bioinformatics, 21 (Suppl. 2), ii243-ii244:
    http://dx.doi.org/10.1093/bioinformatics/bti1140
"""
DoFast=True

# Refine the normalization parameters for each image?
""" This variant of the algorithm deals with normalization errors. For more info see (and please cite) Scheres et. al. (2009) J. Struc. Biol., in press
"""
DoNorm=False

# {expert} Restart after iteration:
""" For previous runs that stopped before convergence, resume the calculations
    after the completely finished iteration. (Use zero to start from the beginning)
"""
RestartIter=0

# {expert} Additional xmipp_ml_align2d parameters:
""" For a complete description see the ml_align2d manual page at:
    http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/MLalign2D
"""
ExtraParamsMLalign2D=""

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of (shared-memory) threads?
""" This option provides shared-memory parallelization on multi-core machines. 
    It does not require any additional software, other than xmipp
"""
NumberOfThreads=1
# Use distributed-memory parallelization (MPI)?
""" This option provides distributed-memory parallelization on multi-node machines. 
    It requires the installation of some MPI flavour, possibly together with a queueing system
"""
DoParallel=False
# Number of MPI processes to use:
NumberOfMpiProcesses=5
# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. 
"""
SystemFlavour=''

#------------------------------------------------------------------------------------------------
# {expert}{file} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_ml2d.py"
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

import os, sys, shutil
from xmipp import MetaData
from protlib_base import *

class ProtML2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(protDict.ml2d.key, scriptname, RunName, project, NumberOfMpiProcesses)

        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        # store script name
        protocolName = sys.argv[0]
        #assume restart
        restart = True
        # This is not a restart
        if (RestartIter < 1):
            # Delete working directory if it exists, make a new one
            if (DoDeleteWorkingDir and DoML2D): 
                if (self.WorkingDir==""):
                   raise RuntimeError,"No working directory given"
                if os.path.exists(WorkingDir):
                    shutil.rmtree(WorkingDir)
            if not os.path.exists(WorkingDir):
                os.makedirs(WorkingDir)


            # Create a selfile with absolute pathname in the WorkingDir
#            mysel = MetaData(InSelFile);
#            mysel.makeAbsPath();
#            InSelFile = os.path.abspath(os.path.join(WorkingDir, InSelFile))
#            mysel.write(InSelFile)

            if (DoMlf and DoCorrectAmplitudes):
                # Copy CTFdat to the workingdir as well
                shutil.copy(InCtfDatFile, os.path.join(WorkingDir, 'my.ctfdat'))

            # Backup script
            log.make_backup_of_script_file(protocolName, os.path.abspath(WorkingDir))
            restart = False

        # Store current directory before moving
        currentDir = os.getcwd()            
        # Execute protocol in the working directory
        os.chdir(WorkingDir)
        self.run(restart)
        # Return to parent dir
        os.chdir(currentDir)

    def run(self, restart=False):
        import os
        import launch_job
        print '*********************************************************************'
        progId = "ml"
        if (DoMlf):
            progId += "f"  
        
        program = "xmipp_%s_align2d" % progId
        action = "Executing"
        if (restart):
            action = "Restarting"
            
        print '*  %s %s program :' % (action, program)
        
        if (restart):
            params= ' --restart ' + utils_xmipp.composeFileName('ml2d_it', RestartIter,'log')
        else: 
            params = ' -i %s --oroot %s2d' % (InSelFile, progId)
            # Number of references will be ignored if -ref is passed as expert option
            if ExtraParamsMLalign2D.find("--ref") == -1:
                params += ' --nref %i' % NumberOfReferences
            params += ' ' + ExtraParamsMLalign2D
            if (DoFast and not DoMlf):
                params += ' --fast'
            if (DoNorm):
                params += ' --norm'
            if (DoMirror):
                params += ' --mirror'
            if (NumberOfThreads > 1  and not DoMlf):
                params += ' --thr %i' % NumberOfThreads
            if (DoMlf):
                if (DoCorrectAmplitudes):
                    params += ' --ctfdat my.ctfdat'
                else:
                    params += ' --no_ctf -pixel_size %f' + PixelSize
                if (not self.ImagesArePhaseFlipped):
                    params += ' --not_phase_flipped'
                if (self.HighResLimit > 0):
                    params += ' --high %f' + HighResLimit
           
        launchJob(program, params, self.log, DoParallel, NumberOfMpiProcesses, NumberOfThreads, SystemFlavour)

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'


if __name__ == '__main__':
    protocolMain(ProtML2D)

