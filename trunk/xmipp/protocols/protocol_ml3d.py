#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict

class ProtML3D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml2d.name, scriptname, project)
        self.Import = 'from protocol_ml3d import *'

        # This is not a restart
        if (RestartIter < 1):

            if (self.DoMlf and self.DoCorrectAmplitudes):
                # Copy CTFdat to the workingdir as well
                shutil.copy(InCtfDatFile,self.WorkingDir+'/my.ctfdat')

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
            md = MetaData(self.SeedsSelfile)
            md_out = MetaData()
            new_fn = FileName()
            for id in md:
                fn = str(md.getValue(MDL_IMAGE))
                new_fn.compose('initial_seed', i, 'vol')
                shutil.copy(fn, new_fn)
                md_out.addObject()
                md_out.setValue(MDL_IMAGE, new_fn)
            md_out.write('ml3d_seeds.sel')

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
        launchJob("xmipp_header_extract",
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
                    ' -compute_neighbors -angular_distance -1'

        launchJob("xmipp_angular_project_library",
                              parameters,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpi,
                              self.NumberOfThreads,
                              self.SystemFlavour)


        print '*********************************************************************'
        print '* Perform projection matching'
        parameters= ' -i '              + dirname + docfile + \
                    ' -o '              + dirname + basename + \
                    ' -ref '            + dirname + refname

        launchJob('xmipp_angular_projection_matching',
                              parameters,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpi,
                              self.NumberOfThreads,
                              self.SystemFlavour)

        print '*********************************************************************'
        print '* Make the class averages '
        parameters =  ' -i '      + dirname + basename + '.doc'  + \
                      ' -lib '    + dirname + refname  + '_angles.doc' + \
                      ' -o '      + dirname + basename 

        launchJob('xmipp_angular_class_average',
                              parameters,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpi,
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
           
        launchJob("xmipp_reconstruct_fourier",
                              parameters,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpi,
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
        launchJob("xmipp_fourier_filter",
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
        launchJob("xmipp_selfile_split",
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
           
        launchJob(program,
                              params,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpi,
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

        launchJob(program,
                              params,
                              self.log,
                              self.DoParallel,
                              self.NumberOfMpi,
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
                    NumberOfMpi,
                    SystemFlavour)
    
    # close 
    ML3D.close()

