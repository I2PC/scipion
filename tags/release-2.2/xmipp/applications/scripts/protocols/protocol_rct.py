#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for random conical tilt reconstruction
#
# Example use:
# ./xmipp_protocol_rct.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {file} Selfile with all untilted images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
UntiltedSelFile="all_images_untilted.sel"
# {file} Selfile with all tilted images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
TiltedSelFile="all_images_tilted.sel"
# Working subdirectory:
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir="RCT/test1"
# Delete working subdirectory if it already exists?
""" Just be careful with this option...
"""
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir="/home/scheres/work/protocols/G40P"
# {expert} Directory name for logfiles:
""" All logfiles will be stored in $ProjectDir/$LogDir
"""
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Previous ML2D classification (WITHOUT INCLUDING MIRRORS!)
#------------------------------------------------------------------------------------------------
# {dir} Directory of previous ML2D-classification on the untilted images
PreviousDirML2D="ML2D/ML3ref"
# {expert} Rootname for ML2D run
""" Only provide something if the ml2d output rootname was different from the default.
    This will never be the case when using the standardized protocols 
"""
PreviousML2DRoot=""
# Which of these classes do you want to reconstruct? (Separate numbers by comma's)
SelectClasses="1,2"
#------------------------------------------------------------------------------------------------
# {section} Prepare images
#------------------------------------------------------------------------------------------------
# Prepare local copies of all images?
""" This will make local copies of all untilted images and generate corresponding selfiles. This has to be done at least once.
"""
DoImagePreparation=True
# Set untilted image headers?
""" This will re-align the untilted particles and set the RCT angles correctly
"""
DoUntiltedHeaders=True
# Set tilted image headers?
""" This will center the tilted particles and set the RCT angles correctly
"""
DoTiltedHeaders=True
# Maximum allowed shift for tilted particles (in pixels):
""" Particles that shift more will be discarded. A value larger than the image size will not discard any particle
"""
CenterMaxShift=999
# {expert} Additional parameters for the align_tilt_pairs program
"""  For example:
    -skip_stretching will skip the cosine-stretching prior to centering
    -skip_centering  will skip the entire centering, so that only the RCT angles will be set.
    -force_x_zero    will force the shift in the X direction to be zero, and will only center in the Y direction
"""
AlignTiltPairsAdditionalParams=""
#------------------------------------------------------------------------------------------------
# {section} Reconstruction for each of the classes
#------------------------------------------------------------------------------------------------
# Perform 3D-reconstructions with ART?
DoArtReconstruct=False
# Relaxation parameter for ART reconstruction:
ArtLambda=0.2
# {expert} Additional ART parameters
""" See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Art
"""
ArtAdditionalParams=""
# Perform 3D-reconstructions with WBP?
DoWbpReconstruct=True
# Threshold parameter for WBP-reconstruction:
""" Higher values give lower resolution reconstruction (but possibly less noisy)
    For RCT reconstructions, values of 0.01-0.1 could be suitable
"""
WbpThreshold=0.02
# {expert} Additional WBP parameters
""" See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Wbp
"""
WbpAdditionalParams=""
#------------------------------------------------------------------------------------------------
# {section} Low-pass filter reconstructions
#------------------------------------------------------------------------------------------------
# Low-pass filter the reconstructed volumes?
""" Filtering may be useful to remove noise, especially when few particles contribute to the reconstruction.
"""
DoLowPassFilter=True
# Resolution of the low-pass filter (in Angstroms):
LowPassFilter=50
# Pixel size (in Angstroms):
PixelSize=5.6
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This variable serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_rct.py"
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class RCT_class:

    #init variables
    def __init__(self,
                 UntiltedSelFile,
                 TiltedSelFile,
                 WorkingDir,
                 DoDeleteWorkingDir,
                 ProjectDir,
                 LogDir,
                 PreviousDirML2D,
                 SelectClasses,
                 DoImagePreparation,
                 DoUntiltedHeaders,
                 DoTiltedHeaders,
                 CenterMaxShift,
                 AlignTiltPairsAdditionalParams,
                 DoArtReconstruct,
                 ArtLambda,
                 ArtAdditionalParams,
                 DoWbpReconstruct,
                 WbpThreshold,
                 WbpAdditionalParams,
                 DoLowPassFilter,
                 LowPassFilter,
                 PixelSize):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log,selfile
        
        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.PreviousDirML2D=os.path.abspath(PreviousDirML2D)
        self.PreviousML2DRoot=PreviousML2DRoot
        self.SelectClasses=SelectClasses
        self.CenterMaxShift=CenterMaxShift
        self.AlignTiltPairsAdditionalParams=AlignTiltPairsAdditionalParams
        self.ArtLambda=ArtLambda
        self.ArtAdditionalParams=ArtAdditionalParams
        self.WbpThreshold=WbpThreshold
        self.WbpAdditionalParams=WbpAdditionalParams
        self.LowPassFilter=LowPassFilter
        self.PixelSize=PixelSize
        
        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Delete working directory if it exists, make a new one, and go there
        if (DoDeleteWorkingDir): 
            if os.path.exists(self.WorkingDir):
                shutil.rmtree(self.WorkingDir)
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Create selfiles with absolute pathnames in the WorkingDir
        mysel=selfile.selfile()
        mysel2=selfile.selfile()
        self.UntiltedSelFile=os.path.abspath(self.WorkingDir+'/'+UntiltedSelFile)
        self.TiltedSelFile=os.path.abspath(self.WorkingDir+'/'+TiltedSelFile)
        mysel.read(UntiltedSelFile)
        mysel2.read(TiltedSelFile)
        newsel=mysel.make_abspath()
        newsel2=mysel2.make_abspath()
        newsel.write(self.UntiltedSelFile)
        newsel2.write(self.TiltedSelFile)

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                       os.path.abspath(self.WorkingDir))
        os.chdir(self.WorkingDir)

        # Always prepare the library with all filenames
        self.prepare_name_library()
        
        if (DoImagePreparation):
            self.make_local_copies()
        if (DoUntiltedHeaders):
            self.set_headers_untilted()
        if (DoTiltedHeaders):
            self.set_headers_tilted()

        if (DoArtReconstruct):
            self.execute_art()

        if (DoWbpReconstruct):
            self.execute_wbp()

        if (DoLowPassFilter):
            self.execute_filter()

        # Return to parent dir
        os.chdir(os.pardir)


    # Make a libray of untilted selfiles and averages for all selected classes
    # And make the corresponding tilted selfile
    def prepare_name_library(self):
        import os,sys,glob
        import selfile
        self.untiltclasslist={}

        print '*********************************************************************'

        # Set self.PreviousML2DRoot to the ML2D Working dir if empty
        if self.PreviousML2DRoot=='':
            self.PreviousML2DRoot='ml2d'
        ml2d_abs_rootname=self.PreviousDirML2D+'/'+self.PreviousML2DRoot
        
        # Check whether the ML2D run has written docfiles already
        docfiles=glob.glob(ml2d_abs_rootname+'_it??????.doc')
        if len(docfiles)==0:
            message='No ML2D selfiles yet. Continue script after ML2D job completion... '
            print '* ',message
            self.log.error(message)
            sys.exit()
        else:
            # Check that no mirror option was used
            lastlogfile=docfiles[-1].replace('.doc','.log')
            fh=open(lastlogfile,'r')
            lines=fh.readlines()
            for line in lines:
                if (line.find('-mirror')>-1):
                    message=" * ERROR: you cannot use the -mirror flag if you want to do RCT\n"
                    message+=" * Repeat the ML2D run without providing this flag\n"
                    print '*',message
                    self.log.error(message)
                    sys.exit()
            # Loop over all classes selected for 3D-reconstruction
            lastitername=docfiles[-1].replace('.doc','')
            refs=self.SelectClasses.split(',')
            import utils_xmipp
            for ref in refs:
                # Copy selfile and average image of ML2DDir to WorkingDir
                unt_selfile=utils_xmipp.composeFileName(ml2d_abs_rootname+'_ref',ref,'sel')
                local_unt_selfile=utils_xmipp.composeFileName('rct_ref',ref,'')
                local_unt_selfile += '_untilted.sel'
                local_til_selfile=utils_xmipp.composeFileName('rct_ref',ref,'')
                local_til_selfile += '_tilted.sel'
                refavg=utils_xmipp.composeFileName(lastitername+'_ref',ref,'xmp')
                local_refavg=utils_xmipp.composeFileName('rct_ref',ref,'')
                local_refavg += '_untilted_avg.xmp'
                self.untiltclasslist[ref]=[local_unt_selfile,]
                self.untiltclasslist[ref].append(local_refavg)
                self.untiltclasslist[ref].append(local_til_selfile)
                self.untiltclasslist[ref].append(unt_selfile)
                self.untiltclasslist[ref].append(refavg)

    def make_local_copies(self):
        import os,shutil
        import selfile
        # Loop over all selected untilted classes
        for ref in self.untiltclasslist:
            local_unt_selfile=self.untiltclasslist[ref][0]
            local_refavg=self.untiltclasslist[ref][1]
            local_til_selfile=self.untiltclasslist[ref][2]
            unt_selfile=self.untiltclasslist[ref][3]
            refavg=self.untiltclasslist[ref][4]
            # Copy selfiles and average of untilted images
            shutil.copy(unt_selfile,local_unt_selfile)
            shutil.copy(refavg,local_refavg)
            # Generate corresponding tilted selfile
            self.make_tilted_selfile(local_unt_selfile,local_til_selfile)
            # Make a local copy of the images
            message='Making a local copy of the images in '+local_unt_selfile+' and '+local_til_selfile
            print '* ',message
            self.log.info(message)
            mysel=selfile.selfile()
            mysel.read(local_unt_selfile)
            newsel=mysel.copy_sel('local_untilted_images')
            newsel.write(local_unt_selfile)
            mysel.read(local_til_selfile)
            newsel=mysel.copy_sel('local_tilted_images')
            newsel.write(local_til_selfile)
            ori_local_til_selfile=local_til_selfile.replace('.sel','_all.sel')
            shutil.copy(local_til_selfile,ori_local_til_selfile)

            
    # This routine makes the corresponding selfile of the subset with tilted images
    # using the subset selfile of untilted images, and the original UntiltedSelFile & TiltedSelFile
    def make_tilted_selfile(self,name_unt_sel,name_til_sel):
        import selfile
        unt=selfile.selfile()
        pat1=selfile.selfile()
        pat2=selfile.selfile()
        pat1.read(self.UntiltedSelFile)
        pat2.read(self.TiltedSelFile)
        unt.read(name_unt_sel)
        til=unt.make_corresponding_subset(pat1,pat2)
        til.write(name_til_sel)
        
    def set_headers_untilted(self):
        import os
        import selfile

        print '*********************************************************************'
        print '*  Re-aligning untilted images of each class to set image headers'
        # Loop over all selected untilted classes
        for ref in self.untiltclasslist:

            # Perform a quick align2d to handle image headers correctly
            selfile=self.untiltclasslist[ref][0]
            reference=self.untiltclasslist[ref][1]
            command='xmipp_align2d -i '+selfile+' -ref '+reference+' -iter 2'
            print '* ',command
            self.log.info(command)
            os.system(command)

    def set_headers_tilted(self):
        import os,shutil
        print '*********************************************************************'
        print '*  Setting image headers of tilted images of each class'
        # Loop over all selected untilted classes
        for ref in self.untiltclasslist:
            unt_selfile=self.untiltclasslist[ref][0]
            til_selfile=self.untiltclasslist[ref][2]
            docfile=til_selfile.replace('.sel','.doc')
            ori_til_selfile=til_selfile.replace('.sel','_all.sel')
            shutil.copy(ori_til_selfile,til_selfile)
            command = 'xmipp_align_tilt_pairs -u ' + unt_selfile + \
                      ' -t '+til_selfile + \
                      ' -doc '+docfile + \
                      ' -max_shift ' + str(self.CenterMaxShift)
            command += self.AlignTiltPairsAdditionalParams
            print '* ',command
            self.log.info(command)
            os.system(command)
                
    def execute_art(self):
        import os
        for ref in self.untiltclasslist:
            til_selfile=self.untiltclasslist[ref][2]
            outname=til_selfile.replace('.sel','')
            outname='art_'+outname
            command='xmipp_reconstruct_art -i ' + str(til_selfile) + \
                     ' -o ' + str(outname) + \
                     ' -l ' + str(self.ArtLambda)
            if not self.ArtAdditionalParams=="":
                command+=' '+str(self.ArtAdditionalParams)

            print '* ',command
            self.log.info(command)
            os.system(command)

    def execute_wbp(self):
        import os
        for ref in self.untiltclasslist:
            til_selfile=self.untiltclasslist[ref][2]
            outname=til_selfile.replace('.sel','.vol')
            outname='wbp_'+outname
            command='xmipp_reconstruct_wbp -i ' + til_selfile + \
                     ' -o ' + outname + \
                     ' -threshold ' + str(self.WbpThreshold)
            if not self.WbpAdditionalParams=="":
                command+=' '+str(self.WbpAdditionalParams)

            print '* ',command
            self.log.info(command)
            os.system(command)

    def execute_filter(self):
        import os
        for ref in self.untiltclasslist:
            til_selfile=self.untiltclasslist[ref][2]
            volname=til_selfile.replace('.sel','.vol')
            wbpname='wbp_'+volname
            artname='art_'+volname
            if os.path.exists(wbpname):
                filname=wbpname.replace('.vol','_filtered.vol')
                command='xmipp_fourier_filter -o ' + filname + \
                 ' -i ' + wbpname  + \
                 ' -sampling ' + str(self.PixelSize) + \
                 ' -low_pass ' + str(self.LowPassFilter)
                print '* ',command
                self.log.info(command)
                os.system(command)
            if os.path.exists(artname):
                filname=wbpname.replace('.vol','_filtered.vol')
                command='xmipp_fourier_filter -o ' + filname + \
                 ' -i ' + artname  + \
                 ' -sampling ' + str(self.PixelSize) + \
                 ' -low_pass ' + str(self.LowPassFilter)
                print '* ',command
                self.log.info(command)
                os.system(command)

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'
#		
# Main
#     
if __name__ == '__main__':

    # create ML2D_class object

    RCT=RCT_class(UntiltedSelFile,
                  TiltedSelFile,
                  WorkingDir,
                  DoDeleteWorkingDir,
                  ProjectDir,
                  LogDir,
                  PreviousDirML2D,
                  SelectClasses,
                  DoImagePreparation,
                  DoUntiltedHeaders,
                  DoTiltedHeaders,
                  CenterMaxShift,
                  AlignTiltPairsAdditionalParams,
                  DoArtReconstruct,
                  ArtLambda,
                  ArtAdditionalParams,
                  DoWbpReconstruct,
                  WbpThreshold,
                  WbpAdditionalParams,
                  DoLowPassFilter,
                  LowPassFilter,
                  PixelSize)

    # close 
    RCT.close()

