#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for random conical tilt reconstruction
#
# Example use:
# ./protocol_rct.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Selfile with all untilted images: (in ProjectDir)
UntiltedSelFile="untSelect.sel"
# Selfile with all tilted images: (in ProjectDir)
TiltedSelFile="tilSelect.sel"
# Working subdirectory:
WorkingDir="test1"
# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=True
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols/G40P"
# {expert} Directory name for logfiles:
""" All logfiles will be stored in $ProjectDir/$LogDir
"""
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} Previous ML2D classification (WITHOUT INCLUDING MIRRORS!)
#------------------------------------------------------------------------------------------------
# Directory of previous ML2D-classification on the untilted images (from ProjectDir)
PreviousDirML2D="ML2D/ML3ref"
# {expert} Rootname for ML2D run (only provide if different from its working directory)
PreviousML2DRoot=""
# Which of these classes do you want to reconstruct? (Separate numbers by comma's)
SelectClasses="1,2"
#------------------------------------------------------------------------------------------------
# {section} Prepare image headers
#------------------------------------------------------------------------------------------------
# Set angles and offsets to headers of untilted images?
DoSetHeadersUntilted=True
# Set angles and offsets to headers of tilted images?
DoSetHeadersTilted=True
# Perform centering of the tilted particles?
DoCenterTilted=True
# Maximum allowed shift for tilted particles (in pixels):
""" Particles that shift more will be discarded
"""
CenterMaxShift=10
# Use cosine stretching in the centering of the tilted particles?
DoUseCosineStretching=True
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
                 DoSetHeadersUntilted,
                 DoSetHeadersTilted,
                 DoCenterTilted,
                 CenterMaxShift,
                 DoUseCosineStretching,
                 DoArtReconstruct,
                 ArtLambda,
                 ArtAdditionalParams,
                 DoWbpReconstruct,
                 WbpThreshold,
                 WbpAdditionalParams):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log
        
        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.UntiltedSelFile=self.ProjectDir+'/'+UntiltedSelFile
        self.TiltedSelFile=self.ProjectDir+'/'+TiltedSelFile
        self.PreviousDirML2D=self.ProjectDir+'/'+PreviousDirML2D
        self.PreviousML2DRoot=PreviousML2DRoot
        self.SelectClasses=SelectClasses
        self.DoCenterTilted=DoCenterTilted
        self.CenterMaxShift=CenterMaxShift
        self.DoUseCosineStretching=DoUseCosineStretching
        self.ArtLambda=ArtLambda
        self.ArtAdditionalParams=ArtAdditionalParams
        self.WbpThreshold=WbpThreshold
        self.WbpAdditionalParams=WbpAdditionalParams

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
        os.chdir(self.WorkingDir)

        self.prepare_classes()
            
        if (DoSetHeadersUntilted):
            self.set_headers_untilted()

        if (DoSetHeadersTilted):
            self.set_headers_tilted()

        if (DoArtReconstruct):
            self.execute_art()

        if (DoWbpReconstruct):
            self.execute_wbp()

        # Return to parent dir
        os.chdir(os.pardir)


    # Make a libray of untilted selfiles and averages for all selected classes
    # And make the corresponding tilted selfile
    def prepare_classes(self):
        import os,sys,glob,shutil
        import selfile
        self.untiltclasslist={}

        print '*********************************************************************'

        # Set self.PreviousML2DRoot to the ML2D Working dir if empty
        if self.PreviousML2DRoot=="":
            head,tail=os.path.split(os.path.normpath(self.PreviousDirML2D))
            self.PreviousML2DRoot=tail
        ml2d_abs_rootname=self.PreviousDirML2D+'/'+self.PreviousML2DRoot
        
        # Check whether the ML2D run has written docfiles already
        docfiles=glob.glob(ml2d_abs_rootname+'_it?????.doc')
        if len(docfiles)==0:
            message='No ML2D selfiles yet. Continue script after ML2D job completion... '
            print '* ',message
            self.log.error(message)
            sys.exit()
        else:
            # Loop over all classes selected for 3D-reconstruction
            lastitername=docfiles[-1].replace('.doc','')
            refs=self.SelectClasses.split(',')
            for ref in refs:
                # Copy selfile and average image of ML2DDir to WorkingDir
                unt_selfile=ml2d_abs_rootname+'_ref'+str(ref).zfill(5)+'.sel'
                local_unt_selfile='rct_ref'+str(ref).zfill(5)+'_untilted.sel'
                refavg=lastitername+'_ref'+str(ref).zfill(5)+'.xmp'
                local_refavg='rct_ref'+str(ref).zfill(5)+'_untilted_avg.xmp'
                shutil.copy(unt_selfile,local_unt_selfile)
                shutil.copy(refavg,local_refavg)
                local_til_selfile=self.make_tilted_selfile(local_unt_selfile)
                self.untiltclasslist[ref]=[local_unt_selfile,]
                self.untiltclasslist[ref].append(local_refavg)
                self.untiltclasslist[ref].append(local_til_selfile)
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
                 
    # This routine makes the corresponding selfile of the subset with tilted images
    # using the subset selfile of untilted images, and the original UntiltedSelFile & TiltedSelFile
    def make_tilted_selfile(self,name_unt_sel):
        import selfile
        unt=selfile.selfile()
        pat1=selfile.selfile()
        pat2=selfile.selfile()
        pat1.read(self.UntiltedSelFile)
        pat2.read(self.TiltedSelFile)
        unt.read(name_unt_sel)
        til=unt.make_corresponding_subset(pat1,pat2)
        name_til_sel=name_unt_sel.replace('_untilted.sel','_tilted.sel')
        til.write(name_til_sel)
        return name_til_sel
        
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
        import os
        print '*********************************************************************'
        print '*  Setting image headers of tilted images of each class'
        # Loop over all selected untilted classes
        for ref in self.untiltclasslist:
            unt_selfile=self.untiltclasslist[ref][0]
            til_selfile=self.untiltclasslist[ref][2]
            docfile=til_selfile.replace('.sel','.doc')
#            command='xmipp_align_tilt_pairs -u '+unt_selfile+\
            command='xmipp_centilt -u '+unt_selfile+\
                     ' -t '+til_selfile+\
                     ' -doc '+docfile+\
                     ' -max_shift '+str(self.CenterMaxShift)
            if not self.DoCenterTilted:
                command+=' -skip_centering'
            if not self.DoUseCosineStretching:
                command+=' -skip_stretching'
            print '* ',command
            self.log.info(command)
            os.system(command)
                
    def execute_art(self):
        import os
        for ref in self.untiltclasslist:
            til_selfile=self.untiltclasslist[ref][2]
            outname=til_selfile.replace('.sel','')
            outname='art_'+outname
#            command='xmipp_reconstruct_art -i ' + til_selfile + \
            command='xmipp_art -i ' + til_selfile + \
                     ' -o ' + outname + \
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
#            command='xmipp_reconstruct_wbp -i ' + til_selfile + \ 
            command='xmipp_wbp -i ' + til_selfile + \
                    ' -o ' + outname + \
                     ' -threshold ' + str(self.WbpThreshold)
            if not self.WbpAdditionalParams=="":
                command+=' '+str(self.WbpAdditionalParams)

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
                  DoSetHeadersUntilted,
                  DoSetHeadersTilted,
                  DoCenterTilted,
                  CenterMaxShift,
                  DoUseCosineStretching,
                  DoArtReconstruct,
                  ArtLambda,
                  ArtAdditionalParams,
                  DoWbpReconstruct,
                  WbpThreshold,
                  WbpAdditionalParams)

    # close 
    RCT.close()

