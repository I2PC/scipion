#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for random conical tilt reconstruction
#
# Example use:
# ./xmipp_protocol_rct.py
#
# Author: Sjors Scheres, March 2007
#
from protlib_base import *

class ProtRCT(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.rct.name, scriptname, project)
        self.Import = 'from protocol_rct import *'

    #init variables
#    def __init__(self,
#                 UntiltedSelFile,
#                 TiltedSelFile,
#                 WorkingDir,
#                 DoDeleteWorkingDir,
#                 ProjectDir,
#                 LogDir,
#                 PreviousDirML2D,
#                 SelectClasses,
#                 DoImagePreparation,
#                 DoUntiltedHeaders,
#                 DoTiltedHeaders,
#                 CenterMaxShift,
#                 AlignTiltPairsAdditionalParams,
#                 DoReconstruct,
#                 ReconstructMethod,
#                 ReconstructAdditionalParams,
#                 DoLowPassFilter,
#                 LowPassFilter,
#                 PixelSize):
#	     
#        import os,sys,shutil
#        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
#        sys.path.append(scriptdir) # add default search path
#        from xmipp import MetaData
#        import log
#        import utils_xmipp
#        
#        self.WorkingDir = WorkingDir
#        self.ProjectDir = ProjectDir
#        self.PreviousDirML2D = os.path.abspath(PreviousDirML2D)
#        self.PreviousML2DRoot = PreviousML2DRoot
#        self.SelectClasses = utils_xmipp.getCommaSeparatedIntegerList(SelectClasses)
#        self.CenterMaxShift = CenterMaxShift
#        self.AlignTiltPairsAdditionalParams = AlignTiltPairsAdditionalParams
#        self.ReconstructMethod = ReconstructMethod
#        self.ReconstructAdditionalParams = ReconstructAdditionalParams
#        self.LowPassFilter = LowPassFilter
#        self.PixelSize = PixelSize
#        
#        # Setup logging
#        self.log=log.init_log_system(self.ProjectDir,
#                                     LogDir,
#                                     sys.argv[0],
#                                     self.WorkingDir)
#                
#        # Delete working directory if it exists, make a new one, and go there
#        if (DoDeleteWorkingDir): 
#            if (self.WorkingDir == ""):
#                raise RuntimeError,"No working directory given"
#            if os.path.exists(self.WorkingDir):
#                shutil.rmtree(self.WorkingDir)
#        if not os.path.exists(self.WorkingDir):
#            os.makedirs(self.WorkingDir)
#
#        # Create metadatas with absolute pathnames in the WorkingDir
#        self.UntiltedSelFile = os.path.abspath(self.WorkingDir+'/'+UntiltedSelFile)
#        self.TiltedSelFile = os.path.abspath(self.WorkingDir+'/'+TiltedSelFile)
#        mysel = MetaData(UntiltedSelFile)
#        mysel2 = MetaData(TiltedSelFile)
#        mysel.write(self.UntiltedSelFile)
#        mysel2.write(self.TiltedSelFile)
#
#        # Backup script
#        log.make_backup_of_script_file(sys.argv[0],
#                                       os.path.abspath(self.WorkingDir))
#        os.chdir(self.WorkingDir)
#
#        # Always prepare the library with all filenames
#        self.prepare_name_library()
#        
#        if (DoImagePreparation):
#            self.make_local_copies()
#            
#        if (DoUntiltedHeaders):
#            self.set_headers_untilted()
#            
#        if (DoTiltedHeaders):
#            self.set_headers_tilted()
#
#        if (DoReconstruct):
#            self.execute_reconstruction()
#
#        if (DoLowPassFilter):
#            self.execute_filter()
#
#        # Return to parent dir
#        os.chdir(os.pardir)


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
        ml2d_abs_rootname = self.PreviousDirML2D+'/'+self.PreviousML2DRoot
        
        # Check whether the ML2D run has written docfiles already
        docfiles = glob.glob(ml2d_abs_rootname+'_it??????.doc')
        lastiter = len(docfiles)
        if (lastiter==0):
            message='No ML2D docfiles yet. Continue script after ML2D job completion... '
            print '* ',message
            self.log.error(message)
            sys.exit()
        else:
            # Check that no mirror option was used
            lastitername=ml2d_abs_rootname+'_it'+str(lastiter).zfill(6)
            self.lastdocfile=lastitername+'.doc'
            lastlogfile=lastitername+'.log'
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
            import utils_xmipp
            for ref in self.SelectClasses:
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
        import launch_job
        from xmipp import MetaData
        # Loop over all selected untilted classes
        for ref in self.untiltclasslist:
            local_unt_selfile = self.untiltclasslist[ref][0]
            local_refavg = self.untiltclasslist[ref][1]
            local_til_selfile = self.untiltclasslist[ref][2]
            unt_selfile = self.untiltclasslist[ref][3]
            refavg = self.untiltclasslist[ref][4]
            # Copy selfiles and average of untilted images
            shutil.copy(unt_selfile, local_unt_selfile)
            shutil.copy(refavg, local_refavg)
            # Generate corresponding tilted selfile
            self.make_tilted_selfile(local_unt_selfile, local_til_selfile)
            # Make a local copy of the images
            message='Making a local copy of the images in '+local_unt_selfile+' and '+local_til_selfile
            print '* ',message
            self.log.info(message)
            # Assign angles to the untilted images
            local_unt_docfile = local_unt_selfile.replace('.sel','.doc')
            command=' -i '+self.lastdocfile+' -sel '+local_unt_selfile+' -o '+local_unt_docfile
            launchJob("xmipp_docfile_select_subset",
                                  command,
                                  self.log,
                                  False,1,1,'')
            mysel = MetaData(local_unt_selfile)
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
        import launch_job

        print '*********************************************************************'
        print '*  Assigning in-plane transformation from ML2D to headers of local untilted images'
        # Loop over all selected untilted classes
        for ref in self.untiltclasslist:
            local_unt_selfile=self.untiltclasslist[ref][0]
            local_unt_docfile = local_unt_selfile.replace('.sel','.doc')
            command=' -i '+local_unt_docfile+' -o '+local_unt_selfile+' -force -columns 0 0 3 4 5'
            launchJob("xmipp_header_assign",
                                  command,
                                  self.log,
                                  False,1,1,'')
            command=' -i '+local_unt_selfile+' -o '+local_unt_docfile
            launchJob("xmipp_header_extract",
                                  command,
                                  self.log,
                                  False,1,1,'')

    def set_headers_tilted(self):
        import os,shutil
        import launch_job
        print '*********************************************************************'
        print '*  Setting image headers of tilted images of each class'
        # Loop over all selected untilted classes
        for ref in self.untiltclasslist:
            unt_selfile=self.untiltclasslist[ref][0]
            til_selfile=self.untiltclasslist[ref][2]
            docfile=til_selfile.replace('.sel','.doc')
            ori_til_selfile=til_selfile.replace('.sel','_all.sel')
            shutil.copy(ori_til_selfile,til_selfile)
            command = ' -u ' + unt_selfile + \
                      ' -t '+til_selfile + \
                      ' -doc '+docfile + \
                      ' -max_shift ' + str(self.CenterMaxShift)
            command += self.AlignTiltPairsAdditionalParams
            launchJob("xmipp_align_tilt_pairs",
                                  command,
                                  self.log,
                                  False,1,1,'')
                
    def execute_reconstruction(self):
        import os
        import launch_job
        for ref in self.untiltclasslist:
            til_selfile=self.untiltclasslist[ref][2]
            outname=til_selfile.replace('.sel','')
            if (self.ReconstructMethod=='art'):
                program='xmipp_reconstruct_art'
                command=' -i ' + til_selfile + \
                        ' -o ' + outname
            elif (self.ReconstructMethod=='wbp'):
                program='xmipp_reconstruct_wbp'
                command=' -i ' + til_selfile + \
                        ' -o ' + outname+'.vol'
            elif (self.ReconstructMethod=='fourier'):
                program='xmipp_reconstruct_fourier'
                command=' -i ' + til_selfile + \
                        ' -o ' + outname+'.vol'
            else:
                message = "Error: unrecognized reconstruction method: ", self.ReconstructMethod
            if not self.ReconstructAdditionalParams=="":
                command += ' ' + self.ReconstructAdditionalParams

            launchJob(program,
                                  command,
                                  self.log,
                                  False,1,1,'')

    def execute_filter(self):
        import os
        import launch_job
        for ref in self.untiltclasslist:
            til_selfile=self.untiltclasslist[ref][2]
            volname=til_selfile.replace('.sel','.vol')
            if os.path.exists(volname):
                filname=volname.replace('.vol','_filtered.vol')
                command=' -o ' + filname + \
                 ' -i ' + volname  + \
                 ' -sampling ' + str(self.PixelSize) + \
                 ' -low_pass ' + str(self.LowPassFilter)
                launchJob("xmipp_fourier_filter",
                                      command,
                                      self.log,
                                      False,1,1,'')

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

