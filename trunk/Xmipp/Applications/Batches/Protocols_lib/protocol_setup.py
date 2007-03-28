#!/usr/bin/env python
#------------------------------------------------------------------------
# General script to setup standard Xmipp protocols
#  - preprocessing (parts A & B) 
#  - 2D image alignment and classification (by ml2d & kerdenSOM)
#  - 2D image analysis by classification of rotational spectra
#  - 3D classification by ml3D
#  - 3D refinement by standard projection matching
#  - 3D refinement by a multi-resolution wavelet approach
#
# Example use:
# ./setup_protocols.py
#
# Author: Sjors Scheres, March 2007
#------------------------------------------------------------------------
# Choose the protocol(s) you want to setup:
#------------------------------------------------------------------------
# {setup-pre} part A: preprocess micrographs
DoSetupPreProcessA=False
# {setup-pre} part B: preprocess particles
DoSetupPreProcessB=False
# {setup-2d} Maximum-likelihood refinement
DoSetupML2D=False
# {setup-2d} kerdenSOM classification 
DoSetupKerdensom=False
# {setup-2d} Rotational spectra classification
DoSetupRotSpectra=False
# {setup-3d} Random Conical Tilt reconstruction
DoSetupRCT=False
#------------------------------------------------------------------------
# {section} Global Parameters
#------------------------------------------------------------------------
# Root directory name for this project:
ProjectDir="/home/scheres/work/protocols"
# {expert} Directory name for logfiles:
LogDir="Logs"
# {expert} Directory name for preprocessing:
PreProcessDir="Preprocessing"
# {expert} Directory name for particle images:
ImagesDir="Images"
# {expert} Directory name for 2D alignment and classification:
ML2DDir="Analyse2d"
# {expert} Directory name for rotational spectra classification:
RotSpectraDir="Rotspectra"
# {expert} Directory name for Random Conical Tilt reconstruction:
RCTDir="RCT"
#
#
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class setup_protocols_class:
       #init variables
        def __init__(self,
                     ProjectDir,
                     LogDir,
                     ImagesDir,
                     PreProcessDir,
                     ML2DDir,
                     RotSpectraDir,
                     RCTDir,
                     DoSetupPreProcessA,
                     DoSetupPreProcessB,
                     DoSetupML2D,
                     DoSetupKerdensom,
                     DoSetupRotSpectra,
                     DoSetupRCT,
                     AutoLaunch):

            import os

            self.ProjectDir=ProjectDir
            self.LogDir=LogDir
            self.ImagesDir=ImagesDir
            self.PreProcessDir=PreProcessDir
            self.ML2DDir=ML2DDir
            self.RotSpectraDir=RotSpectraDir
            self.RCTDir=RCTDir

            self.DoSetupPreProcessA=DoSetupPreProcessA
            self.DoSetupPreProcessB=DoSetupPreProcessB
            self.DoSetupML2D=DoSetupML2D
            self.DoSetupKerdensom=DoSetupKerdensom
            self.DoSetupRotSpectra=DoSetupRotSpectra
            self.DoSetupRCT=DoSetupRCT

            self.AutoLaunch=AutoLaunch
            self.SYSTEMSCRIPTDIR="/home/scheres/Xmipp/Applications/Batches/Protocols_lib"

            # Which scripts and which directories to use
            self.library={}
            self.library['DoSetupPreProcessA']=[self.DoSetupPreProcessA,
                                                self.PreProcessDir,
                                                'protocol_preprocess_A.py']
            self.library['DoSetupPreProcessB']=[self.DoSetupPreProcessB,
                                                self.PreProcessDir,
                                                'protocol_preprocess_B.py']
            self.library['DoSetupML2D']=[self.DoSetupML2D,
                                         self.ML2DDir,
                                         'protocol_ML2D.py']
            self.library['DoSetupKerdensom']=[self.DoSetupKerdensom,
                                              self.ML2DDir,
                                              'protocol_kerdensom.py']
            self.library['DoSetupRotSpectra']=[self.DoSetupRotSpectra,
                                               self.RotSpectraDir,
                                               'protocol_rotspectra.py']
            self.library['DoSetupRCT']=[self.DoSetupRCT,
                                        self.RCTDir,
                                        'protocol_RCT.py']

            # For automated editing of default directories in protocols
            self.DEFAULTDIRS={"ProjectDir":self.ProjectDir,
                              "LogDir":self.LogDir,
                              "PreProcessDir":self.PreProcessDir,
                              "ImagesDir":self.ImagesDir,
                              "ML2DDir":self.ML2DDir,
                              "RotSpectraDir":self.RotSpectraDir,
                              "RCTDir":self.RCTDir}
            

            # Perform the actual setup:
            if (self.AutoLaunch!=""):
                # A. Setup from GUI (Autolaunch)
                # This will copy the (modified) protocol script to the corresponding directory
                # and will automatically launch the GUI for this protocol
                self.setup_protocol(self.library[self.AutoLaunch][1],
                                    self.library[self.AutoLaunch][2])
            else:
                # B. Setup from this script:
                # This will only copy the (modified) protocol script to the corresponding directory
                for var in self.library:
                    if (self.library[var][0]):
                        self.setup_protocol(self.library[var][1],
                                            self.library[var][2])

        def modify_script_header(self,src):

            import sys
            # Read script header and body
            fh=open(src,'r')
            header_lines=[]
            body_lines=[]
            isheader=False
            while (isheader!=True):
                line=fh.readline()
                if line=="":
                    print "Error, this file does not have a {end-of-header} label"
                    sys.exit()
                header_lines.append(line)
                if "{end-of-header}" in line:
                    isheader=True

            body_lines=fh.readlines()
            fh.close()

            # Loop over all project-related directories and preset default directories
            for dir in self.DEFAULTDIRS:
                newheader_lines=[]
                for line in header_lines:
                    if ((not line[0]=="#" and
                         not line[0]==" " and
                         not line[0]=="\t" and
                         not line[0]=="\n" and
                         not line[0]=="\"") and (dir in line)):
                        args=line.split("=")
                        lineb=str(args[0])+'=\"'+self.DEFAULTDIRS[dir]+'\"\n'
                    else:
                        lineb=line
                    newheader_lines.append(lineb)
                header_lines=newheader_lines
            return header_lines+body_lines

        def setup_protocol(self,dir,script):
            import os
            import shutil

            dir=self.ProjectDir+'/'+dir
            if not os.path.exists(dir):
                os.makedirs(dir)
            src=str(self.SYSTEMSCRIPTDIR)+"/"+str(script)
            dst=str(dir)+"/"+str(script)
            if os.path.exists(dst):
                print "File "+dst+" already existed (now updated)"

            text=self.modify_script_header(src)
            fh=open(dst,'w')
            fh.writelines(text)
            fh.close()

            if (self.AutoLaunch!=""):
                os.chdir(dir)
                command='python '+str(self.SYSTEMSCRIPTDIR)+'/protocol_gui.py '+dst+' &'
                os.system(command)
                os.chdir(self.ProjectDir)


	def close(self):
                message=" Done setting up!"
		print '* ',message
		print '*********************************************************************'
#
# Main
#
if __name__ == '__main__':

    AutoLaunch=""
    setup=setup_protocols_class(ProjectDir,
                                LogDir,
                                ImagesDir,
                                PreProcessDir,
                                ML2DDir,
                                RotSpectraDir,
                                RCTDir,
                                DoSetupPreProcessA,
                                DoSetupPreProcessB,
                                DoSetupML2D,
                                DoSetupKerdensom,
                                DoSetupRotSpectra,
                                DoSetupRCT,
                                AutoLaunch)
    setup.close()
