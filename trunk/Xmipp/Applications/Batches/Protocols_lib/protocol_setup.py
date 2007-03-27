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
# {is-setup} Parameters
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
# Micrograph preprocessing (part A)
DoSetupPreProcessA=False
# Micrograph preprocessing (part B)
DoSetupPreProcessB=False
# ML2D multi-reference alignment
DoSetupML2D=False
# 2D kerdenSOM classification 
DoSetupKerdensom=False
# Rotational spectra classification
DoSetupRotSpectra=False
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
                     DoSetupPreProcessA,
                     DoSetupPreProcessB,
                     PreProcessDir,
                     ImagesDir,
                     DoSetupML2D,
                     DoSetupKerdensom,
                     ML2DDir,
                     DoSetupRotSpectra,
                     RotSpectraDir,
                     AutoLaunch):

            import os
            import sys
            import string

            self.ProjectDir=ProjectDir
            self.LogDir=LogDir
            self.DoSetupPreProcessA=DoSetupPreProcessA
            self.DoSetupPreProcessB=DoSetupPreProcessB
            self.PreProcessDir=PreProcessDir
            self.ImagesDir=ImagesDir
            self.DoSetupML2D=DoSetupML2D
            self.DoSetupKerdensom=DoSetupKerdensom
            self.ML2DDir=ML2DDir
            self.DoSetupRotSpectra=DoSetupRotSpectra
            self.RotSpectraDir=RotSpectraDir
            self.AutoLaunch=AutoLaunch
            
            self.SYSTEMSCRIPTDIR="/home/scheres/Xmipp/Applications/Batches/Protocols_lib"
            self.library={}
            self.library['DoSetupPreProcessA']=[self.ProjectDir+'/'+self.PreProcessDir,
                                              'protocol_preprocess_A.py']
            self.library['DoSetupPreProcessB']=[self.ProjectDir+'/'+self.PreProcessDir,
                                              'protocol_preprocess_A.py']
            self.library['DoSetupML2D']=[self.ProjectDir+'/'+self.ML2DDir,
                                              'protocol_ML2D.py']
            self.library['DoSetupKerdensom']=[self.ProjectDir+'/'+self.ML2DDir,
                                              'protocol_kerdensom.py']
            self.library['DoSetupRotSpectra']=[self.ProjectDir+'/'+self.RotSpectraDir,
                                              'protocol_rotspectra.py']
            # For automated editing of default directories
            self.DEFAULTDIRS={"ProjectDir":self.ProjectDir,
                              "LogDir":self.LogDir,
                              "PreProcessDir":self.PreProcessDir,
                              "ImagesDir":self.ImagesDir,
                              "ML2DDir":self.ML2DDir,
                              "RotSpectraDir":self.RotSpectraDir}

            # Execute program
            self.perform_setup()

        def perform_setup(self):

            import os

            if (self.AutoLaunch!=""):
                self.setup_protocol(self.library[self.AutoLaunch][0],
                                    self.library[self.AutoLaunch][1])

            if (self.DoSetupPreProcessA):
                self.setup_protocol(self.library[DoSetupPreProcessA][0],
                                    self.library[DoSetupPreProcessA][1])
            if (self.DoSetupPreProcessB):
                self.setup_protocol(self.library[DoSetupPreProcessB][0],
                                    self.library[DoSetupPreProcessB][1])
            if (self.DoSetupML2D):
                self.setup_protocol(self.library[DoSetupML2D][0],
                                    self.library[DoSetupML2D][1])
            if (self.DoSetupKerdensom):
                self.setup_protocol(self.library[DoSetupKerdensom][0],
                                    self.library[DoSetupKerdensom][1])
            if (self.DoSetupRotSpectra):
                self.setup_protocol(self.library[DoSetupRotSpectra][0],
                                    self.library[DoSetupRotSpectra][1])


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

            os.chdir(dir)
            command='python '+str(self.SYSTEMSCRIPTDIR)+'/protocol_gui.py '+dst+' &'
            os.system(command)


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
                                DoSetupPreProcessA,
                                DoSetupPreProcessB,
                                PreProcessDir,
                                ImagesDir,
                                DoSetupML2D,
                                DoSetupKerdensom,
                                ML2DDir,
                                DoSetupRotSpectra,
                                RotSpectraDir,
                                AutoLaunch)
    setup.close()
