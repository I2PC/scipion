#!/usr/bin/env python
#------------------------------------------------------------------------
# General script to setup standard Xmipp protocols
#  - Preprocessing of micrographs
#  - Manual particle picking
#  - Preprocessing of extracted particles
#  - 2D image alignment and classification (by ML2D & kerdenSOM)
#  - 2D image analysis by classification of rotational spectra
#  - 3D classification by ml3D
#  - 3D projection matching refinement
#  - 3D high resolution refinement
#
# Example use:
# ./setup_protocols.py
#
# {please cite} The Xmipp team (2007) Nature Protocols, xx,xx-xx
#
# Author: Sjors Scheres, March 2007
#------------------------------------------------------------------------
# Choose the protocol(s) you want to setup:
#------------------------------------------------------------------------
# {setup-pre} Preprocess micrographs
SetupPreProcessMicrographs=False
# {setup-pre} Manual particle selection
SetupParticlePick=False
# {setup-pre} Preprocess particles
SetupPreProcessParticles=False
# {setup-2d} ML2D classification
SetupML2D=False
# {setup-2d} kerdenSOM classification 
SetupKerdensom=False
# {setup-2d} Rotational spectra classification
SetupRotSpectra=False
# {setup-3d} Random conical tilt
SetupRCT=False
# {setup-3d} ML3D classification
SetupML3D=False
# {setup-3d} Projection matching refinement
SetupProjMatch=False
# {setup-3d} High resolution refinement
SetupHighRes3d=False
#------------------------------------------------------------------------
# {section} Global Parameters
#------------------------------------------------------------------------
# Absolute path to the root directory name for this project:
ProjectDir='/home/scheres/xmipp/applications/scripts/protocols'
# Directory name for logfiles:
LogDir='Logs'
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
                     SetupPreProcessMicrographs,
                     SetupParticlePick,
                     SetupPreProcessParticles,
                     SetupML2D,
                     SetupKerdensom,
                     SetupRotSpectra,
                     SetupRCT,
                     SetupML3D,
                     SetupProjMatch,
                     SetupHighRes3d,
                     ProjectDir,
                     LogDir,
                     AutoLaunch):

            import os,sys

            self.SetupPreProcessMicrographs=SetupPreProcessMicrographs
            self.SetupParticlePick=SetupParticlePick
            self.SetupPreProcessParticles=SetupPreProcessParticles
            self.SetupML2D=SetupML2D
            self.SetupKerdensom=SetupKerdensom
            self.SetupRotSpectra=SetupRotSpectra
            self.SetupRCT=SetupRCT
            self.SetupML3D=SetupML3D
            self.SetupProjMatch=SetupProjMatch
            self.SetupHighRes3d=SetupHighRes3d

            self.ProjectDir=ProjectDir
            self.LogDir=LogDir
            self.AutoLaunch=AutoLaunch

            scriptdir=os.path.expanduser('~')+'/scripts/'
            sys.path.append(scriptdir) # add default search path
            self.SYSTEMSCRIPTDIR=scriptdir

            # Which scripts and which directories to use
            self.library={}
            self.library['SetupPreProcessMicrographs']=[self.SetupPreProcessMicrographs,
                                                ['protocol_preprocess_micrographs.py','visualize_preprocess_micrographs.py']]
            self.library['SetupParticlePick']=[self.SetupParticlePick,
                                                ['protocol_particle_pick.py']]
            self.library['SetupPreProcessParticles']=[self.SetupPreProcessParticles,
                                                ['protocol_preprocess_particles.py','visualize_preprocess_particles.py']]
            self.library['SetupML2D']=[self.SetupML2D,
                                         ['protocol_ml2d.py','visualize_ml2d.py']]
            self.library['SetupKerdensom']=[self.SetupKerdensom,
                                              ['protocol_kerdensom.py','visualize_kerdensom.py']]
            self.library['SetupRotSpectra']=[self.SetupRotSpectra,
                                               ['protocol_rotspectra.py','visualize_rotspectra.py']]
            self.library['SetupRCT']=[self.SetupRCT,
                                        ['protocol_rct.py','visualize_rct.py']]
            self.library['SetupML3D']=[self.SetupML3D,
                                        ['protocol_ml3d.py','visualize_ml3d.py']]
            self.library['SetupProjMatch']=[self.SetupProjMatch,
                                        ['protocol_projmatch.py','visualize_projmatch.py']]
            self.library['SetupHighRes3d']=[self.SetupHighRes3d,
                                        ['protocol_highres3d.py','visualize_highres3d.py']]

            # For automated editing of default directories in protocols
            self.DEFAULTDIRS={"ProjectDir":self.ProjectDir,
                              "LogDir":self.LogDir
                              }
            

            # Perform the actual setup:
            if (self.AutoLaunch!=""):
                # A. Setup from GUI (Autolaunch)
                # This will copy the (modified) protocol script to the corresponding directory
                # and will automatically launch the GUI for this protocol
                self.setup_protocol(self.library[self.AutoLaunch][1])
            else:
                # B. Setup from this script:
                # This will only copy the (modified) protocol script to the corresponding directory
                for var in self.library:
                    if (self.library[var][0]):
                        self.setup_protocol(self.library[var][1])

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
                if (line.find('{end-of-header}') > -1):
                    isheader=True

            body_lines=fh.readlines()
            fh.close()

            # Loop over all project-related directories and preset default directories
            for dir in self.DEFAULTDIRS:
                newheader_lines=[]
                for line in header_lines:
                    if ( (not line[0]=="#" and
                          not line[0]==" " and
                          not line[0]=="\t" and
                          not line[0]=="\n" and
                          not line[0]=="\"") and (line.find(dir) > -1) ):
                        args=line.split("=")
                        lineb=str(args[0])+'=\"'+self.DEFAULTDIRS[dir]+'\"\n'
                    else:
                        lineb=line
                    newheader_lines.append(lineb)
                header_lines=newheader_lines
            return header_lines+body_lines

        def setup_protocol(self,scripts):
            import os
            import shutil

            os.chdir(self.ProjectDir)
                
            for script in scripts:
                src=str(self.SYSTEMSCRIPTDIR)+"/"+str(script)
                dst=str(script)
                if os.path.exists(dst):
                    src=dst
                    print "* File "+dst+" already existed (now updated)"

                text=self.modify_script_header(src)
                fh=open(dst,'w')
                fh.writelines(text)
                fh.close()
                os.chmod(dst,0755)

            # Only auto-launch the first script
            script=scripts[0]
            if (self.AutoLaunch!=""):
                command='python '+str(self.SYSTEMSCRIPTDIR)+'/protocol_gui.py '+str(script)+' &'
                os.system(command)


#
# Main
#
if __name__ == '__main__':

    import sys
    if (len(sys.argv) < 2):
        AutoLaunch=""
    else:
        AutoLaunch=sys.argv[1]

    setup=setup_protocols_class(SetupPreProcessMicrographs,
                                SetupParticlePick,
                                SetupPreProcessParticles,
                                SetupML2D,
                                SetupKerdensom,
                                SetupRotSpectra,
                                SetupRCT,
                                SetupML3D,
                                SetupProjMatch,
                                SetupHighRes3d,
                                ProjectDir,
                                LogDir,
                                AutoLaunch)

