#!/usr/bin/env python
#-----------------------------------------------------------------------
# General script to setup standard Xmipp protocols
#  - Preprocessing of micrographs
#  - Manual particle picking
#  - Preprocessing of extracted particles
#  - 2D image alignment and classification (by ML2D, CL2D & kerdenSOM)
#  - 2D image analysis by classification of rotational spectra
#  - 3D classification by ml3D
#  - 3D projection matching refinement
#
# Example use:
# ./setup_protocols.py
#
# {please cite} The Xmipp team (2007) Nature Protocols, 3, 977-990
#
# Author: Sjors Scheres, March 2007
#------------------------------------------------------------------------
# Choose the protocol(s) you want to setup:
#------------------------------------------------------------------------
# {setup-} This is need to recognize this script as setup
#------------------------------------------------------------------------
# {section} Global Parameters
#------------------------------------------------------------------------


# Absolute path to the root directory name for this project:
ProjectDir='/home/scheres/test'
# System flavour for parallelization
SystemFlavour=''
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
                     ProjectDir,
                     SystemFlavour,
                     LogDir,
                     AutoLaunch):

            import os,sys

            self.ProjectDir=ProjectDir
            self.SystemFlavour=SystemFlavour
            self.LogDir=LogDir
            self.AutoLaunch=AutoLaunch

            scriptdir = os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
            sys.path.append(scriptdir) # add default search path
            self.SYSTEMSCRIPTDIR=scriptdir

            # Which protocols will appears in main windows
            self.LaunchSections = ["Preprocessing", "2D Analysis", "3D Analysis"]
            self.LaunchButtons = {};
            
            ####### Preprocessing
            self.LaunchButtons['SetupPreProcessMicrographs'] = {'title': 'Preprocess micrographs',
                                       'script':  'xmipp_protocol_preprocess_micrographs.py',
                                       'section': 'Preprocessing'}
            self.LaunchButtons['SetupParticlePick'] = {'title': 'Particle selection',
                                       'script':  'xmipp_protocol_particle_pick.py',
                                       'section': 'Preprocessing'}
            self.LaunchButtons['SetupPreProcessParticles'] = {'title': 'Preprocess particles',
                                       'script':  'xmipp_protocol_preprocess_particles.py',
                                       'section': 'Preprocessing'}           
            
            ######## 2D Analysis
            self.LaunchButtons['SetupAlign'] = {'title': 'Align',
                                       'section': '2D Analysis',
                                       'childs': 'SetupML2D, SetupCL2D'}
            self.LaunchButtons['SetupAlignClassify'] = {'title': 'Align + Classify',
                                       'section': '2D Analysis',
                                       'childs': 'SetupML2D, SetupCL2D'}
            self.LaunchButtons['SetupClassify'] = {'title': 'Classify',
                                       'section': '2D Analysis',
                                       'childs': 'SetupKerDenSOM, SetupRotSpectra'} 
            self.LaunchButtons['SetupML2D'] = {'title': 'ML2D',
                                               'script': 'xmipp_protocol_ml2d.py'}     
            self.LaunchButtons['SetupCL2D'] = {'title': 'CL2D',
                                               'script': 'xmipp_protocol_cl2d.py'}
            self.LaunchButtons['SetupKerDenSOM'] = {'title': 'KerDenSOM',
                                               'script': 'xmipp_protocol_kerdensom.py'}     
            self.LaunchButtons['SetupRotSpectra'] = {'title': 'Rotational Spectra',
                                               'script': 'xmipp_protocol_rotspectra.py'}
            
            ####### 3D Analysis
            self.LaunchButtons['SetupInitialModel'] = {'title': 'Initial Model',
                                       'section': '3D Analysis',
                                       'childs': 'SetupCommonLines, SetupRCT'}
            self.LaunchButtons['SetupCommonLines'] = {'title': 'Common lines',
                                       'script': 'xmipp_protocol_commonlines.py'}
            self.LaunchButtons['SetupRCT'] = {'title': 'Random Conical Tilt',
                                       'script': 'xmipp_protocol_rct.py'} 
            self.LaunchButtons['Setup3DClassify'] = {'title': 'ML3D Classification',
                                       'section': '3D Analysis',
                                       'script': 'xmipp_protocol_ml3d.py'}
            self.LaunchButtons['SetupModelRefine'] = {'title': 'Model refinement',
                                       'section': '3D Analysis',
                                       'script': 'xmipp_protocol_projmatch.py'} 

            # For automated editing of default directories in protocols
            self.DEFAULTDIRS={"ProjectDir":self.ProjectDir,
                              "LogDir":self.LogDir,
                              "SystemFlavour":self.SystemFlavour
                              }
            

            # Perform the actual setup:
            if (self.AutoLaunch != ""):
                # A. Setup from GUI (Autolaunch)
                # This will copy the (modified) protocol script to the corresponding directory
                # and will automatically launch the GUI for this protocol
                self.setup_protocol(self.AutoLaunch)
            else:
                # B. Setup from this script:
                # This will only copy the (modified) protocol script to the corresponding directory
                for buttonData in self.LaunchButtons:
                    if 'script' in buttonData:
                        self.setup_protocol(buttonData['script'])
                        
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
                # Only change DEFAULTDIRS if entry is not empty in the setup protocol
                if (len(self.DEFAULTDIRS[dir]) > 0):
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

        def setup_protocol(self,protocol):
            import os
            import shutil

            os.chdir(self.ProjectDir)
                
            visualize=protocol.replace('xmipp_protocol','visualize')
            myscripts=[protocol, visualize]

            for script in myscripts:
                src=str(self.SYSTEMSCRIPTDIR)+"/"+str(script)
                dst=str(script)
                if os.path.exists(dst):
                   src=dst
                   print "* File "+dst+" already existed (now updated)"
                elif (not os.path.exists(src)):
                   src=str(self.SYSTEMSCRIPTDIR)+"/not_implemented.py"

                text=self.modify_script_header(src)
                fh=open(dst,'w')
                fh.writelines(text)
                fh.close()
                os.chmod(dst,0755)

            # Only auto-launch the first script
            if (self.AutoLaunch!=""):
                command='python '+str(self.SYSTEMSCRIPTDIR)+'/xmipp_protocol_gui.py '+str(protocol)+' &'
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

    setup = setup_protocols_class(ProjectDir, SystemFlavour, LogDir, AutoLaunch)

