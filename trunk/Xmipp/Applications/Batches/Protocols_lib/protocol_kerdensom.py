#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for image classification using self-organizing maps,
# after a preliminary alignment (and classification) using MLalign2D
#
# Example use:
# ./protocol_kerdensom.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Working directory:
WorkingDir="SOM_ML1ref_ref1"
# Batch submission command (use "" to launch without batch submission):
""" This will depend on your queueing system., ask your system administrator...
    If you dont use a queueing system, type: LaunchParallelScript=""
"""
LaunchJobCommand="" 
# {expert} Root directory name for this project:
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles:
""" All logfiles will be stored in $ProjectDir/$LogDir
"""
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} MLalign2D parameters
#------------------------------------------------------------------------------------------------
# Directory where you have ran MLalign2D:
ML2DWorkingDir="ML2ref"
# The number of the class to use:
ML2DReferenceNr=1
#------------------------------------------------------------------------------------------------
# {section} Mask parameters
#------------------------------------------------------------------------------------------------
# Design your mask graphically? (Save as name below!)
DoXmask=True
# Name of the mask:
MaskFileName="mask.msk"
#------------------------------------------------------------------------------------------------
# {section} kerdenSOM parameters
#------------------------------------------------------------------------------------------------
# Perform self-organizing map calculation?
DoSOM=False
# Name of Output SOM:
""" Existing files with this name will be deleted!
"""
SomName="som"
# X-dimension of the self-organizing map:
SomXdim=7
# Y-dimension of the self-organizing map:
SomYdim=7
# Initial regularization factor:
""" The kerdenSOM algorithm anneals from an initial high regularization factor
    to a final lower one, in a user-defined number of steps.
    If the output map is too smooth, lower the regularization factors
    If the output map is not organized, higher the regularization factors
"""
SomReg0=1000
# Final regularization factor:
SomReg1=200
# Number of steps to lower the regularization factor:
SomSteps=5
# {expert} Additional kerdenSOM parameters:
""" For a complete description see http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraParams=""
#------------------------------------------------------------------------------------------------
# {section} Analysis of results
#------------------------------------------------------------------------------------------------
# Visualize the SOM?
DoVisualizeSOM=False
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class kerdensom_class:

    #init variables
    def __init__(self,
                 WorkingDir,
                 LaunchJobCommand,
                 ProjectDir,
                 LogDir,
                 ML2DWorkingDir,
                 ML2DReferenceNr,
                 DoXmask,
                 MaskFileName,
                 DoSOM,
                 SomName,
                 SomXdim,
                 SomYdim,
                 SomReg0,
                 SomReg1,
                 SomSteps,
                 KerdensomExtraParams,
                 DoVisualizeSOM ):

        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.LaunchJobCommand=LaunchJobCommand
        self.ML2DWorkingDir='../'+ML2DWorkingDir
        self.ML2DReferenceNr=ML2DReferenceNr
        self.DoXmask=DoXmask
        self.MaskFileName=MaskFileName
        self.DoSOM=DoSOM
        self.SomName=SomName
        self.SomXdim=SomXdim
        self.SomYdim=SomYdim
        self.SomReg0=SomReg0
        self.SomReg1=SomReg1
        self.SomSteps=SomSteps
        self.KerdensomExtraParams=KerdensomExtraParams
        self.DoVisualizeSOM=DoVisualizeSOM
        
        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Delete working directory if it exists, make a new one, and go there
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Execute MLalign2D in the working directory
        os.chdir(self.WorkingDir)
        self.execute_whole_protocol()

        # Return to parent dir
        os.chdir(os.pardir)


    def assign_header(self):
        import os
        import glob

        docfiles=glob.glob(self.ML2DWorkingDir+'/ml2d_it?????.doc')
        docfile=docfiles[-1]
        command='xmipp_headerinfo -assign -i '+docfile+' -mirror \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def reset_header(self):
        import os
        import glob
        
        selfiles=glob.glob(self.ML2DWorkingDir+'/ml2d_ref?????.sel')
        allimages=[]
        for selfile in selfiles:
            fh=open(selfile,'r')
            allimages+=fh.readlines()
            fh.close()
        fh=open('tmp.sel','w')
        fh.writelines(allimages)
        fh.close()
        command='xmipp_headerinfo -reset -i tmp.sel'
        print '* ',command
        self.log.info(command)
        os.system(command)
        command='rm -f tmp.sel'
        print '* ',command
        self.log.info(command)
        os.system(command)
     
        
    def execute_xmask(self,selfile):
        import os

        command='xmipp_xmask -sel '+str(selfile) + ' \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def execute_img2data(self,selfile,maskname,outname):
        import os
        command='xmipp_img2data -sel '+ str(selfile) + \
                 ' -mname '+ str(maskname) + \
                 ' -fname '+ str(outname) + ' \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def execute_data2img(self,iname,maskname):
        import os
        command='xmipp_data2img -iname '+ str(iname) + \
                 ' -mname '+ str(maskname) + ' \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def delete_existing_som(self,somname):
        import os
        # delete existing files with this somname
        if (os.path.exists(somname+'.sel')):
            print 'Deleting existing som...'
            command= 'xmipp_rmsel '+somname+'.sel \n'
            print '* ',command
            self.log.info(command)
            os.system(command)
            command= 'rm -f '+somname+'.* '+somname+'_* \n'
            print '* ',command
            self.log.info(command)
            os.system(command)

    def execute_kerdensom(self,datafile,somname):
        import os
        
        print '*********************************************************************'
        print '*  Executing kerdenSOM program :' 
        command= 'xmipp_kerdensom'+ ' -din ' + str(datafile) + \
                 ' -cout ' + str(somname)  + \
                 ' -xdim ' + str(self.SomXdim) + \
                 ' -ydim ' + str(self.SomYdim) + \
                 ' -reg0 ' + str(self.SomReg0) + \
                 ' -reg1 ' + str(self.SomReg1) + \
                 ' -steps ' + str(self.SomSteps) + \
                 ' '  + str(self.KerdensomExtraParams) + ' -verb 1 \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

                                               
    def visualize_SOM(self,somname):
         import os

         command='xmipp_show -som ' + str(somname)
         print '* ',command
         self.log.info(command)
         os.system(command)

    def execute_whole_protocol(self):
        
        import os

        self.assign_header()
      
        classselfile=self.ML2DWorkingDir+'/ml2d_ref'+ str(self.ML2DReferenceNr).zfill(5) +'.sel'

        if (self.DoXmask):
            self.execute_xmask(classselfile)

        if (self.DoSOM):
            self.delete_existing_som(self.SomName)
            self.execute_img2data(classselfile,self.MaskFileName,'data.dat')
            self.execute_kerdensom('data.dat',self.SomName)
            self.execute_data2img(str(self.SomName)+'.cod',self.MaskFileName)

        if (self.DoVisualizeSOM):
            self.visualize_SOM(self.SomName)
            
        self.reset_header()

          
    def close(self):
        print '*********************************************************************'
#		
# Main
#     
if __name__ == '__main__':

    # create kerdensom_class object

    kerdensom=kerdensom_class(WorkingDir,
                              LaunchJobCommand,
                              ProjectDir,
                              LogDir,
                              ML2DWorkingDir,
                              ML2DReferenceNr,
                              DoXmask,
                              MaskFileName,
                              DoSOM,
                              SomName,
                              SomXdim,
                              SomYdim,
                              SomReg0,
                              SomReg1,
                              SomSteps,
                              KerdensomExtraParams,
                              DoVisualizeSOM)
    # close 
    kerdensom.close()

