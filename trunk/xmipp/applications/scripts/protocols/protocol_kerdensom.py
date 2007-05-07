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
# Working subdirectory:
WorkingDir="SOM/ML2ref_ref1"
# Delete working subdirectory if it already exists?
DoDeleteWorkingDir=False
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project
"""
ProjectDir="/home2/bioinfo/scheres/work/protocols"
# {expert} Directory name for logfiles (relative path from ProjectDir):
""" All logfiles will be stored here
"""
LogDir="Logs"
#------------------------------------------------------------------------------------------------
# {section} ml_align2d parameters
#------------------------------------------------------------------------------------------------
# {dir} Directory where you have previously ran ML2D classification:
ML2DWorkingDir="ML2D/ML2ref"
# The number of the class to use:
ML2DReferenceNr=1
#------------------------------------------------------------------------------------------------
# {section} Mask parameters
#------------------------------------------------------------------------------------------------
# Design your mask graphically?
DoXmask=True
""" This will launch a graphical program to design your own mask.
    Be careful NOT to submit your job via a queueing system!
"""
# {file} OR provide an already existing mask:
MaskFileName=""
#------------------------------------------------------------------------------------------------
# {section} classify_kerdensom parameters
#------------------------------------------------------------------------------------------------
# Perform self-organizing map calculation?
DoSOM=True
# Name of output map:
""" Existing files with this name will be delete!
"""
SomName="som"
# X-dimension of the self-organizing map:
SomXdim=10
# Y-dimension of the self-organizing map:
SomYdim=5
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
# {expert} Additional classify_kerdensom parameters:
""" For a complete description see http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/KerDenSOM
"""
KerdensomExtraParams=""
#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript="visualize_kerdensom.py"
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class kerdensom_class:

    #init variables
    def __init__(self,
                 _WorkingDir,
                 _DoDeleteWorkingDir,
                 _ProjectDir,
                 _LogDir,
                 _ML2DWorkingDir,
                 _ML2DReferenceNr,
                 _DoXmask,
                 _MaskFileName,
                 _DoSOM,
                 _SomName,
                 _SomXdim,
                 _SomYdim,
                 _SomReg0,
                 _SomReg1,
                 _SomSteps,
                 _KerdensomExtraParams):

        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=_WorkingDir
        self.ProjectDir=_ProjectDir
        self.ML2DWorkingDir=os.path.abspath(_ML2DWorkingDir)
        self.ML2DReferenceNr=_ML2DReferenceNr
        self.DoXmask=_DoXmask
        self.MaskFileName=_MaskFileName
        self.DoSOM=_DoSOM
        self.SomName=_SomName
        self.SomXdim=_SomXdim
        self.SomYdim=_SomYdim
        self.SomReg0=_SomReg0
        self.SomReg1=_SomReg1
        self.SomSteps=_SomSteps
        self.KerdensomExtraParams=_KerdensomExtraParams
        
        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Delete working directory if it exists, make a new one
        if (_DoDeleteWorkingDir): 
            if os.path.exists(self.WorkingDir):
                shutil.rmtree(self.WorkingDir)
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
                                       os.path.abspath(self.WorkingDir))
        os.chdir(self.WorkingDir)

        # Check whether ML2D subdirectory exists
        if not os.path.exists(self.ML2DWorkingDir):
            message=" Error: ML2D directory "+self.ML2DWorkingDir+" does not exists!"
            print '* '+message
            self.log.error(message)
            sys.exit()
            
        # Execute kerdensom protocol in the working directory
        self.execute_whole_protocol()

        # Return to parent dir
        os.chdir(os.pardir)


    def make_local_copy_of_images(self):

        import sys,os,glob
        import selfile

        # Make a directory for the local copies of all relevant images
        if not os.path.exists('local_images'):
            os.makedirs('local_images')

        # Check whether the ML2D run has written docfiles already
        docfiles=glob.glob(self.ML2DWorkingDir+'/*_it?????.doc')
        if len(docfiles)==0:
            message='No ML2D docfiles yet. Continue script after ML2D job completion... '
            print '* ',message
            self.log.error(message)
            sys.exit()
        else:
            # Copy relevant selfile and images of ML2DDir to WorkingDir/local_images
            lastitername=docfiles[-1].replace('.doc','')
            splits=lastitername.split('_it')
            ml2d_abs_rootname=splits[0]
            self.classselfile=ml2d_abs_rootname+'_ref'+str(self.ML2DReferenceNr).zfill(5)+'.sel'
            # Make a local copy of the images
            message='Making a local copy of the images in '+str(self.classselfile)
            print '* ',message
            self.log.info(message)
            mysel=selfile.selfile()
            mysel.read(self.classselfile)
            newsel=mysel.copy_sel('local_images')
            self.classselfile=os.path.basename(self.classselfile)
            newsel.write(self.classselfile)
            # Make a local copy of the (partial) docfile
            self.make_local_copy_docfile(mysel,newsel)

    def make_local_copy_docfile(self,oldsel,newsel):

        import os,glob
        docfiles=glob.glob(self.ML2DWorkingDir+'/*_it?????.doc')
        docfile=docfiles[-1]
        fh=open(docfile,'r')
        doclines=fh.readlines()
        newdoc=[]
        newdoc.append(doclines[0])
        
        for name,state in oldsel.sellines:
            for i in range(len(doclines)):
                if name in doclines[i]:
                    splits=doclines[i].split()
                    name=os.path.basename(splits[1])
                    newname=' ; local_images/'+name+'\n'
                    newdoc.append(newname)
                    newdoc.append(doclines[i+1])
                    i=i+1
                    break

        self.docfilename='ml2d_ref'+str(self.ML2DReferenceNr).zfill(5)+'.doc'
        fh=open(self.docfilename,'w')
        fh.writelines(newdoc)
        fh.close()

    def assign_header(self):
        import os
        import glob

        command='xmipp_header_assign -i '+self.docfilename+' -mirror \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def execute_xmask(self,selfile):
        import os

        self.MaskFileName='mask_design.msk'
        command='xmipp_mask_design -sel '+str(selfile) + ' -save_as mask_design.msk \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def execute_img2data(self,selfile,maskname,outname):
        import os
        if (maskname==''):
            command='xmipp_convert_img2data -i '+ str(selfile) + \
                     ' -nomask '+ \
                     ' -o ' + str(outname) + ' \n'
        else:
            command='xmipp_convert_img2data -i '+ str(selfile) + \
                     ' -mask '+ str(maskname) + \
                     ' -o ' + str(outname) + ' \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def execute_data2img(self,iname,maskname):
        import os
        if (maskname==''):
            command='xmipp_convert_data2img -i '+ str(iname) + \
                     ' -nomask \n'
        else:
            command='xmipp_convert_data2img -i '+ str(iname) + \
                     ' -mask '+ str(maskname) + ' \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def delete_existing_som(self,somname):
        import os
        # delete existing files with this somname
        if (os.path.exists(somname+'.sel')):
            print 'Deleting existing som...'
            command= 'xmipp_selfile_delete '+somname+'.sel \n'
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
        command= 'xmipp_classify_kerdensom'+ ' -i ' + str(datafile) + \
                 ' -o '    + str(somname)  + \
                 ' -xdim ' + str(self.SomXdim) + \
                 ' -ydim ' + str(self.SomYdim) + \
                 ' -reg0 ' + str(self.SomReg0) + \
                 ' -reg1 ' + str(self.SomReg1) + \
                 ' -steps ' + str(self.SomSteps) + \
                 ' '  + str(self.KerdensomExtraParams) + ' -verb 1 \n'
        print '* ',command
        self.log.info(command)
        os.system(command)

    def execute_whole_protocol(self):
        
        import os,glob

        self.make_local_copy_of_images()
        self.assign_header()
      
        if (self.DoXmask):
            self.execute_xmask(self.classselfile)

        if (self.DoSOM):
            self.delete_existing_som(self.SomName)
            self.execute_img2data(self.classselfile,self.MaskFileName,'data.dat')
            self.execute_kerdensom('data.dat',self.SomName)
            self.execute_data2img(str(self.SomName)+'.cod',self.MaskFileName)

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'
        
#		
# Main
#     
if __name__ == '__main__':

    # create kerdensom_class object

    kerdensom=kerdensom_class(WorkingDir,
                              DoDeleteWorkingDir,
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
                              KerdensomExtraParams)
    # close 
    kerdensom.close()

