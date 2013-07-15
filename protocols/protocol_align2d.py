#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment,
#
# Example use:
# ./xmipp_protocol_align2d.py
#
# Author: Carlos Oscar Sanchez Sorzano, January 2010
#
# {begin_of_header}
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# {file} Input images:
""" This selfile points to the spider single-file format images that make up your data set. The filenames can have relative or absolute paths, but it is strictly necessary that you put this selfile IN THE PROJECTDIR. 
"""
InSelFile='all_images.sel'
# Working subdirectory:
""" This directory will be created if it doesn't exist, and will be used to store all output from this run. Don't use the same directory for multiple different runs, instead use a structure like run1, run2 etc. 
"""
WorkingDir='Align2D/run1'
# {expert} Root directory name for this project:
""" Absolute path to the root directory for this project. Often, each data set of a given sample has its own ProjectDir.
"""
ProjectDir='/gpfs/fs1/home/bioinfo/coss/analu/Mcm467_cdt1_sinOli-7AporPix'
# {expert} Directory name for logfiles:
""" All logfiles will be stored here
"""
LogDir='Logs'

#------------------------------------------------------------------------------------------------
# {section} Alignment parameters
#------------------------------------------------------------------------------------------------
# {file} Reference image (optional)
""" If not given, the reference is the average of all the input images. 
"""
ReferenceImage=""

# {expert} Number of iterations
""" Maximum number of iterations within each level
"""
NumberOfIterations=15

# {expert} Band pass filter
""" Apply a band pass filter before clustering """
DoFilter =True

# {expert} Highpass cutoff frequency
""" In (Angstroms/Pixel). Set to 0 if not desired """
Highpass =0.02

# {expert} Lowpass cutoff frequency
""" In (Angstroms/Pixel). Set to 0 if not desired """
Lowpass =0.4

#------------------------------------------------------------------------------------------------
# {section} Parallelization issues
#------------------------------------------------------------------------------------------------
# Number of MPI processes to use:
NumberOfMpi=40

# MPI system Flavour 
""" Depending on your queuing system and your mpi implementation, different mpirun-like commands have to be given.
    Ask the person who installed your xmipp version, which option to use. 
    Or read: http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/ParallelPage. 
"""
SystemFlavour='TORQUE-OPENMPI'

# {hidden} This protocol only works in parallel
DoParallel=True

#------------------------------------------------------------------------------------------------
# {expert} Analysis of results
""" This script serves only for GUI-assisted visualization of the results
"""
AnalysisScript='visualize_align2d.py'
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end_of_header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

import os,sys,shutil,time
scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols', 'r').read()))[0] + '/protocols'
sys.path.append(scriptdir)

def stepPerformed(step,filename):
    import re
    f = open(filename, 'r')
    lines=f.readlines()
    f.close()
    expr = re.compile(step)
    return len(filter(expr.search,lines))>0

class Align2D_class:
    def saveAndCompareParameters(self, listOfParameters):
        fnOut=self.WorkingDir + "/protocolParameters.txt"
        linesNew=[];
        for prm in listOfParameters:
            eval("linesNew.append('"+prm +"='+str("+prm+")+'\\n')")
        if os.path.exists(fnOut):
            f = open(fnOut, 'r')
            linesOld=f.readlines()
            f.close()
            same=True;
            if len(linesOld)==len(linesNew):
                for i in range(len(linesNew)):
                    if not linesNew[i]==linesOld[i]:
                        same=False
                        break;
            else:
                same=False
            if not same:
                print("Deleting")
                self.log.info("Deleting working directory since it is run with different parameters")
                shutil.rmtree(self.WorkingDir)
                os.makedirs(self.WorkingDir)
        f = open(fnOut, 'w')
        f.writelines(linesNew)
        f.close()

    #init variables
    def __init__(self,
                 InSelFile,
                 WorkingDir,
                 ProjectDir,
                 LogDir,
                 ReferenceImage,
                 NumberOfIterations,
                 DoFilter,
                 Highpass,
                 Lowpass,
                 NumberOfMpi,
                 SystemFlavour):
	     
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log

        self.WorkingDir=WorkingDir
        self.ProjectDir=ProjectDir
        self.InSelFile=InSelFile
        self.ReferenceImage=ReferenceImage
        self.NumberOfIterations=NumberOfIterations
        self.DoFilter=DoFilter
        self.Highpass=Highpass
        self.Lowpass=Lowpass
        self.NumberOfMpi=NumberOfMpi
        self.SystemFlavour=SystemFlavour
   
        # Setup logging
        self.log=log.init_log_system(self.ProjectDir,
                                     LogDir,
                                     sys.argv[0],
                                     self.WorkingDir)
                
        # Create directory if does not exist
        if not os.path.exists(self.WorkingDir):
            os.makedirs(self.WorkingDir)
        self.doStep1=True
        self.doStep2=True
                
        # Save parameters and compare to possible previous runs
        self.saveAndCompareParameters([
                 "InSelFile",
                 "ReferenceImage",
                 "NumberOfIterations",
                 "DoFilter",
                 "Highpass",
                 "Lowpass"]);

        # Backup script
        log.make_backup_of_script_file(sys.argv[0],
            os.path.abspath(self.WorkingDir))

        # Update status
        fh=open(self.WorkingDir + "/status.txt", "a")
        fh.write("Step 0: Process started at " + time.asctime() + "\n")
        fh.close()

        # Run
        self.preprocess()
        self.execute_CLalign2D()
        self.postprocess()
        
        fh=open(self.WorkingDir + "/status.txt", "a")
        fh.write("Step F: Process finished at " + time.asctime() + "\n")
        fh.close()

    def preprocess(self):
        import launch_job,xmipp
        if self.DoFilter:
            slope=max(self.Highpass/2,0.01)
            fnOut=self.WorkingDir+'/preprocessedImages.stk'
            self.selFileToUse=fnOut
            if not self.doStep1 and os.path.exists(fnOut):
                numberOfImages=xmipp.getImageSize(fnOut)[3]
                MD=xmipp.MetaData(self.InSelFile)
                if MD.size()!=numberOfImages:
                    self.doStep1=True
            if not stepPerformed("Step 1",self.WorkingDir + "/status.txt"):
                self.doStep1=True
            if self.doStep1:
                params= '-i '+str(self.InSelFile)+\
                        ' --fourier_mask raised_cosine '+str(slope)+\
                        ' -o '+fnOut
                if self.Highpass>0 and self.Lowpass>0:
                    params+=" --band_pass "+str(self.Highpass)+" "+str(self.Lowpass)
                elif self.Highpass>0:
                    params+=" --high_pass "+str(self.Highpass)
                elif self.Lowpass>0:
                    params+=" --low_pass "+str(self.Lowpass)
                launchJob("xmipp_fourier_filter",
                                      params,
                                      self.log,
                                      False,
                                      1,
                                      1,
                                      self.SystemFlavour)
                # Update status    
                if os.path.exists(self.WorkingDir+'/preprocessedImages.stk'):
                    fh=open(self.WorkingDir + "/status.txt", "a")
                    fh.write("Step 1: Preprocessing finished at " + time.asctime() + "\n")
                    fh.close()
        else:
            self.selFileToUse=self.InSelFile

    def execute_CLalign2D(self):
        import launch_job,xmipp
        if self.DoFilter and not stepPerformed("Step 1",self.WorkingDir + "/status.txt"):
            return
        if not stepPerformed("Step 2",self.WorkingDir + "/status.txt"):
            self.doStep2=True
        if not self.doStep2:
            if not os.path.exists(WorkingDir+'/class.sel'):
                self.doStep2=True
            else:
                return
        params= '-i '+str(self.selFileToUse)+' -o '+WorkingDir+'/class '+\
                ' -codes 1'+\
                ' -iter '+str(self.NumberOfIterations)
        if self.ReferenceImage!="":
            MD=xmipp.MetaData();
            MD.addObject()
            MD.setValue(xmipp.MDL_IMAGE,self.WorkingDir+"/reference.xmp")
            MD.write(self.WorkingDir+"/reference.sel")
            os.system("xmipp_convert_image -i "+self.ReferenceImage+" -o "+self.WorkingDir+"/reference.xmp")
            params+=" -codesSel0 "+self.WorkingDir+"/reference.sel"
        else:
            params+=' -codes0 1'

        launchJob("xmipp_classify_CL2D",
                              params,
                              self.log,
                              True,
                              self.NumberOfMpi,
                              1,
                              self.SystemFlavour)
        
        # If the images have been filtered, translate the output of cl2d in
        # terms of the original images
        # Construct dictionary
        fnClasses=WorkingDir+'/class.sel'
        if os.path.exists(fnClasses) and self.DoFilter:
            import xmipp
            xmipp.substituteOriginalImages(fnClasses,self.InSelFile,
                                           self.WorkingDir+"/intermediate.sel",
                                           xmipp.MDL_IMAGE_ORIGINAL,False)
            os.system("mv -f "+self.WorkingDir+"/intermediate.sel "+fnClasses)
        if self.ReferenceImage!="":
            os.remove(self.WorkingDir+"/reference.xmp")
            os.remove(self.WorkingDir+"/reference.sel")
        if os.path.exists(fnClasses) and self.DoFilter:
            fh=open(self.WorkingDir + "/status.txt", "a")
            fh.write("Step 2: CL2D finished at " + time.asctime() + "\n")
            fh.close()

    def postprocess(self):
        if not stepPerformed("Step 2",self.WorkingDir + "/status.txt"):
            return
        if not self.DoFilter or not os.path.exists(self.WorkingDir+'/preprocessedImages.stk'):
            return
        
        print "Analyzing alignment ----------------------------"
        
        # Remove useless files
        import xmipp,glob
        fnClasses=glob.glob(WorkingDir+"/class_level_??.*")
        for file in fnClasses:
            os.remove(file)
        os.remove(self.WorkingDir+'/preprocessedImages.stk')
        os.system('mv -f '+self.WorkingDir+'/class_aligned.stk '+\
                  self.WorkingDir+'/images_aligned.stk')
        
        # Process the selfile
        MD=xmipp.MetaData()
        MD.readBlock(self.WorkingDir+'/class.sel','class_%06d'%0)
        MDout=xmipp.MetaData()
        for id in MD:
            fnImg=MD.getValue(xmipp.MDL_IMAGE)
            fnImgOrig=MD.getValue(xmipp.MDL_IMAGE_ORIGINAL)
            fnImg=fnImg.replace("class_aligned","images_aligned")
            MDout.addObject()
            MDout.setValue(xmipp.MDL_IMAGE,fnImg)
            MDout.setValue(xmipp.MDL_IMAGE_ORIGINAL,fnImgOrig)
        os.remove(self.WorkingDir+'/class.sel')
        MDout.write(self.WorkingDir+'/alignment.sel')
        
        # Calculate the alignment statistics (mean, stddev, PCAbasis)
        # Falta calcular el PCA de las imagenes alineadas
        command="xmipp_classify_analyze_cluster"+\
                " -i "+self.WorkingDir+'/alignment.sel'+\
                " --ref 0@"+self.WorkingDir+'/class.stk'+\
                " -o "+self.WorkingDir+'/pca_zscore.sel'+\
                " --basis "+self.WorkingDir+'/alignment_stats.stk'+\
                " --maxDist -1 --quiet"
        os.system(command)
        os.remove(self.WorkingDir+'/class.stk')

        # Now compute the multivariate zscore
        command="xmipp_sort_by_statistics"+\
                " -i "+self.WorkingDir+'/alignment.sel'+\
                " -o "+self.WorkingDir+'/multivariate_zscore'+\
                " --multivariate --quiet"
        os.system(command)

        # Merge all results
        MD=xmipp.MetaData(self.WorkingDir+'/pca_zscore.sel')
        PCAzscore={}
        for id in MD:
            PCAzscore[MD.getValue(xmipp.MDL_IMAGE)]=MD.getValue(xmipp.MDL_ZSCORE)
        MD=xmipp.MetaData(self.WorkingDir+'/multivariate_zscore.sel')
        Multivariatezscore={}
        for id in MD:
            Multivariatezscore[MD.getValue(xmipp.MDL_IMAGE)]=MD.getValue(xmipp.MDL_ZSCORE)
        MD=xmipp.MetaData(self.WorkingDir+'/alignment.sel')
        for id in MD:
            fnImg=MD.getValue(xmipp.MDL_IMAGE)
            MD.setValue(xmipp.MDL_ZSCORE,(PCAzscore[fnImg]+Multivariatezscore[fnImg])/2)
        MD.write(self.WorkingDir+'/alignment.sel')
        os.remove(self.WorkingDir+'/pca_zscore.sel')
        os.remove(self.WorkingDir+'/multivariate_zscore.sel')
        
        fh=open(self.WorkingDir + "/status.txt", "a")
        fh.write("Step 4: Postprocessing finished at " + time.asctime() + "\n")
        fh.close()

# Preconditions
def checkErrors():
    errors = []
    # Check if there is workingdir
    if WorkingDir == "":
        errors.append("No working directory given")
    # Check that there are any micrograph to process
    if not os.path.exists(InSelFile):
        errors.append("The input selfile is not valid")
    # Check filter parameters
    if DoFilter and (Highpass<0 or Highpass>0.5 or Lowpass<0 or Lowpass>0.5):
        errors.append("The filter frequencies must be between 0 and 0.5")

    return errors

#		
# Main
#     
if __name__ == '__main__':
    Align2D=Align2D_class(
                 InSelFile,
                 WorkingDir,
                 ProjectDir,
                 LogDir,
                 ReferenceImage,
                 NumberOfIterations,
                 DoFilter,
                 Highpass,
                 Lowpass,
                 NumberOfMpi,
                 SystemFlavour)
