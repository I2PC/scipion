#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using hierarchical clustering principles
#
# Example use:
# ./xmipp_protocol_cl2d.py
#
# Author: Carlos Oscar Sanchez Sorzano, July 2009
#

import os,sys,shutil,time
from protlib_base import *
from config_protocols import protDict

def getParameter(prm,filename):
    f = open(filename, 'r')
    lines=f.readlines()
    f.close()
    for line in lines:
        tokens=line.split('=')
        if tokens[0]==prm:
            return tokens[1].strip()
    return ""

def stepPerformed(step,filename):
    import re
    f = open(filename, 'r')
    lines=f.readlines()
    f.close()
    expr = re.compile(step)
    return len(filter(expr.search,lines))>0

class ProtCL2D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.cl2d.key, scriptname, project)
        self.Import = 'from protocol_cl2d import *'    
        
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
#    def __init__(self,
#                 InSelFile,
#                 WorkingDir,
#                 ProjectDir,
#                 LogDir,
#                 NumberOfReferences,
#                 NumberOfReferences0,
#                 NumberOfIterations,
#                 DoFilter,
#                 Highpass,
#                 Lowpass,
#                 ComparisonMethod,
#                 ClusteringMethod,
#                 AdditionalParameters,
#                 thGoodClass,
#                 thZscore,
#                 thPCAZscore,
#                 NumberOfMpiProcesses,
#                 SystemFlavour):
#	     
#        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
#        sys.path.append(scriptdir) # add default search path
#        import log
#
#        self.WorkingDir=WorkingDir
#        self.ProjectDir=ProjectDir
#        self.InSelFile=InSelFile
#        self.NumberOfReferences=NumberOfReferences
#        self.NumberOfReferences0=NumberOfReferences0
#        self.NumberOfIterations=NumberOfIterations
#        self.DoFilter=DoFilter
#        self.Highpass=Highpass
#        self.Lowpass=Lowpass
#        self.ComparisonMethod=ComparisonMethod
#        self.ClusteringMethod=ClusteringMethod
#        self.AdditionalParameters=AdditionalParameters
#        self.thGoodClass=thGoodClass
#        self.thZscore=thZscore
#        self.thPCAZscore=thPCAZscore
#        self.NumberOfMpiProcesses=NumberOfMpiProcesses
#        self.SystemFlavour=SystemFlavour
#   
#        # Setup logging
#        self.log=log.init_log_system(self.ProjectDir,
#                                     LogDir,
#                                     sys.argv[0],
#                                     self.WorkingDir)
#                
#        # Create directory if does not exist
#        if not os.path.exists(self.WorkingDir):
#            os.makedirs(self.WorkingDir)
#            self.doStep1=True
#            self.doStep2=True
#            self.doStep3=True
#        else:
#            fnParam=self.WorkingDir + "/protocolParameters.txt"
#            oldInSelFile=getParameter("InSelFile",fnParam)
#            oldDoFilter=getParameter("DoFilter",fnParam)
#            oldHighpass=getParameter("Highpass",fnParam)
#            oldLowpass=getParameter("Lowpass",fnParam)
#
#            oldNumberOfReferences=getParameter("NumberOfReferences",fnParam)
#            oldNumberOfReferences0=getParameter("NumberOfReferences0",fnParam)
#            oldNumberOfIterations=getParameter("NumberOfIterations",fnParam)
#            oldComparisonMethod=getParameter("ComparisonMethod",fnParam)
#            oldClusteringMethod=getParameter("ClusteringMethod",fnParam)
#            oldAdditionalParameters=getParameter("AdditionalParameters",fnParam)
#            
#            oldthGoodClass=getParameter("thGoodClass",fnParam)
#            oldthZscore=getParameter("thZscore",fnParam)
#            oldthPCAZscore=getParameter("thPCAZscore",fnParam)
#            self.doStep1=False
#            self.doStep2=False
#            self.doStep3=False
#            if oldInSelFile!=InSelFile or str(DoFilter)!=oldDoFilter or \
#               oldHighpass!=str(Highpass) or oldLowpass!=str(Lowpass):
#                self.doStep1=True
#                self.doStep2=True
#                self.doStep3=True
#            elif oldNumberOfReferences!=str(NumberOfReferences) or \
#               oldNumberOfReferences0!=str(NumberOfReferences0) or \
#               oldNumberOfIterations!=str(NumberOfIterations) or \
#               oldComparisonMethod!=ComparisonMethod or \
#               oldClusteringMethod!=ClusteringMethod or \
#               oldAdditionalParameters!=AdditionalParameters:
#                self.doStep2=True
#                self.doStep3=True
#            elif oldthGoodClass!=str(thGoodClass) or oldthZscore!=str(thZscore) or \
#               oldthPCAZscore!=str(thPCAZscore):
#                self.doStep3=True
#                
#        # Save parameters and compare to possible previous runs
#        self.saveAndCompareParameters([
#                 "InSelFile",
#                 "NumberOfReferences",
#                 "NumberOfReferences0",
#                 "NumberOfIterations",
#                 "DoFilter",
#                 "Highpass",
#                 "Lowpass",
#                 "ComparisonMethod",
#                 "ClusteringMethod",
#                 "AdditionalParameters",
#                 "thGoodClass",
#                 "thZscore",
#                 "thPCAZscore"]);
#
#        # Backup script
#        log.make_backup_of_script_file(sys.argv[0],
#            os.path.abspath(self.WorkingDir))
#
#        # Update status
#        fh=open(self.WorkingDir + "/status.txt", "a")
#        fh.write("Step 0: Process started at " + time.asctime() + "\n")
#        fh.close()
#
#        # Run
#        self.preprocess()
#        self.execute_CLalign2D()
#        self.execute_core_analysis()
#        self.postprocess()
#        
#        fh=open(self.WorkingDir + "/status.txt", "a")
#        fh.write("Step F: Process finished at " + time.asctime() + "\n")
#        fh.close()

    def preprocess(self):
        import launch_job,xmipp
        if self.DoFilter:
            slope=max(self.Highpass/2,0.01)
            fnOut=self.WorkingDir+'/preprocessedImages.stk'
            self.selFileToUse=fnOut
            if not self.doStep1 and os.path.exists(fnOut):
                numberOfImages=xmipp.SingleImgSize(fnOut)[3]
                MD=xmipp.MetaData(self.InSelFile)
                if MD.size()!=numberOfImages:
                    self.doStep1=True
            if not stepPerformed("Step 1",self.WorkingDir + "/status.txt"):
                self.doStep1=True
            if self.doStep1:
                params= '-i '+str(self.InSelFile)+\
                        ' -o '+fnOut+\
                        ' --fourier '
                if self.Highpass>0 and self.Lowpass>0:
                    params+=" band_pass "+str(self.Highpass)+" "+str(self.Lowpass)
                elif self.Highpass>0:
                    params+=" high_pass "+str(self.Highpass)
                elif self.Lowpass>0:
                    params+=" low_pass "+str(self.Lowpass)
                params+=' raised_cosine '+str(slope)
                launchJob("xmipp_transform_filter",
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
        import launch_job
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
                ' -codes '+str(self.NumberOfReferences)+\
                ' -codes0 '+str(self.NumberOfReferences0)+\
                ' -iter '+str(self.NumberOfIterations)
        params+=' '+self.AdditionalParameters
        if (self.ComparisonMethod=='correlation'):
            params+= ' -useCorrelation '
        if (self.ClusteringMethod=='classical'):
            params+= ' -classicalMultiref '

        launchJob("xmipp_classify_CL2D",
                              params,
                              self.log,
                              True,
                              self.NumberOfMpiProcesses,
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
        if os.path.exists(fnClasses) and self.DoFilter:
            fh=open(self.WorkingDir + "/status.txt", "a")
            fh.write("Step 2: CL2D finished at " + time.asctime() + "\n")
            fh.close()

    def execute_core_analysis(self):
        import launch_job
        if not stepPerformed("Step 2",self.WorkingDir + "/status.txt"):
            return
        if not stepPerformed("Step 3",self.WorkingDir + "/status.txt"):
            self.doStep3=True
        if not self.doStep3:
            if not os.path.exists(self.WorkingDir+"/class_core_sorted.sel"):
                self.doStep3=True
            else:
                return
        params= WorkingDir+'/class '+\
                str(self.thGoodClass)+' '+\
                str(self.thZscore)+' '+\
                str(self.thPCAZscore)+' '+\
                str(self.NumberOfMpiProcesses)+' '+\
                self.SystemFlavour

        launchJob("xmipp_classify_CL2D_core_analysis",
                              params,
                              self.log,
                              False,
                              1,
                              self.NumberOfMpiProcesses,
                              self.SystemFlavour)
        if os.path.exists(self.WorkingDir+"/class_core_sorted.sel"):
            fh=open(self.WorkingDir + "/status.txt", "a")
            fh.write("Step 3: Core analysis finished at " + time.asctime() + "\n")
            fh.close()
        
    def postprocess(self):
        if not stepPerformed("Step 3",self.WorkingDir + "/status.txt"):
            return
        if not self.DoFilter or not os.path.exists(self.WorkingDir+'/preprocessedImages.stk'):
            return
        import xmipp,glob
        fnClasses=glob.glob(WorkingDir+"/class_level_??.sel")
        os.system('xmipp_metadata_selfile_create -s -q'+\
                  ' -p '+self.WorkingDir+'/class_aligned.stk '+\
                  ' -o '+self.WorkingDir+'/class_aligned.sel')
        for fnClass in fnClasses:
            xmipp.substituteOriginalImages(fnClass,self.WorkingDir+'/class_aligned.sel',
                                           fnClass+"_intermediate.sel",xmipp.MDL_IMAGE,True)
            os.system("mv -f "+fnClass+"_intermediate.sel "+fnClass)
        os.remove(self.WorkingDir+'/class_aligned.sel')
        os.remove(self.WorkingDir+'/preprocessedImages.stk')
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
    # Check that the number of classes is correct
    if NumberOfReferences0<=0:
        errors.append("The number of initial classes must be positive")
    # Check that the number of classes is correct
    if NumberOfReferences0>NumberOfReferences:
        errors.append("The number of initial classes cannot be larger than the number of final classes")
    # Check filter parameters
    if DoFilter and (Highpass<0 or Highpass>0.5 or Lowpass<0 or Lowpass>0.5):
        errors.append("The filter frequencies must be between 0 and 0.5")
    # Check core parameters
    if thGoodClass<0 or thGoodClass>100:
        errors.append("The good class threshold must be between 0 and 100")
    
    return errors
#		
# Main
#     
if __name__ == '__main__':
    CL2D=CL2D_class(
                 InSelFile,
                 WorkingDir,
                 ProjectDir,
                 LogDir,
                 NumberOfReferences,
                 NumberOfReferences0,
                 NumberOfIterations,
                 DoFilter,
                 Highpass,
                 Lowpass,
                 ComparisonMethod,
                 ClusteringMethod,
                 AdditionalParameters,
                 thGoodClass,
                 thZscore,
                 thPCAZscore,
                 NumberOfMpiProcesses,
                 SystemFlavour)
