#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_CL2D.py
#
# Example use:
# python visualize_cl2d.py protocol_cl2d.py
#
# This script requires that protocol_cl2d.py is in the current directory
#
# Author: Carlos Oscar Sanchez Sorzano
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Show results for levels
""" You may specify more than one iteration here 
    This can be done by a sequence of numbers (for instance, "0 2" 
    specifies iteration 0 and 2. The first level is always level 0.
    You can use negative numbers: -1 is the last iteration, -2 is the
    one before last, ...
"""
DisplayLevelsNo='1'

# Show Raw CL2D results
""" Show plain CL2D results as they come out without postprocessing """
DoShowRawCL2D=False

# Show Filtered and sorted CL2D results
""" Show plain CL2D results as they come out without postprocessing """
DoShowFilteredSortedCL2D=True

# Show Analysis Summary
""" Show plain CL2D results as they come out without postprocessing """
DoShowAnalysisSummary=True

#------------------------------------------------------------------------------------------------
# {section} Filter 1: General Z-score
#------------------------------------------------------------------------------------------------
# {expert} Show general image sorting by quality
""" This is the basis for the first filter. F1 selfiles are constructed
    by intersecting the classes with the images whose average image features
    Z-score when considering the whole dataset. If you think there are two
    few images passing this filter, increase the junk Zscore.
"""
DoShowF1GeneralSort=False

# {expert} Show classes after this filter
DoShowF1=False

#------------------------------------------------------------------------------------------------
# {section} Filter 2: Class PCA Z-score
#------------------------------------------------------------------------------------------------
# {expert} Show general image sorting by quality
""" After filtering with F1, images remaining in the classes are
    Z-scored within each class using a PCA projection.
    Those images whose Z-score is larger than Junk Zscore are discarded.
    If you feel too many images have been discarded, increase this threshold.
    For a Gaussian distribution 99.5% of the values are within a Z-score
    of 3.
"""
DoShowF2GeneralSort=False

# {expert} Show classes after this filter
DoShowF2=False

#------------------------------------------------------------------------------------------------
# {section} Filter 3: Class cores
#------------------------------------------------------------------------------------------------
# {expert} Show classes after this filter
""" Filter 3 removes from each class all those images that have not been
    in the same class during the whole hierarchical classification. This is
    called the core. Moreover, if the core of a class is not larger than
    a certain threshold, the whole class is removed. If you think too
    many classes have been discarded, decrease the Good class core size.
"""
DoShowF3=False

#------------------------------------------------------------------------------------------------
# {section} Filter 4: Class Z-score
#------------------------------------------------------------------------------------------------
# {expert} Show classes after this filter
""" After filtering with F3, images remaining in the classes are
    Z-scored again within each class using general image features.
    Those images whose Z-score is larger than Junk Zscore are discarded.
    If you feel too many images have been discarded, increase this threshold.
    For a Gaussian distribution 99.5% of the values are within a Z-score
    of 3.
"""
DoShowF4=False

#------------------------------------------------------------------------------------------------
# {section} Filter 5: Class PCA Z-score
#------------------------------------------------------------------------------------------------
# {expert} Show classes after this filter
""" After filtering with F4, images remaining in the classes are Z-scored
    again within each class using a PCA projection. Those images whose
    Z-score is larger than Junk Zscore are discarded. If you feel too
    many images have been discarded, increase this threshold. For a
    Gaussian distribution 99.5% of the values are within a Z-score of 3.
"""
DoShowF5=False

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_CL2D_class:

    #init variables
    def __init__(self,
                 DisplayLevelsNo,
                 DoShowRawCL2D,
                 DoShowFilteredSortedCL2D,
                 DoShowAnalysisSummary,
                 DoShowF1GeneralSort,
                 DoShowF1,
                 DoShowF2GeneralSort,
                 DoShowF2,
                 DoShowF3,
                 DoShowF4,
                 DoShowF5,
                 ProtocolName):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import log

        # Import the corresponding protocol and get WorkingDir
        shutil.copy(ProtocolName,'protocol.py')
        import protocol
        self.WorkingDir=protocol.WorkingDir
        self.DisplayLevelsNo=DisplayLevelsNo

        # Show
        if DoShowRawCL2D:
            self.showRawCL2D()
        if DoShowFilteredSortedCL2D:
            self.showFilteredSorted()
        if DoShowAnalysisSummary:
            self.showAnalysisSummary()
        if DoShowF1GeneralSort:
            self.showF1GeneralSort()
        if DoShowF1:
            self.showSuffix('class_level_??_.sel',"_F1")
        if DoShowF2GeneralSort:
            self.showF2GeneralSort()
        if DoShowF2:
            self.showSuffix('class_level_??_.sel',"_F12")
        if DoShowF3:
            self.showSuffix('class_level_??_.sel',"_F123")
        if DoShowF4:
            self.showSuffix('class_level_??_.sel',"_F1234")
        if DoShowF5:
            self.showSuffix('class_level_??_F12345.sel',"")

        # Return to parent dir and remove protocol.py(c)
        if (os.path.exists('protocol.py')):
            os.remove('protocol.py')
        if (os.path.exists('protocol.pyc')):
            os.remove('protocol.pyc')
        
        # Finish
        self.close()

    def showRawCL2D(self):
        import os,glob
        
        selfiles=glob.glob(self.WorkingDir+'/class_level_??_.sel')
        selfiles.sort()
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            levels=self.DisplayLevelsNo.split(' ')
            for level in levels:
                selfile=selfiles[int(level)]
                command='xmipp_show -cl2d '+selfile.replace('.sel','')+ ' &'
                print '* ',command
                os.system(command)

    def showFilteredSorted(self):
        import os,glob
        
        selfiles=glob.glob(self.WorkingDir+'/class_level_??_F12345_sorted.sel')
        selfiles.sort()
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            levels=self.DisplayLevelsNo.split(' ')
            for level in levels:
                selfile=selfiles[int(level)]
                command='xmipp_show -cl2d '+selfile.replace('.sel','')+ ' &'
                print '* ',command
                os.system(command)

    def showAnalysisSummary(self):
        import os
        fnSummary=self.WorkingDir+'/class_analysis_summary.txt'
        if os.path.exists(fnSummary):
            command='xmipp_edit -i '+fnSummary
            print '* ',command
            os.system(command)

    def showF1GeneralSort(self):
        import os
        command='xmipp_show -sel '+self.WorkingDir+'/class_F1_sorted.sel '+\
            self.WorkingDir+'/class_F1_good.sel &'
        print '* ',command
        os.system(command)
        command='xmipp_edit -i '+self.WorkingDir+'/class_F1.sumZ &'
        print '* ',command
        os.system(command)

    def showF2GeneralSort(self):
        import os
        command='xmipp_show -sel '+self.WorkingDir+'/class_F12_sorted.sel &'
        print '* ',command
        os.system(command)
        command='xmipp_edit -i '+self.WorkingDir+'/class_F12_sorted_score.txt &'
        print '* ',command
        os.system(command)

    def showSuffix(self,pattern,suffix):
        import os,glob
        
        selfiles=glob.glob(self.WorkingDir+'/'+pattern)
        selfiles.sort()
        if len(selfiles)==0:
            print "No selfiles yet. Visualize after job completion..."
        else:
            levels=self.DisplayLevelsNo.split(' ')
            for level in levels:
                selfile=selfiles[int(level)]
                command='xmipp_show -cl2d '+selfile.replace('.sel','')+\
                    ' -filterSuffix '+suffix+' &'
                print '* ',command
                os.system(command)

    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    import sys
    ProtocolName=sys.argv[1]

    # create ML2D_class object
    visualize_CL2D=visualize_CL2D_class(DisplayLevelsNo,
                                        DoShowRawCL2D,
                                        DoShowFilteredSortedCL2D,
                                        DoShowAnalysisSummary,
                                        DoShowF1GeneralSort,
                                        DoShowF1,
                                        DoShowF2GeneralSort,
                                        DoShowF2,
                                        DoShowF3,
                                        DoShowF4,
                                        DoShowF5,
                                        ProtocolName)
