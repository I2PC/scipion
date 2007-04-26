#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_projmatch.py
#
# Example use:
# python visualize_projmatch.py
#
# This script requires that protocol_projmatch.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
#show results for iteration
DisplayIterNo=1
#------------------------------------------------------------------------------------------------
# {section} Volume visualization
#------------------------------------------------------------------------------------------------
# Visualize volumes in slices along Z?
VisualizeVolZ=False
# Visualize volumes in slices along X?
VisualizeVolX=False
# Visualize volumes in slices along Y?
VisualizeVolY=False
# Visualize volumes in UCSF Chimera?
""" For this to work, you need to have chimera installed!
"""
VisualizeVolChimera=True
# {expert} Width of Matrix-views (even value!):
MatrixWidth=10
#------------------------------------------------------------------------------------------------
# {section}Reference volume 
#------------------------------------------------------------------------------------------------
#show masked volume
DisplayMask=False
#Show projection maching library and alignes classes
DisplayProjectionMatching=False
#display angular distribution
DisplayAngularDistribution=False
#display reconstructed volume
DisplayReconstruction=False



#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_projmatch_class:

    #init variables
    def __init__(self,
                _VisualizeVolZ,
                _VisualizeVolX,
                _VisualizeVolY,
                _VisualizeVolChimera,
                _DisplayIterNo,
                _DisplayMask,
                _DisplayProjectionMatching,
                _DisplayReconstruction,
                _MatrixWidth
                ):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log
        import visualization

        # import corresponding protocol
        import protocol_projmatch
        self._MatrixWidth=_MatrixWidth
        self._WorkDirectory=protocol_projmatch.WorkDirectory
        self._LogDir=protocol_projmatch.LogDir
        self._ProjectDir=protocol_projmatch.ProjectDir
        self._ReferenceVolume=protocol_projmatch.ReferenceVolume
        self._multi_align2d_sel=protocol_projmatch.multi_align2d_sel
        self._SelFileName=self._ProjectDir+'/'+str(protocol_projmatch.SelFileName)
        self._ReferenceVolume=protocol_projmatch.ReferenceVolume
        self.mylog=log.init_log_system(self._ProjectDir,
                                       self._LogDir,
                                       sys.argv[0],
                                       self._WorkDirectory)
        self._iteration_number=_DisplayIterNo
        Iteration_Working_Directory=self._WorkDirectory+'/Iter_'+\
                                   str(self._iteration_number)
        self.ShowVolumes=[] 
        self.ShowVolumes.append(self._ReferenceVolume)
        os.chdir(Iteration_Working_Directory)
        if (_DisplayMask):
           visualization.visualize_volumes(self.ShowVolumes,
                                                _VisualizeVolZ,
                                                _VisualizeVolX,
                                                _VisualizeVolY,
                                                _VisualizeVolChimera)
        self.ShowSelfiles=[] 
        self.ShowSelfiles.append(self._multi_align2d_sel)
        if (_DisplayProjectionMatching):
            visualization.visualize_images(self.ShowSelfiles,
                             True,
                             self._MatrixWidth)
        #if (_DisplayReconstruction):
        #    self.visualize_Reconstruction(self._SomName,self._SpectraName)
            
        # Return to parent dir
        os.chdir(os.pardir)


    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

#		
# Main
#     
if __name__ == '__main__':

    # create projmatch_class object
    visualize_projmatch=visualize_projmatch_class(VisualizeVolZ,
                                                  VisualizeVolX,
                                                  VisualizeVolY,
                                                  VisualizeVolChimera,
                                                  DisplayIterNo,
                                                  DisplayMask,
                                                  DisplayProjectionMatching,
                                                  DisplayReconstruction,
                                                  MatrixWidth)
    # close 
    visualize_projmatch.close()

