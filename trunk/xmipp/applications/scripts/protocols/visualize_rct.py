#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_rct.py
#
# Example use:
# ./visualize_rct.py xmipp_protocol_rct.py
#
# This script requires that xmipp_protocol_rct.py is in the current directory
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# For which classes do you want to perform the visualization?
SelectClasses="1,2"
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
# {expert} Width of selfile visualizations (preferably an even value):
MatrixWidth=10
#------------------------------------------------------------------------------------------------
# {section} Step-by-step visualization
#------------------------------------------------------------------------------------------------
# Visualize untilted average images?
VisualizeUntiltedAverages=True
# Visualize aligned untilted images?
VisualizeUntiltedImages=True
# Visualize aligned tilted images?
VisualizeTiltedImages=True
# Visualize ART reconstructions?
VisualizeArtVols=True
# Visualize WBP reconstructions?
VisualizeWbpVols=True
# Visualize filtered ART reconstructions?
VisualizeFilteredArtVols=True
# Visualize filtered WBP reconstructions?
VisualizeFilteredWbpVols=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
#
class visualize_RCT_class:

    #init variables
    def __init__(self,
                 SelectClasses,
                 VisualizeVolZ,
                 VisualizeVolX,
                 VisualizeVolY,
                 VisualizeVolChimera,
                 MatrixWidth,
                 VisualizeUntiltedAverages,
                 VisualizeUntiltedImages,
                 VisualizeTiltedImages,
                 VisualizeArtVols,
                 VisualizeWbpVols,
                 VisualizeFilteredArtVols,
                 VisualizeFilteredWbpVols,
                 ProtocolName):
	     
        import os,sys,shutil
        scriptdir=os.path.split(os.path.dirname(os.popen('which xmipp_protocols','r').read()))[0]+'/protocols'
        sys.path.append(scriptdir) # add default search path
        import visualization, utils_xmipp
     
        # Import the corresponding protocol, get WorkingDir and go there
        pardir=os.path.abspath(os.getcwd())
        shutil.copy(ProtocolName,'protocol.py')
        import protocol
        self.WorkingDir=protocol.WorkingDir
        os.chdir(self.WorkingDir)
        
        ShowVolumes=[]
        ShowSelfiles=[]
        ShowImages=[]

        refs=SelectClasses.split(',')
        for ref in refs:

            basename=composeFileName('rct_ref',ref,'')
            if VisualizeUntiltedAverages:
                ShowImages.append(basename+'_untilted.med.xmp')

            if VisualizeUntiltedImages:
                ShowSelfiles.append(basename+'_untilted.sel')

            if VisualizeTiltedImages:
                ShowSelfiles.append(basename+'_tilted.sel')

            if VisualizeArtVols:
                ShowVolumes.append('art_'+basename+'_tilted.vol')

            if VisualizeWbpVols:
                ShowVolumes.append('wbp_'+basename+'_tilted.vol')

            if VisualizeFilteredArtVols:
                ShowVolumes.append('art_'+basename+'_tilted_filtered.vol')

            if VisualizeFilteredWbpVols:
                ShowVolumes.append('wbp_'+basename+'_tilted_filtered.vol')

        visualization.visualize_volumes(ShowVolumes,
                                        VisualizeVolZ,
                                        VisualizeVolX,
                                        VisualizeVolY,
                                        VisualizeVolChimera)

        visualization.visualize_images(ShowSelfiles,True,MatrixWidth)

        visualization.visualize_images(ShowImages)
 
        # Return to parent dir and remove protocol.py(c)
        os.chdir(pardir)
        if (os.path.exists('protocol.py')):
            os.remove('protocol.py')
        if (os.path.exists('protocol.pyc')):
            os.remove('protocol.pyc')

    def getname(self,pattern):
        import glob,sys
        names=glob.glob(pattern)
        if len(names)<1:
            print 'No file '+pattern+' exists, ignore...'
            names=[""]
        elif len(names)>1:
            print 'Multiple files '+pattern+\
                  ' exist, taking only first one into account'
        return str(names[0])

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

    visualize_RCT=visualize_RCT_class(SelectClasses,
                                      VisualizeVolZ,
                                      VisualizeVolX,
                                      VisualizeVolY,
                                      VisualizeVolChimera,
                                      MatrixWidth,
                                      VisualizeUntiltedAverages,
                                      VisualizeUntiltedImages,
                                      VisualizeTiltedImages,
                                      VisualizeArtVols,
                                      VisualizeWbpVols,
                                      VisualizeFilteredArtVols,
                                      VisualizeFilteredWbpVols,
                                      ProtocolName)
    visualize_RCT.close()
