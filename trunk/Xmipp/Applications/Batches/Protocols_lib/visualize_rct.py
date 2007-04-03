#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# Protocol for visualization of the results of protocol_rct.py
#
# Example use:
# ./visualize_rct.py
#
# Author: Sjors Scheres, March 2007
#
#------------------------------------------------------------------------------------------------
# {section} Global parameters
#------------------------------------------------------------------------------------------------
# Working subdirectory:
WorkingDir="test1"
# For which classes do you want to perform the visualization?
SelectClasses="1,2"
#------------------------------------------------------------------------------------------------
# {section} Step-by-step visualization
#------------------------------------------------------------------------------------------------
# Visualize untilted average images?
VisualizeUntiltedAverages=True
# Visualize aligned untilted images?
VisualizeUntiltedImages=True
# Visualize aligned tilted images?
VisualizeTiltedImages=True
# Visualize ART reconstructions in slices along Z?
VisualizeArtVolZ=True
# Visualize ART reconstructions in slices along X?
VisualizeArtVolX=True
# Visualize WBP reconstructions in slices along Z?
VisualizeWbpVolZ=True
# Visualize WBP reconstructions in slices along X?
VisualizeWbpVolX=True
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# {end-of-header} USUALLY YOU DO NOT NEED TO MODIFY ANYTHING BELOW THIS LINE ...
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#
class visualize_RCT_class:

    #init variables
    def __init__(self,
                 WorkingDir,
                 SelectClasses,
                 VisualizeUntiltedAverages,
                 VisualizeUntiltedImages,
                 VisualizeTiltedImages,
                 VisualizeArtVolZ,
                 VisualizeArtVolX,
                 VisualizeWbpVolZ,
                 VisualizeWbpVolX):
	     
        import os,sys
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        
        print '*  cd '+str(WorkingDir)
        os.chdir(WorkingDir)

        refs=SelectClasses.split(',')
        for ref in refs:

            basename='*_ref'+str(ref).zfill(5)
            if VisualizeUntiltedAverages:
                imgname=self.getname(basename+'.med.xmp')
                command='xmipp_show -img '+imgname+' &'
                print '* ',command
                os.system(command)

            if VisualizeUntiltedImages:
                selname=self.getname(basename+'.sel')
                command='xmipp_show -sel '+selname+' &'
                print '* ',command
                os.system(command)

            if VisualizeTiltedImages:
                selname=self.getname(basename+'_tilted.sel')
                command='xmipp_show -sel '+selname+' &'
                print '* ',command
                os.system(command)

            if VisualizeArtVolZ:
                volname=self.getname('art_'+basename+'_tilted.vol')
                command='xmipp_show -vol '+volname+' &'
                print '* ',command
                os.system(command)

            if VisualizeArtVolX:
                volname=self.getname('art_'+basename+'_tilted.vol')
                command='xmipp_show -vol '+volname+'x &'
                print '* ',command
                os.system(command)

            if VisualizeWbpVolZ:
                volname=self.getname('wbp_'+basename+'_tilted.vol')
                command='xmipp_show -vol '+volname+' &'
                print '* ',command
                os.system(command)

            if VisualizeWbpVolX:
                volname=self.getname('wbp_'+basename+'_tilted.vol')
                command='xmipp_show -vol '+volname+'x &'
                print '* ',command
                os.system(command)

                
        # Return to parent dir
        print '*  cd ..'
        os.chdir(os.pardir)

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

    visualize_RCT=visualize_RCT_class(WorkingDir,
                                      SelectClasses,
                                      VisualizeUntiltedAverages,
                                      VisualizeUntiltedImages,
                                      VisualizeTiltedImages,
                                      VisualizeArtVolZ,
                                      VisualizeArtVolX,
                                      VisualizeWbpVolZ,
                                      VisualizeWbpVolX)
    visualize_RCT.close()
