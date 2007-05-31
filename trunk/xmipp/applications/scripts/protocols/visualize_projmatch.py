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
#show results for iterations
""" You may specify more than one iteration here 
    This can be done by a sequence of numbers (for instance, "2 8" 
    specifies iteration 2 and 8 (but not 3, 4, 5, 6 and 7)
"""
DisplayIterationsNo="1"
#------------------------------------------------------------------------------------------------
# {section} Volume visualization
#------------------------------------------------------------------------------------------------
# Visualize volumes in slices along Z?
VisualizeVolZ=True
# Visualize volumes in slices along X?
VisualizeVolX=False
# Visualize volumes in slices along Y?
VisualizeVolY=False
# Visualize volumes in UCSF Chimera?
""" For this to work, you need to have chimera installed!
"""
VisualizeVolChimera=True
# {expert} Width of xmipp_show (multiple of 3!):
MatrixWidth=9
#------------------------------------------------------------------------------------------------
# {section}Reference volume 
#------------------------------------------------------------------------------------------------
#show reference volume
""" Volume after filtration and masking
"""
DisplayReference=False
#display angular distribution after projection matching
DisplayAngularDistribution=False
#Show projection maching library and aligned classes
DisplayProjectionMatchingAlign2d=False
#display angular distribution after align2d
DisplayAngularDistributionAlign2d=True
#display reconstructed volume
""" Volume as given by the reconstruction algorithm
"""
DisplayReconstruction=False
#display reconstructed volume
""" Volume after filtration
"""
DisplayFilteredReconstruction=False
#display resolution plots
DisplayResolutionPlots=True


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
                _DisplayIterationsNo,
                _DisplayReference,
                _DisplayProjectionMatchingAlign2d,
                _DisplayReconstruction,
                _DisplayFilteredReconstruction,
                _DisplayResolutionPlots,
                _MatrixWidth,
                _DisplayAngularDistribution,
                _DisplayAngularDistributionAlign2d,
                _WorkingDir,
                ):
	     
        import os,sys,shutil
        scriptdir=os.path.expanduser('~')+'/scripts/'
        sys.path.append(scriptdir) # add default search path
        import log,logging,arg
        import visualization

        # import corresponding protocol
        sys.path.append(_WorkingDir)
        import protocol_projmatch_backup
        #os.chdir(os.pardir)

        self._MatrixWidth=_MatrixWidth
        self._WorkDirectory=protocol_projmatch_backup.WorkDirectory
        self._LogDir=protocol_projmatch_backup.LogDir
        self._ProjectDir=protocol_projmatch_backup.ProjectDir
        self._multi_align2d_sel=protocol_projmatch_backup.multi_align2d_sel
        self._align2d_sel=protocol_projmatch_backup.align2d_sel
        self._align2d_doc=protocol_projmatch_backup.align2d_doc
        self._align2d_doc=protocol_projmatch_backup.align2d_doc
        self._SelFileName=self._ProjectDir+'/'+str(protocol_projmatch_backup.SelFileName)
        self._Filtered_Image=protocol_projmatch_backup.Filtered_Image
        self._Proj_Maching_Output_Root_Name=protocol_projmatch_backup.Proj_Maching_Output_Root_Name
        self._Proj_Maching_Output_Root_Name + '.doc'
        self._mylog=log.init_log_system(self._ProjectDir,
                                       self._LogDir,
                                       sys.argv[0],
                                       self._WorkDirectory)
        self._mylog.setLevel(logging.DEBUG)
        self._DisplayIterationsNo=arg.getListFromVector(_DisplayIterationsNo)
        self._ReconstrucedVolume=protocol_projmatch_backup.ReconstrucedVolume
        _user_suplied_Filtered_Image=os.path.abspath(
                                 protocol_projmatch_backup.ReferenceFileName)
        self._Filtered_Image=[]
        protocol_projmatch_backup.fill_name_vector(
                          "",
                          self._Filtered_Image,
                          protocol_projmatch_backup.NumberofIterations,
                          protocol_projmatch_backup.Filtered_Image)
                          
        self._Reference_volume=[]
        protocol_projmatch_backup.fill_name_vector(
                        "",
                        self._Reference_volume,
                        protocol_projmatch_backup.NumberofIterations,
                        protocol_projmatch_backup.ReferenceVolume)
        
        self._DisplayReference_list=[]
        self._ShowPlotsList=[]
        self._TitleList=[]
        self._ShowSelfilesProjMachingAlign2d=[] 
        self._ShowPlotsafteralig2dList=[]
        self._Titleafteralig2dList=[]
        self._DisplayReconstruction_list=[]
         
        for self._iteration_number in _DisplayIterationsNo:
           if self._iteration_number==' ':
              continue

           if (_DisplayReference):
              self._DisplayReference_list.append(
                           self._Reference_volume[int(self._iteration_number)])
           if (_DisplayAngularDistribution):
              self._ShowPlotsList.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/'+\
                                     protocol_projmatch_backup.Proj_Maching_Output_Root_Name+\
                                     '_classes.doc')
              self._TitleList.append('Angular distribution after *projection matching* for iteration '+\
                       str(self._iteration_number))

           if (_DisplayProjectionMatchingAlign2d):
              self._ShowSelfilesProjMachingAlign2d.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/'+\
                                       self._multi_align2d_sel)

           if (_DisplayAngularDistributionAlign2d):
              print "here"
              self._ShowPlotsafteralig2dList.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+
                                      '/'+\
                                       self._align2d_doc)
              self._Titleafteralig2dList.append('Angular distribution after *alig2d* for iteration '+\
                       str(self._iteration_number))

           if (_DisplayReconstruction):
              self._DisplayReconstruction_list.append('..'+'/Iter_'+\
                                      str(self._iteration_number)+\
                                      '/Iter_'+str(self._iteration_number)+\
                                      '_'+self._ReconstrucedVolume+'.vol')


        #NAMES ARE RELATIVE TO iTER DIRECTORY  
        self._Iteration_Working_Directory=_WorkingDir+'Iter_1'
        self._mylog.debug ("cd " + self._Iteration_Working_Directory)
        os.chdir(self._Iteration_Working_Directory)

        if (_DisplayReference):
           self._mylog.debug ( "_DisplayReference_list "+str(self._DisplayReference_list))
           visualization.visualize_volumes(self._DisplayReference_list,
                                           _VisualizeVolZ,
                                           _VisualizeVolX,
                                           _VisualizeVolY,
                                           _VisualizeVolChimera)


        if (_DisplayAngularDistribution):
           self._mylog.debug ( "_ShowPlotsList "+str(self._ShowPlotsList))
           self._mylog.debug ( "_TitleList "+str(self._TitleList))
           show_ang_distribution(self._ShowPlotsList,
                                 self._iteration_number,
                                 self._TitleList,
                                 self._mylog)

        if (_DisplayProjectionMatchingAlign2d):
           self._mylog.debug ( "_ShowSelfilesProjMachingAlign2d "+\
                              str(self._ShowSelfilesProjMachingAlign2d))
           visualization.visualize_images(self._ShowSelfilesProjMachingAlign2d,
                                          True,
                                          self._MatrixWidth,
                                          "",
                                          True)

        if (_DisplayAngularDistributionAlign2d):
           self._mylog.debug ( "_ShowPlotsafteralig2dList "+str(self._ShowPlotsafteralig2dList))
           self._mylog.debug ( "_Titleafteralig2dList "+str(self._Titleafteralig2dList))
           show_ang_distribution(self._ShowPlotsafteralig2dList,
                                 self._iteration_number,
                                 self._Titleafteralig2dList,
                                 self._mylog)

        if (_DisplayReconstruction):
           self._mylog.debug ( "_DisplayReconstruction_list "+str(self._DisplayReconstruction_list))
           visualization.visualize_volumes(self._DisplayReconstruction_list,
                                           _VisualizeVolZ,
                                           _VisualizeVolX,
                                           _VisualizeVolY,
                                           _VisualizeVolChimera)

##           if (_DisplayResolutionPlots):
##              plot=visualization.gnuplot()
##              plot_name='..'+'/Iter_'+\
##                         str(self._iteration_number)+\
##                         '/'+'split_sel_2.vol.frc'
##              plot.plot_xy_file(plot_name,
##                          Title="Resolution",
##                          X_Label="Armstrong^-1",
##                          Y_Label="y",
##                          X_col=1,
##                          Y_col=2)
##
##
    def close(self):
        message='Done!'
        print '*',message
        print '*********************************************************************'

def show_ang_distribution(_ShowPlots,_iteration_number,_title,_mylog):
        import os
        import docfiles
        import visualization
        for i in range(len(_ShowPlots)):
#        for plots in _ShowPlots: 
            doc=docfiles.docfile(_ShowPlots[i])
            doc.check_angle_range()
            mini=doc.minimum_of_column(7)
            maxi=doc.maximum_of_column(7)
            _mylog.debug("mini "+ str(mini) +" maxi "+ str(maxi))
            if mini==0:
               mini=1
            if maxi<mini:
               maxi=mini   
            doc.write_several(_ShowPlots[i],
                              10,
                              7,
                              mini,
                              maxi,
                              )
            plot=visualization.gnuplot()
            _title[i] =_title[i]+'\\n min= '+str(mini)+', max= '+str(maxi) 
            plot.plot_xy1y2_several_angular_doc_files(_ShowPlots[i],
                                                      _title[i],
                                                      'degrees',
                                                      'degrees',
                                                      3,
                                                      4)


#		
# Main
#     
if __name__ == '__main__':

    import sys
    WorkingDir=sys.argv[1]

    # create projmatch_class object
    visualize_projmatch=visualize_projmatch_class(VisualizeVolZ,
                                                  VisualizeVolX,
                                                  VisualizeVolY,
                                                  VisualizeVolChimera,
                                                  DisplayIterationsNo,
                                                  DisplayReference,
                                                  DisplayProjectionMatchingAlign2d,
                                                  DisplayReconstruction,
                                                  DisplayFilteredReconstruction,
                                                  DisplayResolutionPlots,
                                                  MatrixWidth,
                                                  DisplayAngularDistribution,
                                                  DisplayAngularDistributionAlign2d,
                                                  WorkingDir)
    # close 
    visualize_projmatch.close()

