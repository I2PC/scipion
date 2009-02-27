# Visualize a list of volumes using Xmipp or Chimera
def visualize_volumes(Names,ShowSliceZ,ShowSliceX,ShowSliceY,ShowChimera):
    import os
    if (len(Names)>0):
        chimera_command='chimera '
        xmipp_command='xmipp_show -vol '
        for name in Names:
            if not os.path.exists(name):
                print '* WARNING: '+name+' does not exist, skipping...'
            else:
                if (ShowSliceZ):
                    xmipp_command+=name+' '
                if (ShowSliceX):
                    xmipp_command+=name+'x '
                if (ShowSliceY):
                    xmipp_command+=name+'y '
                if (ShowChimera):
                    chimera_command+='spider:'+name+' '
        
        if (ShowSliceZ or ShowSliceX or ShowSliceY): 
            xmipp_command+=' &'
            print '* ',xmipp_command
            os.system(xmipp_command)

        if (ShowChimera):
            chimera_command+=' &'
            print '* ',chimera_command
            os.system(chimera_command)

def angDistributionChimera(_ShowPlots,
                           _displayVolList,
                           _mylog):
    import string,spider_header,os
    print "***********************************************************\n* "
    for i in range(len(_ShowPlots)):
        doc=_ShowPlots[i]
        suffix = ".doc"
        docfile_without_ext=doc[: -len(suffix)]
        bildfile = docfile_without_ext + ".bild"
        vol=_displayVolList[i]
        myheader=spider_header.spiderheader(vol)
        nslices=myheader.nz
        radius = str(int(nslices/2.))
        command  = "xmipp_angular_distribution_show "
        command += "-ang " + doc + " "
        command += "-R " + radius + " "
        command += "-bild " + bildfile + " "
        command += "-shift_center " + radius + " "
        command += "-wcol 5 -r 5"
        print "*" , command, "\n"
        _mylog.debug(command)
        os.system(command)
        chimera_command='chimera '  + 'spider:' + vol + ' ' + bildfile 
        chimera_command+=' &'
        _mylog.debug(chimera_command)
        print "*" , chimera_command, "\n"
        os.system(chimera_command)
        
        #call chimera
    print "***********************************************************\n* "


def visualize_images(Names,
                     Areselfile=False,
                     SelFileWidth="",
                     SelFileHeight="",
                     ShowAll=False):
    import os
    if (len(Names)>0):
        if (Areselfile):
            if(ShowAll):
               command='xmipp_show -showall '
            else:
               command='xmipp_show '
            if not (SelFileWidth==""):
                command+=' -w '+str(SelFileWidth)
            if not (SelFileHeight==""):
                command+=' -h '+str(SelFileHeight)
            command+=' -sel '
        else:
            command='xmipp_show -img '
        
        for name in Names:
            if not os.path.exists(name):
                print '* WARNING: '+name+' does not exist, skipping...'
            else:
                command+=name+' '
        command+=' &'
        print '* ',command
        os.system(command)
        
#use gnuplot from python
import os
class gnuplot:
    def __init__(self,persistent=True):
        gnuplot_base="gnuplot"
        if persistent:
           gnuplot_base="gnuplot -persist"
        self.session = os.popen(gnuplot_base,"w")
    def __del__(self):
        print "closing gnuplot session..."
        self.session.close()
    def send(self, cmd):
        self.session.write(cmd+'\n')
        self.session.flush() 

    def prepare_empty_plot(self,
                           Title="",
                           X_Label="x",
                           Y_Label="y"):
       self.send(" set title '"+Title+"'")   
       self.send(" set xlabel '"+X_Label+"'")   
       self.send(" set ylabel '"+Y_Label+"'")   

    def plot_xy_file(self,DataFile,
                          Title="",
                          X_Label="x",
                          Y_Label="y",
                          X_col=1,
                          Y_col=2):
       """ plots a file using gnuplot
       """
       self.prepare_empty_plot(Title,X_Label,Y_Label)
       self.send(" plot '" + DataFile + "' using "+str(X_col)+":"+str(Y_col)+" with lines")   

    def plot_xy1y2_file(self,DataFile,
                             Title="",
                             titley1="",
                             titley2="",
                             X_Label="x",
                             Y_Label="y",
                             X_col=1,
                             Y1_col=2,
                             Y2_col=3):
       """ plots a file using gnuplot
       """
       self.prepare_empty_plot(Title,X_Label,Y_Label)
       self.send(" plot '" + DataFile + "'using "+str(X_col)+":"+str(Y1_col)+" title '"+\
                                        titley1+"' with lines")   
       self.send(" replot '" + DataFile + "'using "+str(X_col)+":"+str(Y2_col)+" title '"+\
                                        titley2+"' with lines")   

    def plot_xy1y2_several_angular_doc_files(self, docfilename,
                                Title="",
                                X_Label="x",
                                Y_Label="y",
                                X_col=3,
                                Y_col=4):
       import glob
       file_patern=docfilename+'[0-9]*.doc'
       self.send(' set title "'+Title+'"')   
       self.send(" set xlabel '"+X_Label+"'")   
       self.send(" set ylabel '"+Y_Label+"'") 
       self.send(" set polar") 
       self.send(" set size square")
       self.send(" set angles degrees") 
       self.send(" set xrange[-95:95]") 
       self.send(" set yrange[-95:95]") 

# The following line is not functional in gnuplot 3.7 (jumilla)
# pre-treat docfiles to remove comment lines!
#       self.send(' set datafile commentschars "#;"') 
       self.send(" set grid polar 30") 
       self.send(" set xtics 10") 
       i = 0
       for _docfilename in glob.glob(file_patern):
           point_size = _docfilename.replace(docfilename,'')
           point_size = point_size.replace('.doc','')
           point_size=(float(point_size)+2)/2
           if (i==0):
               plot='plot'
           else:
               plot='replot'
           i=i+1
           plotcommand=plot+" '" + _docfilename+ "' using "+str(X_col)+":"+str(Y_col)+ \
                     " title '' with points pt 7 ps "+ str(point_size)
           self.send(plotcommand)
#           print plotcommand

