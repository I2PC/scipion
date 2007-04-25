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

def visualize_images(Names,Areselfile=False,SelFileWidth="",SelFileHeight=""):
    import os
    if (len(Names)>0):
        if (Areselfile):
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
    def plot_xy_file(self,DataFile,
                          Title="",
                          X_Label="x",
                          Y_Label="y",
                          X_col=1,
                          Y_col=2):
       """ plots a file using gnuplot
           data file may have comments (line starting with #)
       """
       print DataFile
       print Title
       print X_Label
       print Y_Label
       self.send(" set title '"+Title+"'")   
       self.send(" set xlabel '"+X_Label+"'")   
       self.send(" set ylabel '"+Y_Label+"'")   
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
           data file may have comments (line starting with #)
       """
       print DataFile
       print Title
       print X_Label
       print Y_Label
       self.send(" set title '"+Title+"'")   
       self.send(" set xlabel '"+X_Label+"'")   
       self.send(" set ylabel '"+Y_Label+"'")   
       self.send(" plot '" + DataFile + "'using "+str(X_col)+":"+str(Y1_col)+" title '"+\
                                        titley1+"' with lines")   
       self.send(" replot '" + DataFile + "'using "+str(X_col)+":"+str(Y2_col)+" title '"+\
                                        titley2+"' with lines")   
