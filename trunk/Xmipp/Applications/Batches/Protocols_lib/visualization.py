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

def visualize_images(Names,AreSelFiles=False,SelFileWidth="",SelFileHeight=""):
    import os
    if (len(Names)>0):
        if (AreSelFiles):
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


