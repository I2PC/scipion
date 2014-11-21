
import sys
from pyworkflow.manager import Manager
from pyworkflow.em.packages.xmipp3.nma.viewer_nma import createDistanceProfilePlot
from pyworkflow.em.packages.xmipp3.protocol_movie_alignment import createPlots, PLOT_POLAR, PLOT_CART, PLOT_POLARCART
from pyworkflow.em.packages.xmipp3.nma.viewer_nma import createVmdView
from  pyworkflow.gui.plotter import Plotter

# Import possible Object commands to be handled
from pyworkflow.em.showj import OBJCMD_NMA_PLOTDIST, OBJCMD_NMA_VMD, OBJCMD_MOVIE_ALIGNPOLAR, OBJCMD_MOVIE_ALIGNCARTESIAN, OBJCMD_MOVIE_ALIGNPOLARCARTESIAN


if __name__ == '__main__':

    cmd = sys.argv[1]
    projectId = sys.argv[2]
    protocolId = int(sys.argv[3])
    objId = int(sys.argv[4])
    project = Manager().loadProject(projectId)
    protocol = project.mapper.selectById(protocolId)
    
    if cmd == OBJCMD_NMA_PLOTDIST:
        Plotter.setInteractive(True) # use 
        plotter = createDistanceProfilePlot(protocol, modeNumber=objId)
        plotter.show(block=True)
        
    elif cmd == OBJCMD_NMA_VMD:
        vmd = createVmdView(protocol, modeNumber=objId)
        vmd.show()

    elif cmd == OBJCMD_MOVIE_ALIGNPOLAR:
        Plotter.setInteractive(True) # use
        plotter = createPlots(PLOT_POLAR, protocol, objId)
        plotter.show(block=True)

    elif cmd == OBJCMD_MOVIE_ALIGNCARTESIAN:
        Plotter.setInteractive(True) # use
        plotter = createPlots(PLOT_CART, protocol, objId)
        plotter.show(block=True)

    elif cmd == OBJCMD_MOVIE_ALIGNPOLARCARTESIAN:
        Plotter.setInteractive(True) # use
        plotter = createPlots(PLOT_POLARCART, protocol, objId)
        plotter.show(block=True)


