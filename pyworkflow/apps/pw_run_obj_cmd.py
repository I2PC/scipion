import sys
from pyworkflow.em.showj import *
from pyworkflow.manager import Manager
from pyworkflow.em.packages.xmipp3.nma.viewer_nma import createDistanceProfilePlot
from  pyworkflow.gui.plotter import Plotter




if __name__ == '__main__':



    cmd = sys.argv[1]
    if cmd == OBJ_NMAPLOT:

        Plotter.setInteractive(True)
        projectId=sys.argv[2]
        protocolId=int(sys.argv[3])
        modeNumber=int(sys.argv[4])
        project = Manager().loadProject(projectId)
        protocol = project.mapper.selectById(protocolId)
        plotter = createDistanceProfilePlot(protocol, modeNumber)

        plotter.show()
        #showPlot()


