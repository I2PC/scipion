import sys
from pyworkflow.em.showj import *
from pyworkflow.manager import Manager
from pyworkflow.em.packages.xmipp3.nma import XmippProtNMA
from pyworkflow.em.packages.xmipp3.nma.viewer_nma import createDistanceProfilePlot


def showPlot():
    print "aaaaaaaaaaaaaaaaa"
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    def update_line(num, data, line):
        line.set_data(data[...,:num])
        return line,

    fig1 = plt.figure()

    data = np.random.rand(2, 25)
    l, = plt.plot([], [], 'r-')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.title('test')
    line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
        interval=50, blit=True)
    #line_ani.save('lines.mp4')

    fig2 = plt.figure()

    x = np.arange(-9, 10)
    y = np.arange(-9, 10).reshape(-1, 1)
    base = np.hypot(x, y)
    ims = []
    for add in np.arange(15):
        ims.append((plt.pcolor(x, y, base + add, norm=plt.Normalize(0, 30)),))

    im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000,
        blit=True)
    #im_ani.save('im.mp4', metadata={'artist':'Guido'})

    plt.show()

if __name__ == '__main__':
    #TODO: REMOVE THIS AFTER DEBUGGING
    print "ARGS: "
    for i, arg in enumerate(sys.argv):
        print "%02d: %s" % (i, arg)


    cmd = sys.argv[1]
    if cmd == OBJ_NMAPLOT:
        print "OBJ_NMAPLOT"

        projectId=sys.argv[2]
        protocolId=int(sys.argv[3])
        modeNumber=int(sys.argv[4])
        project = Manager().loadProject(projectId)
        protocol = project.mapper.selectById(protocolId)
        plotter = createDistanceProfilePlot(protocol, modeNumber)
        plotter.show()
        #showPlot()


