from views_util import *

############## VIEWER SPIDER WARD ############
def viewerWARD(request, protocol, protocolViewer):
    ioDict = {}
    
    if protocolViewer.doShowClasses:
        typeUrl, url = doVisualizeClasses(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowDendrogram:
        typeUrl, url = doVisualizeDendrogram(request, protocol, protocolViewer)
        ioDict[typeUrl]= url

    return ioDict

def doVisualizeClasses(request, protocol, protocolViewer):
    pass

def doVisualizeDendrogram(request, protocol, protocolViewer):
    pass

