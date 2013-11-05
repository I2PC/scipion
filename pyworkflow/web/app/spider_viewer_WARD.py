from views_util import *

############## VIEWER SPIDER WARD ############
def viewerWARD(request, protocolViewer):
    ioDict = {}
    
    if protocolViewer.doShowClasses:
        typeUrl, url = doVisualizeClasses(request, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowDendrogram:
        typeUrl, url = doVisualizeDendrogram(request, protocolViewer)
        ioDict[typeUrl]= url

    return ioDict

def doVisualizeClasses(request, protocolViewer):
    pass

def doVisualizeDendrogram(request, protocolViewer):
    pass

