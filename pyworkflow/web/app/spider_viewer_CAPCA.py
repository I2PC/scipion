from views_util import *

############## VIEWER SPIDER CAPCA ############
def viewerCAPCA(request, protocol, protocolViewer):
    ioDict = {}
    # SHOWJS
    if protocolViewer.doShowEigenImages and protocolViewer.doShowReconsImages:
        typeUrl, url = doShowImagesCAPCA(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
    else:   
        if protocolViewer.doShowEigenImages:
            typeUrl, url = doShowEigenImages(request, protocol, protocolViewer)
            ioDict[typeUrl]= url
        elif protocolViewer.doShowReconsImages:
            typeUrl, url = doShowReconsImages(request, protocol, protocolViewer)
            ioDict[typeUrl]= url
    # PLOTS
    if protocolViewer.doShowHistogram and protocolViewer.doShowFactorMaps:
        typeUrl, url = doPlotsCAPCA(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
    else:
        if protocolViewer.doShowHistogram:
            typeUrl, url = doPlotsCAPCA(request, protocol, protocolViewer)
            ioDict[typeUrl]= url
        elif protocolViewer.doShowFactorMaps:
            typeUrl, url = doPlotFactorMaps(request, protocol, protocolViewer)
            ioDict[typeUrl]= url
    # FILE VIEWER
    if protocolViewer.doShowPcaFile:
        typeUrl, url = doShowPcaFile(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
        
    return ioDict

def doShowImagesCAPCA(request, protocol, protocolViewer):
    _, eigenUrl = doShowEigenImages(request, protocol, protocolViewer)
    _, reconsUrl = doShowReconsImages(request, protocol, protocolViewer)
    return "showjs", [str(eigenUrl) , str(reconsUrl)]

def doShowEigenImages(request, protocol, protocolViewer):
    return "showj", "/visualize_object/?path="+ protocol._getFileName('eigenimages')

def doShowReconsImages(request, protocol, protocolViewer):
    return "showj", "/visualize_object/?path="+ protocol._getFileName('reconstituted')

def doPlotsCAPCA(request, protocol, protocolViewer):
    _, histogram = doPlotHistogram(request, protocol, protocolViewer)
    _, factorMaps = doPlotFactorMaps(request, protocol, protocolViewer)
    return "plots", [str(histogram) , str(factorMaps)]

def doPlotHistogram(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=plotHistogram&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def plotHistogram(request, protocol, protocolViewer):
    fn = protocol._getFileName('eigFile')
    xplotter = protocolViewer.prepPlotHistogram(fn)
    return xplotter

def doPlotFactorMaps(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=plotFactorMaps&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def plotFactorMaps(request, protocol, protocolViewer):
    fn = protocol._getFileName('imcFile')
    xplotter = protocolViewer.prepPlotFactorMaps(fn)
    return xplotter

def doShowPcaFile(request, protocol, protocolViewer):
    html = textfileViewer("PCA files", [protocol.imcFile.filename.get()])
    return "html", html
