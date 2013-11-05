from views_util import *

############## VIEWER SPIDER CAPCA ############
def viewerCAPCA(request, protocolViewer):
    ioDict = {}
    # SHOWJS
    if protocolViewer.doShowEigenImages and protocolViewer.doShowReconsImages:
        typeUrl, url = doShowImagesCAPCA(request, protocolViewer)
        ioDict[typeUrl]= url
    else:   
        if protocolViewer.doShowEigenImages:
            typeUrl, url = doShowEigenImages(request, protocolViewer)
            ioDict[typeUrl]= url
        elif protocolViewer.doShowReconsImages:
            typeUrl, url = doShowReconsImages(request, protocolViewer)
            ioDict[typeUrl]= url
    # PLOTS
    if protocolViewer.doShowHistogram and protocolViewer.doShowFactorMaps:
        typeUrl, url = doPlotsCAPCA(request, protocolViewer)
        ioDict[typeUrl]= url
    else:
        if protocolViewer.doShowHistogram:
            typeUrl, url = doPlotsCAPCA(request, protocolViewer)
            ioDict[typeUrl]= url
        elif protocolViewer.doShowFactorMaps:
            typeUrl, url = doPlotFactorMaps(request, protocolViewer)
            ioDict[typeUrl]= url
    # FILE VIEWER
    if protocolViewer.doShowPcaFile:
        typeUrl, url = doShowPcaFile(request, protocolViewer)
        ioDict[typeUrl]= url
        
    return ioDict

def doShowImagesCAPCA(request, protocolViewer):
    _, eigenUrl = doShowEigenImages(request, protocolViewer)
    _, reconsUrl = doShowReconsImages(request, protocolViewer)
    return "showjs", [str(eigenUrl) , str(reconsUrl)]

def doShowEigenImages(request, protocolViewer):
    return "showj", "/visualize_object/?path="+ protocolViewer.protocol._getFileName('eigenimages')

def doShowReconsImages(request, protocolViewer):
    return "showj", "/visualize_object/?path="+ protocolViewer.protocol._getFileName('reconstituted')

def doPlotsCAPCA(request, protocolViewer):
    _, histogram = doPlotHistogram(request, protocolViewer)
    _, factorMaps = doPlotFactorMaps(request, protocolViewer)
    return "plots", [str(histogram) , str(factorMaps)]

def doPlotHistogram(request, protocolViewer):
    return "plot","/view_plots/?function=plotHistogram&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocolViewer.protocol.getObjId())

def plotHistogram(request, protocolViewer):
    xplotter = protocolViewer.prepPlotHistogram()
    return xplotter

def doPlotFactorMaps(request, protocolViewer):
    return "plot","/view_plots/?function=plotFactorMaps&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocolViewer.protocol.getObjId())

def plotFactorMaps(request, protocolViewer):
    print "first factor:", protocolViewer.firstFactor.get()
    print "second factor:", protocolViewer.secondFactor.get()
    xplotter = protocolViewer.prepPlotFactorMaps()
    return xplotter

def doShowPcaFile(request, protocolViewer):
    html = textfileViewer("PCA files", [protocolViewer.protocol.imcFile.filename.get()])
    return "html", html
