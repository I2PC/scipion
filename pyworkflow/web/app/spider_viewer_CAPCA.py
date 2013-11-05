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
    path = protocolViewer.protocol._getFileName('eigenimages')
    return "showj", "/visualize_object/?path="+ path

def doShowReconsImages(request, protocolViewer):
    path = protocolViewer.protocol._getFileName('reconstituted')
    return "showj", "/visualize_object/?path="+ path

def doPlotsCAPCA(request, protocolViewer):
    _, histogram = doPlotHistogram(request, protocolViewer)
    _, factorMaps = doPlotFactorMaps(request, protocolViewer)
    return "plots", [str(histogram) , str(factorMaps)]

def doPlotHistogram(request, protocolViewer):
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=plotHistogram&protViewerClass="+ protViewerClass + "&protId="+ protId

def plotHistogram(request, protocolViewer):
    xplotter = protocolViewer._plotHistogram()
    return xplotter

def doPlotFactorMaps(request, protocolViewer):
    first = str(protocolViewer.firstFactor.get())
    second = str(protocolViewer.secondFactor.get())
    protViewerClass = str(protocolViewer.getClassName())
    protId = str(protocolViewer.protocol.getObjId())
    return "plot","/view_plots/?function=plotFactorMaps&protViewerClass="+ protViewerClass + "&protId="+ protId +"&first="+ first +"&second="+second 

def plotFactorMaps(request, protocolViewer):
    protocolViewer.firstFactor.set(request.GET.get('first', None))
    protocolViewer.secondFactor.set(request.GET.get('second', None))
    xplotter = protocolViewer._plotFactorMaps()
    return xplotter

def doShowPcaFile(request, protocolViewer):
    html = textfileViewer("PCA files", [protocolViewer.protocol.imcFile.filename.get()])
    return "html", html
