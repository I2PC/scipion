from pyworkflow.viewer import createPlots

############## VIEWER XMIPP ML2D ##############
def viewerML2D(request, protocol, protocolViewer):
    ioDict = {}

    if protocolViewer.doShowClasses:
        typeUrl, url = doShowClasses(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
    if protocolViewer.doShowPlots:
        typeUrl, url = doAllPlotsML2D(request, protocol, protocolViewer)
        ioDict[typeUrl]= url
    else:
        typeUrl, url = doSomePlotsML2D(protocol, protocolViewer)
        if url != None:
            ioDict[typeUrl]= url
        
    return ioDict

def doShowClasses(request, protocol, protocolViewer):
    objId = protocol.outputClasses.getObjId()
    return "showj","/visualize_object/?objectId="+str(objId)
            
def doAllPlotsML2D(request, protocol, protocolViewer):
    return "plots","/view_plots/?function=allPlotsML2D&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def doSomePlotsML2D(protocol, protocolViewer):
    plots=""
    for p in protocolViewer._plotVars:
        if protocolViewer.getAttributeValue(p):
            plots = plots + p + "-"
    
    if plots != "":
        return "plots","/view_plots/?function=somePlotsML2D&plots="+str(plots)+"&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())
    else:
        return "", None

def doShowLL(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=somePlotsML2D&plots=doShowLL&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def doShowPmax(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=somePlotsML2D&plots=doShowPmax&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def doShowSignalChange(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=somePlotsML2D&plots=doShowSignalChange&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())
    
def doShowMirror(request, protocol, protocolViewer):
    return "plot","/view_plots/?function=somePlotsML2D&plots=doShowMirror&protViewerClass="+ str(protocolViewer.getClassName())+ "&protId="+ str(protocol.getObjId())

def allPlotsML2D(request, protocol, protocolViewer):
    xplotter = createPlots(protocol, protocolViewer._plotVars)
    return xplotter

def somePlotsML2D(request, protocol, protocolViewer):
    plots = request.GET.get('plots', None)
    plots = str(plots).split("-")
    if len(plots) > 1:
        plots.remove(plots[len(plots)-1])
    
    xplotter = createPlots(protocol, plots)
    return xplotter