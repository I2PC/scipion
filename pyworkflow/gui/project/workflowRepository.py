import os, sys
import poster.streaminghttp
import urllib2
import webbrowser
import json
import tempfile

sys.path.append(os.path.join(os.environ['SCIPION_HOME'], 'config'))
from config.config import WORKFLOW_REPOSITORY_SERVER, \
    WORKFLOW_PROG_STEP1, WORKFLOW_PROG_STEP2

def searchWorkflow():
    """ open browser with the workflow repository home page"""
    search_url = WORKFLOW_REPOSITORY_SERVER
    webbrowser.open(search_url)

def exportUploadProtocols(project, protocols):
    """ upload workflow to workflow repository"""

    def uploadWorkflow(upload_file_url, upload_metadata_url, jsonFileName):
        """ Upload workflow 'jsonFileName'to server upload_file_url
            First the file is uploaded, then the metadata is uploaded.
            The script uploads the file and then opens a browser for the metadata
            Note that the two steps are needed since noinitial value can be passed
            to a file field. poster module is needed. Poster is pure python
            so it may be added to the directory rather than installed if needed.

            The server is django a uses filefield and csrf_exempt.
            csrf_exempt disable csrf checking. filefield
        """
        # json file to upload
        params = dict(json=open(jsonFileName, 'rb'))

        # we are going to upload a file so this is a multipart
        # connection
        datagen, headers = poster.encode.multipart_encode(params)
        opener = poster.streaminghttp.register_openers()

        # create request and connect to server
        response = opener.open(urllib2.Request(upload_file_url,
                                               datagen, headers))
        # server returns dictionary as json with remote name of the saved file
        _dict = json.loads(response.read())

        # server has stored file in file named _dict['jsonFileName']
        # TODO: This version thing is a horrible hack but it is not
        # my fault. Scipion does not offer a better way to know the
        # version
        import imp
        currentDir = os.path.dirname(os.path.abspath(__file__))
        currentFile = os.path.join(currentDir,'../../../scipion')
        mod = imp.load_source('scipion', currentFile)
        version = mod.getVersion(False)
        # version hack end

        fileNameUrl = "?jsonFileName=%s&versionInit=%s"%\
                       (_dict['jsonFileName'], version) # use GET
        # open browser to fill metadata, fileNAme will be saved as session
        # variable. Note that I cannot save the file nave in the
        # session in the first conection because the browser changes
        # from urlib2 to an actual browser
        # so sessions are different
        webbrowser.open(upload_metadata_url+fileNameUrl)
        #delete temporary file
        os.remove(jsonFileName) if os.path.exists(jsonFileName) \
                                else None

    jsonFileName = os.path.join(tempfile.mkdtemp(), 'workflow.json')
    upload_file_url = WORKFLOW_REPOSITORY_SERVER + WORKFLOW_PROG_STEP1
    upload_metadata_url = WORKFLOW_REPOSITORY_SERVER + WORKFLOW_PROG_STEP2
    try:
        project.exportProtocols(protocols, jsonFileName)
        uploadWorkflow(upload_file_url, upload_metadata_url, jsonFileName)
        return 0
    except Exception as ex:
        #self.windows.showError(str(ex))
        return(ex)
