import os
import urllib2
import webbrowser
import json

import config


class WorkflowRepository(object):
    """ Manager to communicate with the workflow repository services.
    It will provide functions to:
    - Search workflows (open the url in a browser).
    - Upload a given workflow json file.
    """
    def __init__(self,
                 repositoryUrl=config.WORKFLOW_REPOSITORY_SERVER,
                 uploadFileSuffix=config.WORKFLOW_PROG_STEP1,
                 uploadMdSuffix=config.WORKFLOW_PROG_STEP2):
        self._url = repositoryUrl
        self._uploadFileUrl = repositoryUrl + uploadFileSuffix
        self._uploadMdUrl = repositoryUrl + uploadMdSuffix

    def search(self):
        """ Open the repository URL in a webbrowser. """
        webbrowser.open(self._url)

    def upload(self, jsonFileName):
        """ Upload a given workflow providing the path ot the json file.

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
        import poster.streaminghttp
        datagen, headers = poster.encode.multipart_encode(params)
        opener = poster.streaminghttp.register_openers()

        # create request and connect to server
        response = opener.open(urllib2.Request(self._uploadFileUrl, datagen,
                                               headers))
        # server returns dictionary as json with remote name of the saved file
        _dict = json.loads(response.read())

        version = os.environ['SCIPION_VERSION']
        # version hack end

        fnUrl = "?jsonFileName=%s&versionInit=%s" % (_dict['jsonFileName'],
                                                     version)  # use GET
        # open browser to fill metadata, fileNAme will be saved as session
        # variable. Note that I cannot save the file nave in the
        # session in the first connection because the browser changes
        # from urlib2 to an actual browser
        # so sessions are different
        webbrowser.open(self._uploadMdUrl + fnUrl)

    # try:
    #     jsonFileName = os.path.join(tempfile.mkdtemp(), 'workflow.json')
    #     project.exportProtocols(protocols, jsonFileName)
    #     uploadWorkflow(upload_file_url, upload_metadata_url, jsonFileName)
    #     return 0
    # except Exception as ex:
    #     #self.windows.showError(str(ex))
    #     return(ex)
