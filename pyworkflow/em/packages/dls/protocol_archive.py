from pyworkflow.em import ProtProcessMovies
import json
import os
import pyworkflow.protocol.params as params


class ProtArchive(ProtProcessMovies):
    """
    Sends a message to the dials archiving queue
    """
    _label = 'DLS Archive'

    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

        group = form.addGroup('Message Queue')
        group.addParam('configuration', params.FileParam,
                       label="Configuration", important=True,
                       help='a file path for a StompTransport configuration configuration file')
        group.addParam('dropFileDir', params.FolderParam,
                       label='Dropfile dir',
                       help='The location of the dropfile dir where archive instructions will be sent')

    def processMovieStep(self, movieDict, hasAlignment):
        SCIPION_HOME = os.environ['SCIPION_HOME']
        LOGFILE = os.path.join(SCIPION_HOME, 'pyworkflow', 'em', 'packages','dls','archive_recipe.json')
        with open(LOGFILE) as recipe:
            recipe_data = json.load(recipe)

        real_path = os.path.realpath(movieDict['_filename'])

        self._setValues(recipe_data, 'filelist', [real_path ])
        self._setValues(recipe_data, 'dropfile-dir', self.dropFileDir.get())

        self._sendMessage(recipe_data)

    def _setValues(self, data, key, value):
        for k in data:
            if k == key:
                data[k] = value
            elif type(data[k]) is dict:
                self._setValues(data[k], key, value)

    def _sendMessage(self, message):
        # type: (str) -> void
        try:
            try:
                from workflows.transport.stomp_transport import StompTransport
            except ImportError as ex:
                print "You need stomp transport to run the archiver, run 'scipion run pip install workflows' first"
                raise ex
            StompTransport.load_configuration_file(self.configuration.get())
            self.stomp.send('archive.filelist', message, headers={ 'workflows-recipe': True })
        except AttributeError:
            self.stomp = StompTransport()
            self.stomp.connect()
            self._sendMessage(message)

    def createOutputStep(self):
        pass
