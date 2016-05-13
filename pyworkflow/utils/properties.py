# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *              Adrian Quintana (aquintana@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module defines the text used in the application.
"""

class Message():
    # Example Usage: 
    # MyMessage = Message()
    # print MyMessage.label

    # Header List
    VIEW_PROJECTS = 'Projects'
    VIEW_PROTOCOLS = 'Protocols'
    VIEW_DATA = 'Data'
    VIEW_UPLOAD = 'Upload'
    
    # Projects Template
    LABEL_PROJECTS = 'Projects'
    LABEL_CREATE_PROJECT = 'Create Project'
    TITLE_CREATE_PROJECT = 'Enter the project name'
    TITLE_CREATE_PROJECT_NAME = 'Project Name: '
    TITLE_EDIT_OBJECT = 'Edit Object properties'
    MESSAGE_DELETE_PROJECT = 'This will *delete* the project and all its *data*. Are you sure?'
    LABEL_DELETE_PROJECT = '[Delete]'
    TITLE_DELETE_PROJECT = 'Confirm project deletion'
    LABEL_RENAME_PROJECT = '[Rename]'
    TITLE_RENAME_PROJECT = 'Confirm project renaming'
    LABEL_CREATED = 'Created: '
    LABEL_MODIFIED = 'Modified: '
    
    # Project Content Template
    LABEL_PROJECT = 'Project '
    #-- Protocol Treeview --
    LABEL_WORKFLOW = 'Workflow View: '
    LABEL_PROTTREE_NONE = 'None'
    #-- Toolbar --
    LABEL_NEW = 'New'
    LABEL_NEW_ACTION = 'New     '
    LABEL_EDIT = 'Edit'
    LABEL_EDIT_ACTION = 'Edit     '
    LABEL_COPY = 'Copy'
    LABEL_COPY_ACTION = 'Copy   '
    LABEL_DELETE = 'Delete'
    LABEL_DELETE_ACTION = 'Delete    '
    LABEL_STEPS = 'Steps'
    LABEL_BROWSE = 'Browse'
    LABEL_BROWSE_ACTION = 'Browse '
    LABEL_DB = 'Db'
    LABEL_STOP = 'Stop'
    LABEL_STOP_ACTION = 'Stop execution'
    LABEL_ANALYZE = 'Analyze Results'
    LABEL_TREE = 'Tree'
    LABEL_SMALLTREE = 'Small Tree'
    LABEL_LIST = 'List'
    LABEL_REFRESH = 'Refresh'
    LABEL_DEFAULT = 'Default'
    LABEL_CONTINUE = 'Continue'
    LABEL_CONTINUE_ACTION = 'Approve continue'
    LABEL_EXPORT = 'Export'
    
    #-- Tabs --
    LABEL_DATA = 'Data'
    LABEL_SUMMARY = 'Summary'
    LABEL_INPUT = 'Input'
    LABEL_OUTPUT = 'Output'
    LABEL_COMMENT = 'Comments'
    
    LABEL_OBJSUMMARY = 'Object Summary'
    LABEL_OBJINFO = 'Info'
    LABEL_OBJCREATED = 'Created'
    LABEL_OBJLABEL = 'Label'
    
    LABEL_METHODS = 'Methods'
    LABEL_LOGS = 'Output Logs'
    LABEL_LOGS_OUTPUT = 'Output Log'
    LABEL_LOGS_ERROR = 'Error Log'
    LABEL_LOGS_SCIPION = 'Scipion Log'
    
    LABEL_RUNNAME = 'Run name'
    LABEL_EXECUTION = 'Run mode'
    LABEL_PARALLEL = 'Parallel'
    
    NO_INFO_SUMMARY = 'No summary information.'
    NO_INFO_METHODS = 'No methods information.'
    NO_INFO_LOGS = 'No logs information.'
    NO_SAVE_SETTINGS = 'Error try to save settings.'
    
    #-------- Protocol Form messages ----------
    LABEL_CITE = 'Cite'
    LABEL_HELP = 'Help'
    TEXT_HELP = 'The file selected will be uploaded to the project folder. If the file was uploaded before, It will be replaced.'
    LABEL_RUNNAME = 'Run name'
    LABEL_RUNMODE = 'Mode' 
    LABEL_HOST = 'Host'
    LABEL_THREADS = 'Threads'
    LABEL_MPI = 'MPI'
    LABEL_QUEUE = 'Use queue?'
    
    LABEL_EXPERT = 'Expert Level'
    LABEL_EXPERT_NORMAL = 'Normal'
    LABEL_EXPERT_ADVANCE = 'Advanced'
    LABEL_EXPERT_EXPERT = 'Expert'
    
    HELP_RUNMODE = """  
Normally, each protocol is composed by several atomic steps.
Each step could be computationally intensive, that's why
the *Continue* execution mode will try to continue from the
last completed step. On the other hand, the *Restart* mode
will clean the whole run directory and start from scratch.    
    """    
    HELP_MPI_THREADS = """  
Define the number of processors to be used in the execution.
*MPI*: This is a number of independent processes
       that communicate through message passing
       over the network (or the same computer).
*Threads*: This refers to different execution threads 
       in the same process that can share memory. They run in
       the same computer.     
    """
    HELP_USEQUEUE = """  
Select *Yes* if you want to submit the job to a Queue system.
The queue commands for launch and stop jobs should be configured
for the current host in the _hosts.conf_ file.    
    """     
    
    
    TITLE_NAME_RUN = ' Protocol Run: '
    TITLE_RUN = 'Run'
    TITLE_LABEL = 'Label'
    LABEL_OPT_COMMENT = 'Describe your run here...'
    TITLE_COMMENT = 'Comment'
    LABEL_RUN_MODE_RESUME = 'resume'
    LABEL_RUN_MODE_RESTART = 'restart'
    TITLE_EXEC = 'Execution'
    TITLE_BROWSE_DATA = 'Protocol data'
    LABEL_QUEUE_YES = 'Yes'
    LABEL_QUEUE_NO = 'No'
    LABEL_PARAM_YES = 'Yes'
    LABEL_PARAM_NO = 'No'
    LABEL_BUTTON_CLOSE = 'Close'
    LABEL_BUTTON_SAVE = 'Save'
    LABEL_BUTTON_EXEC = 'Execute'
    LABEL_BUTTON_VIS = 'Visualize'
    LABEL_BUTTON_WIZ = 'Wizard'
    LABEL_BUTTON_HELP = 'Help'
    LABEL_BUTTON_RETURN = 'Save'
    # VARS
    VAR_EXEC_HOST = 'hostName'
    VAR_EXPERT = 'expertLevel'
    VAR_MPI = 'numberOfMpi'
    VAR_QUEUE = '_useQueue'
    VAR_RUN_NAME = 'runName'
    VAR_COMMENT = 'comment'
    VAR_RUN_MODE = 'runMode'
    VAR_THREADS = 'numberOfThreads'
    
    
    LABEL_PATTERN = 'Pattern'
    TEXT_PATTERN = """\
Pattern (that can include wildcards) of the files to import.
For example:
  *data/particles/***.spi*
  *~/Micrographs/mic/***.mrc*"""
    ERROR_PATTERN_EMPTY = 'The *pattern* cannot be empty.'
    ERROR_PATTERN_FILES = 'There are no files matching the *pattern*'
    LABEL_CHECKSTACK = 'Check stack files?'
    LABEL_COPYFILES = 'Copy files?'
    LABEL_VOLTAGE = 'Microscope voltage (kV)'
    TEXT_VOLTAGE = "Microscope voltage"
    LABEL_SPH_ABERRATION = 'Spherical aberration (mm)'
    TEXT_SPH_ABERRATION = """\
Optical effect due to the increased refraction of light rays when they
strike the lens near its edge, in comparison with those that strike near
the center."""
    LABEL_AMPLITUDE = 'Amplitude Contrast'
    TEXT_AMPLITUDE = """\
Produced by the loss of amplitude (i.e. electrons) from the beam.

For a weak phase and weak amplitude object, the amplitude contrast ratio Qo
is automatically computed. It should be a positive number, typically between
0.05 and 0.3."""
    LABEL_PATTERNU = 'Pattern untilted'
    LABEL_PATTERNT = 'Pattern tilted'

    
    LABEL_SAMP_MODE = 'Sampling rate mode'
    TEXT_SAMP_MODE = """\
You can specify the sampling rate (pixel size) directly from the image
(A/pixel, Ts) or by specifying the magnification rate (M) and the scanner
pixel size (microns/pixel, Tm).

They are related by  Ts = Tm / M"""
    LABEL_SAMP_MODE_1 = 'From image'
    LABEL_SAMP_MODE_2 = 'From scanner'
    LABEL_SAMP_RATE = 'Pixel size ("sampling rate") (A/px)'
    TEXT_SAMP_RATE = "Pixel size"
    LABEL_MAGNI_RATE = 'Magnification rate'
    TEXT_MAGNI_RATE = """\
Electron optical magnification (M). It can be used to compute the Image Pixel
Size ("Sampling Rate") (Ts) using the Scanner Pixel Size (Tm), Ts = Tm / M.

It is used by a few programs like Ctffind or Frealign."""
    LABEL_SCANNED = 'Scanned pixel size (microns/px)'


    ERROR_IMPORT_VOL = 'importVolumes:There is not filePaths matching pattern'
    
    LABEL_CTF_ESTI = 'CTF Estimation'
    LABEL_INPUT = 'Input'
    LABEL_INPUT_MIC = 'Input Micrographs'
    LABEL_INPUT_PART = 'Input Particles'
    LABEL_INPUT_VOLS = 'Input Volumes'
    LABEL_INPUT_MOVS = 'Input Movies'
    LABEL_ALIG_PART = 'Write aligned particles?'
    TEXT_ALIG_PART = 'If set to *Yes*, the alignment will be applied to \n'+'input particles and a new aligned set will be created.'

    TEXT_NO_INPUT_MIC = 'No *Input Micrographs* selected.'
    TEXT_NO_CTF_READY = 'CTF of *Input Micrographs* not ready yet.'
    TEXT_NO_OUTPUT_CO = 'Output coordinates not ready yet.'
    ERROR_NO_EST_CTF = '_estimateCTF should be implemented'
    
    
    TITLE_LAUNCHED = 'Success'
    LABEL_LAUNCHED = 'The protocol was launched successfully.'
    LABEL_FOUND_ERROR = 'Errors found'
    TITLE_SAVED_FORM = 'Success'
    LABEL_SAVED_FORM = 'The protocol was saved successfully.'
    TITLE_DELETE_FORM = 'Confirm DELETE'
    TITLE_RESTART_FORM = 'Confirm RESTART'
    LABEL_DELETE_FORM = """
You are going to *DELETE* the run(s): 
  - %s
*ALL DATA* related will be permanently removed.

Do you really want to continue?'
"""
    LABEL_RESTART_FORM = """
You are going to *RESTART* the run: 
  - %s
*ALL DATA* related will be permanently removed.

Do you really want to continue?'
"""    
    TITLE_STOP_FORM = 'Confirm STOP'
    LABEL_STOP_FORM = 'Do you really want to *STOP* this run?'
    
    NO_VIEWER_FOUND = 'There is not viewer for protocol:' 
    
    TITLE_SAVE_OUTPUT = 'Save protocol output'
    LABEL_SAVE_OUTPUT = 'Do you wish to save protocol output?'
    
    #SHOWJ_WEB
    SHOWJ_TITLE = 'Showj'
    
    LABEL_RESLICE = 'Reslice'
    RESLICE_Z_NEGATIVE = 'Z Negative (Front)'
    RESLICE_Y_NEGATIVE = 'Y Negative (Top)'
    RESLICE_Y_POSITIVE = 'Y Positive (Bottom)'
    RESLICE_X_NEGATIVE = 'X Negative (Left)'
    RESLICE_X_POSITIVE = 'X Positive (Right)'
    
    LABEL_COLS = 'Cols'
    LABEL_ROWS = 'Rows'
    
    LABEL_MIRRORY = 'Invert Y Axis'
    LABEL_APPLY_TRANSFORM = 'Apply Transform Matrix'
    LABEL_ONLY_SHIFTS = 'Only Shifts'
    LABEL_WRAP = 'Wrap'
    
    LABEL_BLOCK_SELECTION = 'Select Block'
    LABEL_LABEL_SELECTION = 'Select Label'
    LABEL_VOLUME_SELECTION = 'Select Volume'
    
    LABEL_ENABLE = 'Enable'
    LABEL_DISABLE = 'Disable'
    LABEL_SELECT_ALL = 'Select All'
    LABEL_SELECT_FROM = 'Select From Here'
    LABEL_SELECT_TO = 'Select To Here'
    
    LABEL_DISPLAY_TABLE_CONFIG = 'Display Table Configuration'
    
    LABEL_LABEL = 'Label'
    LABEL_VISIBLE = 'Visible'
    LABEL_RENDER = 'Render'
    
    LABEL_BUTTON_OK = 'Ok'
    LABEL_BUTTON_CANCEL = 'Cancel'
    
    LABEL_THRESHOLD = 'Threshold:'
    
    ERROR_DIMENSIONS = 'Incorrect table width or height: '
    ERROR_WEBGL = 'Your web browser does not support or is not configured for WebGL. See [[http://get.webgl.org/][WebGL Support]] for more information.'
    
    TOOLTIP_SEARCH = 'Search a given world in the text. '
    TOOLTIP_REFRESH = 'Reload the content of the files in the viewer. '
    TOOLTIP_EXTERNAL = 'Open the viewer in an external window. '

    TITLE_PICK_GAUSS = 'Automatic gaussian picking'
    LABEL_PICK_GAUSS = 'Do you wish to perform an automatic gaussian picking for the remaining micrographs?'


# To get font awesome icons into png use: http://fa2png.io/
class Icon():
    # Project Content Template
    RUNS_TREE = 'fa-sitemap.png'
    RUNS_LIST = 'fa-bars.png'
    ACTION_NEW = 'fa-plus-circle.png'
    ACTION_EDIT = 'fa-pencil.png'
    ACTION_COPY = 'fa-files-o.png'
    ACTION_DELETE = 'fa-trash-o.png'
    ACTION_REFRESH = 'fa-refresh.png'
    # TODO: change action_steps icon - fa-codefork?
    ACTION_STEPS = 'fa-list-ul.png'
    ACTION_BROWSE = 'fa-folder-open.png'
    ACTION_DB = 'fa-database.png'
    ACTION_TREE = None
    ACTION_LIST = 'fa-bars.png'
    ACTION_STOP = 'fa-stop.png'
    ACTION_CONTINUE = 'fa-play-circle-o.png'
    ACTION_RESULTS = 'fa-eye.png'
    ACTION_CLOSE = 'fa-times.png'
    ACTION_SAVE = 'fa-save.png'
    ACTION_VISUALIZE = 'fa-eye.png'
    ACTION_WIZ = 'fa-magic.png'
    ACTION_HELP = 'fa-question-circle.png'
    ACTION_REFERENCES = 'fa-external-link.png'
    ACTION_EXPORT = 'fa-external-link.png'
    
    ACTION_SEARCH = 'fa-search.png'
    ACTION_EXECUTE = 'fa-cogs.png'
    ACTION_IN = 'fa-sign-in.png'
    ACTION_OUT = 'fa-sign-out.png'
    #Host template
    BUTTON_SELECT = 'fa-check.png'
    BUTTON_CLOSE = 'fa-times.png'
    BUTTON_CANCEL = 'fa-ban.png'
    BUTTON_SAVE = ACTION_SAVE
    BUTTON_PC = 'fa-laptop.png'
    
    ARROW_UP = 'fa-arrow-up.png'
    ARROW_LEFT = 'fa-arrow-left.png'
    BRUSH = 'fa-paint-brush.png'
    TAGS = 'fa-tags.png'
    HOME = 'fa-home.png'
    LIGHTBULB = 'fa-lightbulb-o.png'
    PLUS_CIRCLE = 'fa-plus-circle.png'


class Color():
    RED_COLOR = 'Firebrick' # Red color for background label
#    RED_COLOR = '#B22222'
    LIGHT_RED_COLOR = '#F3CBCB' # Very light red for row selection
    LIGHT_BLUE_COLOR = '#EAEBFF' # Very light blue for even rows
    LIGHT_GREY_COLOR = '#EAEBEC' # Light grey for background color in form, protocol, table header and west container
    LIGHT_GREY_COLOR_2 = '#F2F2F2' # Very light grey for odd rows, input background, etc
    DARK_GREY_COLOR = '#6E6E6E' # Very dark grey for project title, tubes, etc
    
    STATUS_SAVED = '#D9F1FA', 
    STATUS_LAUNCHED = '#D9F1FA', 
    STATUS_RUNNING = '#FCCE62', 
    STATUS_FINISHED = '#D2F5CB', 
    STATUS_FAILED = '#F5CCCB', 
    STATUS_INTERACTIVE = '#F3F5CB',
    STATUS_ABORTED = '#F5CCCB',
    #STATUS_SAVED = '#124EB0',

class colorText:
    """printing in colors, bold, etc,
       example: print colorText.BOLD + 'Hello World !' + color.END
    """
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

