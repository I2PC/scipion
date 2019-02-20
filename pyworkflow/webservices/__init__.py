
# Export all configuration variables in case they will be needed outside this module
# Nonetheless, the recommend way is to import the classes or functions provided
# by the module
from config import *

from notifier import ProjectWorkflowNotifier
from repository import WorkflowRepository