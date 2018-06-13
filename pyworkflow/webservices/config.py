# Variable related to the online services (workflow stats and repository) that Scipion reports to or communicate.
# These variables may change between Scipion versions but do not depend on the user or system

# Web that gathers protocol usage
SCIPION_STATS_SERVER = 'http://scipion.i2pc.es'
SCIPION_STATS_WORKFLOW_APP = SCIPION_STATS_SERVER + '/report_protocols/api/workflow/workflow/'

# Web that handles workflows
WORKFLOW_REPOSITORY_SERVER = 'http://workflows.scipion.i2pc.es/'
WORKFLOW_PROG_STEP1 = 'workflowProgStep1_add/'
WORKFLOW_PROG_STEP2 = 'workflowProgStep2_add/'