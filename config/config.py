# variables that may change between python version
# but do not depend on the user or system

# web that gathers protocol usage
SCIPION_STATS_SERVER = 'http://scipion.i2pc.es'
SCIPION_STATS_WORKFLOW_APP = SCIPION_STATS_SERVER + '/report_protocols/api/workflow/workflow/'

# web that handles workflows
WORKFLOW_REPOSITORY_SERVER='http://workflows.scipion.i2pc.es/'
WORKFLOW_PROG_STEP1='workflowProgStep1_add/'
WORKFLOW_PROG_STEP2='workflowProgStep2_add/'