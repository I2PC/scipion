#----------------- SOFTWARE PACKAGES -----------------------
# We define a common root to install all software packages
# if the paths are relative, they will be taken from SCIPION_HOME

EM_ROOT = 'software/em/'

# Xmipp 3
XMIPP_HOME = EM_ROOT + 'xmipp'

# Spider 
SPIDER_DIR = EM_ROOT + 'spider'

# Eman2
EMAN2DIR = EM_ROOT + 'EMAN2'

# frealign 8.11
FREALIGN_HOME = EM_ROOT + 'brandeis/frealign_v9.07'

# Ctffind 3
CTFFIND_HOME = EM_ROOT + 'brandeis/ctffind'

# Relion
RELION_HOME = EM_ROOT + 'relion-1.2'

# xmipp_optical_alignment
OPT_ALIGN_HOME = EM_ROOT + 'opt_flow_alignment'

# BSOFT
BSOFT_HOME = EM_ROOT + 'bsoft'


# TODO: maybe it would be nicer to allow the user choose where to put
# the external packages, ask for it at installation time or something.
#
# Also, when a package is going to be used but cannot be found, say
# where the system expects it (and why).
