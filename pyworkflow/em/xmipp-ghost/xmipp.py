

# FIXME: This is a bypass untill xmipp is replaced to xmippLib
# FIXME:   in all imports related to the binding

print('import xmipp is deprecated for the xmipp binding.\n'
      '  > Please change it to: import xmippLib ')

import traceback
traceback.print_stack()
from xmippLib import *