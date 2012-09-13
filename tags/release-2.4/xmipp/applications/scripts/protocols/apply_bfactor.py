#!/usr/bin/env python
#apply bfactor to a vector of volumes
#
""" This utility boost up the high frequencies. Do not use the automated
    mode [default] for maps with resolutions lower than 12-15 Angstroms.
    It does not make sense to apply the Bfactor to the firsts iterations
    see http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Correct_bfactor
"""
def apply_bfactor(_DisplayReference_list,\
        bFactorExtension,\
        _SamplingRate,\
        _MaxRes,\
        _CorrectBfactorExtraCommand,\
        volExtension,\
        _mylog\
        ):
    import os
    if len(_CorrectBfactorExtraCommand)<1:
        _CorrectBfactorExtraCommand=' -auto '
    for name in _DisplayReference_list:
       xmipp_command='xmipp_correct_bfactor '
       aux_name = name.replace(bFactorExtension,'')
       if not os.path.exists(aux_name):
            print '* WARNING: '+ aux_name +' does not exist, skipping...'
       else:
            argument = ' -i ' + name.replace(bFactorExtension,'') +\
                       ' -o ' + name +\
                       ' -sampling ' + str(_SamplingRate)+\
                       ' -maxres '   + str (_MaxRes) +\
                       ' '
            xmipp_command = xmipp_command + argument
            xmipp_command = xmipp_command + ' ' + _CorrectBfactorExtraCommand
            _mylog.debug (xmipp_command)
            print "*************************************************"
            print "* " + xmipp_command
            os.system(xmipp_command)

