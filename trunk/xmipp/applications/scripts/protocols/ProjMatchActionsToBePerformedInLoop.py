import launch_job
import os, shutil

def execute_mask(_log, dict):
    _mylog = _log
    _mylog.debug("execute_mask")
    if (dict['DoMask']):
        print '*********************************************************************'
        print '* Mask the reference volume'
        command = ' -i ' + dict['reconstructedFileName'] + \
                  ' -o ' + dict['maskedFileName']

        if (dict['DoSphericalMask']):
            command += ' --mask circular -' + str(dict['maskRadius'])
        else:
            command += ' --mask ' + dict[userSuppliedMask]

        launch_job.launch_job("xmipp_mask",
                              command,
                              _mylog,
                              False, 1, 1, '')

    else:
        shutil.copy(dict['reconstructedFileName'], dict['maskedFileName'])
        _mylog.info("Skipped Mask")
        _mylog.info("cp " + dict['reconstructedFileName'] + " " + dict['maskedFileName'])
        print '*********************************************************************'
        print '* Skipped Mask'
