'''
/***************************************************************************
 * Authors:     RobertoMarabini (roberto@cnb.csic.es)
 *              Jose Miguel de la Rosa
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''
from emx_data_model import *
from xmipp import FileName, MetaData
MICFILE   = 'micrograph.xmd'
PARTFILE  = 'images.xmd'
CTFENDING = '_ctf.param'

def ctfMicXmippFromEmx(emxData,emxFileName,oroot=''):
        
    for micrograph in emxData.listMicrographs:
        index     = micrograph.get(INDEX)
        _fileName = fileName = micrograph.get(FILENAME)
        if fileName is None:
            if index is None:
                raise Exception("ctfMicXmippFromEmx: Micrograph has neither filename not index")
            else: #only index? for xmipp index should behave as filename
                _fileName = str(index)
                index    = None
        elif index is None:
            _fileName = micrograph.get('fileName')
        #micrograph is a stack. Unlikely but not impossible
        else:
            _fileName.compose(index,micrograph.get('fileName'))

        mdMic = MetaData()
        mdMicId = MDmic.addObject()
        mdMic.setValue(xmipp.MDL_MICROGRAPH, _fileName, mdMicId)
        ctfModelFileName = fileName.withoutExtension() + CTFENDING
        mdMic.setValue(xmipp.MDL_CTF_MODEL, mdMicId)

        acceleratingVoltage = micrograph.get('acceleratingVoltage')
        amplitudeContrast   = micrograph.get('amplitudeContrast')
        cs                  = micrograph.get('cs')
        defocusU            = micrograph.get('defocusU')
        defocusV            = micrograph.get('defocusV')
        if defocusV is None:
            defocusV = defocusU
            defocusUAngle = 0.
        else:
            defocusUAngle   = micrograph.get('defocusUAngle')
        pixelSpacing__XX    = micrograph.get('pixelSpacing__XX')
        pixelSpacing__YY    = micrograph.get('pixelSpacing__YY')
        if not pixelSpacing__YY is None:
            if pixelSpacing__XX != pixelSpacing__YY:
                raise Exception ('pixelSpacing__XX != pixelSpacing__YY. xmipp does not support it') 

        MD = MetaData()
        objId = MD.addObject()
        MD.setValue(xmipp.MDL_CTF_SAMPLING_RATE, float(pixelSpacing__XX), objId)
        MD.setValue(xmipp.MDL_CTF_VOLTAGE,       float(acceleratingVoltage), objId)
        MD.setValue(xmipp.MDL_CTF_CS,            float(Cs), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUSU,      float(defocusU), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUSV,      float(defocusV), objId)
        MD.setValue(xmipp.MDL_CTF_DEFOCUS_ANGLE, float(defocusUAngle), objId)
        MD.setValue(xmipp.MDL_CTF_Q0,            float(amplitudeContrast), objId)
        MD.setValue(xmipp.MDL_CTF_K,             1.0, objId)
        MD.write(ctfModelFileName)
    mdMic.write(MICFILE)