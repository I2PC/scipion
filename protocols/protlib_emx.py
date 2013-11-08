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
#from emx_data_model import *
from xmipp import *
from emx.emx import *
from emx.emxmapper import *
from numpy import eye
from protlib_filesystem import join, dirname, abspath, replaceBasenameExt
from protlib_xmipp import RowMetaData

BINENDING     = '.mrc'
CTFENDING     = '_ctf.param'
FAMILY        = 'DefaultFamily'
FIRSTIMAGE    =  1
MICFILE       = 'micrographs.xmd'
PARTFILE      = 'images.xmd'
POSENDING     = ".pos"
STACKENDING   = '.stk'


class CTF(object):
    # Dictionary for matching between EMX var and Xmipp labes
    ctfVarLabels = {'acceleratingVoltage': MDL_CTF_VOLTAGE,
                    'amplitudeContrast': MDL_CTF_Q0,
                    'cs': MDL_CTF_CS,
                    'defocusU': MDL_CTF_DEFOCUSU,
                    'defocusV': MDL_CTF_DEFOCUSV,
                    'defocusUAngle': MDL_CTF_DEFOCUS_ANGLE,
                    'pixelSpacing__X': MDL_CTF_SAMPLING_RATE,
                    }

#*******************************************************************
#*            EMX to Xmipp convertion functions                    *
#*******************************************************************

def emxMicsToXmipp(emxData, outputFileName=MICFILE, filesPrefix=None, ctfRoot=None):
    """ This function will iterate over the EMX micrographs and create
    the equivalent Xmipp 'micrographs.xmd' with the list of micrographs.
    If CTF information is found, for each micrograph, a file .ctfparam
    will be generated with the CTF parameters as expected by Xmipp.
    """
    #iterate though emxData
    mdMic = MetaData()
    if ctfRoot is None:
        ctfRoot = dirname(outputFileName)

    samplingRate = 0.
    voltage      = 0.
    cs           = 0.
    oldSamplingRate = -1.
    oldVoltage      = -1.
    oldCs           = -1.
    
    for micrograph in emxData.iterClasses(MICROGRAPH):
        micIndex = micrograph.get(INDEX)
        micFileName = micrograph.get(FILENAME)

        if micFileName is None:
            raise Exception("emxMicsToXmipp: Xmipp doesn't support Micrograph without filename")
        
        if filesPrefix is not None:
            micFileName = join(filesPrefix, micFileName)
        # micIndex is ignored now, in a more general solution
        # Xmipp should be able to handle micrographs in a stack
        # where the index has sense....

        mdMicId = mdMic.addObject()
        mdMic.setValue(MDL_MICROGRAPH, micFileName, mdMicId)
        ctfModelFileName = join(ctfRoot, replaceBasenameExt(micFileName, CTFENDING))
        mdMic.setValue(MDL_CTF_MODEL, ctfModelFileName, mdMicId)

        ctf = CTF()
        ctf.pixelSpacing__Y = micrograph.get('pixelSpacing__Y')
        
        # Set the variables in the dict
        for var in CTF.ctfVarLabels.keys():
            setattr(ctf, var, micrograph.get(var))

        # Now do some adjusments to vars        
        if ctf.defocusU is not None:
            ctf.defocusU *= 10.
        
        if ctf.defocusV is None:
            ctf.defocusV = ctf.defocusU
            ctf.defocusUAngle = 0.
        else:
            ctf.defocusV *= 10.
            
        if ctf.pixelSpacing__Y is not None:
            if ctf.pixelSpacing__X != ctf.pixelSpacing__Y:
                raise Exception ('pixelSpacingX != pixelSpacingY. Xmipp does not support it') 
            samplingRate    = ctf.pixelSpacing__X
            if (oldSamplingRate > -1):
                if oldSamplingRate != samplingRate:
                    raise Exception ('Xmipp emx import cannot import emx files with different samplingRate') 
            oldSamplingRate = samplingRate

        if ctf.voltage is not None:
            voltage         = ctf.voltage
            oldVoltage      = voltage
            if (oldVoltage > -1):
                if oldVoltage != voltage:
                    raise Exception ('Xmipp emx import cannot import emx files with different voltage') 

        if ctf.cs is not None:
            cs              = ctf.cs
            oldCs           = cs
            if (oldCs > -1):
                if oldCs != cs:
                    raise Exception ('Xmipp emx import cannot import emx files with different cs') 
            
        # Create the .ctfparam, replacing the micrograph name
        mdCtf = RowMetaData()
        
        for var, label in CTF.ctfVarLabels.iteritems():
            v = getattr(ctf, var)
            if v is not None:
                mdCtf.setValue(label, float(v))
        
        mdCtf.setValue(MDL_CTF_K, 1.0)
        mdCtf.write(ctfModelFileName)
        
    # Sort metadata by micrograph name
    mdMic.sort(MDL_MICROGRAPH)
    # Write micrographs metadata
    mdMic.write('Micrographs@' + outputFileName)
    return voltage, cs, samplingRate

def emxCoordsToXmipp(emxData, filesRoot):
    """ This function will iterate for each particle and 
    create several .pos metadata (per micrograph) where 
    the coordinates will written.
    """
    mdDict = {}
    # Create a dictionary for each micrograph
    for mic in emxData.iterClasses(MICROGRAPH):
        micFileName = mic.get(FILENAME)

        if micFileName is None:
            raise Exception("emxCoordsToXmipp: Xmipp doesn't support Micrograph without filename")
        
        mdDict[micFileName] = MetaData()
                  
    for part in emxData.iterClasses(PARTICLE):
        mic = part.getForeignObject(MICROGRAPH)
        micFileName = mic.get(FILENAME)
        md = mdDict[micFileName]
        objId = md.addObject()
        md.setValue(MDL_XCOOR, int(part.get('centerCoord__X')), objId)
        md.setValue(MDL_YCOOR, int(part.get('centerCoord__Y')), objId)        
    
    for mic in emxData.iterClasses(MICROGRAPH):
        micFileName = mic.get(FILENAME)
        md = mdDict[micFileName]
        mdFn = join(filesRoot, replaceBasenameExt(micFileName, POSENDING))
        md.write('coordinates@' + mdFn)
        
    part = emxData.iterClasses(PARTICLE)[0]
    if part.has('boxSize__X'):
        boxSize = part.get('boxSize__X')
        md = MetaData()
        objId = md.addObject()
        md.setValue(MDL_PICKING_PARTICLE_SIZE, int(boxSize), objId)
        md.write('properties@' + join(filesRoot, 'config.xmd'))

    
def emxParticlesToXmipp(emxData, outputFileName=PARTFILE, filesPrefix=None, ctfRoot=None):
    """ This function will iterate over the EMX particles and create
    the equivalent Xmipp 'images.xmd' with the list of particles.
    If CTF information is found, for each particle will contains information about the CTF
    """
    #iterate though emxData
    md    = MetaData()
    mdMic = MetaData()
    if ctfRoot is None:
        ctfRoot = dirname(outputFileName)
    
    samplingRate = 0.
    
    for particle in emxData.iterClasses(PARTICLE):
        pIndex = particle.get(INDEX)
        pFileName = particle.get(FILENAME)

        if pFileName is None:
            raise Exception("emxParticlesToXmipp: Xmipp doesn't support Particles without filename")
        
        if filesPrefix is not None:
            pFileName = join(filesPrefix, pFileName)
            
        if pIndex is not None:
            pFileName = '%06d@%s' % (pIndex, pFileName)
        # pIndex is ignored now, in a more general solution
        # Xmipp should be able to handle particles in a stack
        # where the index has sense....

        objId = md.addObject()
        md.setValue(MDL_IMAGE, pFileName, objId)
        
        mic = particle.getForeignObject(MICROGRAPH)

        if mic is not None:
            micFileName = mic.get(FILENAME)
            if filesPrefix is not None:
                micFileName = join(filesPrefix, micFileName)
            md.setValue(MDL_MICROGRAPH, micFileName, objId)



        if particle.has('centerCoord__X'):
            md.setValue(MDL_XCOOR, int(particle.get('centerCoord__X')), objId)
            md.setValue(MDL_YCOOR, int(particle.get('centerCoord__Y')), objId)

        if particle.has('centerCoord__X'):
            emxTransformToXmipp(md, objId, particle)
        
    # Sort metadata by particle name
    md.sort(MDL_IMAGE)
    # Write particles metadata
    md.write('Particles@' + outputFileName)   
    
    
def emxTransformToXmipp(md, objId, particle):
    """ Transform the particle transformation matrix
    in the euler angles and shift as expected by Xmipp.
    """
    #eulerangle2matrix
    _array = eye(3)# unitary matrix

    for i, j, label in iterTransformationMatrix():
        if particle.has(label):
            _array[i][j] = particle.get(label)
            
    rot, tilt, psi = Euler_matrix2angles(_array)
    
    _X = particle.get('transformationMatrix__t14', 0.)
    _Y = particle.get('transformationMatrix__t24', 0.)
    _Z = particle.get('transformationMatrix__t34', 0.)

    md.setValue(MDL_ANGLE_ROT , rot , objId)
    md.setValue(MDL_ANGLE_TILT, tilt, objId)
    md.setValue(MDL_ANGLE_PSI , psi , objId)

    md.setValue(MDL_SHIFT_X, _X, objId)
    md.setValue(MDL_SHIFT_Y, _Y, objId)
    md.setValue(MDL_SHIFT_Z, _Z, objId)
    

#*******************************************************************
#*            Xmipp to EMX convertion functions                    *
#*******************************************************************
    
def hasGeoLabels(md):
    """ Return true if there is geometric information. """
    for label in [MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_ANGLE_PSI, MDL_SHIFT_X, MDL_SHIFT_Y, MDL_SHIFT_Z]:
        if md.containsLabel(label):
            return True
    return False


def xmippCtfToEmx(md, objId, micrograph):
    """ Read a CTF file ctfparam from Xmipp 
    and set the attributes in micrograph.
    """
    ctfModel = md.getValue(MDL_CTF_MODEL, objId)
    mdCTF = RowMetaData(ctfModel)
    ctf = CTF()
    # Set the variables in the dict
    for var, label in CTF.ctfVarLabels.iteritems():
        setattr(ctf, var, mdCTF.getValue(label))
        
    ctf.defocusU /= 10.
    ctf.defocusV /= 10.

    while ctf.defocusUAngle < 0:
        ctf.defocusUAngle += 180.
    
    while ctf.defocusUAngle > 180.:
        ctf.defocusUAngle -= 180.; 
        
    for var in CTF.ctfVarLabels.keys():
        micrograph.set(var, getattr(ctf, var))           

    micrograph.set('pixelSpacing__Y', ctf.pixelSpacing__X)

    
def _writeEmxData(emxData, filename):
    xmlMapper = XmlMapper(emxData)
    xmlMapper.writeEMXFile(filename)      
    
    
def xmippMicrographsToEmx(micMd, emxData, emxDir):
    """ Export micrographs from xmipp metadata to EMX.
    """
    from protlib_particles import readPosCoordinates
    md = MetaData(micMd)
    micFn = 'mic%06d.mrc'
    index = 0
    pIndex = 0
    img = Image()
    hasCtf = md.containsLabel(MDL_CTF_MODEL)
    filesRoot = dirname(micMd)

    for objId in md:
        fnIn = md.getValue(MDL_MICROGRAPH, objId)
        index += 1
        fnOut = micFn % index
        img.read(fnIn)
        img.write(join(emxDir, fnOut))
        
        micrograph = Emxmicrograph(fnOut)
        # Set CTF parameters if present
        if hasCtf:
            xmippCtfToEmx(md, objId, micrograph)
            
        posFile = join(filesRoot, 'extra', replaceBasenameExt(fnIn, POSENDING))
        
        if exists(posFile):
            mdPos = readPosCoordinates(posFile)
            print "mdPos.size: ", mdPos.size()
            
            for pId in mdPos:
                pIndex += 1
                # We are using here a dummy filename, since not make sense
                # for particles with just coordinates
                particle = Emxparticle(str(pIndex), pIndex)
                particle.set('centerCoord__X', mdPos.getValue(MDL_XCOOR, pId))
                particle.set('centerCoord__Y', mdPos.getValue(MDL_YCOOR, pId))
                particle.setForeignKey(micrograph)
                emxData.addObject(particle)
                
        emxData.addObject(micrograph)
    # Write EMX particles
    _writeEmxData(emxData, join(emxDir, 'micrographs.emx'))
 

def xmippParticlesToEmx(imagesMd, emxData, emxDir):
    """ Export particles from xmipp metadata to EMX.
    imagesMd: the filename of the images metadata
    emxData: the emxData object to be populated.
    """
    md = MetaData(imagesMd)
    md.removeDisabled()
    
    imgFn = 'particles.mrc'
    mrcFn = join(emxDir, imgFn)
    index = 0
    img = Image()
    hasCoords = md.containsLabel(MDL_XCOOR)
    hasGeo = hasGeoLabels(md)
    hasMicrograph = md.containsLabel(MDL_MICROGRAPH)
    hasCtf = md.containsLabel(MDL_CTF_MODEL)
    micsDict = {}
    
    
    for objId in md:            
        # Read image from filename
        fn = FileName(md.getValue(MDL_IMAGE, objId))
        img.read(fn)
        
        # Write to mrc stack
        index += 1
        img.write('%d@%s' % (index, mrcFn))
        particle = Emxparticle(imgFn, index)
        # Create the micrograph object if exists
        if hasMicrograph:
            fn = md.getValue(MDL_MICROGRAPH, objId)
            if fn in micsDict:
                micrograph = micsDict[fn]
            else:
                idx, path = FileName(fn).decompose()
                idx = idx or None
                micrograph = Emxmicrograph(path, idx)
                emxData.addObject(micrograph)
                micsDict[fn] = micrograph 
            
            particle.setForeignKey(micrograph)
            
            if hasCtf:
                xmippCtfToEmx(md, objId, micrograph)
                particle.set('defocusU', micrograph.get('defocusU'))
                particle.set('defocusV', micrograph.get('defocusV'))

        # Check if there are coordinates
        if hasCoords:
            particle.set('centerCoord__X', md.getValue(MDL_XCOOR, objId))
            particle.set('centerCoord__Y', md.getValue(MDL_YCOOR, objId))
        
        # If there is geometrical info, convert and set
        if hasGeo:
            xmippTransformToEmx(md, objId, particle)
        
        # Add particle to emxData
        emxData.addObject(particle)
        
    # Write EMX particles
    _writeEmxData(emxData, join(emxDir, 'particles.emx'))
       
        
def xmippTransformToEmx(md, objId, particle):
    rot = md.getValue(MDL_ANGLE_ROT, objId) or 0.
    tilt = md.getValue(MDL_ANGLE_TILT, objId) or 0.
    psi = md.getValue(MDL_ANGLE_PSI , objId) or 0.
    
    tMatrix = Euler_angles2matrix(rot, tilt, psi)

    for i, j, label in iterTransformationMatrix():
        particle.set(label, tMatrix[i][j])

    x = md.getValue(MDL_SHIFT_X, objId) or 0.
    particle.set('transformationMatrix__t14', x)
    y = md.getValue(MDL_SHIFT_Y, objId) or 0.
    particle.set('transformationMatrix__t24', y)
    z = md.getValue(MDL_SHIFT_Z, objId) or 0.
    particle.set('transformationMatrix__t34', z)


def iterTransformationMatrix():
    """ This function will return an iterator over the indexes 
    and the EMX label for transformation matrix.
    """
    r = range(3)
    for i in r:
        for j in r:
            label = "transformationMatrix__t%d%d" % (i+1, j+1)
            yield (i, j, label)
    
    
    
    
################ Old Roberto convertion functions ############


def ctfMicXmippToEmx(emxData,xmdFileName):
    
    md    = MetaData()
    mdCTF = MetaData()
    md.read(xmdFileName)
    for objId in md:
        micrographName = md.getValue(MDL_MICROGRAPH, objId)
        ctfModel       = md.getValue(MDL_CTF_MODEL,objId)
        m1             = Emxmicrograph(fileName=micrographName)

        mdCTF.read(ctfModel)
        objId2 = mdCTF.firstObject()
        pixelSpacing        = mdCTF.getValue(MDL_CTF_SAMPLING_RATE, objId2)####
        acceleratingVoltage = mdCTF.getValue(MDL_CTF_VOLTAGE, objId2)
        cs                  = mdCTF.getValue(MDL_CTF_CS, objId2)
        defocusU            = mdCTF.getValue(MDL_CTF_DEFOCUSU, objId2)/10.
        defocusV            = mdCTF.getValue(MDL_CTF_DEFOCUSV, objId2)/10.
        defocusUAngle       = mdCTF.getValue(MDL_CTF_DEFOCUS_ANGLE, objId2)
        amplitudeContrast   = mdCTF.getValue(MDL_CTF_Q0, objId2)
        while defocusUAngle < 0:
            defocusUAngle += 180.
        while defocusUAngle > 180.:
            defocusUAngle -= 180.;            
        
        m1.set('acceleratingVoltage',acceleratingVoltage)
        m1.set('amplitudeContrast',amplitudeContrast)
        m1.set('cs',cs)
        m1.set('defocusU',defocusU)
        m1.set('defocusV',defocusV)
        m1.set('defocusUAngle',defocusUAngle)
        m1.set('pixelSpacing__X',pixelSpacing)
        m1.set('pixelSpacing__Y',pixelSpacing)
        emxData.addObject(m1)

def ctfMicXmippToEmxChallenge(emxData,xmdFileName):
    
    md    = MetaData()
    mdCTF = MetaData()
    md.read(xmdFileName)
    for objId in md:
        micrographName = md.getValue(MDL_MICROGRAPH, objId)
        ctfModel       = md.getValue(MDL_CTF_MODEL,objId)
        m1             = Emxmicrograph(fileName=micrographName)

        mdCTF.read(ctfModel)
        objId2 = mdCTF.firstObject()
        defocusU            = mdCTF.getValue(MDL_CTF_DEFOCUSU, objId2)/10.
        defocusV            = mdCTF.getValue(MDL_CTF_DEFOCUSV, objId2)/10.
        defocusUAngle       = mdCTF.getValue(MDL_CTF_DEFOCUS_ANGLE, objId2)
        amplitudeContrast   = mdCTF.getValue(MDL_CTF_Q0, objId2)
        while defocusUAngle < 0:
            defocusUAngle += 180.
        while defocusUAngle > 180.:
            defocusUAngle -= 180.;            
            
        m1.set('defocusU',defocusU)
        m1.set('defocusV',defocusV)
        m1.set('defocusUAngle',defocusUAngle)
        emxData.addObject(m1)
        

def alignXmippToEmx(emxData,xmdFileName):
    ''' convert a set of particles including geometric information '''
    md       = MetaData(xmdFileName)
    for objId in md:
        #get fileName
        fileName = FileName(md.getValue(MDL_IMAGE, objId))
        (index,fileName)=fileName.decompose()
        if index==0:
            index= None
        p1=Emxparticle(fileName=fileName, index=index)
        rot  = md.getValue(MDL_ANGLE_ROT ,objId)
        tilt = md.getValue(MDL_ANGLE_TILT,objId)
        psi  = md.getValue(MDL_ANGLE_PSI ,objId)
        if rot is  None:
            rot = 0
        if tilt is  None:
            tilt = 0
        if psi is  None:
            psi = 0
        tMatrix=Euler_angles2matrix(rot,tilt,psi)

        p1.set('transformationMatrix__t11', tMatrix[0][0])
        p1.set('transformationMatrix__t12', tMatrix[0][1])
        p1.set('transformationMatrix__t13', tMatrix[0][2])
        x = md.getValue(MDL_SHIFT_X, objId)
        p1.set('transformationMatrix__t14',x if x is not None else 0.)

        p1.set('transformationMatrix__t21', tMatrix[1][0])
        p1.set('transformationMatrix__t22', tMatrix[1][1])
        p1.set('transformationMatrix__t23', tMatrix[1][2])
        y = md.getValue(MDL_SHIFT_Y, objId)
        p1.set('transformationMatrix__t24',y if y is not None else 0.)

        p1.set('transformationMatrix__t31', tMatrix[2][0])
        p1.set('transformationMatrix__t32', tMatrix[2][1])
        p1.set('transformationMatrix__t33', tMatrix[2][2])
        z = md.getValue(MDL_SHIFT_Z, objId)
        p1.set('transformationMatrix__t34',z if z is not None else 0.)

        emxData.addObject(p1)
        
def alignEMXToXmipp(emxData,mode,xmdFileName):
    """import align related information
    DEPRECATED: Should be used the more general emxParticlesToXmipp
    """
    #iterate though emxData
    mdParticle     = MetaData()
    for particle in emxData.objLists[mode]:
        mdPartId   = mdParticle.addObject()

        partIndex     = particle.get(INDEX)
        partFileName  = particle.get(FILENAME)
        if partFileName is None:
            if partIndex is None:
                raise Exception("emxCoordsToXmipp: Particle has neither filename not index")
            else: #only index? for xmipp index should behave as filename
                partFileName = FileName(str(partIndex).zfill(FILENAMENUMBERLENGTH))
                partIndex    = None
                fileName = FileName(partFileName)
        elif partIndex is None:
            fileName = FileName(partFileName)
        #particle is a stack. Unlikely but not impossible
        else:
            fileName=FileName()
            fileName.compose(int(partIndex),partFileName)
        mdParticle.setValue(MDL_IMAGE , fileName , mdPartId)
        #eulerangle2matrix
        _array = eye(3)# unitary matrix

        t = particle.get('transformationMatrix__t11')
        _array[0][0]   = t if t is not None else 1.
        t = particle.get('transformationMatrix__t12')
        _array[0][1]   = t if t is not None else 0.
        t = particle.get('transformationMatrix__t13')
        _array[0][2]   = t if t is not None else 0.
        
        t = particle.get('transformationMatrix__t21')
        _array[1][0]   = t if t is not None else 0.
        t = particle.get('transformationMatrix__t22')
        _array[1][1]   = t if t is not None else 1.
        t = particle.get('transformationMatrix__t23')
        _array[1][2]   = t if t is not None else 0.

        t = particle.get('transformationMatrix__t31')
        _array[2][0]   = t if t is not None else 0.
        t = particle.get('transformationMatrix__t32')
        _array[2][1]   = t if t is not None else 0.
        t = particle.get('transformationMatrix__t33')
        _array[2][2]   = t if t is not None else 1.

        rot,tilt,psi = Euler_matrix2angles(_array)
        
        t = particle.get('transformationMatrix__t14')
        _X   = t if t is not None else 0.
        t = particle.get('transformationMatrix__t24')
        _Y   = t if t is not None else 0.
        t = particle.get('transformationMatrix__t34')
        _Z   = t if t is not None else 0.

        mdParticle.setValue(MDL_ANGLE_ROT , rot , mdPartId)
        mdParticle.setValue(MDL_ANGLE_TILT, tilt, mdPartId)
        mdParticle.setValue(MDL_ANGLE_PSI , psi , mdPartId)

        mdParticle.setValue(MDL_SHIFT_X, _X, mdPartId)
        mdParticle.setValue(MDL_SHIFT_Y, _Y, mdPartId)
        mdParticle.setValue(MDL_SHIFT_Z, _Z, mdPartId)
    mdParticle.write(xmdFileName)
    
    
def coorrXmippToEmx(emxData,xmdFileName):
    ''' convert a single file '''
    md    = MetaData()
    md.read(xmdFileName)
    xmdFileNameNoExt = FileName(xmdFileName.withoutExtension())
    xmdFileNameNoExtNoBlock = xmdFileNameNoExt.removeBlockName()
    micrographName = xmdFileNameNoExtNoBlock + BINENDING
    particleName   = xmdFileNameNoExtNoBlock + STACKENDING
    m1 = Emxmicrograph(fileName=micrographName)
    emxData.addObject(m1)
    counter = FIRSTIMAGE
    #check if MDL_XCOOR column exists if not 
    #print error no coordenates found. did you rememebr to include the block name in the filename
    for objId in md:

        coorX = md.getValue(MDL_XCOOR, objId)
        coorY = md.getValue(MDL_YCOOR, objId)
        p1    = Emxparticle(fileName=particleName, index=counter)
        p1.set('centerCoord__X',coorX)
        p1.set('centerCoord__Y',coorY)
        p1.setForeignKey(m1)
        emxData.addObject(p1)

        counter += 1
