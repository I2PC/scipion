#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# {begin_of_header}

#{eval} expandCommentRun()

# {file}(*.descr) Input phantom motif
"""Input phantom motif that will be used to create the micrograph"""
InputMotif = "thinSlice.descr"

# Output root for files generation
Prefix = "micrograph"

# Micrograph size
MicrographSize = 3000

# Number of Micrographs
NumberOfMicrographs = 1

# {list}(random, from_file) How to create coordinates?
Coordinates = "random"

#{condition}(Coordinates=="random")  Number of Particles
""" Number of copies of the input motif"""
NumberOfElements = 100

# {expert}{condition}(Coordinates=="random")  Minimun distance between particles:
MinDist = 150

# {expert}{condition}(Coordinates=="random")  Minimun distance to borders:
BorderMinDist = 100

# {file}(*.txt){condition}(Coordinates=="from_file") Coordinates text (x y) file:
CoordinatesFile = "coords.txt"

# Rotate in-plane the motif?
DoRotate = True

# Include flipped particles?
""" Randomly select to particles to be flipped"""
DoFlip = False

# {condition}(DoFlip) Rate of flipped particles
"""
1 - all particles will be flipped
0 - no particle flippled
0.5 - half of the particles 
...and so on
"""
FlipRate = 0.4

# Add gaussian noise to coordinates?
AddNoiseU = True

# Add gaussian noise to tilted coordinates?
AddNoiseT = True

#{expert}{condition}(AddNoiseU || AddNoiseT) Stdev for noise?
NoiseSigma = 3.0

# Alpha untilted:
AlphaUListStr = "0"

# Alpha tilted:
AlphaTListStr = "0"

# Tilt angle:
TiltAngle = 45

#{hidden} Prefix to runs
CustomPrefix = "create_phantom_micrograph"

#{hidden} Usage of the program
Usage = """
Given an phantom description file, create a micrograph 
with several copies of the phantom, rotated and shifted
(can also be flipped)
"""
#------------------------------------------------------------------------------------------------
# {end_of_header} 

from os.path import join
from numpy import zeros, deg2rad, random, sin, cos, array, dot
from numpy.linalg import inv
from protlib_transformations import em_euler_matrix, translation_matrix, translation_from_matrix
from protlib_base import CustomProtocol, protocolMain
from protlib_utils import *
from protlib_filesystem import *
from xmipp import *

AlphaTList=[float(s) for s in AlphaTListStr.split()]
AlphaUList=[float(s) for s in AlphaUListStr.split()]

class CustomCreatePhantomMic(CustomProtocol):
    def __init__(self, scriptname, project):
        CustomProtocol.__init__(self, scriptname, project)
        self.ORoot = self.workingDirPath(self.Prefix)
       
    def defineSteps(self):
        # Insert the step to create the stack from a custom Python function(see function below)
        for i in range(self.NumberOfMicrographs):
            self.insertStep('createPhantomMicrograph', verifyfiles=[], WorkingDir=self.WorkingDir,
                            ORoot= "%s%03d" % (self.ORoot, i+1), micNo=(i+1))
        self.insertStep('createTiltPairsMd', verifyfiles=['tilted_pairs.xmd'],
                        WorkingDir=self.WorkingDir, ORoot=self.ORoot)

    def visualize(self):
        root = self.ORoot
        runShowJ("%(root)sU.mrc %(root)sT.mrc" % locals())
        
def createPhantomMicrograph(log, WorkingDir, ORoot, micNo):
    '''Function to create micrographs and particles metadata'''
    md = MetaData("block1@" + InputMotif)
    half = MicrographSize / 2
    v_half = array([half, half, 0, 0])
    MU = em_euler_matrix(0,0, AlphaUList[micNo-1], 'rzyz')
    M = em_euler_matrix(AlphaTList[micNo-1], TiltAngle, 0, 'rzyz')
    
    # Used later for projection
    dims = [float(MicrographSize) for i in range(3)]
    md.setValue(MDL_DIMENSIONS_3D, dims, md.firstObject())
    
    mdAll = MetaData()
    mdPos = MetaData() # Metadata with the coordinates of each particle
    mdPosT = MetaData() # same for tilted coordinates
    mdPos_geo = MetaData() # Metadata with the coordinates of each particle
    mdPosT_geo = MetaData() # same for tilted coordinates
    mdAlignedU = MetaData() # 2D alignment file untilted
    mdAlignedT = MetaData() # 2D alignment file tilted
    fnAlignedU = os.path.join(WorkingDir,"untilted.xmd")
    fnAlignedT = os.path.join(WorkingDir,"tilted.xmd")
    if os.path.exists(fnAlignedU):
        mdAlignedU.read(fnAlignedU)
    if os.path.exists(fnAlignedT):
        mdAlignedT.read(fnAlignedT)
    program = ""
    
    if DoRotate:
        # rotations = random.random_integers(0, 360, NumberOfElements)
        rotations = []
        step=360.0/(NumberOfElements-1)
        for i in range(NumberOfElements):
            rotations.append(i*step)
    else:
        rotations = zeros(NumberOfElements)
    
    if Coordinates == "random":
        xcoords, ycoords = generateRandomPositions(MicrographSize, NumberOfElements)
    elif Coordinates == "from_file":
        xcoords, ycoords = readPositionFromFile(CoordinatesFile)
        
    def runProgram(args):
        runJob(log, program, args)
        
    def addPos(md, x, y, nx, ny):
        xy=dot(MU,array([x,y,0,1]))
        xy[0]+=nx
        xy[1]+=ny
        objId = md.addObject()
        md.setValue(MDL_XCOOR, int(round(half + xy[0])), objId)
        md.setValue(MDL_YCOOR, int(round(half + xy[1])), objId) 
        return objId 
    
    def addPosT(md, x, y, nx, ny):
        xy = dot(M, array([x, y, 0, 1]))
        xy[0]+=nx
        xy[1]+=ny
        objId = md.addObject()
        md.setValue(MDL_XCOOR, int(round(half + xy[0])), objId)
        md.setValue(MDL_YCOOR, int(round(half + xy[1])), objId) 

    def addGeo(md, micNo, imgNo, newx, newy, x, y, nx, ny, rot=0, flip=False):
        objId = md.addObject()
        imgFn = "%06d@Images/Extracted/run_001/extra/micrograph%03dU.stk" % (imgNo, micNo)
        md.setValue(MDL_IMAGE, imgFn, objId)
        md.setValue(MDL_XCOOR, half + int(x), objId)
        md.setValue(MDL_YCOOR, half + int(y), objId)
        md.setValue(MDL_XCOOR_TILT, half + int(x), objId)
        md.setValue(MDL_YCOOR_TILT, half + int(y), objId)

        if flip:
           rot=-rot

        r = float(rot)
        rr = deg2rad(-r)
        cr = cos(rr)
        sr = sin(rr)

        shiftX = nx * cr + ny * sr
        shiftY = ny * cr - nx * sr
        if flip:
           shiftX=-shiftX
        
        md.setValue(MDL_SHIFT_X, shiftX, objId)
        md.setValue(MDL_SHIFT_Y, shiftY, objId)    
        md.setValue(MDL_ANGLE_PSI, -r, objId)    
        md.setValue(MDL_FLIP, flip, objId) 
        return objId 
    
    def addAlignmentU(md, micNo, imgNo, nx, ny, psi, flip):
        T=translation_matrix([nx,ny,0])
        if flip:
            psi=-psi
        angle=-float(psi)-AlphaUList[micNo-1]
        psi+=AlphaUList[micNo-1]
        R=em_euler_matrix(0,0,psi,'rzyz')
        M=T.dot(R)
        Rinv=inv(R)
        newT=translation_from_matrix(Rinv.dot(M))
        if flip:
            newT[0]=-newT[0]
        
        objId = md.addObject()
        imgFn = "%06d@Images/Extracted/run_001/extra/micrograph%03dU.stk" % (imgNo, micNo)
        md.setValue(MDL_IMAGE, imgFn, objId)
        md.setValue(MDL_SHIFT_X, float(newT[0]), objId)
        md.setValue(MDL_SHIFT_Y, float(newT[1]), objId)
        md.setValue(MDL_ANGLE_PSI, angle, objId)
        md.setValue(MDL_FLIP, flip, objId)

    def addAlignmentT(md, micNo, imgNo, nx, ny, psi, flip):
        #T=translation_matrix([nx,ny,0])
        #R=em_euler_matrix(AlphaT,0,0,'rzyz')
        #M=T.dot(R)
        #Rinv=inv(R)
        #newT=translation_from_matrix(Rinv.dot(M))
        newT=[nx,ny]
        #print("nx,ny="+str(nx)+","+str(ny)+" newT="+str(newT[0])+","+str(newT[1]))
        objId = md.addObject()
        imgFn = "%06d@Images/Extracted/run_001/extra/micrograph%03dT.stk" % (imgNo, micNo)
        md.setValue(MDL_IMAGE, imgFn, objId)
        md.setValue(MDL_SHIFT_X, float(newT[0]), objId)
        md.setValue(MDL_SHIFT_Y, float(newT[1]), objId)
        md.setValue(MDL_ANGLE_ROT, float(psi), objId)
        angle=TiltAngle
        if flip:
            angle+=180
        md.setValue(MDL_ANGLE_TILT, float(angle), objId)
        md.setValue(MDL_ANGLE_PSI, float(AlphaTList[micNo-1]), objId)

    def addGeoT(md, x, y, nx, ny):
        posT = dot(M, array([x, y, 0, 1]))
        PosT = dot(M, array([nx, ny, 0, 1]))
        addGeo(md, posT[0], posT[1], nPosT[0], nPosT[1])
    
    tmpDescr = ORoot + '_tmp.descr'
    program = "xmipp_phantom_transform"
    
    imgNo = 0
    for x, y, rot in zip(xcoords, ycoords, rotations):
        # Make a copy of the phantom file
        copyFile(log, InputMotif, tmpDescr)
        # Rotate the phantom
        runProgram(" -i %(tmpDescr)s --operation rotate_axis 0 0 1 %(rot)d" % locals())
        # Flip particles if necessary
        flip = DoFlip and random.rand() <= FlipRate
        if flip:            
            runProgram(" -i %(tmpDescr)s --operation rotate_axis 0 1 0 180" % locals())
        # Apply shift to motif phantom
        runProgram(" -i %(tmpDescr)s --operation shift %(x)d %(y)d 0" % locals())

        # Read the metadata with features to group all
        md2 = MetaData('block2@' + tmpDescr)
        mdAll.unionAll(md2)
        
        nx, ny = (0, 0)
        if AddNoiseU: # Add some noise to exact coordinates
          nx, ny = random.normal(0., NoiseSigma, 2).round()
          nx, ny = (2., 4.)
    
        addPos(mdPos, x, y, nx, ny)
        imgNo = imgNo + 1
        addGeo(mdPos_geo, micNo, imgNo, x+nx, y+ny, x, y, nx, ny, rot, flip)
        addAlignmentU(mdAlignedU, micNo, imgNo, nx, ny, rot, flip)
    
        nx, ny = (0, 0)
        if AddNoiseT: # Add some noise to exact coordinates
          nx, ny = random.normal(0., NoiseSigma, 2).round()
        addPosT(mdPosT, x, y, nx, ny)

        #addGeoT(mdPosT_geo, x, y, nx, ny)
        addAlignmentT(mdAlignedT, micNo, imgNo, nx, ny, rot, flip)
        
    # Write phantom .descr file
    outDescr = ORoot + '.descr'

    # Write the phantom metadata with all phantom particles
    md.write('block1@' + outDescr)
    mdAll.write('block2@' + outDescr, MD_APPEND)
    mdAlignedU.write("images@"+fnAlignedU)
    mdAlignedT.write("images@"+fnAlignedT)
    
    # Write .pos with coordinates
    def writePos(mdPos, outPos):
        family = "DefaultFamily"
#        mdFamily = MetaData()
#        objId = mdFamily.addObject()
#        mdFamily.setValue(MDL_PICKING_FAMILY, family, objId)
#        mdFamily.setValue(MDL_PICKING_MICROGRAPH_FAMILY_STATE, 'Manual', objId)
#        mdFamily.write('families@%(outPos)s' % locals())
        mdPos.write("%(family)s@%(outPos)s" % locals())
        
    writePos(mdPos, ORoot + 'U.pos')
    #print mdPos
    writePos(mdPosT, ORoot + 'T.pos')
    mdPos_geo.write(ORoot + 'U_geo.xmd')
    mdPosT_geo.write(ORoot + 'T_geo.xmd')
    
    program = "xmipp_phantom_project"
    outMicU = ORoot + "U.mrc"
    alpha= AlphaUList[micNo-1]
    runProgram(" -i %(outDescr)s -o %(outMicU)s --angles %(alpha)d 0 0" % locals())
    outMicT = ORoot + "T.mrc"
    
    tilt = TiltAngle # Register in locals() dictionary
    alpha = AlphaTList[micNo-1]
    runProgram(" -i %(outDescr)s -o %(outMicT)s --angles 0 %(tilt)d %(alpha)d" % locals())

def createTiltPairsMd(log, WorkingDir, ORoot):
    md = MetaData()
    for i in range(NumberOfMicrographs):
        f= "./%s%03dU.mrc" % (ORoot, i+1)

        objId = md.addObject()
        md.setValue(MDL_MICROGRAPH, f, objId)
        md.setValue(MDL_MICROGRAPH_TILTED, f.replace('U.mrc', 'T.mrc'), objId)
        md.setValue(MDL_ANGLE_Y, float(AlphaUList[i]), objId)
        md.setValue(MDL_ANGLE_Y2, float(AlphaTList[i]), objId)
        md.setValue(MDL_ANGLE_TILT, float(TiltAngle), objId)
    md.write('tilted_pairs.xmd')

    mdGeo = MetaData(os.path.join(WorkingDir,"untilted.xmd"))
    fn = 'classes.xmd'
    md.clear()    
    objId = md.addObject()
    md.setValue(MDL_REF, 1, objId)
    md.setValue(MDL_CLASS_COUNT, mdGeo.size(), objId)
    md.write('classes@%s' % fn)
    mdGeo.write('class%06d_images@%s' % (1, fn), MD_APPEND)
            
def generateRandomPositions(MicrographSize, NumberOfElements):
    xcoords = []
    ycoords = []
    half = MicrographSize / 2
       
    for i in range(NumberOfElements):
        # Find some random coordinates (x, y) that doesn't overlap
        badPos = True
        while badPos:
            x = random.random_integers(-half + BorderMinDist, half - BorderMinDist)
            y = random.random_integers(-half + BorderMinDist, half - BorderMinDist)
            badPos = False
        #Check overlapping with existing coordinates
            for xx, yy in zip(xcoords, ycoords):
                if abs(x - xx) < MinDist and abs(y - yy) < MinDist:
                    badPos = True
                    break
        xcoords.append(x)
        ycoords.append(y)
    
    return xcoords, ycoords

def readPositionFromFile(FileName):
    xcoords = []
    ycoords = []
    f = open(FileName)
    
    for l in f:
        parts = l.split()
        xcoords.append( int(float(parts[0])) )
        ycoords.append( int(float(parts[1])) )
        
    return xcoords, ycoords
        

if __name__ == '__main__':
    protocolMain(CustomCreatePhantomMic)
