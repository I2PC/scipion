from frm import *
from pytom.basic.structures import PyTomClass

class AngularConstraint(PyTomClass):
    """docstring for AngularConstraint"""
    FIXED_ANGLE = 'Fixed Angle'
    FIXED_AXIS = 'Fixed Axis'
    def __init__(self, type=FIXED_ANGLE):
        super(AngularConstraint, self).__init__()
        self.type = type
        self._cv = None
        self._bw = None

    def isViolated(self, angle):
        raise NotImplementedError("This virtual method must be overridden.")

    def getConstraintVolume(self, bw):
        raise NotImplementedError("This virtual method must be overridden.")

    def toXML(self):
        raise NotImplementedError("This virtual method must be overridden.")

    def fromXML(self, xml_obj):
        try:
            type = xml_obj.get('Type')
            if type == str(AngularConstraint.FIXED_ANGLE):
                c = FixedAngleConstraint(0,0,0,0)
            elif type == str(AngularConstraint.FIXED_AXIS):
                c = FixedAxisConstraint(0,0,0,0)
            else:
                raise Exception('Constraint type wrong!')

            c.fromXML(xml_obj)
            return c
        except Exception, e:
            raise e


class FixedAngleConstraint(AngularConstraint):
    """docstring for FixedAngleConstraint"""
    def __init__(self, phi, psi, the, nearby):
        """Consistent with PyTom notation: Z1 Z2 X
        """
        super(FixedAngleConstraint, self).__init__(AngularConstraint.FIXED_ANGLE)

        assert phi >= 0 and phi < 360 and psi >= 0 and psi < 360 and the >=0 and the < 180

        self.phi = phi
        self.psi = psi
        self.the = the
        self.nearby = nearby

    def isViolated(self, angle):
        pass

    def getConstraintVolume(self, bw):
        if self._cv is None or self._bw != bw:
            from math import pi
            cv = np.zeros((8*bw**3,), dtype='double')

            # the naming is inconsistent with the low-level c, but the result is right. To be changed.
            swig_frm.get_constraint_vol(cv, bw, self.psi*pi/180, self.phi*pi/180, self.the*pi/180, self.nearby*pi/180)

            self._cv = cv.reshape(2*bw, 2*bw, 2*bw)
            self._bw = bw

        return self._cv

    def toXML(self):
        from lxml import etree

        xml_obj = etree.Element('AngularConstraint', Type=str(self.type), Phi=str(self.phi), Psi=str(self.psi), Theta=str(self.the), Nearby=str(self.nearby))
        
        return xml_obj

    def fromXML(self, xml_obj):
        try:
            self.phi = float(xml_obj.get('Phi'))
            self.psi = float(xml_obj.get('Psi'))
            self.the = float(xml_obj.get('Theta'))
            self.nearby = float(xml_obj.get('Nearby'))
        except Exception, e:
            raise e


class FixedAxisConstraint(AngularConstraint):
    """docstring for FixedAxisConstraint"""
    def __init__(self, x, y, z, nearby=0):
        super(FixedAxisConstraint, self).__init__(AngularConstraint.FIXED_AXIS)
        self.x = x/np.linalg.norm([x,y,z])
        self.y = y/np.linalg.norm([x,y,z])
        self.z = z/np.linalg.norm([x,y,z])
        self.nearby = nearby
    
    def isViolated(self, angle):
        pass

    def getConstraintVolume(self, bw):
        if self._cv is None or self._bw != bw:
            cv = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')

            if self.nearby > 0: # get an array of axes nearby
                the = np.array(range(1, int(self.nearby*bw/90.)+1))*(90./bw)
                the = np.append(the, self.nearby)
                phi = range(0, 360, int(round(180./bw)))
                the, phi = np.meshgrid(the, phi)
                the = the.reshape((the.size,))*np.pi/180
                phi = phi.reshape((phi.size,))*np.pi/180
                xx = np.cos(phi)*np.sin(the)
                yy = np.sin(phi)*np.sin(the)
                zz = np.cos(the)

                # rotate all the vectors
                from pytom.angles.angleFnc import vector2euler, zxzToMat
                m = zxzToMat(vector2euler([self.x, self.y, self.z]))
                mat = np.matrix(np.zeros((3,3)))
                mat[0,0] = m[0,0]; mat[0,1] = m[0,1]; mat[0,2] = m[0,2]
                mat[1,0] = m[1,0]; mat[1,1] = m[1,1]; mat[1,2] = m[1,2]
                mat[2,0] = m[2,0]; mat[2,1] = m[2,1]; mat[2,2] = m[2,2]
                coord = mat*np.matrix(np.vstack((xx,yy,zz)))
                axes = coord.swapaxes(0,1)
            else:
                axes = np.matrix([[self.x, self.y, self.z]])

            # here we make the angular sampling as 2 degree
            from pytom.angles.angleFnc import axisAngleToZXZ
            for axis in axes:
                axis = [axis[0,0], axis[0,1], axis[0,2]]
                for ang in xrange(0, 360, 2):
                    euler_ang = axisAngleToZXZ(axis, ang)
                    i,j,k = frm_angle2idx(bw, euler_ang.getPhi(), euler_ang.getPsi(), euler_ang.getTheta())
                    # set the cv
                    cv[i,j,k] = 1

            self._cv = cv
            self._bw = bw

        return self._cv

    def toXML(self):
        from lxml import etree

        xml_obj = etree.Element('AngularConstraint', Type=str(self.type), X=str(self.x), Y=str(self.y), Z=str(self.z), Nearby=str(self.nearby))
        
        return xml_obj

    def fromXML(self, xml_obj):
        try:
            self.x = float(xml_obj.get('X'))
            self.y = float(xml_obj.get('Y'))
            self.z = float(xml_obj.get('Z'))
            self.nearby = float(xml_obj.get('Nearby'))
        except Exception, e:
            raise e


def frm_find_best_constrained_angle_interp(corr, b=None, constraint=None):
    if not b:
        b = corr.shape[0]/2
    if constraint is not None:
        corr = corr * constraint.getConstraintVolume(b)
    return frm_find_best_angle_interp(corr, b)


def frm_find_topn_constrained_angles_interp(corr, n=5, dist=3.0, constraint=None):
    b = corr.shape[0]/2

    # when the angular constraint conflicts the dist
    if constraint.__class__ == FixedAngleConstraint and float(dist+1)/b*180 > constraint.nearby:
        print 'Warning: angular distance cut is overwritten by angular constraint.'
        from math import floor
        dist = floor(constraint.nearby*b/180.)-1
        if dist < 1.:
            raise Exception("The angular constraint is too small. Please re-specify!")

    if constraint is not None:
        corr = corr * constraint.getConstraintVolume(b)
    return frm_find_topn_angles_interp(corr, n, dist)


def frm_constrained_align(vf, wf, vg, wg, b, max_freq, peak_offset=None, mask=None, constraint=None, weights=None, position=None, num_seeds=5):
    """Find the best alignment (translation & rotation) of volume vg (reference) to match vf.
    For details, please check the paper.

    Parameters
    ----------
    vf: Volume Nr. 1
        pytom_volume.vol

    wf: Mask of vf in Fourier space.
        pytom.basic.structures.Wedge. If none, no missing wedge.

    vg: Volume Nr. 2 / Reference
        pytom_volume.vol

    wg: Mask of vg in Fourier space.
        pytom.basic.structures.Wedge. If none, no missing wedge.

    b: Bandwidth range of spherical harmonics.
       None -> [4, 64]
       List -> [b_min, b_max]
       Integer -> [b, b]

    max_freq: Maximal frequency involved in calculation.
              Integer.

    peak_offset: The maximal offset which allows the peak of the score to be.
                 Or simply speaking, the maximal distance allowed to shift vg to match vf.
                 This parameter is needed to prevent shifting the reference volume out of the frame.
                 pytom_volume.vol / Integer. By default is half of the volume radius.

    mask: Mask volume for vg in real space.
          pytom_volume.vol

    constraint: Angular constraint
                sh_alignment.constrained_frm.AngularConstraint

    weights: Obsolete.

    position: If the translation is already known or not. If provided, no translational search will be conducted.
              List: [x,y,z], default None.

    num_seeds: Number of threads for the expectation maximization procedure. The more the better, yet slower.
               Integer, default is 5.

    Returns
    -------
    (The best translation and rotation (Euler angle, ZXZ convention [Phi, Psi, Theta]) to transform vg to match vf.
    (best_translation, best_rotation, correlation_score)
    """
    from pytom_volume import vol, rotateSpline, peak
    from pytom.basic.transformations import shift
    from pytom.basic.correlation import FLCF
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import Mask, SingleTiltWedge
    from pytom_volume import initSphere
    from pytom_numpy import vol2npy

    if vf.sizeX()!=vg.sizeX() or vf.sizeY()!=vg.sizeY() or vf.sizeZ()!=vg.sizeZ():
        raise RuntimeError('Two volumes must have the same size!')

    if wf is None:
        wf = SingleTiltWedge(0)
    if wg is None:
        wg = SingleTiltWedge(0)

    if peak_offset is None:
        peak_offset = vol(vf.sizeX(), vf.sizeY(), vf.sizeZ())
        initSphere(peak_offset, vf.sizeX()/4, 0,0, vf.sizeX()/2,vf.sizeY()/2,vf.sizeZ()/2)
    elif peak_offset.__class__ == int:
        peak_radius = peak_offset
        peak_offset = vol(vf.sizeX(), vf.sizeY(), vf.sizeZ())
        initSphere(peak_offset, peak_radius, 0,0, vf.sizeX()/2,vf.sizeY()/2,vf.sizeZ()/2)
    elif peak_offset.__class__ == vol:
        pass
    else:
        raise RuntimeError('Peak offset is given wrong!')

    # cut out the outer part which normally contains nonsense
    m = vol(vf.sizeX(), vf.sizeY(), vf.sizeZ())
    initSphere(m, vf.sizeX()/2, 0,0, vf.sizeX()/2,vf.sizeY()/2,vf.sizeZ()/2)
    vf = vf*m
    vg = vg*m
    if mask is None:
        mask = m
    else:
        vg = vg*mask

    if position is None: # if position is not given, we have to find it ourself
        # first roughtly determine the orientation (only according to the energy info)
        # get multiple candidate orientations
        numerator, denominator1, denominator2 = frm_correlate(vf, wf, vg, wg, b, max_freq, weights, True, None, None, False)
        score = numerator/(denominator1 * denominator2)**0.5
        res = frm_find_topn_constrained_angles_interp(score, num_seeds, get_adaptive_bw(max_freq, b)/16., constraint)
    else:
        # the position is given by the user
        vf2 = shift(vf, -position[0]+vf.sizeX()/2, -position[1]+vf.sizeY()/2, -position[2]+vf.sizeZ()/2, 'fourier')
        score = frm_correlate(vf2, wf, vg, wg, b, max_freq, weights, ps=False)
        orientation, max_value = frm_find_best_constrained_angle_interp(score, constraint=constraint)

        return position, orientation, max_value

    # iteratively refine the position & orientation
    from pytom.tools.maths import euclidianDistance
    max_iter = 10 # maximal number of iterations
    mask2 = vol(mask.sizeX(), mask.sizeY(), mask.sizeZ()) # store the rotated mask
    vg2 = vol(vg.sizeX(), vg.sizeY(), vg.sizeZ())
    lowpass_vf = lowpassFilter(vf, max_freq, max_freq/10.)[0]

    max_position = None
    max_orientation = None
    max_value = -1.0
    for i in xrange(num_seeds):
        old_pos = [-1, -1, -1]
        lm_pos = [-1, -1, -1]
        lm_ang = None
        lm_value = -1.0
        orientation = res[i][0] # initial orientation
        for j in xrange(max_iter):
            rotateSpline(vg, vg2, orientation[0], orientation[1], orientation[2]) # first rotate
            rotateSpline(mask, mask2, orientation[0], orientation[1], orientation[2]) # rotate the mask as well
            vg2 = wf.apply(vg2) # then apply the wedge
            vg2 = lowpassFilter(vg2, max_freq, max_freq/10.)[0]
            score = FLCF(lowpass_vf, vg2, mask2) # find the position
            pos = peak(score, peak_offset)
            pos, val = find_subpixel_peak_position(vol2npy(score), pos)
            if val > lm_value:
                lm_pos = pos
                lm_ang = orientation
                lm_value = val
        
            if euclidianDistance(lm_pos, old_pos) <= 1.0:
                # terminate this thread
                if lm_value > max_value:
                    max_position = lm_pos
                    max_orientation = lm_ang
                    max_value = lm_value

                break
            else:
                old_pos = lm_pos

            # here we shift the target, not the reference
            # if you dont want the shift to change the energy landscape, use fourier shift
            vf2 = shift(vf, -lm_pos[0]+vf.sizeX()/2, -lm_pos[1]+vf.sizeY()/2, -lm_pos[2]+vf.sizeZ()/2, 'fourier')
            score = frm_correlate(vf2, wf, vg, wg, b, max_freq, weights, False, denominator1, denominator2, True)
            orientation, val = frm_find_best_constrained_angle_interp(score, constraint=constraint)
            
        else: # no converge after the specified iteration, still get the best result as we can
            if lm_value > max_value:
                max_position = lm_pos
                max_orientation = lm_ang
                max_value = lm_value

        # print max_value # for show the convergence of the algorithm

    return max_position, max_orientation, max_value
    
