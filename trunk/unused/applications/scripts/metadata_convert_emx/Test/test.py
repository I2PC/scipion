from transformations import *
#shear = (0,0,0)
#persp = numpy.random.random(4) - 0.5
pi = math.pi

#first image
scale = (2.,1.,1.)
#scale = (1.,1.,1.)
rot, tilt, psi = pi/4.4,pi/2.2,pi/3.3
print rot*180./pi, tilt*180./pi, psi*180./pi
print rot, tilt, psi
#angles = (psi,tilt,rot)
#psi first rotation
angles = (rot,tilt,psi)
translate = (1,2,3)
MD = compose_matrix(scale=scale,
                    angles=angles,
		    translate=translate)
print "FIRST"
print MD
scale, shear, angles, trans, persp = decompose_matrix(MD)
print "scale", scale
#print "shear", shear
print "angles",angles
print "trans", trans
#print "persp", persp

print "%0.6f %0.6f %0.6f %0.6f"% (MD[0][0],MD[0][1],MD[0][2],MD[0][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[1][0],MD[1][1],MD[1][2],MD[1][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[2][0],MD[2][1],MD[2][2],MD[2][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[3][0],MD[3][1],MD[3][2],MD[3][3])

#second image
scale = (1.0,1.0,1.0)
#invert sign
print "SECOND"
rot, tilt, psi = -pi/5.,0,0
print rot*180./pi, tilt*180./pi, psi*180./pi
print rot, tilt, psi
angles = (rot,tilt,psi)
translate = (1.1,2.2,3.3)
MD = compose_matrix(scale=scale,
                    angles=angles,
		    translate=translate)
print MD
scale, shear, angles, trans, persp = decompose_matrix(MD)
print "scale", scale
print "shear", shear
print "angles",angles
print "trans", trans
print "persp", persp

print "%0.6f %0.6f %0.6f %0.6f"% (MD[0][0],MD[0][1],MD[0][2],MD[0][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[1][0],MD[1][1],MD[1][2],MD[1][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[2][0],MD[2][1],MD[2][2],MD[2][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[3][0],MD[3][1],MD[3][2],MD[3][3])

#third image
scale = (1.0,1.1,1.2)
rot, tilt, psi = pi/4.4,pi/2.2,pi/3.3
print "THIRD"
print rot*180./pi, tilt*180./pi, psi*180./pi
print rot, tilt, psi
angles = (rot,tilt,psi)
translate = (-1,+1,-1.1)
MD = compose_matrix(scale=scale,
                    angles=angles,
		    translate=translate)
print MD
scale, shear, angles, trans, persp = decompose_matrix(MD)
print "scale", scale
print "shear", shear
print "angles",angles
print "trans", trans
print "persp", persp

print "%0.6f %0.6f %0.6f %0.6f"% (MD[0][0],MD[0][1],MD[0][2],MD[0][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[1][0],MD[1][1],MD[1][2],MD[1][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[2][0],MD[2][1],MD[2][2],MD[2][3]),
print "%0.6f %0.6f %0.6f %0.6f"% (MD[3][0],MD[3][1],MD[3][2],MD[3][3])
