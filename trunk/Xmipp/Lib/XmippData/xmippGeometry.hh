/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

/* ------------------------------------------------------------------------- */
/* GEOMETRY FUNCTIONS                                                        */
/* ------------------------------------------------------------------------- */
#ifndef _GEOMETRY_HH
#  define _GEOMETRY_HH

// structures to store a point and a plane
// w is a weight, quality assigned to that point
// if you do not want to use it set it to 1

#include "xmippMatrices3D.hh"
#ifndef FLT_EPSILON
#   define FLT_EPSILON 1.19209e-07
#endif
#include <vector>
#include <iostream>
struct fit_point{
   double x;
   double y;
   double z;
   double w;} ;
/*   
ostream &operator<<(ostream &os, const fit_point &s)
   {
   os << "(" << s.x << "," << s.y <<"," << s.z << ") " << s.w << endl;      
   }
*/
/**@name Geometry
   \begin{description}
   \item[Geometrical Operations]
       Project a vector to a plane
   \item[Euler Operations]
       Matrix <---> Angles, Euler angle corrections and transformations,
       Euler rotation
   \item[Intersections]
       ray with unit sphere, unit cylinder
   \end{description}
*/
//@{

/* ------------------------------------------------------------------------- */
/* Geometrical Operations                                                    */
/* ------------------------------------------------------------------------- */
/**@name Geometrical operations */
//@{
/** Project a point to a plane (direction vector, distance).
    Given the coordinates for a vector in R3, and a plane (defined by its
    direction vector and the minimum distance from the plane to the
    coordinate system origin). This function computes the perpendicular
    projection from the vector to the plane. This function has tried
    to be optimized in speed as it is used in core routines within huge
    loops, this is why the result is given as an argument, and why
    no check about the dimensionality of the vectors is performed. The
    routine performs faster if the result vector is already in R3.
    
    The following example projects the point P=(1,1,1) to the XY-plane
    storing the result in Pp (belonging to R3), the result is obviously
    Pp=(1,1,0).
    \\Ex:
    \begin{verbatim}
    matrix1D<double> Z=vector_R3(0,0,1), P=vector_R3(1,1,1), Pp(3);
    Uproject_to_plane(P,Z,0,Pp);
    cout << "After projecting: Pp=" << Pp.transpose() << endl;
    \end{verbatim}
    
    The starting U in the function name stands for the fact that the plane,
    and the point are in the same reference system (called by default
    Universal).
    
    The result and point vectors can be the same one.
    */
void Uproject_to_plane(const matrix1D<double> &point,
   const matrix1D<double> &direction, double distance,
   matrix1D<double> &result);

/** Project a vector to a plane (Euler angles).
    These planes are restricted to have 0 distance to the universal
    coordinate system. In this special case, a plane can be defined by 3
    Euler angles (this is specially suited for projections, where the
    projection plane is defined by its 3 Euler angles). Then, the coordinate
    system associated to the 3 Euler angles (let's call its vectors X',Y',Z')
    defines a projection plane where Z' is the direction vector, and X'Y'
    are in-plane vectors. Actually, X' coincides with the X vector in the
    \URL[matrix2D definition]{../../../Extra_Docs/Conventions.html#Physical}
    and Y' with the Y vector.
     
    The resulting vector is in R3, and the function has been optimized for
    speed, so the result is passed as a parameter. This function is based
    in the one which projects a point given the Euler matrix of the
    projection plane.
    If you are to project several points to the same plane, you should
    better use the project to plane function where you give the Euler matrix.
    
    The following example projects the point P=(1,1,1) to the XY-plane
    storing the result in Pp (belonging to R3), the result is obviously
    Pp=(1,1,0).
    \\Ex:
    \begin{verbatim}
    matrix1D<double> P=vector_R3(1,1,1), Pp(3);
    Uproject_to_plane(P,0,0,0,Pp);
    cout << "After projecting: Pp=" << Pp.transpose() << endl;
    \end{verbatim}
    
    The starting U in the function name stands for the fact that the plane,
    and the point are in the same reference system (called by default
    Universal).
    
    The result and point vectors can be the same one.
    */
void Uproject_to_plane(const matrix1D<double> &r,
   double rot, double tilt, double psi, matrix1D<double> &result);

/** Project a vector to a plane (Euler matrix).
    These planes are restricted to have 0 distance to the universal
    coordinate system. In this special case, a plane can be defined by 3
    Euler angles (this is specially suited for projections, where the
    projection plane is defined by its 3 Euler angles). Then, the coordinate
    system associated to the 3 Euler angles (let's call its vectors X',Y',Z')
    defines a projection plane where Z' is the direction vector, and X'Y'
    are in-plane vectors. Actually, X' coincides with the X vector in the
    \URL[matrix2D definition]{../../../Extra_Docs/Conventions.html#Physical}
    and Y' with the Y vector.
     
    The resulting vector is in R3, and the function has been optimized for
    speed, if the result vector is already 3 dimensional when entering the
    function, no resize is performed; and the result is passed as a parameter.
    If you are to project several points to the same plane, you should
    better use the project to plane function where you give the Euler matrix.
    
    The following example projects the point P=(1,1,1) to the XY-plane
    storing the result in Pp (belonging to R3), the result is obviously
    Pp=(1,1,0).
    \\Ex:
    \begin{verbatim}
    matrix1D<double> P=vector_R3(1,1,1), Pp(3);
    matrix2D<double> euler=Euler_angles2matrix(0,0,0);
    Uproject_to_plane(P,euler,Pp);
    cout << "After projecting: Pp=" << Pp.transpose() << endl;
    \end{verbatim}
    
    The starting U in the function name stands for the fact that the plane,
    and the point are in the same reference system (called by default
    Universal).

    The result and point vectors can be the same one.
    */
void Uproject_to_plane(const matrix1D<double> &r,
   const matrix2D<double> &euler, matrix1D<double> &result);

/** Spherical distance.
    This function returns the distance over a sphere, not the straight line
    but the line which goes from one point to the other going over the
    surface of a sphere, supposing that both points lie on the same sphere. */
double spherical_distance(const matrix1D<double> &r1, const matrix1D<double> &r2);

/** Point to line distance in 3D.
    Let a line in 3-D be specified by the point a and the vector v,
    this fuction returns the minimum distance of this line to the point p.
     */
double point_line_distance_3D(const matrix1D<double> &p, 
                              const matrix1D<double> &a,
			      const matrix1D<double> &v);

/** Point to plane distance in 3D.
    Let a plane in 3-D be specified by the point a and the perpendicular vector v,
    this fuction returns the minimum distance of this plane to the point p.
     */
double point_plane_distance_3D(const matrix1D<double> &p, 
                               const matrix1D<double> &a,
			       const matrix1D<double> &v);
/** Least-squares-fit a plane to an arbitrary number of (x,y,z) points
    PLane described as Ax + By + C = z
    Returns -1  if  A²+B²+C² <<1
    
    Points are defined using the stuct
        \\Ex:
    \begin{verbatim}
    struct fit_point{
       double x;
       double y;
       double z;
       double w;};
    \end{verbatim}
     where w is a weighting factor. Set it to 1 if you do not want to use it
     */
void least_squares_plane_fit(  const vector<fit_point> &IN_points, 
                               double &plane_A,
			       double &plane_B,
			       double &plane_C);
         
/** Rectangle which encloses a deformed rectangle.
    Given a rectangle characterized by the top-left corner and the right-bottom
    corner, and given a matrix after which the rectangle is deformed. Which
    is the minimum rectangle which encloses the preceeding one? This function
    is useful for stablishing for loops which will cover for sure the deformed
    rectangle. All vectors are supposed to be 2x1 and the deformation matrix
    is 2x2. The corner (x0,y0) goes to V*(x0,y0)' and (xF,yF) to V*(xf,yF)'.
    After that you can make a loop from corner1 to corner2.
    
    The v0 and vF vectors can be reused as outputs.*/
void rectangle_enclosing(const matrix1D<double> &v0, const matrix1D<double> &vF,
    const matrix2D<double> &V, matrix1D<double> &corner1,
    matrix1D<double> &corner2);

/** Box which encloses a deformed box.
    Given a box characterized by the top-left corner (most negative)
    and the right-bottom (most positive)
    corner, and given a matrix after which the box is deformed. Which
    is the minimum box which encloses the preceeding one? This function
    is useful for stablishing for loops which will cover for sure the deformed
    box. All vectors are supposed to be 3x1 and the deformation matrix
    is 3x3. The corner (x0,y0,z0) goes to V*(x0,y0,z0)' and (xF,yF,zF) to
    V*(xf,yF,zF)'.
    After that you can make a loop from corner1 to corner2.
    
    The v0 and vF vectors can be reused as outputs.*/
void box_enclosing(const matrix1D<double> &v0, const matrix1D<double> &vF,
    const matrix2D<double> &V, matrix1D<double> &corner1,
    matrix1D<double> &corner2);

/**  Line Plane Intersection 
Let ax+by+cz+D=0 be the equation of your plane
(if your plane is defined by a normal vector N + one point M, then
(a,b,c) are the coordinates of the normal N, and d is calculated by using
the coordinates of M in the above equation).

Let your line be defined by one point P(d,e,f) and a vector V(u,v,w), the
points on your line are those which verify

x=d+lu
y=e+lv
z=f+lw
where l takes all real values.

for this point to be on the plane, you have to have

ax+by+cz+D=0, so,

a(d+lu)+b(e+lv)+c(f+lw)+D=0
that is

l(au+bv+cw)=-(ad+be+cf+D)

note that, if au+bv+cw=0, then your line is either in the plane, or
parallel to it... otherwise you get the value of l, and the intersection
has coordinates:
x=d+lu
y=e+lv
z=f+lw
where

l = -(ad+be+cf+D)/(au+bv+cw)

a= XX(normal_plane);
b= YY(normal_plane);
c= ZZ(normal_plane);
D= intersection_point;

d=XX(point_line)
e=YY(point_line)
f=ZZ(point_line)

u=XX(vector_line);
v=YY(vector_line);
w=ZZ(vector_line);

XX(intersection_point)=x;
YY(intersection_point)=y;
ZZ(intersection_point)=z;

return 0 if sucessful
return -1 if line paralell to plane
return +1 if line in the plane

*/
int line_plane_intersection(const matrix1D<double> normal_plane,
                         const matrix1D<double> vector_line, 
			 matrix1D<double> &intersection_point,
                         const matrix1D<double> point_line,
                         double point_plane_at_x_y_zero=0.);
//@}

/* ------------------------------------------------------------------------- */
/* Euler Operations                                                          */
/* ------------------------------------------------------------------------- */
/**@name Euler Operations */
//@{
/** Getting the Euler angles to a range.
    Adjust angles until they are in the range
    \begin{verbatim}
    0 < rot  < 360
    0 < tilt < 360
    0 < psi  < 360
    \end{verbatim}
    No direction equivalence is applied, ie, there is no correction of the
    direction of projection making use that a view from the top is
    the same as a view from the bottom but reversed ... Just a wrapping
    of the angles is done until the angles fall in the specified ranges.
    The angles given must be variables and they are modified with the new
    values.
    \\Ex: EULER_CLIPPING(rot,tilt,psi); */
#define EULER_CLIPPING(rot,tilt,psi) \
    rot=realWRAP(rot,0,360); \
    tilt=realWRAP(tilt,0,360); \
    psi=realWRAP(psi,0,360);

/** Getting the Euler angles to a range (0-2*PI).
    The same as before but the angles are expressed in radians. */
#define EULER_CLIPPING_RAD(rot,tilt,psi) \
    rot=realWRAP(rot,0,2.*PI); \
    tilt=realWRAP(tilt,0,2.*PI); \
    psi=realWRAP(psi,0,2.*PI);

/** Euler angles --> "Euler" matrix.
    This function returns the transformation matrix associated to the
    3 given Euler angles (in degrees). The returned matrix is the
    following one:
    \begin{verbatim}
    alpha, first rotation around Z (right-hand thumb along axis)
    beta, second rotation around Y (right-hand thumb along axis)
    gamma, third rotation around Z (right-hand thumb too)

    ca = cos(alpha); cb = cos(beta); cg = cos(gamma);
    sa = sin(alpha); sb = sin(beta); sg = sin(gamma);
    cc = cb*ca; cs = cb*sa;
    sc = sb*ca; ss = sb*sa;

    A(0,0) =  cg*cc-sg*sa; A(0,1) =  cg*cs+sg*ca; A(0,2) = -cg*sb;
    A(1,0) = -sg*cc-cg*sa; A(1,1) = -sg*cs+cg*ca; A(1,2) = sg*sb;
    A(2,0) =  sc;          A(2,1) =  ss;          A(2,2) = cb;
    \end{verbatim}
    In the Euler matrix, the rows define a coordinate system for the
    rotated volume, inside the volume this coordinate system is still
    (1,0,0), (0,1,0) and (0,0,1), but in the universal coordinate system
    the directions of these 3 vectors are defined by the rows of the
    mentioned Euler matrix respectively. Notice that these 3 vectors are
    the result of rotating the original universal coordinate system
    as a block following the instructions given by the 3 Euler angles.

    This matrix can be used to relate both coordinate systems, the
    universal one and another defined within an object which is rotated
    after these angles
    \begin{verbatim}
       Rv=Euler*Ru
    \end{verbatim}
    that means that a universal vector (Ru) is seen in the volume coordinate
    system as Rv.
    
    As an implementation note you might like to know that this function
    calls always to matrix2D::resize */
void Euler_angles2matrix(double alpha, double beta, double gamma,
   matrix2D<double> &A);
   
/** Let be two volumes f and g related by g(x,y,z) = f(D(x,y,z)) (where D is a
  lineal transformation) then the projection direction parallel to the vector w
  in f is going to be related with the projection direction parallel to the
  vector w_prime in g. Given the w Euler angles this routine provide the
  w_prime angles */
void Euler_Angles_after_compresion(const double rot, double tilt, double psi,
   double &new_rot, double &new_tilt, double &new_psi, matrix2D<double> &D);

/** Euler direction.
    This function returns  a vector parallel to the  projection direction.
    Resizes v if needed */
void Euler_direction(double alpha, double beta, double gamma,
    matrix1D<double> &v);

/** Euler direction2angles.
    This function returns the 3 Euler angles associated to the direction
    given by the vector v. The 3rd Euler angle is set always to 0 */
void Euler_direction2angles(matrix1D<double> &v,
   double &alpha, double &beta, double &gamma);

/** "Euler" matrix --> angles.
    This function compute a set of Euler angles which result in an "Euler"
    matrix as the one given. See \Ref{Euler_angles2matrix} to know more
    about how this matrix is computed and what each row means. The result
    angles are in degrees. Alpha, beta and gamma are respectively the
    first, second and third rotation angles. If the input matrix is not
    3x3 then an exception is thrown, the function doesn't check that
    the Euler matrix is truly representing a coordinate system.
    \\Ex: Euler_matrix2angles(Euler, alpha, beta, gamma); */
void Euler_matrix2angles(matrix2D<double> &A, double &alpha, double &beta,
   double &gamma);

/** Up-Down projection equivalence.
    As you know a projection view from a point has got its homologous from
    its diametrized point in the projection sphere. This function takes
    a projection defined by its 3 Euler angles and computes an equivalent
    set of Euler angles from which the view is exactly the same but in the
    other part of the sphere (if the projection is taken from the bottom
    then the new projection from the top, and viceversa). The defined
    projections are exactly the same except for a flip over X axis, ie,
    an up-down inversion. Exactly the correction performed is:
    \begin{verbatim}
    newrot  = rot;
    newtilt = tilt+180;
    newpsi  = -(180+psi);
    \end{verbatim}
    Ex: Euler_up_down(rot,tilt,psi,newrot, newtilt, newpsi); */
void Euler_up_down(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi);

/** The same view but differently expressed.
    As you know a projection view from a point can be expressed with
    different sets of Euler angles. This function gives you another
    expression of the Euler angles for this point of view. Exactly
    the operation performed is:
    \begin{verbatim}
    newrot  = rot+180;
    newtilt = -tilt;
    newpsi  = -180+psi;
    \end{verbatim}
    Ex: Euler_another_set(rot,tilt,psi,newrot, newtilt, newpsi); */
void Euler_another_set(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi);

/** Mirror over Y axis.
    Given a set of Euler angles this function returns a new set which define
    a mirrored (over Y axis) version of the former projection.
    \begin{verbatim}
     -----> X               X<------
    |                              |
    |                              |
    |               ======>        |
    v                              v
    Y                              Y
    \end{verbatim}
    The operation performed is
    \begin{verbatim}
    newrot  = rot;
    newtilt = tilt+180;
    newpsi  = -psi;
    \end{verbatim}
    Ex: Euler_mirrorY(rot,tilt,psi,newrot, newtilt, newpsi); */
void Euler_mirrorY(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi);

/** Mirror over X axis.
    Given a set of Euler angles this function returns a new set which define
    a mirrored (over X axis) version of the former projection.
    \begin{verbatim}
     -----> X               Y  
    |                       ^
    |                       |
    |               ======> |
    v                       |
    Y                        -----> X
    \end{verbatim}
    The operation performed is
    \begin{verbatim}
    newrot  = rot;
    newtilt = tilt+180;
    newpsi  = 180-psi;
    \end{verbatim}
    Ex: Euler_mirrorX(rot,tilt,psi,newrot, newtilt, newpsi); */
void Euler_mirrorX(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi);

/** Mirror over X and Y axes.
    Given a set of Euler angles this function returns a new set which define
    a mirrored (over X and Y axes at the same time) version of the former
    projection.
    \begin{verbatim}
     -----> X                       Y  
    |                               ^
    |                               |
    |               ======>         |
    v                               |
    Y                        X<-----
    \end{verbatim}
    The operation performed is
    \begin{verbatim}
    newrot  = rot;
    newtilt = tilt;
    newpsi  = 180+psi;
    \end{verbatim}
    Ex: Euler_mirrorX(rot,tilt,psi,newrot, newtilt, newpsi); */
void Euler_mirrorXY(double rot, double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi);

/** Apply a geometrical transformation.
    The idea behind this function is the following. 3 Euler angles define
    a point of view for a projection, but also a coordinate system. You
    might apply a goemetrical transformation to this system, and then
    compute back what the Euler angles for the new system are. This
    could be used to "mirror" points of view, rotate them and all the
    stuff. The transformation matrix must be 3x3 but it must transform
    R3 vectors into R3 vectors (that is a normal 3D transformation matrix
    when vector coordinates are not homogeneous) and it will be applied
    in the sense:
    \begin{verbatim}
    New Euler matrix = L * Old Euler matrix * R
    \end{verbatim}
    where you know that the Euler matrix rows represent the different
    system axes. See \Ref{Euler_angles2matrix} for more information
    about the Euler coordinate system.
    \\Ex:
    \begin{verbatim}
    matrix2D<double> R60=rot3D_matrix(60,'Z');
    R60.resize(3,3);  // Get rid of homogeneous part
    matrix2D<double> I(3,3); I.init_identity();
    Euler_apply_transf(L,rot,tilt,psi,newrot,newtilt,newpsi);
    \end{verbatim} */
void Euler_apply_transf (const matrix2D<double> &L, const matrix2D<double> &R,
   double rot,     double tilt, double psi,
   double &newrot, double &newtilt, double &newpsi);

/** 3D Rotation matrix after 3 Euler angles.
    Creates a rotational matrix (4x4) for volumes around the combination
    of the 3 rotations around ZYZ. All angles are in degrees.
    You must use it with IS_NOT_INV in \Ref{apply_geom}.
    \\Ex: matrix2D<float> euler=Euler_rot3D_matrix(60,30,60);*/
matrix2D<double> Euler_rot3D_matrix(double rot, double tilt, double psi);

/** Rotate a volume after 3 Euler angles.
    The following prototype of the function is faster. */
matrix3D<double> Euler_rotate(const matrix3D<double> &V,
   double rot, double tilt, double psi);

/** Rotate a volume after 3 Euler angles.
    Input and output volumes cannot be the same one.*/
void Euler_rotate(const matrix3D<double> &V,
   double rot, double tilt, double psi, matrix3D<double> &result);
//@}

/**@name Intersections */
//@{
/** Intersection of a ray with a unit sphere.
    The sphere is centered at (0,0,0) and has got unit radius. The ray
    is defined by its direction (u) and a passing point (r). The function
    returns the length of the intersection. If the ray is tangent to the
    sphere the length is 0. See \Ref{Ellipsoid} to know
    how you can intersect any ray with any ellipsoid/sphere. */
double intersection_unit_sphere(const matrix1D<double> &u, 
   const matrix1D<double> &r);

/** Intersection of a ray with a unit cylinder.
    The cylinder is centered at (0,0,0), has got unit radius on the plane XY,
    and in Z goes from -h/2 to +h/2. The ray is defined by its
    direction (u) and a passing point (r). If the ray belongs to the lateral
    circular wall of the cylinder the length returned is h (h is computed as
    1/ZZ(u), this is so because it is supposed that this intersection
    is computed after a coordinate transformation process from any cylinder
    to a unit one).
    See \Ref{Cylinder} to know how you can intersect any ray with any cylinder.*/
double intersection_unit_cylinder(const matrix1D<double> &u, 
   const matrix1D<double> &r);

/** Intersection of a ray with a unit cube.
    The cube is centered at (0,0,0) and has got unit size length in all
    directions, i.e., the cube goes from (-0.5,-0.5,-0.5) to (0.5,0.5,0.5).
    The ray is defined by its direction (u) and a passing
    point (r). See \Ref{Cube} to know how you can intersect any ray with
    any cube.*/
double intersection_unit_cube(const matrix1D<double> &u, 
   const matrix1D<double> &r);
//@}
//@}
#endif
