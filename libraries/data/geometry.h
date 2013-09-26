/***************************************************************************
*
* Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
*  e-mail address 'xmipp@cnb.csic.es'
***************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "multidim_array.h"
#include "multidim_array_generic.h"

#ifndef FLT_EPSILON
#define FLT_EPSILON 1.19209e-07
#endif

#include <vector>
#include <iostream>

/// @defgroup Geometry Geometry
/// @ingroup DataLibrary
//@{
/// @name Geometrical operations
/// @{

/** Project a point to a plane (direction vector, distance)
 *
 * Given the coordinates for a vector in R3, and a plane (defined by its
 * direction vector and the minimum distance from the plane to the coordinate
 * system origin). This function computes the perpendicular projection from the
 * vector to the plane. This function has tried to be optimized in speed as it
 * is used in core routines within huge loops, this is why the result is given
 * as an argument, and why no check about the dimensionality of the vectors is
 * performed. The routine performs faster if the result vector is already in R3.
 *
 * The following example projects the point P=(1,1,1) to the XY-plane storing
 * the result in Pp (belonging to R3), the result is obviously Pp=(1,1,0).
 *
 * @code
 * Matrix1D< double > Z = vectorR3(0, 0, 1), P = vector_R3(1, 1, 1), Pp(3);
 * Uproject_to_plane(P,Z,0,Pp);
 * std::cout << "After projecting: Pp=" << Pp.transpose() << std::endl;
 * @endcode
 *
 * The starting U in the function name stands for the fact that the plane and
 * the point are in the same reference system (called by default Universal).
 *
 * The result and point vectors can be the same one.
 */
void Uproject_to_plane(const Matrix1D< double >& point,
                       const Matrix1D< double >& direction,
                       double distance,
                       Matrix1D< double >& result);

/** Project a vector to a plane (Euler angles)
 *
 * These planes are restricted to have 0 distance to the universal
 * coordinate system. In this special case, a plane can be defined by 3 Euler
 * angles (this is specially suited for projections, where the projection plane
 * is defined by its 3 Euler angles). Then, the coordinate system associated to
 * the 3 Euler angles (let's call its vectors X',Y',Z') defines a projection
 * plane where Z' is the direction vector, and X'Y' are in-plane vectors.
 * Actually, X' coincides with the X vector in the Matrix2D definition and Y'
 * with the Y vector.
 *
 * The resulting vector is in R3, and the function has been optimized for speed,
 * so the result is passed as a parameter. This function is based in the one
 * which projects a point given the Euler matrix of the projection plane. If you
 * are to project several points to the same plane, you should better use the
 * project to plane function where you give the Euler matrix.
 *
 * The following example projects the point P=(1,1,1) to the XY-plane storing
 * the result in Pp (belonging to R3), the result is obviously Pp=(1,1,0).
 *
 * @code
 * Matrix1D< double > P = vectorR3(1, 1, 1), Pp(3);
 * Uproject_to_plane(P, 0, 0, 0, Pp);
 * std::cout << "After projecting: Pp=" << Pp.transpose() << std::endl;
 * @endcode
 *
 * The starting U in the function name stands for the fact that the plane, and
 * the point are in the same reference system (called by default Universal).
 *
 * The result and point vectors can be the same one.
 */
void Uproject_to_plane(const Matrix1D< double >& r,
                       double rot,
                       double tilt,
                       double psi,
                       Matrix1D< double >& result);

/** Project a vector to a plane (Euler matrix)
 *
 * These planes are restricted to have 0 distance to the universal coordinate
 * system. In this special case, a plane can be defined by 3 Euler angles (this
 * is specially suited for projections, where the projection plane is defined by
 * its 3 Euler angles). Then, the coordinate system associated to the 3 Euler
 * angles (let's call its vectors X',Y',Z') defines a projection plane where Z'
 * is the direction vector, and X'Y' are in-plane vectors. Actually, X'
 * coincides with the X vector in the Matrix2D definition and Y' with the Y
 * vector.
 *
 * The resulting vector is in R3, and the function has been optimized for speed,
 * if the result vector is already 3 dimensional when entering the function, no
 * resize is performed; and the result is passed as a parameter. If you are to
 * project several points to the same plane, you should better use the project
 * to plane function where you give the Euler matrix.
 *
 * The following example projects the point P=(1,1,1) to the XY-plane storing
 * the result in Pp (belonging to R3), the result is obviously Pp=(1,1,0).
 *
 * @code
 * Matrix1D< double > P = vectorR3(1, 1, 1), Pp(3);
 * Matrix2D< double > euler;
 * Euler_angles2matrix(0, 0, 0, euler);
 * Uproject_to_plane(P, euler, Pp);
 * std::cout << "After projecting: Pp=" << Pp.transpose() << std::endl;
 * @endcode
 *
 * The starting U in the function name stands for the fact that the plane, and
 * the point are in the same reference system (called by default Universal).
 *
 * The result and point vectors can be the same one.
 */
void Uproject_to_plane(const Matrix1D< double >& r,
                       const Matrix2D< double >& euler,
                       Matrix1D< double >& result);

/** Spherical distance
 *
 * This function returns the distance over a sphere, not the straight line but
 * the line which goes from one point to the other going over the surface of a
 * sphere, supposing that both points lie on the same sphere.
 */
double spherical_distance(const Matrix1D< double >& r1,
                          const Matrix1D< double >& r2);

/** Point to line distance in 3D
 *
 * Let a line in 3-D be specified by the point a and the vector v, this fuction
 * returns the minimum distance of this line to the point p.
 */
double point_line_distance_3D(const Matrix1D< double >& p,
                              const Matrix1D< double >& a,
                              const Matrix1D< double >& v);

/** Point to plane distance in 3D
 *
 * Let a plane in 3-D be specified by the point a and the perpendicular vector
 * v, this fuction returns the minimum distance of this plane to the point p.
 */
inline double point_plane_distance_3D(const Matrix1D< double >& p,
                                      const Matrix1D< double >& a,
                                      const Matrix1D< double >& v)
{
    static Matrix1D< double > p_a(3);

    V3_MINUS_V3(p_a, p, a);
    return (dotProduct(p_a, v) / v.module());
}

/** Structure of the points to do model fitting
 */
struct FitPoint
{
    /// x coordinate
    double x;
    /// y coordinate
    double y;
    /// z coordinate, assumed to be a function of x and y
    double z;
    /// Weight of the point in the Least-Squares problem
    double w;
};

/** Least-squares-fit a plane to an arbitrary number of (x,y,z) points
 *
 * Plane described as Ax + By + C = z
 *
 * Points are defined using the struct
 *
 * @code
 * struct fit_point
 * {
 *     double x;
 *     double y;
 *     double z;
 *     double w;
 * };
 * @endcode
 *
 * where w is a weighting factor. Set it to 1 if you do not want to use it
 */
void least_squares_plane_fit(FitPoint *IN_points,
        					 int Npoints,
        					 double& plane_A,
                             double& plane_B,
                             double& plane_C);

/** Structure of the points to do model fitting
 */
struct fit_point2D
{
    /// x coordinate
    double x;
    /// y coordinate (assumed to be a fucntion of x)
    double y;
    /// Weight of the point in the Least-Squares problem
    double w;
};

/** Least-squares-fit a line to an arbitrary number of (x,y) points
 *
 * Plane described as Ax + B = y
 *
 * Points are defined using the struct
 *
 * @code
 * struct fit_point2D
 * {
 *     double x;
 *     double y;
 *     double w;
 * };
 * @endcode
 *
 * where w is a weighting factor. Set it to 1 if you do not want to use it
 */
void least_squares_line_fit(const std::vector< fit_point2D >& IN_points,
                            double& line_A,
                            double& line_B);

/** Bspline model class
 *
 * When you fit a Bspline model this is the type returned. You can use it to
 * evaluate it anywhere.
 *
 * The model is f(x,y)=sum_{l=l0}^{lF} {sum_{m=m0}^{mF}
 * {c_{ml}Beta_n((x-x0)/h_x-l) Beta_n((y-y0)/h_y-m) } }.
 *
 * The parameter n is the Bspline degree. l0, lF, m0 and mF are the Bspline
 * indexes. hx and hy are related to the extent of the Bspline.
 */
class Bspline_model
{
public:
    /// l0
    int l0;
    /// lF;
    int lF;
    /// m0
    int m0;
    /// mF;
    int mF;

    /// x0
    double x0;
    /// y0
    double y0;

    /// Order of the Bspline
    int SplineDegree;

    /// Scale X
    double h_x;
    /// Scale Y
    double h_y;

    /** Bspline coefficients, c_{ml}
     *
     * The logical indexes of this matrix go from Y=[m0...mF] and X=[l0...lF]
     */
    MultidimArray< double > c_ml;

    /// Evaluate the model at the point (x,y)
    inline double evaluate(double x, double y) const
    {
        int SplineDegree_1 = SplineDegree - 1;
        double x_arg = (x - x0) / h_x;
        double y_arg = (y - y0) / h_y;

        int l1 = CLIP(CEIL(x_arg - SplineDegree_1), l0, lF);
        int l2 = CLIP(l1 + SplineDegree, l0, lF);
        int m1 = CLIP(CEIL(y_arg - SplineDegree_1), m0, mF);
        int m2 = CLIP(m1 + SplineDegree, m0, mF);
        double columns = 0.0;
        for (int m = m1; m <= m2; m++)
        {
            double rows = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x_arg - (double)l;
                double Coeff = c_ml(m, l);
                switch (SplineDegree)
                {
                case 2:
                    rows += Coeff * Bspline02(xminusl);
                    break;
                case 3:
                    rows += Coeff * Bspline03(xminusl);
                    break;
                case 4:
                    rows += Coeff * Bspline04(xminusl);
                    break;
                case 5:
                    rows += Coeff * Bspline05(xminusl);
                    break;
                case 6:
                    rows += Coeff * Bspline06(xminusl);
                    break;
                case 7:
                    rows += Coeff * Bspline07(xminusl);
                    break;
                case 8:
                    rows += Coeff * Bspline08(xminusl);
                    break;
                case 9:
                    rows += Coeff * Bspline09(xminusl);
                    break;
                }
            }

            double yminusm = y_arg - (double)m;
            switch (SplineDegree)
            {
            case 2:
                columns += rows * Bspline02(yminusm);
                break;
            case 3:
                columns += rows * Bspline03(yminusm);
                break;
            case 4:
                columns += rows * Bspline04(yminusm);
                break;
            case 5:
                columns += rows * Bspline05(yminusm);
                break;
            case 6:
                columns += rows * Bspline06(yminusm);
                break;
            case 7:
                columns += rows * Bspline07(yminusm);
                break;
            case 8:
                columns += rows * Bspline08(yminusm);
                break;
            case 9:
                columns += rows * Bspline09(yminusm);
                break;
            }
        }
        return columns;
    }
};

/** Least-squares fit of a B-spline 2D model
 *
 * For fitting a set of values that are distributed between (x0,y0) and (xF,yF)
 * with a cubic Bspline centered on each corner, the right call is
 *
 * @code
 * Bspline_model model;
 * Bspline_model_fitting(list_of_points, 3, -1, 2, -1, 2, xF-x0, yF-y0, x0,
 *     y0, model);
 * @endcode
 *
 * Once the model is returned you can evaluate it at any point simply by
 *
 * @code
 * model.evaluate(x,y);
 * @endcode
 */
void Bspline_model_fitting(const std::vector< FitPoint >& IN_points,
                           int SplineDegree,
                           int l0,
                           int lF,
                           int m0,
                           int mF,
                           double h_x,
                           double h_y,
                           double x0,
                           double y0,
                           Bspline_model& result);

/** Rectangle which encloses a deformed rectangle
 *
 * Given a rectangle characterized by the top-left corner and the right-bottom
 * corner, and given a matrix after which the rectangle is deformed. Which is
 * the minimum rectangle which encloses the preceeding one? This function is
 * useful for stablishing for loops which will cover for sure the deformed
 * rectangle. All vectors are supposed to be 2x1 and the deformation matrix is
 * 2x2. The corner (x0,y0) goes to V*(x0,y0)' and (xF,yF) to V*(xf,yF)'. After
 * that you can make a loop from corner1 to corner2.
 *
 * The v0 and vF vectors can be reused as outputs.
 */
void rectangle_enclosing(const Matrix1D< double >& v0,
                         const Matrix1D< double >& vF,
                         const Matrix2D< double >& V,
                         Matrix1D< double >& corner1,
                         Matrix1D< double >& corner2);

/** Box which encloses a deformed box
 *
 * Given a box characterized by the top-left corner (most negative) and the
 * right-bottom (most positive) corner, and given a matrix after which the box
 * is deformed. Which is the minimum box which encloses the preceeding one? This
 * function is useful for stablishing for loops which will cover for sure the
 * deformed box. All vectors are supposed to be 3x1 and the deformation matrix
 * is 3x3. The corner (x0,y0,z0) goes to V*(x0,y0,z0)' and (xF,yF,zF) to
 * V*(xf,yF,zF)'. After that you can make a loop from corner1 to corner2.
 *
 * The v0 and vF vectors can be reused as outputs.
 */
void box_enclosing(const Matrix1D< double >& v0,
                   const Matrix1D< double >& vF,
                   const Matrix2D< double >& V,
                   Matrix1D< double >& corner1,
                   Matrix1D< double >& corner2);

/** Point inside polygon
 *
 * Given a polygon described by a list of points (the last one and the first one
 * ust be the same), determine whether another point is inside the polygon or
 * not.
 */
bool point_inside_polygon(const std::vector< Matrix1D< double > > & polygon,
                          const Matrix1D< double >& point);

/** Line Plane Intersection
 *
 * Let ax+by+cz+D=0 be the equation of your plane (if your plane is defined by a
 * normal vector N + one point M, then (a,b,c) are the coordinates of the normal
 * N, and d is calculated by using the coordinates of M in the above equation).
 *
 * Let your line be defined by one point P(d,e,f) and a vector V(u,v,w), the
 * points on your line are those which verify
 *
 * x = d + lu
 * y = e + lv
 * z = f + lw
 *
 * where l takes all real values.
 *
 * for this point to be on the plane, you have to have
 *
 * ax + by + cz + D = 0, so,
 *
 * a(d + lu) + b(e + lv) + c(f + lw) + D = 0
 *
 * that is
 *
 * l(au + bv + cw) = -(ad + be + cf + D)
 *
 * note that, if au + bv + cw = 0, then your line is either in the plane, or
 * parallel to it... otherwise you get the value of l, and the intersection has
 * coordinates:
 *
 * x = d + lu
 * y = e + lv
 * z = f + lw
 *
 * where
 *
 * l = -(ad + be + cf + D) / (au + bv + cw)
 *
 * a = XX(normal_plane);
 * b = YY(normal_plane);
 * c = ZZ(normal_plane);
 * D = intersection_point;
 *
 * d = XX(point_line)
 * e = YY(point_line)
 * f = ZZ(point_line)
 *
 * u = XX(vector_line);
 * v = YY(vector_line);
 * w = ZZ(vector_line);
 *
 * XX(intersection_point) = x;
 * YY(intersection_point) = y;
 * ZZ(intersection_point) = z;
 *
 * return 0 if sucessful
 * return -1 if line paralell to plane
 * return +1 if line in the plane
 */
int line_plane_intersection(const Matrix1D< double > normal_plane,
                            const Matrix1D< double > vector_line,
                            Matrix1D< double >& intersection_point,
                            const Matrix1D< double > point_line,
                            double point_plane_at_x_y_zero = 0.);
//@}

/// @name Euler operations
/// @{

/** Getting the Euler angles to a range (0-360).
 * No direction equivalence is applied, ie, there is no correction of the
 * direction of projection making use that a view from the top is the same as a
 * view from the bottom but reversed ... Just a wrapping of the angles is done
 * until the angles fall in the specified ranges. The angles given must be
 * variables and they are modified with the new values.
 *
 */
#define EULER_CLIPPING(rot,tilt,psi) \
    rot = realWRAP(rot, 0, 360); \
    tilt = realWRAP(tilt, 0, 360); \
    psi = realWRAP(psi, 0, 360);

/** Getting the Euler angles to a range (0-2*PI).
 *
 * The same as before but the angles are expressed in radians.
 */
#define EULER_CLIPPING_RAD(rot, tilt, psi) \
    rot = realWRAP(rot, 0, 2.0*PI); \
    tilt = realWRAP(tilt, 0, 2.0*PI); \
    psi = realWRAP(psi, 0, 2.0*PI);

/** Euler angles --> Euler matrix.
 *
 * This function returns the transformation matrix associated to the 3 given
 * Euler angles (in degrees).
 *
 * As an implementation note you might like to know that this function calls
 * always to Matrix2D::resize
 *
 * See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/EulerAngles for a
 * description of the Euler angles.
 */
void Euler_angles2matrix(double a, double b, double g, Matrix2D< double >& A,
                         bool homogeneous=false);

/** Euler angles --> Euler matrix.
 *
 * This function returns the transformation matrix associated to the 3 given
 * Euler angles (in degrees).
 */
void Euler_anglesZXZ2matrix(double a, double b, double g, Matrix2D< double >& A,
                         bool homogeneous=false);

/** Distance between two Euler matrices.
 *
 * The distance is defined as 1/3*(X1.X2 + Y1.Y2 + Z1.Z2)
 */
double Euler_distanceBetweenMatrices(const Matrix2D<double> &E1,
                                     const Matrix2D<double> &E2);

/** Average distance between two angle sets.
 * If the only_projdir is set, then only the projection direction is considered.
 */
double Euler_distanceBetweenAngleSets(double rot1, double tilt1, double psi1,
                                      double rot2, double tilt2, double psi2,
                                      bool only_projdir);

/** Average distance between two angle sets.
 * E1 must contain the Euler matrix corresponding to set1, E2 is used as
 * an auxiliary variable for storing the second Euler matrix.
 */
double Euler_distanceBetweenAngleSets_fast(const Matrix2D<double> &E1,
        								   double rot2, double tilt2, double psi2,
        								   bool only_projdir, Matrix2D<double> &E2);

/** Angles after compresion
 *
 * Let be two volumes f and g related by g(x,y,z) = f(D(x,y,z)) (where D is a
 * lineal transformation) then the projection direction parallel to the vector w
 * in f is going to be related with the projection direction parallel to the
 * vector w_prime in g. Given the w Euler angles this routine provide the
 * w_prime angles
 */
void Euler_Angles_after_compresion(const double rot,
                                   double tilt,
                                   double psi,
                                   double& new_rot,
                                   double& new_tilt,
                                   double& new_psi,
                                   Matrix2D< double >& D);

/** Euler direction
 *
 * This function returns  a vector parallel to the  projection direction.
 * Resizes v if needed
 */
void Euler_direction(double alpha,
                     double beta,
                     double gamma,
                     Matrix1D< double >& v);

/** Euler direction2angles
 *
 * This function returns the 3 Euler angles associated to the direction given by
 * the vector v. The 3rd Euler angle is set always to 0
 */
void Euler_direction2angles(Matrix1D< double >& v,
                            double& alpha,
                            double& beta,
                            double& gamma);

/** "Euler" matrix --> angles
 *
 * This function compute a set of Euler angles which result in an "Euler" matrix
 * as the one given. See \ref Euler_angles2matrix to know more about how this
 * matrix is computed and what each row means. The result angles are in degrees.
 * Alpha, beta and gamma are respectively the first, second and third rotation
 * angles. If the input matrix is not 3x3 then an exception is thrown, the
 * function doesn't check that the Euler matrix is truly representing a
 * coordinate system.
 *
 * @code
 * Euler_matrix2angles(Euler, alpha, beta, gamma);
 * @endcode
 */
void Euler_matrix2angles(const Matrix2D< double >& A,
                         double& alpha,
                         double& beta,
                         double& gamma);

/** Up-Down projection equivalence
 *
 * As you know a projection view from a point has got its homologous from its
 * diametrized point in the projection sphere. This function takes a projection
 * defined by its 3 Euler angles and computes an equivalent set of Euler angles
 * from which the view is exactly the same but in the other part of the sphere
 * (if the projection is taken from the bottom then the new projection from the
 * top, and viceversa). The defined projections are exactly the same except for
 * a flip over X axis, ie, an up-down inversion. Exactly the correction
 * performed is:
 *
 * @code
 * newrot = rot;
 * newtilt = tilt + 180;
 * newpsi = -(180 + psi);
 * @endcode
 *
 * @code
 * Euler_up_down(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_up_down(double rot,
                   double tilt,
                   double psi,
                   double& newrot,
                   double& newtilt,
                   double& newpsi);

/** The same view but differently expressed
 *
 * As you know a projection view from a point can be expressed with different
 * sets of Euler angles. This function gives you another expression of the Euler
 * angles for this point of view. Exactly the operation performed is:
 *
 * @code
 * newrot = rot + 180;
 * newtilt = -tilt;
 * newpsi = -180 + psi;
 * @endcode
 *
 * @code
 * Euler_another_set(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_another_set(double rot,
                       double tilt,
                       double psi,
                       double& newrot,
                       double& newtilt,
                       double& newpsi);

/** Mirror over Y axis
 *
 * Given a set of Euler angles this function returns a new set which define a
 * mirrored (over Y axis) version of the former projection.
 *
 * @code
 *  -----> X               X<------
 *  |                              |
 *  |                              |
 *  |               ======>        |
 *  v                              v
 *  Y                             Y
 * @endcode
 *
 * The operation performed is
 *
 * @code
 * newrot = rot;
 * newtilt = tilt + 180;
 * newpsi = -psi;
 * @endcode
 *
 * @code
 * Euler_mirrorY(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_mirrorY(double rot,
                   double tilt,
                   double psi,
                   double& newrot,
                   double& newtilt,
                   double& newpsi);

/** Mirror over X axis
 *
 * Given a set of Euler angles this function returns a new set which define a
 * mirrored (over X axis) version of the former projection.
 *
 * @code
 *  -----> X               Y
 *  |                       ^
 *  |                       |
 *  |               ======> |
 *  v                       |
 *  Y                        -----> X
 * @endcode
 *
 * The operation performed is
 *
 * @code
 * newrot = rot;
 * newtilt = tilt + 180;
 * newpsi = 180 - psi;
 * @endcode
 *
 * @code
 * Euler_mirrorX(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_mirrorX(double rot,
                   double tilt,
                   double psi,
                   double& newrot,
                   double& newtilt,
                   double& newpsi);

/** Mirror over X and Y axes
 *
 * Given a set of Euler angles this function returns a new set which define a
 * mirrored (over X and Y axes at the same time) version of the former
 * projection.
 *
 * @code
 *  -----> X                       Y
 *  |                               ^
 *  |                               |
 *  |               ======>         |
 *  v                               |
 *  Y                        X<-----
 * @endcode
 *
 * The operation performed is
 *
 * @code
 * newrot = rot;
 * newtilt = tilt;
 * newpsi = 180 + psi;
 * @endcode
 *
 * @code
 * Euler_mirrorX(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_mirrorXY(double rot,
                    double tilt,
                    double psi,
                    double& newrot,
                    double& newtilt,
                    double& newpsi);

/** Apply a geometrical transformation
 *
 * The idea behind this function is the following. 3 Euler angles define a point
 * of view for a projection, but also a coordinate system. You might apply a
 * geometrical transformation to this system, and then compute back what the
 * Euler angles for the new system are. This could be used to "mirror" points of
 * view, rotate them and all the stuff. The transformation matrix must be 3x3
 * but it must transform R3 vectors into R3 vectors (that is a normal 3D
 * transformation matrix when vector coordinates are not homogeneous) and it
 * will be applied in the sense:
 *
 * @code
 * New Euler matrix = L * Old Euler matrix * R
 * @endcode
 *
 * where you know that the Euler matrix rows represent the different system
 * axes. See Euler_angles2matrix for more information about the Euler coordinate
 * system.
 *
 * @code
 * Matrix2D< double > R60 = rotation3DMatrix(60, 'Z');
 * R60.resize(3, 3); // Get rid of homogeneous part
 * Matrix2D< double > I(3, 3);
 * I.initIdentity();
 * Euler_apply_transf(I, R60, rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_apply_transf(const Matrix2D< double >& L,
                        const Matrix2D< double >& R,
                        double rot,
                        double tilt,
                        double psi,
                        double& newrot,
                        double& newtilt,
                        double& newpsi);

/** Rotate a volume after 3 Euler angles
 *
 * Input and output volumes cannot be the same one.
 */
void Euler_rotate(const MultidimArray< double >& V,
                  double rot,
                  double tilt,
                  double psi,
                  MultidimArray< double >& result);

/** Rotate a volume after 3 Euler angles
 *
 * Input and output volumes cannot be the same one.
 */
void Euler_rotate(const MultidimArrayGeneric &V,
                  double rot,
                  double tilt,
                  double psi,
                  MultidimArray<double> &result);
/** Compute circle around Euler matrix
 *
 * Given an input Euler matrix, this function returns a set of Euler
 * angles such that they sample a circle around the original projection
 * direction (a sample every angStep). The projection directions in the
 * circle are separated by angCircle.
 *
 * The output is in outputEulerAngles whose structure is
 * (newrot1,newtilt1,newpsi1,newrot2,newtilt2,newpsi2,...)
 */
void computeCircleAroundE(const Matrix2D<double> &E,
                          double angCircle, double angStep, std::vector<double> &outputEulerAngles);
//@}

/// @name Intersections
/// @{

/** Intersection of a ray with a unit sphere
 *
 * The sphere is centered at (0,0,0) and has got unit radius. The ray is defined
 * by its direction (u) and a passing point (r). The function returns the length
 * of the intersection. If the ray is tangent to the sphere the length is 0. See
 * Ellipsoid to know how you can intersect any ray with any ellipsoid/sphere.
 */
double intersection_unit_sphere(const Matrix1D< double >& u,
                                const Matrix1D< double >& r);

/** Intersection of a ray with a unit cylinder
 *
 * The cylinder is centered at (0,0,0), has got unit radius on the plane XY, and
 * in Z goes from -h/2 to +h/2. The ray is defined by its direction (u) and a
 * passing point (r). If the ray belongs to the lateral circular wall of the
 * cylinder the length returned is h (h is computed as 1/ZZ(u), this is so
 * because it is supposed that this intersection is computed after a coordinate
 * transformation process from any cylinder to a unit one).
 *
 * See Cylinder to know how you can intersect any ray with any cylinder.
 */
double intersection_unit_cylinder(const Matrix1D< double >& u,
                                  const Matrix1D< double >& r);

/** Intersection of a ray with a unit cube
 *
 * The cube is centered at (0,0,0) and has got unit size length in all
 * directions, i.e., the cube goes from (-0.5, -0.5, -0.5) to (0.5, 0.5, 0.5).
 * The ray is defined by its direction (u) and a passing point (r). See Cube to
 * know how you can intersect any ray with any cube.
 */
double intersection_unit_cube(const Matrix1D< double >& u,
                              const Matrix1D< double >& r);
//@}
//@}
#endif
