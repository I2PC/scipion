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
/* ------------------------------------------------------------------------- */
/* PHANTOMS                                                                  */
/* ------------------------------------------------------------------------- */

#ifndef _PHANTOM_HH
#define _PHANTOM_HH

#include <vector>

#include "multidim_array.h"
#include "blobs.h"
#include "projection.h"

/**@defgroup Phantoms Phantoms
 * @ingroup DataLibrary
 *  Phantoms are mathematical description of volumes such that in the
 *  reconstruction process, we can know exactly which was the original
 *  volume in a mathematical way. In this package phantoms are considered
 *  to be a collection of features plus some more information about the
 *  background and the phantom size. A feature is a cone, a box, a cylinder,
 *  or any other geometrical figure described by its parameters and not
 *  as a volume with voxels at a given value.
 *
 *  The file format generated and accepted by this library is the following:
 *  @code
 *  "#Phantom Xdim Ydim Zdim Background density [scale factor]\n"
 *  "          64   64   64            0             1\n"
 *  "#Type +/= Density X_Center Y_Center Z_Center\n"
 *  " sph   +     1      <x0>     <y0>     <z0>    <radius>\n"
 *  " blo   +     1      <x0>     <y0>     <z0>    <radius>   <alpha> <order>\n"
 *  " gau   +     1      <x0>     <y0>     <z0>    <sigma>\n"
 *  " cyl   +     1      <x0>     <y0>     <z0>    <xradius> <yradius> <height>               <rot> <tilt> <psi>\n"
 *  " dcy   +     1      <x0>     <y0>     <z0>    <radius>            <height>  <separation> <rot> <tilt> <psi>\n"
 *  " cub   =     1      <x0>     <y0>     <z0>    <xdim>     <ydim>    <zdim>                <rot> <tilt> <psi>\n"
 *  " ell   =     1      <x0>     <y0>     <z0>    <xradius> <yradius> <zradius>              <rot> <tilt> <psi>\n"
 *  " con   +     1      <x0>     <y0>     <z0>    <radius>            <height>               <rot> <tilt> <psi>\n"
 *  @endcode
 *  where spheres, blobs, gaussians, cylinders, double cylinders, cubes, ellipsoids and cones
 *  are defined (in this order). The '+' sign means that this feature will
 *  be added at those positions of the volume ocuupied by it. '=' instead
 *  means that those voxels will be set to the value of the density of this
 *  feature (if two features overlap the density of the last one is kept
 *  in the overlapping voxels). The density is the grey level of voxels
 *  affected by that feature. In the preceeding example the final volume
 *  is 64x64x64 and has got a background density of 0. The center of the
 *  features might be negative, and they represent mathematical positions
 *  in R3. The phantom dimension, instead, define the phantom in this case
 *  to go from -32 to 31, in this R3 space.
 *
 *  If the scale factor, which by default is 1, is not unity then the whole
 *  phantom is scaled (0.5 means to its half and 1.5 enlarged by one half)
 *  just after reading it.
*/
//@{
/* FEATURE ================================================================= */
/** Feature superclass.
    This is a superclass from which all features (cones, cylinders, ...)
    inherit. It contains all general information common to all feature
    classes, then the subclasses give the specific parameters for the
    feature.
*/
class Feature
{

    /* Structure --------------------------------------------------------------- */
public:
    /** Feature type.
        A three letters string telling what kind of feature this object is.
        For example, "cyl", "con", "cub", ... See the specific classes to
        know exactly which is each label. */
    std::string       Type;

    /** Feature behaviour.
        This flag indicates how the feature behaves inside the voxel volume.
        If this flag is set to '+' then the voxels occupied by this feature
        are incremented with the feature value. If the flag is set to '='
        then the voxels are set to to the same value of the feature. */
    char              Add_Assign;

    /** Feature Density.
        Density of the feature in grey level values. It needn't be between
        0 and 1, or 0 and 255. What is more, it can be even negative, and so
        you can build holes inside volumes. */
    double             Density;

    /** Center of the feature.
        The center of the feature is understood differently according to the
        specific class, see them to know exactly how this value is interpreted. */
    Matrix1D<double>   Center;

    /** Maximum distance from the center.
        This value is a precalculated and tells the maximum distance from any
        point belonging to the feature to its center. This is used to speed up
        some functions not considering voxels which we know are beyond the
        scope of this feature. */
    double             max_distance;

public:
    virtual ~Feature()
    {}

    /** Prepare feature for work.
        This function computes the maximum distance and possibly the Euler
        and inverse Euler matrices. */
    virtual void prepare() = 0;

    /// Assignment
    Feature & operator = (const Feature &F);

    /** Another function for assigmnet.*/
    void assign(const Feature &F);

    /** Rotate the whole feature.
        Rotate this feature using a rotation matrix. The center as well as
        the feature itself is rotated. */
    virtual void rotate(const Matrix2D<double> &E);

    /** Rotate only the center.
        Rotate the center of this feature only around the phantom center*/
    virtual void rotate_center(const Matrix2D<double> &E);

    /* Inside ------------------------------------------------------------------ */
    /** Speeded up point inside a feature, VIRTUAL!!.
        This function MUST be implemented for each subclass and tells you if
        a point is inside the feature or not. If the point is inside returns 1
        and if not, returns 0. The point is given in the vector r, and aux
        is an auxiliar vector with dimension 3 (must be externally resized)
        used to do some computations. This auxiliar vector must be supplied in
        order to gain speed as no memory allocating and freeing is needed. */
    virtual int point_inside(const Matrix1D<double> &r,
                             Matrix1D<double> &aux) const = 0;

    /** Point inside a feature.
        This function is based in the previous one. It makes the same but
        you needn't supply the auxiliar vector. */
    int point_inside(const Matrix1D<double> &r) const
    {
        Matrix1D<double> aux(3);
        return point_inside(r, aux);
    }

    /** Speeded up density inside a feature, VIRTUAL!!.
        This function MUST be implemented for each subclass with NON constant
        density as blobs and tells you if
        a point is inside the feature or not plus give you information about
         the density of the feature at that point. If the point is inside returns 1
        multiplied by the density of the NORMALIZED feature and if not, returns 0.
        The point is given in the vector r, and aux
        is an auxiliar vector with dimension 3 (must be externally resized)
        used to do some computations. This auxiliar vector must be supplied in
        order to gain speed as no memory allocating and freeing is needed. */
    virtual double density_inside(const Matrix1D<double> &r,
                                  Matrix1D<double> &aux) const = 0;


    /** Speeded up voxel inside a feature.
        A voxel is compound of 8 subvoxels. This function returns the number
        of subvoxel centers falling inside the feature. The voxel size is
        supposed to be 1, and r is the center of the voxel in R3.
        This speeded up function needs
        two vectors with dimension 3 externally resized. */
    int voxel_inside(const Matrix1D<double> &r, Matrix1D<double> &aux1,
                     Matrix1D<double> &aux2) const;

    /** Voxel inside a feature.
        This function is based in the previous one. It makes the same but
        you needn't supply the auxiliar vectors. */
    int voxel_inside(const Matrix1D<double> &r) const
    {
        Matrix1D<double> aux1(3), aux2(3);
        return voxel_inside(r, aux1, aux2);
    }

    /** A voxel is compound of 8 subvoxels. This function returns the number
        of subvoxel centers falling inside the feature multiplied by the normalized
        feature density. That is, the value of the density is always 1 at the origin
        The voxel size is supposed to be 1, and r is the center of the voxel
        in R3. This speeded up function needs
        two vectors with dimension 3 externally resized. */

    double voxel_inside_by_normalized_density(const Matrix1D<double> &r, Matrix1D<double> &aux1,
            Matrix1D<double> &aux2) const;


    /** Speeded up sphere intersecting feature.
        This function returns TRUE if a sphere of a given radius has any voxel
        inside this feature. r is the center of the sphere in R3.
        This speeded up function needs
        two vectors with dimension 3 externally resized. */
    int intersects_sphere(const Matrix1D<double> &r, double radius,
                          Matrix1D<double> &aux1, Matrix1D<double> &aux2, Matrix1D<double> &aux3)
    const;

    /** Sphere intersecting feature.
        This function is based in the previous one. It makes the same but
        you needn't supply the auxiliar vectors. */
    int intersects_sphere(const Matrix1D<double> &r, double radius) const
    {
        Matrix1D<double> aux1(3), aux2(3), aux3(3);
        return intersects_sphere(r, radius, aux1, aux2, aux3);
    }

    /** Produce a sphere envolving the feature.
        This function returns a pointer to a feature (a sphere) with
        the same center, density, and feature behaviour as the given feature.
        The radius could be given in voxel units or if you leave it 0 then
        the radius is computed as 1.5 times the \ref max_distance of this
        feature */
    Feature *encircle(double radius = 0) const;

    /** Return a scaled version of this feature, VIRTUAL!!!.
        This function returns a pointer to a feature (of the same type as the
        feature for which the function was called, ie, if you scale a sphere
        the result is a sphere, if you scale a cone, the result is a cone, ...)
        that is a scaled version of the actual feature (see each specific
        implementation to see which is the relatioship between the two features).
        This function is useful to define backgrounds with the same shape
        of the given feature. */
    virtual Feature *scale(double factor) const = 0;

#define ENLARGE_MODE 1
#define SPHERE_MODE  2
    /** Return a pointer to a feature which is the background of the actual one.
        You can specify two ways of computing the background, either as a sphere
        surrounding the feature or as an enlarged version of the feature. This
        is chosen giving the modes ENLARGE_MODE or SPHERE_MODE. Depending on
        the mode used the background parameter is understood as the sphere
        radius or as the scaling factor. */
    Feature *background(int back_mode, double back_param) const;

    /** Speeded Up intersection of a feature with a ray, VIRTUAL!!!.
        This function returns the length of the intersection between the ray
        defined by its direction and a passing point and the actual feature
        (whose center and dimensions are known by itself).

        r and u are auxiliar vectors of dimension 3, they must be supplied in
        order to gain speed. r is the passing point expressed in the
        feature coordinate system, and u is the direction in the same
        coordinate system. */
    virtual double intersection(const Matrix1D<double> &direction,
                                const Matrix1D<double> &passing_point, Matrix1D<double> &r,
                                Matrix1D<double> &u) const = 0;

    /** Intersection of a feature with a ray.
        This function does the same as the previous one but you needn't provide
        extra auxiliar vectors. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point) const
    {
        Matrix1D<double> r(3), u(3);
        return intersection(direction, passing_point, r, u);
    }

    /** Volume of the feature in (voxel units)³ , VIRTUAL!!!.
        This function returns the volume of each feature supposing that a voxel
        is of size 1x1x1. */
    virtual double volume() const = 0;

    /** Mean and variance in a given plane.
        Given a plane z=z0, this function returns the mean and variance values
        of the volume (*V) in the voxels inside this feature. A voxel is
        considered to belong to the feature if it is totally inside
        the feature. If the plane is outside the volume
        scope the result is mean=variance=0 */
    void mean_variance_in_plane(Volume *V, double z, double &mean, double &var);

    /** Project feature onto a projection plane.
        Projection is a class itself which has got inside the direction of
        projection. The projection plane is not cleaned (set all values to 0)
        when entering this function, but the projection of this feature is added
        to the already stored one. The projection plane is supposed to have its
        logical origin in the right place (normally, the middle of the image) before
        entering the function. An oversampling technique is used in order to
        produce better projections, the projection in one pixel is computed
        from the value of 4 points near the pixel, you can modify the subsampling
        value (actually 2x2=4) by changing a constant at the beginning of the
        function.

        The matrix VP (3x3) is used to define how to relate a point in the 3D
        universal coordinate to the projection plane. So a point ru in the
        universal coordinate system, projects to VP*r. PV must be the inverse
        of the forwarding projection matrix. */
    void project_to(Projection &P, const Matrix2D<double> &VP,
                    const Matrix2D<double> &PV) const;

    /** Define 3D corners for a feature.
        This function returns two Z3 points where the feature is confined.
        The volume borders are taken into account and you might make a for
        using these two values like this:
        \\Ex:
        @code
        F.corners(V,corner1,corner2);
        for (int k=ZZ(corner1); k<=ZZ(corner2); k++)
            for (int i=YY(corner1); i<=YY(corner2); i++)
                for (int j=XX(corner1); j<=XX(corner2); j++) {
                    ...
        }
        @endcode*/
    void corners(const Volume *V, Matrix1D<double> &corner1,
                 Matrix1D<double> &corner2);

#define INTERNAL 0
#define EXTERNAL 1
    /** Draw a feature in a Volume.
        This function draws a feature in a volume (*V). The word "colour" is
        a little misleading as the colour really refers to a grey-level global
        density. Ie, the default mode (INTERNAL) draws the feature with its own
        density and discards the "colour" given in the function call. However,
        you may change this colour, set a new global density and draw the
        volume with this given density. This option is useful to generate easily
        a labelled volume.

        The Add-Assign behaviour is the one of the feature if the colour is
        internally defined, or Assign if the colour is externally given.
        In the assign behaviour a voxel value is assigned to the volume if the
        voxel there is smaller than the value we pretend to assign.

        A voxel is drawn with a colour proportional to the number of vertices
        inside the features.
        The volume is not cleaned at the beginning.
        \\ Ex:
        @code
        f.draw_in(&V)              --> Internal density
        f.draw_in(&V,INTERNAL,0.5) --> Internal density, 0.5 is discarded
        f.draw_in(&V,EXTERNAL,0.5) --> External density=0.5
        @endcode
    */
    void draw_in(Volume *V, int color_mode = INTERNAL, double colour = -1);

    /** Draw the surface of the feature.
        This function draws the surface of the feature at the given volume.
        A voxel is said to belong to the surface if the number of corners
        inside the feature meets 1<=n<=7. The default gray_level with which
        voxels will be drawn is 2 (normally higher than the usual volume grey
        levels, and the behaviour of the voxels drawn is Assign. */
    void sketch_in(Volume *V, double colour = 2);

    /** Shift.
        Shift the feature a given amount of voxels. The new center is the old
        center plus the given shift, and that's all*/
    void shift(double shiftX, double shiftY, double shiftZ);

    /** Apply a general geometric transformation.
        The transformation must be
        defined by a 4x4 matrix that can be generated using the geometric functions
        in xmippGeometry or xmippMatrix2D. The matrix must be the desired
        transformation (i.e., new coordinate=A*old_coordinate.
        Don't worry because the selfApplyGeometry of Phantom take
        cares of passing to this function the apropriate matrix. No check is done
        about the size of A.

        Only the center is transformed, the feature will keep the same size.*/
    void selfApplyGeometry(const Matrix2D<double> &A);

    /** Print the feature in the Feature format, VIRTUAL!!!.
        This function prints the feature in the Standard Feature format (readable
        by this library). Notice that the standard format is different for each
        specific feature. See \ref Phantoms for more information about the
        supported file format. */
    virtual void feat_printf(FILE *fh) const = 0;

    /** Read common part of the feature description.
        The common part is the feature type, the behaviour, density and center.
        The description is passed as a line. Exceptions are thrown if the
        description doesn't conform the standard specification.*/
    void read_common(char *line);

    /** Read a feature from a file, VIRTUAL!!!.
        The format must be the one given in \ref Phantoms, and each subclass
        must implement its own I/O routines. These routines must fill only the
        non common part of the feature description, but they receive the whole
        line with the description. */
    virtual void read_specific(char *line) = 0;

    /** Show feature not in the standard format but more informatively.
        This function is based on the std::cout << ... of each subclass. First
        shows the common part of the feature and then its specific part. Be
        careful that you must show a pointer to the feature!!
        \\ Ex: Sphere sphere; std::cout << (Feature *) \&sphere; */
    friend std::ostream& operator << (std::ostream &o, const Feature *F);
};

/* ORIENTED FEATURE ======================================================== */
/** Oriented Features.
    The oriented features are defined in a canonical position (usually
    along Z axis) and then they are rotated after the 3 Euler angles.
    The corresponding Euler matrix is stored in the field "euler" while
    its inverse is in "eulert".
*/
class Oriented_Feature: public Feature
{
public:
    /// First Euler angle
    double             rot;

    /// Second Euler angle
    double             tilt;

    /// Third Euler angle
    double             psi;

    /// Euler matrix
    Matrix2D<double>   euler;

    /// Inverse Euler matrix
    Matrix2D<double>   eulert;
public:
    virtual ~Oriented_Feature()
    {}

    /** Compute Euler and inverse Euler matrices from the Euler angles. */
    void prepare_Euler()
    {
        Euler_angles2matrix(rot, tilt, psi, euler);
        eulert = euler.transpose();
    }

/// Assigment
    Oriented_Feature & operator = (const Oriented_Feature & OF);

    /** Another function for assigment.*/
    void assign(const Oriented_Feature & OF);

    /** Rotate.
        Rotate this feature. The center as well as the feature itself are
        rotated. */
    virtual void rotate(const Matrix2D<double> &E);
};

/* SPHERES ================================================================= */
/** Spheres.
    The spheres are defined by its center and a radius. The label is "sph".

    A point (r) in the universal coordinate system, where spheres are defined
    in general, can be expressed (rp) with respect to a system where the sphere
    is centered at the origin and its radius is unity by
    @code
    V3_MINUS_V3(rp,r,sph.Center);
    V3_BY_CT(rp, rp, 1/radius);
    @endcode
    Directions (free vectors) are transformed in the same fashio except
    that for the first vector substraction (referring to the origin).
*/
class Sphere: public Feature
{
public:
    /// Sphere radius
    double             radius;
public:
    /** Prepare sphere for work.
        Computes the maximum distance, in this case, equal to "radius" */
    void prepare();

/// Assignment
    Sphere & operator = (const Sphere &F);

    /** Another function for assignment.*/
    void assign(const Sphere &F);

    /** Speeded up point inside a sphere.
        This function tells you if a point is inside the sphere or not.
        See \ref Feature::point_inside */
    int point_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Density inside an Sphere.
        This function tells you the density of the sphere at point r
        for constant valued features is trivial*/
    double density_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const
    {
        return (1.);
    }

    /** Return a scaled sphere.
        The center, density, and behaviour of the new sphere is exactly the
        same as the actual one. The radius is multiplied by the scale factor
        and the maximum distance is recalculated for the new sphere. */
    Feature *scale(double factor) const;

    /** Another function for return a scaled sphere.*/
    void scale(double factor, Feature *_f) const;

    /** Intersection of a ray with a sphere.
        See \ref Feature::intersection to know more about the parameters
        meaning. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point,
                        Matrix1D<double> &r, Matrix1D<double> &u) const;

    /** Volume of a sphere.
        This function returns 4/3*PI*radius*radius*radius.
        See \ref Feature::volume */
    double volume() const
    {
        return 4 / 3*PI*radius*radius*radius;
    }

    /** Read specific description for a sphere.
        An exception is thrown if the line doesn't conform the standard
        specification. See \ref Feature::read_specific */
    void read_specific(char *line);

    /** Print sphere in the standard feature format.
        \ See {Feature::feat_printf}, \ref Phantoms */
    void feat_printf(FILE *fh) const;

    /** Show feature not in the standard format but more informatively.
        This function only shows the non common part of the feature. Use
        the << operator of Feature to show the whole feature. */
    friend std::ostream& operator << (std::ostream &o, const Sphere &f);

    /** Generate a random sphere.
        A sphere is generated randomly within the range of the parameters
        given. Notice that most of them are set by default. The exact
        parameters are picked uniformly from the range. The maximum distance
        is adjusted properly according to the randomly generated feature. 'den'
        stands for density and (x0,y0,z0) is the center of the feature. The
        behaviour of the new sphere is always Add */
    void init_rnd(
        double minradius,   double maxradius,
        double minden = 1.0,  double maxden = 1.0,
        double minx0 = 0,     double maxx0 = 0,
        double miny0 = 0,     double maxy0 = 0,
        double minz0 = 0,     double maxz0 = 0);
};

/* BLOB ================================================================= */
/** Blobs.
    The blobs are defined by its center, radius, alpha and m (bessel function
    order). The label is "blo".

    A point (r) in the universal coordinate system, where blobs are defined
    in general, can be expressed (rp) with respect to a system where the blob
    is centered at the origin and its radius is unity by
    @code
    V3_MINUS_V3(rp,r,sph.Center);
    V3_BY_CT(rp, rp, 1/radius);
    @endcode
    Directions (free vectors) are transformed in the same fashion except
    that for the first vector substraction (referring to the origin).
*/

class Blob: public Feature
{
public:

    /* I do not use the structure "blobtype" because common code with sphere becames
    more difficult*/
    /// Blob radius
    double             radius;
    /// Blob alpha
    double             alpha;
    ///
    int         m;
public:
    /** Prepare blob for work.
        Computes the maximum distance, in this case, equal to "radius" */
    void prepare();

    /// Assignment
    Blob & operator = (const Blob &F);

    /** Another function for assignment.*/
    void assign(const Blob &F);

    /** Speeded up point inside a blob.
        This function tells you if a point is inside the blob or not.
        See \ref Feature::point_inside */
    int point_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Density inside a blob.
        This function tells you the density of the blob at point r  */
    double density_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Return a scaled Blob.
        The center, density, and behaviour of the new blob is exactly the
        same as the actual one. The radius is multiplied by the scale factor
        and the maximum distance is recalculated for the new sphere. Alpha
        is kept constant*/
    Feature *scale(double factor) const;

    /** Another function for return a scaled Blob.*/
    void scale(double factor, Feature *_f) const;

    /** Intersection of a ray with a blob.
        See \ref Feature::intersection to know more about the parameters
        meaning. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point,
                        Matrix1D<double> &r, Matrix1D<double> &u) const;

    /** Mass of a Blob.
        This function returns mass inside a blob. 3 is the dimension
        See \ref Feature::volume */
    double volume() const
    {
        return basvolume(radius, alpha, m, 3);
    }

    /** Read specific description for a blob.
        An exception is thrown if the line doesn't conform the standard
        specification. See \ref Feature::read_specific */
    void read_specific(char *line);

    /** Print blob in the standard feature format.
        \ See {Feature::feat_printf}, \ref Phantoms */
    void feat_printf(FILE *fh) const;

    /** Show feature not in the standard format but more informatively.
        This function only shows the non common part of the feature. Use
        the << operator of Feature to show the whole feature. */
    friend std::ostream& operator << (std::ostream &o, const Blob &f);

    /** Generate a random blob.
        A sphere is generated randomly within the range of the parameters
        given. Notice that most of them are set by default. The exact
        parameters are picked uniformly from the range. The maximum distance
        is adjusted properly according to the randomly generated feature. 'den'
        stands for density and (x0,y0,z0) is the center of the feature. The
        behaviour of the new sphere is always Add */
    void init_rnd(
        double minradius,   double maxradius,
        double minalpha,    double maxalpha,
        double minorder,    double maxorder,
        double minden = 1.0,  double maxden = 1.0,
        double minx0 = 0,     double maxx0 = 0,
        double miny0 = 0,     double maxy0 = 0,
        double minz0 = 0,     double maxz0 = 0);
};

/* GAUSSIAN ================================================================ */
/** Gaussians.
    The Gaussians are defined by its center, and sigma. The label is "gau".

    A point (r) in the universal coordinate system, where Gaussians are defined
    in general, can be expressed (rp) with respect to a system where the Gaussian
    is centered at the origin and its sigma is unity by
    @code
    V3_MINUS_V3(rp,r,gau.Center);
    V3_BY_CT(rp, rp, 1/sigma);
    @endcode
    Directions (free vectors) are transformed in the same fashion except
    that for the first vector substraction (referring to the origin).
*/

class Gaussian: public Feature
{
public:

    /// Sigma
    double             sigma;
public:
    /** Prepare Gaussian for work.
        Computes the maximum distance, in this case, equal to "4sigma" */
    void prepare();

    /// Assignment
    Gaussian & operator = (const Gaussian &F);

    /** Another function for assignment.*/
    void assign(const Gaussian &F);

    /** Speeded up point inside a Gaussian.
        This function tells you if a point is inside the Gaussian or not.
        See \ref Feature::point_inside */
    int point_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Density inside a Gaussian.
        This function tells you the density of the Gaussian at point r  */
    double density_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Return a scaled Gaussian.
        The center, density, and behaviour of the new Gaussian is exactly the
        same as the actual one. The radius is multiplied by the scale factor
        and the maximum distance is recalculated for the new sphere. Alpha
        is kept constant*/
    Feature *scale(double factor) const;

    /** Another function for return a scaled Blob.*/
    void scale(double factor, Feature *_f) const;

    /** Intersection of a ray with a Gaussian.
        See \ref Feature::intersection to know more about the parameters
        meaning. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point,
                        Matrix1D<double> &r, Matrix1D<double> &u) const;

    /** Mass of a Gaussian.
        This function returns mass inside a Gaussian. 3 is the dimension
        See \ref Feature::volume */
    double volume() const
    {
        return 1;
    }

    /** Read specific description for a Gaussian.
        An exception is thrown if the line doesn't conform the standard
        specification. See \ref Feature::read_specific */
    void read_specific(char *line);

    /** Print Gaussian in the standard feature format.
        \ See {Feature::feat_printf}, \ref Phantoms */
    void feat_printf(FILE *fh) const;

    /** Show feature not in the standard format but more informatively.
        This function only shows the non common part of the feature. Use
        the << operator of Feature to show the whole feature. */
    friend std::ostream& operator << (std::ostream &o, const Gaussian &f);

    /** Generate a random Gaussian.
        A sphere is generated randomly within the range of the parameters
        given. Notice that most of them are set by default. The exact
        parameters are picked uniformly from the range. The maximum distance
        is adjusted properly according to the randomly generated feature. 'den'
        stands for density and (x0,y0,z0) is the center of the feature. The
        behaviour of the new sphere is always Add */
    void init_rnd(
        double minsigma,   double maxsigma,
        double minden = 1.0,  double maxden = 1.0,
        double minx0 = 0,     double maxx0 = 0,
        double miny0 = 0,     double maxy0 = 0,
        double minz0 = 0,     double maxz0 = 0);
};

/* CYLINDERS =============================================================== */
/** Cylinder.
    A cylinder is defined by its center, a radius, a height and a set of
    Euler angles, label="cyl". The center of the cylinder is at the
    geometrical center,
    ie, first the cylinder is placed with its base parallel to XY plane
    (the base is of radius "radius"), and the cylinder is defined between
    "-height/2" to "+height/2". Then the cylinder is rotated after the
    three Euler angles.

    A point (r) in the universal coordinate system, where cylinders are defined
    in general, can be expressed (rp) with respect to a system where the cylinder
    is centered at the origin, its radius and height are unity and its base is
    parallel to the XY plane by
    @code
    V3_MINUS_V3(rp,r,cyl.Center);
    M3x3_BY_V3x1(rp,cyl.euler,rp);
    XX(rp) /= cyl.radius;
    YY(rp) /= cyl.radius;
    ZZ(rp) /= cyl.height;
    @endcode
    Directions (free vectors) are transformed in the same fashion except
    that for the first vector substraction (referring the origin).
*/
class Cylinder: public Oriented_Feature
{
public:
    /// Cylinder X radius
    double             xradius;

    /// Cylinder Y radius
    double             yradius;

    /// Cylinder height
    double             height;
public:
    /** Prepare cylinder for work.
        Computes the maximum distance, in this case, equal to
        "sqrt(height*height/4+radius*radius)", and computes the Euler and inverse
        Euler matrices as a function of the Euler angles */
    void prepare();

/// Assignment
    Cylinder & operator = (const Cylinder &F);

    /** Another function for assignment.*/
    void assign(const Cylinder &F);

    /** Speeded up point inside a cylinder.
        This function tells you if a point is inside the cylinder or not.
        See \ref Feature::point_inside */
    int point_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Density inside an cylinder.
        This function tells you the density of the cylinder at point r
        for constant valued features is trivial*/
    double density_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const
    {
        return (1.);
    }

    /** Return a scaled cylinder.
        The center, density, angles and behaviour of the new cylinder is exactly the
        same as the actual one. The radius and height are multiplied by the
        scale factor and the maximum distance is recalculated for the new cylinder.*/
    Feature *scale(double factor) const;

    /** Another function for return a scaled cylinder.*/
    void scale(double factor, Feature *_f) const;

    /** Intersection of a ray with a cylinder.
        See \ref Feature::intersection to know more about the parameters
        meaning. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point,
                        Matrix1D<double> &r, Matrix1D<double> &u) const;

    /** Volume of a cylinder.
        This function returns 4/3*PI*radius*radius*height.
        See \ref Feature::volume */
    double volume() const
    {
        return 4 / 3*PI*xradius*yradius*height;
    }

    /** Read specific description for a cylinder.
        An exception is thrown if the line doesn't conform the standard
        specification. See \ref Feature::read_specific */
    void read_specific(char *line);

    /** Print cylinder in the standard feature format.
        \ See {Feature::feat_printf}, \ref Phantoms */
    void  feat_printf(FILE *fh) const;

    /** Show feature not in the standard format but more informatively.
        This function only shows the non common part of the feature. Use
        the << operator of Feature to show the whole feature. */
    friend std::ostream& operator << (std::ostream &o, const Cylinder &f);

    /** Generate a random cylinder.
        A cylinder is generated randomly within the range of the parameters
        given. Notice that most of them are set by default. The exact
        parameters are picked uniformly from the range. The maximum distance
        is adjusted properly according to the randomly generated feature. 'den'
        stands for density and (x0,y0,z0) is the center of the feature. The
        behaviour of the new sphere is always Add */
    void  init_rnd(
        double minxradius,  double maxxradius,
        double minyradius,  double maxyradius,
        double minheight,   double maxheight,
        double minden = 1.0,  double maxden = 1.0,
        double minx0 = 0,     double maxx0 = 0,
        double miny0 = 0,     double maxy0 = 0,
        double minz0 = 0,     double maxz0 = 0,
        double minrot = 0,    double maxrot = 360,
        double mintilt = 0,   double maxtilt = 180,
        double minpsi = 0,    double maxpsi = 360);
};

/* DOUBLE CYLINDERS ======================================================== */
/** Double cylinders.
    A double cylinder is defined by its center, a radius, a height and a set of
    Euler angles, label="dcy". The double cylinder are two cylinders whose
    centers are
    (before rotating) "+separation/2+height/2" and "-separation/2-height/2"
    respectively. Ie, two cylinders along Z axis of height "height" and
    separated "separation" units. These two cylinders are
    placed with its base parallel to XY plane
    (the base is of radius "radius"), and the cylinder is defined between
    "-height/2" to "+height/2" starting from their own centers.
    Then the cylinders, as a block, are rotated after the
    three Euler angles.

    A point (r) in the universal coordinate system, where cylinders are defined
    in general, can be expressed (rp) with respect to a system where the FIRST
    cylinder is centered at the origin, its radius is unity and its base is
    parallel to the XY plane (notice that the height is not transformed) by
    @code
    V3_MINUS_V3(rp,r,dcy.Center+dcy.);
    ZZ(r) -= (separation/2+height/2);
    M3x3_BY_V3x1(rp,dcy.euler,rp);
    XX(rp) /= dcy.radius;
    YY(rp) /= dcy.radius;
    @endcode
    and for the SECOND dcyinder
    @code
    V3_MINUS_V3(rp,r,dcy.Center+dcy.);
    ZZ(r) += (separation/2+height/2);     // This is the only one line changing
    M3x3_BY_V3x1(rp,dcy.euler,rp);
    XX(rp) /= dcy.radius;
    YY(rp) /= dcy.radius;
    @endcode
    Directions (free vectors) are transformed in the same fashion except
    that for the first two vector substraction (referring the origin).
*/
class DCylinder: public Oriented_Feature
{
public:
    /// Cylinder radius
    double             radius;

    /// Each cylinder height
    double             height;

    /// Separation between cylinders
    double             separation;
public:
    /** Prepare double cylinder for work.
        Computes the maximum distance, in this case, equal to
        "sqrt((height+separation)*(height+separation)/4+radius*radius)",
        and computes the Euler and inverse
        Euler matrices as a function of the Euler angles */
    void prepare();

/// Assignment
    DCylinder & operator = (const DCylinder &F);

    /** Another function for assignment.*/
    void assign(const DCylinder &F);

    /** Speeded up point inside a double cylinder.
        This function tells you if a point is inside any of the cylinders or not.
        See \ref Feature::point_inside */
    int point_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Density inside a double cylinder.
        This function tells you the density of the double cylinder at point r
        for constant valued features is trivial*/
    double density_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const
    {
        return (1.);
    }

    /** Return a scaled double cylinder.
        The center, density, angles and behaviour of the new double cylinder
        is exactly the
        same as the actual one. The radius and height are multiplied by the
        scale factor and the maximum distance is recalculated for the new
        double cylinder.*/
    Feature *scale(double factor) const;

    /** Another function for return a scaled double cylinder.*/
    void scale(double factor, Feature *_f) const;

    /** Intersection of a ray with a double cylinder.
        See \ref Feature::intersection to know more about the parameters
        meaning. The ray intersection consider both cylinders, of course. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point,
                        Matrix1D<double> &r, Matrix1D<double> &u) const;

    /** Volume of a double cylinder.
        This function returns 2* 4/3*PI*radius*radius*height.
        See \ref Feature::volume */
    double volume() const
    {
        return 2* 4 / 3*PI*radius*radius*height;
    }

    /** Read specific description for a double cylinder.
        An exception is thrown if the line doesn't conform the standard
        specification. See \ref Feature::read_specific */
    void read_specific(char *line);

    /** Print double cylinder in the standard feature format.
        \ See {Feature::feat_printf}, \ref Phantoms */
    void feat_printf(FILE *fh) const;

    /** Show feature not in the standard format but more informatively.
        This function only shows the non common part of the feature. Use
        the << operator of Feature to show the whole feature. */
    friend std::ostream& operator << (std::ostream &o, const DCylinder &f);

    /** Generate a random double cylinder.
        A double cylinder is generated randomly within the range of the parameters
        given. Notice that most of them are set by default. The exact
        parameters are picked uniformly from the range. The maximum distance
        is adjusted properly according to the randomly generated feature. 'den'
        stands for density and (x0,y0,z0) is the center of the feature. The
        behaviour of the new sphere is always Add. 'Sep' is the range for the
        separation */
    void init_rnd(
        double minradius,   double maxradius,
        double minheight,   double maxheight,
        double minsep,      double maxsep,
        double minden = 1.0,  double maxden = 1.0,
        double minx0 = 0,     double maxx0 = 0,
        double miny0 = 0,     double maxy0 = 0,
        double minz0 = 0,     double maxz0 = 0,
        double minrot = 0,    double maxrot = 360,
        double mintilt = 0,   double maxtilt = 180,
        double minpsi = 0,    double maxpsi = 360);
};

/* CUBE ==================================================================== */
/** Cube.
    A cube is defined by its center, and the length in the 3 dimensions
    (X, Y and Z, before rotating), label="cub". The center of the cube is
    at the geometrical center, and the dimensions define the cube from
    "-Xdim/2" to "Xdim/2", ... around the center. Then the cube is rotated
    after the
    three Euler angles. The corresponding Euler matrix is stored in the
    field "euler" while its inverse is in "eulert".

    A point (r) in the universal coordinate system, where cubes are defined
    in general, can be expressed (rp) with respect to a system where the cube
    is centered at the origin, all its dimensions are unity and its axes
    are aligned with XYZ by
    @code
    V3_MINUS_V3(rp,r,cub.Center);
    M3x3_BY_V3x1(rp,cub.euler,rp);
    XX(rp) /= cub.xdim;
    YY(rp) /= cub.ydim;
    ZZ(rp) /= cub.zdim;
    @endcode
    Directions (free vectors) are transformed in the same fashion except
    that for the first vector substraction (referring the origin).
*/
class Cube: public Oriented_Feature
{
public:
    /// X dimension before rotating
    double             xdim;

    /// Y dimension before rotating
    double             ydim;

    /// Z dimension before rotating
    double             zdim;
public:
    /** Prepare cube for work.
        Computes the maximum distance, in this case, equal to
        "sqrt(xdim*xdim+ydim*ydim+zdim*zdim)", and computes the Euler and inverse
        Euler matrices as a function of the Euler angles */
    void prepare();

/// Assignment
    Cube & operator = (const Cube &F);

    /** Another function for assignment.*/
    void assign(const Cube &F);

    /** Speeded up point inside a cube.
        This function tells you if a point is inside the cube or not.
        See \ref Feature::point_inside */
    int point_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Density inside an Cube.
        This function tells you the density of the cube at point r
        for constant valued features is trivial*/
    double density_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const
    {
        return (1.);
    }

    /** Return a scaled cube.
        The center, density, angles and behaviour of the new cube is exactly the
        same as the actual one. The dimensions are multiplied by the
        scale factor and the maximum distance is recalculated for the new cube.*/
    Feature *scale(double factor) const;

    /** Another function for return a scaled cube.*/
    void scale(double factor, Feature *_f) const;

    /** Intersection of a ray with a cube, NOT IMPLEMENTED!!!.
        See \ref Feature::intersection to know more about the parameters
        meaning. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point,
                        Matrix1D<double> &r, Matrix1D<double> &u) const;

    /** Volume of a cube.
        This function returns xdim*ydim*zdim.
        See \ref Feature::volume */
    double volume() const
    {
        return xdim*ydim*zdim;
    }

    /** Read specific description for a cube.
        An exception is thrown if the line doesn't conform the standard
        specification. See \ref Feature::read_specific */
    void read_specific(char *line);

    /** Print cube in the standard feature format.
        \ See {Feature::feat_printf}, \ref Phantoms */
    void feat_printf(FILE *fh) const;

    /** Show feature not in the standard format but more informatively.
        This function only shows the non common part of the feature. Use
        the << operator of Feature to show the whole feature. */
    friend std::ostream& operator << (std::ostream &o, const Cube &f);

    /** Generate a random cube.
        A cube is generated randomly within the range of the parameters
        given. Notice that most of them are set by default. The exact
        parameters are picked uniformly from the range. The maximum distance
        is adjusted properly according to the randomly generated feature. 'den'
        stands for density and (x0,y0,z0) is the center of the feature. The
        behaviour of the new sphere is always Add. If Ydim and Zdim ranges are
        restricted to 0, then it's applied the same range as for Xdim. */
    void  init_rnd(
        double minXdim,   double maxXdim,
        double minYdim = 0, double maxYdim = 0,
        double minZdim = 0, double maxZdim = 0,
        double minden = 1.0,  double maxden = 1.0,
        double minx0 = 0,   double maxx0 = 0,
        double miny0 = 0,   double maxy0 = 0,
        double minz0 = 0,   double maxz0 = 0,
        double minrot = 0,  double maxrot = 360,
        double mintilt = 0, double maxtilt = 180,
        double minpsi = 0,  double maxpsi = 360);
};

/* ELLIPSOID =============================================================== */
/** Ellipsoid.
    An elliposid is defined by its center, and the radius in the 3 dimensions
    (X, Y and Z, before rotating), label="ell". The center of the ellipsoid is
    at the geometrical center, and the dimensions define the elliposid from
    "-Xradius/2" to "Xradius/2", ... around the center. Then the ellipsoid
    is rotated after the
    three Euler angles. The corresponding Euler matrix is stored in the
    field "euler" while its inverse is in "eulert".

    A point (r) in the universal coordinate system, where elliposids are defined
    in general, can be expressed (rp) with respect to a system where the ellipsoid
    is centered at the origin, all its radii are unity and its axes
    are aligned with XYZ by
    @code
    V3_MINUS_V3(rp,r,ell.Center);
    M3x3_BY_V3x1(rp,ell.euler,rp);
    XX(rp) /= ell.xradius;
    YY(rp) /= ell.yradius;
    ZZ(rp) /= ell.zradius;
    @endcode
    Directions (free vectors) are transformed in the same fashion except
    that for the first vector substraction (referring the origin).
*/
class Ellipsoid: public Oriented_Feature
{
public:
    /// X radius before rotating
    double             xradius;

    /// Y radius before rotating
    double             yradius;

    /// Z radius before rotating
    double             zradius;
public:
    /** Prepare ellipsoid for work.
        Computes the maximum distance, in this case, equal to
        "MAX(MAX(xradius,yradius),zradius)", and computes the Euler and inverse
        Euler matrices as a function of the Euler angles */
    void prepare();

/// Assignment
    Ellipsoid & operator = (const Ellipsoid &F);

    /** Another function for assignment.*/
    void assign(const Ellipsoid &F);

    /** Speeded up point inside an ellipsoid.
        This function tells you if a point is inside the elliposoid or not.
        See \ref Feature::point_inside */
    int point_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Density inside a Ellipsoid.
        This function tells you the density of the ellipsoid at point r
        for constant valued features is trivial*/
    double density_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const
    {
        return (1.);
    }

    /** Return a scaled elliposoid.
        The center, density, angles and behaviour of the new ellipsoid is exactly the
        same as the actual one. The dimensions are multiplied by the
        scale factor and the maximum distance is recalculated for the new ellipsoid.*/
    Feature *scale(double factor) const;

    /** Another function for return a scaled elliposoid.*/
    void scale(double factor, Feature *_f) const;

    /** Intersection of a ray with an ellipsoid.
        See \ref Feature::intersection to know more about the parameters
        meaning. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point,
                        Matrix1D<double> &r, Matrix1D<double> &u) const;

    /** Volume of an ellipsoid.
        This function returns 4/3*PI*xradius*yradius*zradius.
        See \ref Feature::volume */
    double volume() const
    {
        return 4 / 3*PI*xradius*yradius*zradius;
    }

    /** Read specific description for an ellipsoid.
        An exception is thrown if the line doesn't conform the standard
        specification. See \ref Feature::read_specific */
    void read_specific(char *line);

    /** Print ellipsoid in the standard feature format.
        \ See {Feature::feat_printf}, \ref Phantoms */
    void feat_printf(FILE *fh) const;

    /** Show feature not in the standard format but more informatively.
        This function only shows the non common part of the feature. Use
        the << operator of Feature to show the whole feature. */
    friend std::ostream& operator << (std::ostream &o, const Ellipsoid &f);

    /** Generate a random ellipsoid.
        An ellipsoid is generated randomly within the range of the parameters
        given. Notice that most of them are set by default. The exact
        parameters are picked uniformly from the range. The maximum distance
        is adjusted properly according to the randomly generated feature. 'den'
        stands for density and (x0,y0,z0) is the center of the feature. The
        behaviour of the new sphere is always Add. */
    void init_rnd(
        double minXradius, double maxXradius,
        double minYradius, double maxYradius,
        double minZradius, double maxZradius,
        double minden = 1.0, double maxden = 1.0,
        double minx0 = 0,    double maxx0 = 0,
        double miny0 = 0,    double maxy0 = 0,
        double minz0 = 0,    double maxz0 = 0,
        double minrot = 0,   double maxrot = 360,
        double mintilt = 0,  double maxtilt = 180,
        double minpsi = 0,   double maxpsi = 360);
};

/* CONE ==================================================================== */
/** Cone.
    A cone is defined by its center, a radius, a height and a set of
    Euler angles, label="con". The center of the cone is at the
    geometrical center,
    ie, first the cone is placed with its base parallel to XY plane
    (the base is of radius "radius"), and the cone is defined between
    "-height/2" to "+height/2". Then the cone is rotated after the
    three Euler angles. The corresponding Euler matrix is stored in the
    field "euler" while its inverse is in "eulert".

    A point (r) in the universal coordinate system, where cones are defined
    in general, can be expressed (rp) with respect to a system where the cone
    is centered at the origin, its radius and height is unity and its base is
    parallel to the XY plane by
    @code
    V3_MINUS_V3(rp,r,con.Center);
    M3x3_BY_V3x1(rp,con.euler,rp);
    XX(rp) /= con.radius;
    YY(rp) /= con.radius;
    ZZ(rp) /= con.height;
    @endcode
    Directions (free vectors) are transformed in the same fashion except
    that for the first vector substraction (referring the origin).
*/
class Cone: public Oriented_Feature
{
public:
    /// Cone base radius
    double             radius;

    /// Cone height
    double             height;

public:
    /** Prepare cone for work.
        Computes the maximum distance, in this case, equal to
        "sqrt(height*height/4+radius*radius)", and computes the Euler and inverse
        Euler matrices as a function of the Euler angles */
    void prepare();

/// Assignment
    Cone & operator = (const Cone &F);

    /** Another function for assignment.*/
    void assign(const Cone &F);

    /** Speeded up point inside a cone.
        This function tells you if a point is inside the cone or not.
        See \ref Feature::point_inside */
    int point_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const;

    /** Density inside a cone.
        This function tells you the density of the cone at point r
        for constant valued features is trivial*/
    double density_inside(const Matrix1D<double> &r, Matrix1D<double> &aux) const
    {
        return (1.);
    }

    /** Return a scaled cone.
        The center, density, angles and behaviour of the new cone is exactly the
        same as the actual one. The radius and height are multiplied by the
        scale factor and the maximum distance is recalculated for the new cone.*/
    Feature *scale(double factor) const;

    /** Another function for return a scaled cone.*/
    void scale(double factor, Feature *_f) const;

    /** Intersection of a ray with a cone, NOT IMPLEMENTED!!!!.
        See \ref Feature::intersection to know more about the parameters
        meaning. */
    double intersection(const Matrix1D<double> &direction,
                        const Matrix1D<double> &passing_point,
                        Matrix1D<double> &r, Matrix1D<double> &u) const;

    /** Volume of a cone.
        This function returns 1/3*PI*radius*radius*height.
        See \ref Feature::volume */
    double volume() const
    {
        return 1 / 3*PI*radius*radius*height;
    }

    /** Read specific description for a cone.
        An exception is thrown if the line doesn't conform the standard
        specification. See \ref Feature::read_specific */
    void read_specific(char *line);

    /** Print cone in the standard feature format.
        \ See {Feature::feat_printf}, \ref Phantoms */
    void feat_printf(FILE *fh) const;

    /** Show feature not in the standard format but more informatively.
        This function only shows the non common part of the feature. Use
        the << operator of Feature to show the whole feature. */
    friend std::ostream& operator << (std::ostream &o, const Cone &f);

    /** Generate a random cone.
        A cone is generated randomly within the range of the parameters
        given. Notice that most of them are set by default. The exact
        parameters are picked uniformly from the range. The maximum distance
        is adjusted properly according to the randomly generated feature. 'den'
        stands for density and (x0,y0,z0) is the center of the feature. The
        behaviour of the new sphere is always Add */
    void init_rnd(
        double minradius,   double maxradius,
        double minheight,   double maxheight,
        double minden = 1.0,  double maxden = 1.0,
        double minx0 = 0,     double maxx0 = 0,
        double miny0 = 0,     double maxy0 = 0,
        double minz0 = 0,     double maxz0 = 0,
        double minrot = 0,    double maxrot = 360,
        double mintilt = 0,   double maxtilt = 180,
        double minpsi = 0,    double maxpsi = 360);
};

/* PHANTOM ================================================================= */
/** Phantom class.
    The phantom class is simply a list (STL vector) of features plus some
    information about the size of the final volume to generate and its
    background density. This is the class that will interact with the
    reconstruction programs as the features classes themselves haven't
    got enough information to generate the final volume. The file format
    to generate the phantom is described in the previous page (\ref Phantoms).

    This class is thought to be filled from a file, and doesn't give
    many facilities to update it from program. This is something
    to do.

    Here goes an example of how to manage loops in the phantom class,
    @code
       // Show all features
       for (int i=1; i<=P.FeatNo(); i++) std::cout << P(i);
    @endcode
*/
class Phantom
{
public:
    /// Filename
    FileName       fn;

    /// Final volume X dimension
    int            xdim;

    /// Final volume Y dimension
    int            ydim;

    /// Final volume Z dimension
    int            zdim;

    /// Final volume background density
    double          Background_Density;

    /// Has been the  scale applied?
    double          current_scale;

    ///  Param file scale
    double          param_file_scale;

    ///  Param file scale
    double          phantom_scale;

    /// List with the features
    std::vector<Feature*> VF;
public:
    /** Empty constructor.
        The empty phantom is 0x0x0, background density=0, no feature is inside
        and no name. */
    Phantom();

    /** Construct from a phantom file.
        Construct the phantom according to the specifications of the given file.
        The file must accomplish the structure given in \ref Phantoms. */
    Phantom(const FileName &fn_phantom)
    {
        read(fn_phantom);
    }

    /** Destructor.
        All features are freed. */
    ~Phantom()
    {
        clear();
    }

    /** Clear the phantom.
        Force the phantom to be empty. All features are freed. */
    void clear();

    /** Returns the number of features in the list. */
    int FeatNo()
    {
        return VF.size();
    }

    /** Access to a feature pointer.
        You can address each one of the features using this operator,
        remember that features are numbered as 1, 2, ... FeatNo() */
    Feature * operator()(int i)
    {
        return VF[i-1];
    }

    /** Constant Access to a feature pointer.
        Same as the previous one but with constant access. */
    const Feature * operator()(int i) const
    {
        return VF[i-1];
    }

    /** Add a feature to the phantom.
        This function allows you to add new features at the end of the phantom.
        \\ Ex: Sphere S; S.radius(); S.prepare(); Phantom.add(&S); */
    void add(Feature *f)
    {
        VF.push_back(f);
    }

    /// Assignment
    Phantom & operator = (const Phantom &P);

    /** Another function for assignment.*/
    void assign(const Phantom &P);

    /// Prepare for work.
    void prepare();

    /// Return the maximum distance of any feature to the volume center
    double max_distance() const;

    /// Return the volume of all the features
    double volume() const;

    /** Read a phantom file.
        The file must accomplish the structure given in \ref Phantoms .

        If you don't apply the scale then all spatial coordinates are
        expressed in the given scale units. I.e., if the scale is 0.25 that
        means that every voxel in the voxel phantom represent 4 length units
        (usually Angstroms). If the scale is applied, then coordinates are
        expressed in voxel edge units. Otherwise, coordinates are expressed
        in length units (usually Angstroms).

        We recommend apply_scale=false only if you plan to modify the
        description file without produstd::cing projections, voxel phantoms, ... */
    void read(const FileName &fn_phantom, bool apply_scale = true);

    /** Show a phantom file.
        The more descriptive format is used instead of the standard Feature format */
    friend std::ostream& operator << (std::ostream &o, const Phantom &f);

    /** Write a phantom file in the standard feature format.
        You may rename the file or not giving a different name in the write call. */
    void write(const FileName &fn_phantom = "");

    /** Speeded up voxel inside any feature.
        This function tells you if the voxel of size 1x1x1 whose center is at
        position r is inside any of the features. A voxel is said to be inside
        a feature if it shares at least 1 corner with the feature.
        This function returns the first feature of the list containing the voxel.
        Features are numbered as 1, 2, ...
        If the voxel is not in any feature the function returns 0.

        aux1 and aux2, are two vectors of dimension 3. They must be supplied
        in order to gain speed in the calculations. This is very useful
        when checking if many voxels are inside any feature. See also
        \ref Feature::voxel_inside to know more. */
    int voxel_inside_any_feat(const Matrix1D<double> &r,
                              Matrix1D<double> &aux1, Matrix1D<double> &aux2) const;

    /** Voxel inside any feature.
        The same as the previous one but you needn't supply the extra auxiliar
        vectors. */
    int voxel_inside_any_feat(const Matrix1D<double> &r) const
    {
        Matrix1D<double> aux1(3), aux2(3);
        return voxel_inside_any_feat(r, aux1, aux2);
    }

    /** Speeded up sphere intersecting any feature.
        This function returns the first feature in the list intersecting
        a sphere with center r in R3 and the given radius. In none, 0 is
        returned.
        This speeded up function needs
        two vectors with dimension 3 externally resized. */
    int any_feature_intersects_sphere(const Matrix1D<double> &r, double radius,
                                      Matrix1D<double> &aux1, Matrix1D<double> &aux2, Matrix1D<double> &aux3)
    const;

    /** Sphere intersecting feature.
        This function is based in the previous one. It makes the same but
        you needn't supply the auxiliar vectors. */
    int any_feature_intersects_sphere(const Matrix1D<double> &r,
                                      double radius) const
    {
        Matrix1D<double> aux1(3), aux2(3), aux3(3);
        return any_feature_intersects_sphere(r, radius, aux1, aux2, aux3);
    }

    /** Draw the phantom in the volume.
        The volume is cleaned, resized to the phantom size and its origin
        is set at the center of the volume. Then every feature is drawn into
        the volume */
    void draw_in(Volume *V);

    /** Label a volume after the phantom.
        The volume is cleaned, resized to the phantom size and its origin
        is set at the center of the volume. Then every feature is drawn into
        the volume using different densities. Background has got density 0,
        the border for first feature density -1,
        the first feature has got density 1, the second 2 and its border -2,
        ... */
    void label(Volume *V);

#define DONT_CLEAN 0
#define CLEAN      1
    /** Sketch the surface of the phantom in the volume.
        This function allows you to draw only the surface of every feature
        (see \ref Feature::sketch_in to see when a voxel is said to belong
        to the feature surface). The input volume might be cleaned (resized
        and the logical origin set at the physical center) or not according
        to the labels CLEAN or DONT_CLEAN (by default). The grey level for
        those voxels can be defined (by default, 2) and the colour is
        assigned */
    void sketch_in(Volume *V, int clean = DONT_CLEAN, double colour = 2);

    /** Shift.
        Shift all features in the phantom a given amount of voxels. See
        \ref Feature::shift. */
    void shift(double shiftX, double shiftY, double shiftZ);

    /** Rotate.
        Rotate a phantom after a given 3D rotation. */
    void rotate(const Matrix2D<double> &E);

    /** Apply a general geometric transformation.
        The transformation must be defined by a 4x4 matrix that can be
        generated using the geometric functions in xmippGeometry or
        xmippMatrix2D. The inv argument is either IS_INV or IS_NOT_INV.
        Matrices coming out from xmippGeometry for rotations, shifts, ...
        are not INV.

        An exception is thrown if the matrix A is not valid.*/
    void selfApplyGeometry(const Matrix2D<double> &A, int inv);

    /** Project phantom from a direction.
        The direction is specified by the 3 Euler angles (as usual, rot=1st,
        tilt=2nd, psi=3rd), the projection size is (Ydim x Xdim) given
        in the function call, the projection logical origin is set at the
        physical center of the image. The projection is cleaned and then
        the projection of the phantom is computed. A matrix A (3x3) can be supplied
        in order to apply a deformation in the projection plane. A must
        be such that
        @code
        deformed position=A*undeformed position
        @endcode*/
    void project_to(Projection &P, int Ydim, int Xdim,
                    double rot, double tilt, double psi, const Matrix2D<double> *A = NULL) const;

    /** Project phantom from a direction.
        The same as before but this time the projection is supposed to be
        already resized and with the right center. The phantom projection
        is added to the already drawn projection. */
    void project_to(Projection &P,
                    double rot, double tilt, double psi, const Matrix2D<double> *A = NULL) const;

    /** Project phantom using a conversion matrix.
        The same as before but this time the projection is supposed to be
        already resized and with the right center. The phantom projection
        is added to the already drawn projeciton. Euler angles of the projection
        are not set and the phantom volume is projected to the projection plane
        making use of the matrix provided (3x3) which is the needed transformation
        from the volume to the projection plane. Inside the projection plane
        only the first two coordinates are valid. If disappearing_th < 1
        then the ith phantom feature is skiped with a provability =
        1-disappearing_th*/

    void project_to(Projection &P, const Matrix2D<double> &VP, double    disappearing_th = 1.0) const;

#define POS_NEG 1
#define NEG_POS 2
    /** Surface of the phantom as in AFM.
        This function produces the surface of a phantom in a given image
        when measuring rays go POS_NEG or NEG_POS (from Z positive to Z negative,
        and viceversa). If rays reach point z0
        and the phantom is not reached the surface for that point is z0,
        otherwise it is the z of the last voxel which is still outside the phantom.

        This function throws an exception if z0 is outside the volume definition.
        If z0 is equal to zdim then it is not used.

        The radius is the probe radius for the surface generation.

        If the output image is not resized, it is resized to the Y and X dimensions
        of the phantom.
        */
    void surface(double z0, double radius, int direction, Image *P)
    const;
};
//@}
#endif
