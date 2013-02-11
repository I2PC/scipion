/***************************************************************************
 *
 * Authors:    Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
/* GRIDS                                                                     */
/* ------------------------------------------------------------------------- */

#ifndef _GRIDS_HH
#define _GRIDS_HH

#include <vector>

#include "xmipp_image.h"
#include "geometry.h"
#include "args.h"

/* Forward declarations ---------------------------------------------------- */
template <class T>
class GridVolumeT;
template <class T>
GridVolumeT<T> operator -(T f, const GridVolumeT<T> &GV);
template <class T>
GridVolumeT<T> operator /(T f, const GridVolumeT<T> &GV);
template <class T>
std::ostream& operator << (std::ostream &o, const GridVolumeT<T> &GV);
class Basis;

/**@defgroup Grids Grids
   @ingroup DataLibrary
    The grids are one of the most basic things in the reconstruction
    process, since the reconstructed volumes are expressed as a linear
    combination of a volume basis function weighted and shifted to all
    positions defined by the grid. Grids in Xmipp may be as complex as
    you liked, a complex grid is supposed to be a superposition of
    simpler grids. Simple grids maintain information about the simple
    grids themselves (orientation, spacing, ...) while complex
    ones are only collections of simple grids. Usual grids as BCC
    (Body Centered Cubic) or FCC (Face Centered Cubic),
    can be expressed as the superposition of two CC (Cubic) grids.

    It is important to notice that the existence of this class makes
    the reconstruction algorithms independent from the underlying
    grid, incrementing the reusability of the code.
*/
//@{
/*****************************************************************************/
/* Simple Grids                                                              */
/*****************************************************************************/
/** Basic grid class.
    A Simple grid is defined as a set of 3 directions (X,Y,Z) which need not
    to coincide with the universal X,Y,Z (ie, e1, e2, e3), an origin,
    a relative size to the voxels in the universal grid, and the lowest
    and highest indexes within the grid (these two last parameters set the
    grid size and the space where the grid is defined). There will be
    volume elements at each point of the grid.

    The concept of a simple grid is the following: you are defining a
    coordinate system (translated, rotated and whatever you like with
    respect to the Universal Coordinate System, and inside this system
    you must define a "measuring unit". This coordinate system will be
    valid only in a "box" of the space. A "box" is quoted as because
    it is a parallelepiped defined by this coordinate system axes, and
    not by the Universal axes.

    So to define the simple grid you must address the following topics:
    <ul>
    <li>Coordinate system:
       The coordinate system is defined by three R3 vectors, called
       X, Y and Z. They must form a true three dimensional coordinate
       system, ie, must not be linearly dependent, but they are not
       restricted to be orthogonal, normal, or whatever. You can specify
       any R3 vectors with the only condition that the range of the
       matrix formed by the 3 is not 0.

       When a basis function is placed at position (-1,1,2) in the
       grid coordinate system, it is really placed at (-X+Y+2Z)+origin
       in the
       universal coordinate system (these X,Y,Z are the grid axis
       vectors and must not be confused with the e1(=x), e2(=y), and
       e3(=z) of the Universal Coordinate System). This is true
       only if the measuring unit is 1, see the Measuring unit section
       to know exactly where the samples is placed in space.

       The origin of the grid is the position in the Universal coordinate
       system of the (0,0,0) sample inside the grid, ie, where is the
       grid sample (0,0,0) in the universal coordinate system?
    <li> Measuring unit:
       The measuring unit is the distance between two samples in the
       grid lattice, and in the direction of the grid axes. Suppose
       we have a sample at position (-1,1,2) but the measuring
       unit instead of being 1 is 2. Then the sample really is at
       2(-X+Y+2Z)+origin in the Universal Coordinate System.

       This measuring unit is controlled in the class by the public
       variable \ref relative_size .
    <li>Size:
       Again the size is given again in the grid system units and it
       expresses which the minimum and maximum indexes will be in the
       grid coordinate system. This means that if the lowest and highest
       indexes are (-1,-1,-1) and (1,1,1) then the only valid indexes
       for each direction are -1, 0 and 1. The grid is defined in the
       Universal Coordinate System between
       (-X-Y-Z)*relative_size to (X+Y+Z)*relative_size.
    </ul>

    This is the default grid after creation:
    @code
       (e1,e2,e3) vectors
       grid relative_size = 1
       origin    = ( 0, 0, 0)
       lowest    = (-5,-5,-5)
       highest   = ( 5, 5, 5)
    @endcode

    It should be convenient that you modify these parameters at your
    convinience. And after setting them PREPARE!! the grid to be used
    with the function \ref SimpleGrid::prepare_grid . There are several public
    variables which you might directly modify, namely, \ref origin ,
    \ref relative_size , \ref lowest , and \ref highest .

*/
class SimpleGrid
{
    /* Structure --------------------------------------------------------------- */
public:
    Matrix2D<double> basis;                // These are the 3 unitary vectors
    // which define the grid. First
    // column is vector X, second
    // column is vector Y, and the third
    // is Z
    Matrix2D<double> inv_basis;            // Inverse matrix of the basis
    // It is used to save computation
    // time
    /** Lowest index inside the grid coordinate system.
        Although it is a float vector, it should be kept integer all time */
    Matrix1D<double> lowest;
    /** Highest index inside the grid coordinate system.
        Although it is a float vector, it should be kept integer all time */
    Matrix1D<double> highest;
    /// Measuring unit in the grid coordinate system
    double relative_size;
    /// Origin of the grid in the Universal coordinate system
    Matrix1D<double> origin;
    /** Reconstructed sphere squared radius.
        The reconstruction is supposed to fit within this sphere centered
        at the origin. If this radius is -1, then this feature is not used
        and the reconstruction is done all over the cube. */
    double R2;

    /* Prototypes --------------------------------------------------------------- */
public:
    /** Default constructor.
        Be careful that this is not an empty constructor as usual but it
        sets a grid by default (see the class documentation to see exactly
        which grid is).
        \\ Ex: SimpleGrid sg; */
    SimpleGrid();

    /** Copy constructor.
        This constructor builds an exact copy of the simple grid.
        \\ Ex: SimpleGrid sg2(sg1); */
    SimpleGrid(const SimpleGrid &SG);

    /** Show a Simple grid.
        Shows all information about the simple grid.
        \\Ex: std::cout << sg; */
    friend std::ostream& operator <<(std::ostream& o, const SimpleGrid &grid);

    /** Assignment.
        \\ Ex: sg2=sg1; */
    SimpleGrid& operator = (const SimpleGrid &SG);

    /** Another function for assigment.*/
    void assign(const SimpleGrid &SG);

    /** Set X vector of the grid.
        \\Ex: Matrix1D<double> X=vectorR3(1,0,0); sg.set_X(X); */
    void set_X(const Matrix1D<double> &v)
    {
        basis.setCol(0, v);
    }

    /** Set Y vector of the grid.
        \\Ex: Matrix1D<double> Y=vectorR3(0,1,0); sg.set_Y(Y); */
    void set_Y(const Matrix1D<double> &v)
    {
        basis.setCol(1, v);
    }

    /** Set Z vector of the grid.
        \\Ex: Matrix1D<double> Z=vectorR3(1,1,1); sg.set_Z(Z); */
    void set_Z(const Matrix1D<double> &v)
    {
        basis.setCol(2, v);
    }

    /** Get X vector of the grid.
        \\Ex: Matrix1D<double> X; sg.get_X(X) << std::endl; */
    void get_X(Matrix1D<double> &v) const
    {
        basis.getCol(0, v);
    }

    /** Get Y vector of the grid.
        \\Ex: Matrix1D<double> Y; sg.get_Y(Y) << std::endl; */
    void get_Y(Matrix1D<double> &v) const
    {
        basis.getCol(1, v);
    }

    /** Get Z vector of the grid.
        \\Ex: Matrix1D<double> Z; sg.get_Z(Z) << std::endl; */
    void get_Z(Matrix1D<double> &v) const
    {
        basis.getCol(2, v);
    }

    /** Get grid number of samples in each direction.
        This function returns the number of samples on each direction,
        this number of samples can be used to resize a volume which might
        hold a grid like this one. The number of samples is computed as
        highest-lowest+1.
        \\ Ex: Build a volume able to hold this grid
        @code
        grid.getSize(Zdim,Ydim,Xdim);
        vol().resize(Zdim,Ydim,Xdim);
        STARTINGX(Vol_aux())=(int) XX(grid.lowest);
        STARTINGY(Vol_aux())=(int) YY(grid.lowest);
        STARTINGZ(Vol_aux())=(int) ZZ(grid.lowest);
        @endcode
    */
    void getSize(int &Zdim, int &Ydim, int &Xdim) const
    {
        Zdim = (int)(ZZ(highest) - ZZ(lowest)) + 1;
        Ydim = (int)(YY(highest) - YY(lowest)) + 1;
        Xdim = (int)(XX(highest) - XX(lowest)) + 1;
    }

    /** Get number of samples.
       This function returns the number of samples within the grid.
       If the grid has a radius of interest, then only those samples
       within that radius are accounted. */
    int get_number_of_samples() const;

    /// Set reconstruction radius
    void set_interest_radius(double _R)
    {
        if (_R == -1)
            R2 = -1;
        else
            R2 = _R * _R;
    }

    /// Get reconstruction radius
    double get_interest_radius() const
    {
        if (R2 == -1)
            return R2;
        else
            return sqrt(R2);
    }

    /// Set relative_size
    void set_relative_size(double size)
    {
        relative_size = size;
    }

    /// Get relative_size
    double get_relative_size() const
    {
        return (relative_size);
    }

    /** Prepare grid for work.
        This function MUST BE CALLED before working with a simple grid
        which is not the default one. What it actually does is to check
        that the grid vectors form a true 3D coordinate system (if not
        it throws an exception), and if it is calculate the pass matrices
        from the grid to the universal system (called basis) and viceversa
        (called inv_basis).
        \\ Ex:
        @code
        SimpleGrid sg;
        sg.set_Z(vectorR3(1,1,1)); --> Change grid vectors
        sg.prepare_grid();          --> Now the grid is ready to work
        @endcode */
    void prepare_grid();

    /** Grid --> Universe.
        This function transforms a vector in the grid coordinate system
        to the universal one. For example the sample at grid coordinate
        (-1,2,1) is translated to its position in the universe. Notice
        that usually this operation is performed with "integer" grid
        vectors giving "float" universal vectors.
        The operation performed exactly is
        @code
        return origin+basis*gv*relative_size;
        @endcode
        or what is the same:
        @code
        return origin+(XX(gv)*X+YY(gv)*Y+ZZ(gv)*Z)*relative_size;
        @endcode
        \\Ex: Matrix1D<double> uv; sg.grid2universe(vectorR3(-1,2,1),uv); */
    void grid2universe(const Matrix1D<double> &gv, Matrix1D<double> &uv) const
    {
        SPEED_UP_temps012;
        uv.resize(3);
        M3x3_BY_V3x1(uv, basis, gv);
        V3_BY_CT(uv, uv, relative_size);
        V3_PLUS_V3(uv, uv, origin);
    }

    /** Universe --> Grid.
        This function transforms a vector in the universal coordinate system
        to the one of the grid. Ie, it gives the position in the grid system
        of a point in the universe. Notice that this time the position in the
        grid needs not to be "integer".

        The operation performed is
        @code
        return inv_basis*(uv-origin)/relative_size;
        @endcode
        where inv_basis is an internal matrix which contains information
        about the X, Y and Z vectors. In fact is the inverse of the
        matrix formed by X, Y and Z. See \ref grid2universe .

        Next example should give (0,0,0) as result.
        \\Ex:
        @code
           SimpleGrid sg;
           sg.origin=vectorR3(10,10,10);
           std::cout << "What is the position within the grid of the origin? "
                << sg.universe2grid(sg.origin);
        @endcode */
    void universe2grid(const Matrix1D<double> &uv, Matrix1D<double> &gv) const
    {
        SPEED_UP_temps012;
        gv.resize(3);
        V3_MINUS_V3(gv, uv, origin);
        V3_BY_CT(gv, gv, 1 / relative_size);
        M3x3_BY_V3x1(gv, inv_basis, gv);
    }

    /** TRUE if the selected universe coordinate is interesting.
        A point is interesting if it is inside the reconstruction radius.
        In case that this radius is -1, then all points are interesting. */
    bool is_interesting(const Matrix1D<double> &uv) const
    {
        if (R2 == -1)
            return true;
        else
            return XX(uv)*XX(uv) + YY(uv)*YY(uv) + ZZ(uv)*ZZ(uv) < R2;
    }

    /** Project a grid vector onto a plane (Euler angles).
        This function projects a vector defined in the Grid coordinate
        system onto a plane defined by its three Euler angles (measured
        with respect to the universal coordinate system). This plane
        is restricted to pass through the Universe origin. This function
        is used to produce the projection of a volume according to
        a certain direction (the 3 Euler angles in the universal system).
        Notice that this function only tells you where a given point in
        the volume is projecting, but it doesn't tell you the value.
        See \ref Uproject_to_plane for more
        information.
        \\ Ex: Matrix1D<double> up=sg.Gproject_to_plane(vectorR3(1,0,0),45,45,60); */
    void Gproject_to_plane(const Matrix1D<double> &gr,
                           double rot, double tilt, double psi, Matrix1D<double> &result) const
    {
        grid2universe(gr, result);
        Uproject_to_plane(result, rot, tilt, psi, result);
    }

    /** Project a grid vector onto a plane (Euler matrix).
        This function does the same as the preceeding one, but it accepts
        an Euler matrix, this might help you to save time when making
        the projection of a whole volume.
        See \ref Uproject_to_plane for more
        information.
        \\ Ex:
        @code
        Matrix2D<double> euler; Euler_angles2matrix(45,45,60,euler);
        Matrix1D<double> up=sg.Gproject_to_plane(vectorR3(1,0,0),euler);
        @endcode */
    void Gproject_to_plane(const Matrix1D<double> &gr,
                           const Matrix2D<double> &euler, Matrix1D<double> &result) const
    {
        grid2universe(gr, result);
        Uproject_to_plane(result, euler, result);
    }

    /** Project a grid direction onto a plane (Euler angles).
        More or less this function goes in the fashion of the preceeding
        ones of projecting a vector. The difference between them is while
        in the first one the grid vector is translated directly to the
        universal coordinate system (it is like a point translation, this
        point in the grid system is this other one in the universal
        system), in this second function only the direction is translated
        into the universal coordinate system. Ie, a point in the grid
        indicates a direction in the grid system, then this same direction
        is expressed in the universal coordinate system. And then this
        direction is projected onto the given plane. See
        \ref Uproject_to_plane for more
        information about the projection process.

        The direction can be seen then as the free vector associated to
        a given vector in the grid. It can also be seen as the translation
        to the universal coordinate system of the grid vector when the
        grid origin is moved to the Universe origin.
        A pseudocode to get the direction of a given grid vector (gr) is the
        following
        @code
        Udirection=grid2universe(gr)-origin;
        @endcode

        This function can be used, for instance, to compute the projections
        of the grid axes onto the projection plane. */
    void Gdir_project_to_plane(const Matrix1D<double> &gr,
                               double rot, double tilt, double psi, Matrix1D<double> &result) const
    {
        grid2universe(gr, result);
        V3_MINUS_V3(result, result, origin);
        Uproject_to_plane(result, rot, tilt, psi, result);
    }

    /** Project a grid direction onto a plane (Euler angles).
        This function does the same as the preceeding one, but it accepts
        an Euler matrix, this might help you to save time when making
        the projection of a whole volume.
        See \ref Uproject_to_plane for more
        information about the projecting process. */
    void Gdir_project_to_plane(const Matrix1D<double> &gr,
                               const Matrix2D<double> &euler, Matrix1D<double> &result) const
    {
        grid2universe(gr, result);
        V3_MINUS_V3(result, result, origin);
        Uproject_to_plane(result, euler, result);
    }
};

/*****************************************************************************/
/* Complex Grids                                                             */
/*****************************************************************************/
/** Complex Grids.
    Grid is the structure where "true" grids like (BCC, FCC, CC, ...) are
    stored. A complex grid is nothing else than a list of Simple Grids, and
    this is so because a complex grid is supposed to be the superposition
    of several simple grids, each one with all the information that a
    \ref SimpleGrid has got.
    When projecting the volume or whatever what you must do is to project
    each simple grid, or apply the function you want to each simple grid.
    For instance, here you have a function to show a complex grid
    @code
       std::ostream& operator << (std::ostream& o, Grid &grid) {
          o << "Complex Grid -------------------------------------\n";
          for (int i=0; i<grid.GridsNo(); i++) o << grid(i);
          return o;
       }
    @endcode

    Initially the complex grid is empty, and you must fill it adding new
    simple grids to it with the function \ref Grid::add_grid . The simple
    grids must be prepared to work before entering into the complex one.
*/
class Grid
{
    /* Structure --------------------------------------------------------------- */
    std::vector<SimpleGrid>   LG;                 // List of grids
public:
    /* Protoypes --------------------------------------------------------------- */
    /** Add a grid to the set.
        The complex grid is a list of simple grids, use this function to
        add a simple grid to the complex one, remember that before using
        this function the simple grid must have been prepared to work
        (see \ref SimpleGrid::prepare_grid ). See class documentation for
        an example of use. */
    void add_grid(const SimpleGrid &SG)
    {
        LG.push_back(SG);
    }

    /** Get a constant "pointer" to a simple grid.
        The complex grid is built upon the STL vectors, using this function
        you can access (constant access, you cannot modify the simple grid)
        any of the simple grids inside the list. Remember that the first
        grid is the number 0. An exception is thrown if you try to access
        a simple grid beyond the number of actual simple grids inside the
        complex one.
        \\ Ex: std::cout << "The first grid in the BCC grid is " << BCC(0); */
    const SimpleGrid & operator()(size_t n) const
    {
        if (n>LG.size())
            REPORT_ERROR(ERR_VALUE_INCORRECT, "The Grid hasn't got so many Simple Grids");
        return LG[n];
    }

    /** Get a "pointer" to a simple grid.
        The complex grid is built upon the STL vectors, using this function
        you can access
        any of the simple grids inside the list. Remember that the first
        grid is the number 0. An exception is thrown if you try to access
        a simple grid beyond the number of actual simple grids inside the
        complex one.
        \\ Ex: BCC(0).origin=vectorR3(1,1,1); */
    SimpleGrid& operator()(size_t n)
    {
        if (n>LG.size())
            REPORT_ERROR(ERR_VALUE_INCORRECT, "The Grid hasn't got so many Simple Grids");
        return LG[n];
    }

    /** Another function for get a "pointer" to a simple grid.*/
    void get_SimpleGrid(int n, SimpleGrid& G) const
    {
        G = (*this)(n);
    }

    /** Number of simple grids inside.
        This function returns the number of simple grids inside the complex
        grid.
        \\ Ex: std::cout << "In BCC there are " << BCC.GridsNo() << " grids\n"; */
    size_t GridsNo() const
    {
        return LG.size();
    }

    /** Show a complex grid.
        Show all simple grids inside the complex one.
        \\ Ex: std::cout << BCC; */
    friend std::ostream& operator << (std::ostream& o, const Grid &grid)
    {
        o << "Complex Grid -------------------------------------\n";
        for (size_t i = 0; i < grid.GridsNo(); i++)
            o << grid(i);
        return o;
    }

    /** Assignment.
        \\ Ex: Grid BCC2=BCC1; */
    Grid& operator = (const Grid &G)
    {
        if (&G != this)
            LG = G.LG;
        return *this;
    }

    /** Another function for assignment.*/
    void assign(const Grid &G)
    {
        *this = G;
    }

    /** Clear. */
    void clear()
    {
        LG.clear();
    }

    /** Minimum size for a voxel volume if it is to hold the whole grid.
        The returned vectors Gcorner1 and Gcorner2 enclose this grid.
        You can supply a deformation matrix (see \ref blobs2voxels )*/
    void voxel_corners(Matrix1D<double> &Gcorner1, Matrix1D<double> &Gcorner2,
                       const Matrix2D<double> *V = NULL) const;
};

/*****************************************************************************/
/* Some useful Grids                                                         */
/*****************************************************************************/
/**@name UsefulGrids Some useful grids
    These are already-built grids to make your task easier. You can
    find here Simple Cubic (CC), Face-Centered (FCC) and Body-Centered
    Cubic (BCC) grids. Most of them are supposed to have its origin in
    the middle of the two defining corners. The defining corners
    (corner1 and corner2) are directly related with lowest and highest
    indexes inside the grid respectively. The difference is that the
    corners are given ALWAYS in the universal coordinate system while
    the lowest and highest are vectors in the grid coordinate system */
//@{
/// CC identifier
#define  CC 0
/// FCC identifier
#define FCC 1
/// BCC identifier
#define BCC 2

/** Create a CC grid as a Simple grid.
    The parameters of the CC grid are:
    @code
    axes:            (1,0,0),(0,1,0),(0,0,1)
    lowest:          FLOORnD(grid.universe2grid(corner1))
    highest:         CEILnD (grid.universe2grid(corner2))
    relative_size:   specified
    origin:          specified
    @endcode
    where the specified means that it has been specified as a parameter of
    the function. The resulting grid is ready to work.
    There are two roundings for the lowest and highest indexes, they are
    such that corner1 and corner2 are guaranteed to be inside the Cubic
    grid. It might be that due to the relative size you are including more
    points in the universal grid that needed. Look at the following example,
    imagine that you have a relative size of 2.75 and that you want a cubic
    grid between (-1,-1,-1) and (1,1,1) (in the universal coord. system),
    with origin at (0,0,0); then the smallest cubic grid that includes the
    corners is also defined between (-1,-1,-1) and (1,1,1) in the grid
    system, but in the universal system these points corresponds to
    (-2.75,-2.75,-2.75) and (2.75,2.75,2.75). You see that you are including
    points (specifically, (-2,-2,-2), (-2,-2,-1), ...) that you didn't
    intend at the beginning, but remember that with this relative size,
    the resulting CC grid is the smallest one which includes the given corners.
*/
SimpleGrid Create_CC_grid(double relative_size,
                          const Matrix1D<double> &corner1,
                          const Matrix1D<double> &corner2,
                          const Matrix1D<double> &origin);

/** Create a CC grid as a Complex grid (two corners).
    This function makes a call to the previous one with
    @code
    corner1:         specified
    corner2:         specified
    relative_size:   specified
    origin:          ROUNDnD((corner1+corner2)/2)
    @endcode

    That is to say, the origin of the grid is at its center, lowest and
    highest indexes are nearly symmetrical (lowest=-highest). The possible
    asymetries come from the ROUNDnD in the origin calculation.
*/
Grid Create_CC_grid(double relative_size,
                    const Matrix1D<double> &corner1,
                    const Matrix1D<double> &corner2);

/** Create a CC grid as a Complex grid (size).
    This function makes a call to the previous one with
    @code
    corner1:         -origin
    corner2:         vectorR3(Xdim,Ydim,Zdim)-origin-1
    relative_size:   specified
    origin:          vectorR3((int)(Xdim/2.0),(int)(Ydim/2.0),(int)(Zdim/2.0))
    @endcode

    That is to say, the origin of the grid is at its center, corner1
    and corner2 are chosen such that from corner1 to corner2 (both included)
    there are exactly (Xdim,Ydim,Zdim) samples in the universal grid. The
    resulting volume is symmetrical (corner1=-corner2) if the sizes are odd
    in all directions, and is not if the sizes are even.
*/
Grid Create_CC_grid(double relative_size,
                    int Zdim, int Ydim, int Xdim);

/** Create a BCC grid as a Complex grid (two corners).
    This function constructs two CC Simple grids with the following
    parameters
    @code
    Grid(0) -------------------------------------
    corner1:         specified
    corner2:         specified
    relative_size:   specified
    origin:          ROUNDnD((corner1+corner2)/2)

    Grid(1) -------------------------------------
    corner1:         specified
    corner2:         specified
    relative_size:   specified
    origin:          ROUNDnD((corner1+corner2)/2)+relative_size/2*vectorR3(1,1,1)
    @endcode

    As you see the BCC grid is the superposition of 2 CC grids one shifted
    with the other by (0.5,0.5,0.5) units in the grid coordinate system.
*/
Grid Create_BCC_grid(double relative_size,
                     const Matrix1D<double> &corner1, const Matrix1D<double> &corner2);

/** Create a FCC grid as a Complex grid (two corners).
    This function constructs four CC Simple grids with the following
    parameters

    @code
    origin:          ROUNDnD((corner1+corner2)/2)

    Grid(0) -------------------------------------
    corner1:         specified
    corner2:         specified
    relative_size:   specified
    origin:          origin

    Grid(1) -------------------------------------
    corner1:         specified
    corner2:         specified
    relative_size:   specified
    origin:          origin+relative_size/2*vectorR3(0,1,1)

    Grid(2) -------------------------------------
    corner1:         specified
    corner2:         specified
    relative_size:   specified
    origin:          origin+relative_size/2*vectorR3(1,0,1)

    Grid(3) -------------------------------------
    corner1:         specified
    corner2:         specified
    relative_size:   specified
    origin:          origin+relative_size/2*vectorR3(1,1,0)
    @endcode
*/
Grid Create_FCC_grid(double relative_size,
                     const Matrix1D<double> &corner1, const Matrix1D<double> &corner2);

/** Create a simple grid fitting the given sphere.
    Not all index combinations are within the sphere, so attention
    must be  paid whether a certain point is within it or not.

    The formula for the grid limits are computed as
    @code
       XX(highest)=CEIL((R/relative_size)/X.module());
       XX(lowest)=-XX(highest);
    @endcode*/
SimpleGrid Create_grid_within_sphere(double relative_size,
                                     const Matrix1D<double> &origin,
                                     const Matrix1D<double> &X, const Matrix1D<double> &Y,
                                     const Matrix1D<double> &Z, double R2);

/** Create a CC grid such that a sphere of radius R and centered
    at the origin is inside.
    @code
    origin=(0,0,0);
    x=(1,0,0);
    y=(0,1,0);
    z=(0,0,1);
    @endcode*/
Grid Create_CC_grid(double relative_size, double R);

/** Create a BCC grid such that a sphere of radius R and centered
    at the origin is inside.
    @code
    origin=(0,0,0);
    x=(0.5,0.5,-0.5);
    y=(0.5,-0.5,0.5);
    z=(-0.5,0.5,0.5);
    @endcode*/
Grid Create_BCC_grid(double relative_size, double R);

/** Create a FCC grid such that a sphere of radius R and centered
    at the origin is inside.
    @code
    origin=(0,0,0);
    x=(0.5,0.5,0);
    y=(0.5,0,0.5);
    z=(0,0.5,0.5);
    @endcode*/
Grid Create_FCC_grid(double relative_size, double R);
//@}

/*****************************************************************************/
/* Grid Volumes                                                              */
/*****************************************************************************/
/** Grid Volume class.
    The grid volumes are a special kind of volumes used in reconstructions
    where the volumes are understood as a 3D matrix of coefficients placed
    at the points defined by a \ref Grid (BCC, FCC, ...) Remember that
    a complex grid like these ones are seen as lists of simpler grids. That
    is why this Grid Volumes are compound of a list of volumes. The first
    volume of the grid volume should contain the coefficients for the first
    grid in the list of grids, and so on.
*/
template <class T>
class GridVolumeT
{
    // Structure ---------------------------------------------------------------
    std::vector<Image<T> * > LV;               // List of volumes
    Grid                 G;                 // Grid associated to this volume

public:
    // Protoypes ---------------------------------------------------------------
    /** Empty constructor.
        The volume list is empty, and the grid is also empty */
    GridVolumeT()
    {
        LV.clear();
        (Grid) G;
    }

    /** Copy constructor. */
    GridVolumeT(const GridVolumeT &_RV)
    {
        *this = _RV;
    }

    /** Create using a grid as pattern and basis.
        Create the grid volume using this grid as pattern. This function calls
        \ref adapt_to_grid , so look at there to know exactly what is
        performed */
    GridVolumeT(const Grid &_grid)
    {
        adapt_to_grid(_grid);
    }

    /** Assignment. */
    GridVolumeT& operator = (const GridVolumeT& RV)
    {
        if (this != &RV)
        {
            clear();
            G = RV.G;
            for (size_t i = 0; i < RV.VolumesNo(); i++)
            {
                Image<T>  *V = new Image<T>;
                *V = RV(i);
                LV.push_back(V);
            }
        }
        return *this;
    }

    /** Another function for assignment.*/
    void assign(const GridVolumeT& RV)
    {
        *this = RV;
    }

    /** Destructor. */
    ~GridVolumeT()
    {
        clear();
    }

    /** Reset current volume and use the given grid as pattern.
        The current grid volume is cleared. Then the given grid is taken
        as the grid for this grid volume. Zero-valued volumes are added
        to the volume list according to the number of grids inside the complex
        grid. The size and starting point of the volumes added are
        fixed by the lowest and highest fields of each \ref SimpleGrid . */
    void adapt_to_grid(const Grid &_grid)
    {
        // Clear old list of volumes
        for (size_t i = 0; i < VolumesNo(); i++)
            if (LV[i]!=NULL)
                delete LV[i];
        LV.clear();

        // Keep the incoming Grid at the same time the old grid is forgotten
        G = _grid;

        // Generate a volume for each subgrid
        int                        Zdim, Ydim, Xdim;
        Image<T> *                 Vol_aux;
        for (size_t i = 0; i < G.GridsNo(); i++)
        {
            SimpleGrid & grid = G(i);
            grid.getSize(Zdim, Ydim, Xdim);
            Vol_aux = new Image<T>;
            (*Vol_aux)().resize(Zdim, Ydim, Xdim);  // Using this function
            // after empty creation the volume
            // is zero-valued.
            STARTINGX((*Vol_aux)()) = (int) XX(grid.lowest); // This values are already
            STARTINGY((*Vol_aux)()) = (int) YY(grid.lowest); // integer although they
            STARTINGZ((*Vol_aux)()) = (int) ZZ(grid.lowest); // are stored as float
            LV.push_back(Vol_aux);
        }
    }

    /** Resize.
        The grid volume is resized such that the space delimited by the
        two given corners is covered. Overlapping grid points are retained
        while non-overlapping ones are set to 0. */
    void resize(const Matrix1D<double> &corner1,
                const Matrix1D<double> &corner2)
    {
        Image<T> *         Vol_aux;
        std::vector<Image<T> * > LV_aux;

        for (size_t n = 0; n < G.GridsNo(); n++)
        {
            SimpleGrid &grid = G(n);

            // Resize grid
            grid.universe2grid(corner1, grid.lowest);
            grid.lowest.selfFLOOR();
            grid.universe2grid(corner2, grid.highest);
            grid.highest.selfCEIL();

            // Resize auxiliary volume
            int Zdim, Ydim, Xdim;
            grid.getSize(Zdim, Ydim, Xdim);
            Vol_aux = new Image<T>;
            (*Vol_aux)().resize(Zdim, Ydim, Xdim);
            STARTINGX((*Vol_aux)()) = (int) XX(grid.lowest); // This values are already
            STARTINGY((*Vol_aux)()) = (int) YY(grid.lowest); // integer although they
            STARTINGZ((*Vol_aux)()) = (int) ZZ(grid.lowest); // are stored as float

            // Copy values in common
            Image<T> * origin = LV[n];
            SPEED_UP_tempsInt;
            FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(VOLMATRIX(*Vol_aux), VOLMATRIX(*origin))
            {
                VOLVOXEL(*Vol_aux, k, i, j) = VOLVOXEL(*origin, k, i, j);
            }

            // Extract old volume and push new one
            delete LV[n];
            LV_aux.push_back(Vol_aux);
        }
        LV = LV_aux;
    }

    /** Resize after a pattern.
        The volume is of the same shape, the grids are copied, and the
        volume coeeficients are set to 0. */
    template <class T1>
    void resize(const GridVolumeT<T1> &GV)
    {
        clear();
        for (size_t n = 0; n < GV.VolumesNo(); n++)
        {
            SimpleGrid grid;
            grid = GV.grid(n);
            G.add_grid(grid);

            Image<T> *Vol_aux;
            Vol_aux = new Image<T>;
            (*Vol_aux)().resize(GV(n)());
            LV.push_back(Vol_aux);
        }
    }

    /** Set to zero with the actual size and origin. */
    void initZeros()
    {
        for (size_t i = 0; i < VolumesNo(); i++)
            (*this)(i)().initZeros();
    }

    /** Clear the volume */
    void clear()
    {
        for (size_t i = 0; i < VolumesNo(); i++)
            delete LV[i];
        LV.clear();
        G.clear();
    }

    /** Access to one of the volumes in the list.
        The first volume is the number 0. */
    Image<T> & operator()(size_t n)
    {
        if (n>LV.size())
            REPORT_ERROR(ERR_VALUE_INCORRECT, "The Grid Volume hasn't got so many Simple Volumes");
        return *(LV[n]);
    }

    /** Another function for access to one of the volumes in the list.*/
    void get_volume(size_t n, Image<T> &V)
    {
        V = (*this)(n);
    }

    /** Constant access to a volume in the list.
        The first volume is the number 0. */
    const Image<T> & operator()(size_t n) const
    {
        if (n>LV.size())
            REPORT_ERROR(ERR_VALUE_INCORRECT, "The Grid Volume hasn't got so many Simple Volumes");
        return *(LV[n]);
    }

    /** Constant access to a simple grid.
        The grid is the \ref SimpleGrid  associated to the volume which
        occupies position n in the volume list. */
    const SimpleGrid & grid(size_t n) const
    {
        return G(n);
    }

    /** Get simple grid. */
    void get_SimpleGrid(size_t n, SimpleGrid &G)
    {
        G = grid(n);
    }

    /** Constant access to the whole grid.
        The grid is the whole \ref Grid  associated to the volume. */
    const Grid & grid() const
    {
        return G;
    }

    /** Get grid. */
    void get_Grid(Grid &G)
    {
        G = grid();
    }

    /** Number of volumes inside structure. */
    size_t VolumesNo() const
    {
        return LV.size();
    }

#define GRIDVOLUME_BY_SCALAR(op) \
    GridVolumeT<T> result; \
    result.G = G; \
    result.LV.reserve(VolumesNo()); \
    for (size_t i=0; i<VolumesNo(); i++) \
        array_by_scalar((*this)(i)(),f,result(i)(),op); \
    return result;

    /** Sum a constant.
        The constant is added to all simple volumes.
        \\Ex: V2=V1+6; */
    GridVolumeT<T> operator + (T f) const
    {
        GRIDVOLUME_BY_SCALAR('+');
    }

    /** Substract a constant.
        The constant is substracted from all simple volumes.
        \\Ex: V2=V1-6; */
    GridVolumeT<T> operator - (T f) const
    {
        GRIDVOLUME_BY_SCALAR('-');
    }

    /** Multiply by a constant.
        The constant is multiplied with all simple volumes.
        \\Ex: V2=V1*6; */
    GridVolumeT<T> operator *(T f) const
    {
        GRIDVOLUME_BY_SCALAR('*');
    }

    /** Divide by a constant.
        The constant divides all simple volumes.
        \\Ex: V2=V1/6; */
    GridVolumeT<T> operator / (T f) const
    {
        GRIDVOLUME_BY_SCALAR('/');
    }

    /** Sum a constant.
        The constant is added to all simple volumes.
        \\Ex: V2=6+V1; */
    GridVolumeT<T> friend operator + (T f, const GridVolumeT<T> &GV)
    {
        return GV + f;
    }

    /** Multiply by a constant.
        The constant is multiplied with all simple volumes.
        \\Ex: V2=6*V1; */
    GridVolumeT<T> friend operator *(T f, const GridVolumeT<T> &GV)
    {
        return GV*f;
    }

#define GRIDVOL_BY_GRIDVOL(op) \
    GridVolumeT<T> result; \
    Image<T> * Vol_aux; \
    \
    if (VolumesNo()!=GV.VolumesNo()) \
        REPORT_ERROR(ERR_GRID_SIZE,(std::string)"GridVolume::"+op+": Different number of subvolumes");\
    \
    result.G = G;\
    result.LV.reserve(VolumesNo());\
    \
    for (size_t i=0; i<VolumesNo(); i++) { \
        try { \
            Vol_aux = new Image<T>; \
            arrayByArray((*this)(i)(),GV(i)(),(*Vol_aux)(),op); \
            result.LV.push_back(Vol_aux); \
        } catch (XmippError XE) {\
            std::cout << XE; \
            REPORT_ERROR(ERR_GRID_SIZE,(std::string)"GridVolume::"+op+": Different shape of volume " +\
                         integerToString(i)); \
        } \
    } \
    \
    return result;

    /** Sum another volume.
        The two volumes must be equally the same in shape (size,
        origin, number of simple volumes, ...) if they aren't an
        exception is thrown.
        \\Ex: V3=V1+V2; */
    GridVolumeT<T> operator + (const GridVolumeT<T> &GV)
    {
        GRIDVOL_BY_GRIDVOL('+');
    }

    /** Substract another volume.
        The two volumes must be equally the same in shape (size,
        origin, number of simple volumes, ...) if they aren't an
        exception is thrown.
        \\Ex: V3=V1-V2; */
    GridVolumeT<T> operator - (const GridVolumeT<T> &GV)
    {
        GRIDVOL_BY_GRIDVOL('-');
    }

    /** Multiply by another volume.
        The two volumes must be equally the same in shape (size,
        origin, number of simple volumes, ...) if they aren't an
        exception is thrown.
        \\Ex: V3=V1*V2; */
    GridVolumeT<T> operator *(const GridVolumeT<T> &GV)
    {
        GRIDVOL_BY_GRIDVOL('*');
    }

    /** Divide by another volume.
        The two volumes must be equally the same in shape (size,
        origin, number of simple volumes, ...) if they aren't an
        exception is thrown.
        \\Ex: V3=V1/V2; */
    GridVolumeT<T> operator / (const GridVolumeT<T> &GV)
    {
        GRIDVOL_BY_GRIDVOL('/');
    }

#define GRIDVOL_BY_GRIDVOLASSIG(op) \
    if (VolumesNo()!=GV.VolumesNo()) \
        REPORT_ERROR(ERR_GRID_SIZE,(std::string)"GridVolume::"+op+"=: Different number of subvolumes");\
    \
    for (size_t i=0; i<VolumesNo(); i++) { \
        try { \
            arrayByArray((*this)(i)(),GV(i)(),(*this)(i)(),op); \
        } catch (XmippError XE) {\
            std::cout << XE; \
            REPORT_ERROR(ERR_GRID_SIZE,(std::string)"GridVolume::"+op+"=: Different shape of volume " +\
                         integerToString(i)); \
        } \
    }

    /** Sum another volume. */
    void operator += (const GridVolumeT<T> &GV)
    {
        GRIDVOL_BY_GRIDVOLASSIG('+');
    }

    /** Substract another volume.*/
    void operator -= (const GridVolumeT<T> &GV)
    {
        GRIDVOL_BY_GRIDVOLASSIG('-');
    }

    /** Multiply by another volume.*/
    void operator *= (const GridVolumeT<T> &GV)
    {
        GRIDVOL_BY_GRIDVOLASSIG('*');
    }

    /** Divide by another volume.*/
    void operator /= (const GridVolumeT<T> &GV)
    {
        GRIDVOL_BY_GRIDVOLASSIG('/');
    }

    /** @name Grids I/O
       The reconstructing volume is stored in the
       Xmipp format. Although it
       is not a true volume some special structure is used to store each
       subgrid and subvolume.
       - a first slice contains control information about the subvolume
             which is behind
       - the subvolume is stored after its control information and it can
             be represented by a volume viewer.
       The structure inside the control slice is the following:
       @code
       * Subgrid information: (19 doubles)
            - basis               9 doubles
            - lowest              3 doubles
            - highest             3 doubles
            - relative_size       1 double
            - origin              3 doubles
       * Subvolume information (6 doubles)
            - Zdim, Ydim, Xdim    3 doubles
            - Zinit, Yinit, Xinit 3 doubles
       @endcode

       With respect to the volume structure, if the volume is smaller
       (Ydim,Xdim) than the size assigned to the file then the data is right
       justified and the extra data is filled with zeros.
    */
    //@{
    /** Write grid volume.
        The Filename is compulsory, an exception is thrown if the volume
        is too small to hold the control information of each layer. */
    void write(const FileName &fn) const
    {
        Image<T>    V;
        float temp_float;
        size_t floatsize;
        const std::type_info &typeinfoT = typeid(T); // We need to know what kind
        // of variable is T
        const std::type_info &typeinfoD = typeid(double);
        const std::type_info &typeinfoI = typeid(int);

        floatsize = (size_t) sizeof(float);

        if (VolumesNo() == 0)
            return;

        // Create the writing volume ............................................
        size_t Zdim = 0, Ydim = 0, Xdim = 0;
        for (size_t v = 0; v < VolumesNo(); v++)
        {
            const Image<T> & this_vol = (*this)(v);
            Zdim += ZSIZE(this_vol());
            Ydim = XMIPP_MAX(Ydim, YSIZE(this_vol()));
            Xdim = XMIPP_MAX(Xdim, XSIZE(this_vol()));
        }

        // Check if there is enough space for the control slice
        if (Xdim*Ydim < 25)
            Ydim = (int) CEIL(25.0f / Xdim);

        // A slice is added for control information for each subvolume
        VOLMATRIX(V).initZeros(Zdim + VolumesNo(), Ydim, Xdim);

        // Write Grid volume ....................................................
#define PACK_DOUBLE(v) \
{jj=pos%Xdim; ii=pos/Xdim; pos++; VOLVOXEL(V,sli,ii,jj)=(T)(v);}
#define PACK_INT(v) \
    {jj=pos%Xdim; ii=pos/Xdim; pos++; \
        temp_float = (float) (v); \
        memcpy( &(VOLVOXEL(V,sli,ii,jj)) , &temp_float, floatsize); \
    }

        int sli = 0;
        for (size_t v = 0; v < VolumesNo(); v++)
        {
            int pos, ii, jj;           // Position inside the control slice
            size_t k, i, j;               // Auxiliar counters

            // Choose grid and volume
            const SimpleGrid & this_grid = grid(v);
            const Image<T> & this_vol = (*this)(v);

            // Store Grid data ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
            pos = 0;

            if (typeinfoT == typeinfoD)
            {
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        PACK_DOUBLE((this_grid.basis)(i, j));
                for (i = 0; i < 3; i++)
                    PACK_DOUBLE((this_grid.lowest)(i));
                for (i = 0; i < 3; i++)
                    PACK_DOUBLE((this_grid.highest)(i));
                PACK_DOUBLE(this_grid.relative_size);
                for (i = 0; i < 3; i++)
                    PACK_DOUBLE((this_grid.origin)(i));
                PACK_DOUBLE(this_grid.R2);

                // Store volume control ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
                PACK_DOUBLE(ZSIZE(this_vol()));
                PACK_DOUBLE(YSIZE(this_vol()));
                PACK_DOUBLE(XSIZE(this_vol()));
                PACK_DOUBLE(STARTINGZ(this_vol()));
                PACK_DOUBLE(STARTINGY(this_vol()));
                PACK_DOUBLE(STARTINGX(this_vol()));
            }
            else if (typeinfoT == typeinfoI)
            {
                // We use a trick to save the grid information in the volume
                // If the following if is true the trick can not be used
                if ((sizeof(float) != sizeof(int)))
                    REPORT_ERROR(ERR_TYPE_INCORRECT,
                                 "GridVolume is integer and (sizeof(float)!= sizeof(int)");

                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        PACK_INT((this_grid.basis)(i, j));
                for (i = 0; i < 3; i++)
                    PACK_INT((this_grid.lowest)(i));
                for (i = 0; i < 3; i++)
                    PACK_INT((this_grid.highest)(i));
                PACK_INT(this_grid.relative_size);
                for (i = 0; i < 3; i++)
                    PACK_INT((this_grid.origin)(i));
                PACK_INT(this_grid.R2);

                // Store volume control ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
                PACK_INT(ZSIZE(this_vol()));
                PACK_INT(YSIZE(this_vol()));
                PACK_INT(XSIZE(this_vol()));
                PACK_INT(STARTINGZ(this_vol()));
                PACK_INT(STARTINGY(this_vol()));
                PACK_INT(STARTINGX(this_vol()));
            }
            else
                REPORT_ERROR(ERR_TYPE_INCORRECT, "GridVolume must be double or int\n");

            sli++;

            // Write the whole volume ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
            for (k = 0; k < ZSIZE(VOLMATRIX(this_vol)); k++)
            {
                for (i = 0; i < YSIZE(VOLMATRIX(this_vol)); i++)
                    for (j = 0; j < XSIZE(VOLMATRIX(this_vol)); j++)
                        DIRECT_VOLVOXEL(V, sli, i, j) = DIRECT_VOLVOXEL(this_vol, k, i, j);
                sli++;
            }
        }
#undef PACK_DOUBLE
#undef PACK_INT

        // Effectively write the volume .........................................
        V.write(fn);
    }

    /** Read grid volume.
        The volume is read from a Xmipp volume with a special structure
        at several slices. */
    void read(const FileName &fn, const std::string &basisName)
    {
        Image<T>       V;
        Image<T>       * sV;
        SimpleGrid     sG;
        size_t         sli = 0;

        float temp_float;
        size_t floatsize;
        const std::type_info &typeinfoT = typeid(T); // We need to know what kind
        // of variable is T
        const std::type_info &typeinfoD = typeid(double);
        const std::type_info &typeinfoI = typeid(int);

        floatsize = (size_t) sizeof(float);
        // We use a trick to save the grid information in the volume
        // If the following if is true the trick can not be used
        if ((typeid(T) == typeid(int)) && (sizeof(float) != sizeof(int)))
            REPORT_ERROR(ERR_TYPE_INCORRECT,"Error: GridVolume is integer and (sizeof(float)!= sizeof(int)");

        // Allocate memory ......................................................
        sG.basis.resize(3, 3);
        sG.lowest.resize(3);
        sG.highest.resize(3);
        sG.origin.resize(3);

        // Read Reconstructing volume from file .................................
        V.read(fn);
        if (basisName=="voxels")
        {
        	V().setXmippOrigin();
        	sV=new Image<double>;
        	*sV=V;
        	LV.push_back(sV);
        	G=Create_CC_grid(1.0,ZSIZE(V()),YSIZE(V()),XSIZE(V()));
        }
        else
        {
#define UNPACK_DOUBLE(v,cast) \
			{jj=pos%VOLMATRIX(V).xdim; ii=pos/VOLMATRIX(V).xdim; pos++; \
			(v)=(cast)VOLVOXEL(V,sli,ii,jj);}
#define UNPACK_INT(v,cast) \
			{jj=pos%VOLMATRIX(V).xdim; ii=pos/VOLMATRIX(V).xdim; pos++; \
			memcpy( &temp_float, &(VOLVOXEL(V,sli,ii,jj)),floatsize);\
			(v)=(cast)temp_float;}

            while (sli < ZSIZE(V()))
            {
                int pos, ii, jj;           // Position inside the control slice
                size_t k, i, j;               // Auxiliary counters
                size_t Zdim, Ydim, Xdim;
                int   Zinit, Yinit, Xinit;

                // Read Grid data ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
                pos = 0;
                if (typeinfoT == typeinfoD)
                {
                    for (i = 0; i < 3; i++)
                        for (j = 0; j < 3; j++)
                            UNPACK_DOUBLE((sG.basis)(i, j), double);
                    for (i = 0; i < 3; i++)
                        UNPACK_DOUBLE((sG.lowest)(i), int);
                    for (i = 0; i < 3; i++)
                        UNPACK_DOUBLE((sG.highest)(i), int);
                    UNPACK_DOUBLE(sG.relative_size, double);
                    for (i = 0; i < 3; i++)
                        UNPACK_DOUBLE((sG.origin)(i), double);
                    UNPACK_DOUBLE(sG.R2, double);
                }
                else if (typeinfoT == typeinfoI)
                {
                    // We use a trick to save the grid information in the volume
                    // If the following if is true the trick can not be used
                    if ((sizeof(float) != sizeof(int)))
                        REPORT_ERROR(ERR_TYPE_INCORRECT,
                                     "GridVolume is integer and (sizeof(float)!= sizeof(int)");

                    for (i = 0; i < 3; i++)
                        for (j = 0; j < 3; j++)
                            UNPACK_INT((sG.basis)(i, j), double);
                    for (i = 0; i < 3; i++)
                        UNPACK_INT((sG.lowest)(i), int);
                    for (i = 0; i < 3; i++)
                        UNPACK_INT((sG.highest)(i), int);
                    UNPACK_INT(sG.relative_size, double);
                    for (i = 0; i < 3; i++)
                        UNPACK_INT((sG.origin)(i), double);
                    UNPACK_INT(sG.R2, double);
                }
                sG.inv_basis = sG.basis.inv();

                // Store Grid in the list of the grid volume
                G.add_grid(sG);

                // Read Volume Control Information ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
                if (typeinfoT == typeinfoD)
                {
                    UNPACK_DOUBLE(Zdim, int);
                    UNPACK_DOUBLE(Ydim, int);
                    UNPACK_DOUBLE(Xdim, int);
                    UNPACK_DOUBLE(Zinit, int);
                    UNPACK_DOUBLE(Yinit, int);
                    UNPACK_DOUBLE(Xinit, int);
                }
                else if (typeinfoT == typeinfoI)
                {
                    UNPACK_INT(Zdim, int);
                    UNPACK_INT(Ydim, int);
                    UNPACK_INT(Xdim, int);
                    UNPACK_INT(Zinit, int);
                    UNPACK_INT(Yinit, int);
                    UNPACK_INT(Xinit, int);
                }

                // Set volume size and origin
                sV = new Image<T>;
                VOLMATRIX(*sV).initZeros(Zdim, Ydim, Xdim);
                STARTINGZ(VOLMATRIX(*sV)) = Zinit;
                STARTINGY(VOLMATRIX(*sV)) = Yinit;
                STARTINGX(VOLMATRIX(*sV)) = Xinit;
#ifdef DEBUG

                std::cout << "The read grid is \n" << sG;
                std::cout << "Volume dimensions: " << Zdim << " x " << Ydim << " x "
                << Xdim << std::endl;
                std::cout << "Volume init: " << Zinit << " x " << Yinit << " x "
                << Xinit << std::endl;
#endif

                sli++;

                // Read volume ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
                for (k = 0; k < ZSIZE(VOLMATRIX(*sV)); k++)
                {
                    for (i = 0; i < YSIZE(VOLMATRIX(*sV)); i++)
                        for (j = 0; j < XSIZE(VOLMATRIX(*sV)); j++)
                        {
#ifdef DEBUG
                            std::cout << "Reading from file position (" << sli << "," << i
                            << "," << j << ") to subvolume position ("
                            << k << "," << i << "," << j << ")\n";
#endif

                            DIRECT_VOLVOXEL(*sV, k, i, j) = DIRECT_VOLVOXEL(V, sli, i, j);
                        }
                    sli++;
                }

                // Store volume in the list
                LV.push_back(sV);
            }
#undef UNPACK_DOUBLE
#undef UNPACK_INT
        }
    }
#undef DEBUG

    /** Show volume.
        \\Ex: std::cout << V; */
    friend std::ostream& operator <<< > (std::ostream &o, const GridVolumeT &GV);
    //@}
};

typedef GridVolumeT<double> GridVolume;

// Show a grid volume ------------------------------------------------------
/** Show a grid. */
template <class T>
std::ostream& operator << (std::ostream &o, const GridVolumeT<T> &GV)
{
    o << "Grid Volume -----------\n";
    o << GV.G;
    o << "Number of volumes= " << GV.VolumesNo() << std::endl;
    for (int i = 0; i < GV.VolumesNo(); i++)
    {
        o << "Volume " << i << "------------" << std::endl;
        o << GV(i)();
    }
    return o;
}
//@}
#endif
