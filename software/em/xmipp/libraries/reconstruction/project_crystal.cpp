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

#include "project_crystal.h"
#include "art_crystal.h"

#include <data/args.h>

/* Empty constructor ======================================================= */
Crystal_Projection_Parameters::Crystal_Projection_Parameters()
{
    crystal_Xdim = 0;
    crystal_Ydim = 0;
    orthogonal = false;
    a.clear();
    b.clear();
    Nshift_avg = 0;
    Nshift_dev = 0;
    disappearing_th = 0;
    DF_shift_bool = false;
    DF_shift.clear();
}

/* Read Crystal Projection Parameters ====================================== */
void Crystal_Projection_Parameters::read(const FileName &fn, double scale)
{
    MetaData MD;
    size_t objId;
    FileName fn_crystal = fn.removeBlockName();
    MD.read(formatString("crystal@%s", fn_crystal.c_str()));
    if (MD.isEmpty())
        REPORT_ERROR(ERR_IO_NOTOPEN,
                     (String)"Prog_Project_Parameters::read: There is a problem "
                     "opening the file " + fn_crystal);

    std::vector <double> ParamVec;
    objId = MD.firstObject();
    MD.getValue(MDL_DIMENSIONS_2D, ParamVec, objId);
    crystal_Xdim = (int)ParamVec[0];
    crystal_Ydim = (int)ParamVec[1];
    crystal_Xdim = ROUND(scale * crystal_Xdim);
    crystal_Ydim = ROUND(scale * crystal_Ydim);
    MD.getValue(MDL_CRYSTAL_LATTICE_A, ParamVec, objId);
    a.resize(3);
    XX(a) = scale * ParamVec[0];
    YY(a) = scale * ParamVec[1];
    ZZ(a) = 0;
    MD.getValue(MDL_CRYSTAL_LATTICE_B, ParamVec, objId);
    b.resize(3);
    XX(b) = scale * ParamVec[0];
    YY(b) = scale * ParamVec[1];
    ZZ(b) = 0;
    MD.getValue(MDL_CRYSTAL_NOISE_SHIFT,ParamVec, objId);
    Nshift_dev = scale *  ParamVec[0];
    if (ParamVec.size() < 2)
        Nshift_avg = 0;
    else
        Nshift_avg = scale * ParamVec[1];
    MD.getValue(MDL_CRYSTAL_DISAPPEAR_THRE,disappearing_th, objId);
    MD.getValue(MDL_CRYSTAL_ORTHO_PRJ,orthogonal, objId);
    if (MD.getValue(MDL_CRYSTAL_SHFILE, fn_shift, objId))
        DF_shift_bool = true;
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "crystal_Xdim"   << crystal_Xdim<<std::endl;
    std::cerr << "crystal_Ydim"   << crystal_Ydim<<std::endl;
    std::cerr << "a vector"  << a <<std::endl;
    std::cerr << "b vector"  << b <<std::endl;
    std::cerr << "Nshift_devr"  << Nshift_dev <<std::endl;
    std::cerr << "Nshift_avg"  << Nshift_avg <<std::endl;
#endif
#undef DEBUG
}

/* Write =================================================================== */
void Crystal_Projection_Parameters::write(const FileName &fn_crystal)
{
    if (fn_crystal.isMetaData())
    {
        MetaData MD1;  //MetaData for crystal projection parameters
        std::vector<double> FCVect(2);  //For the center of feature
        size_t id;
        std::vector<double> TVector;  //For general use
        MD1.setColumnFormat(false);
        id = MD1.addObject();
        TVector[0] = crystal_Xdim;
        TVector[1] = crystal_Ydim;
        MD1.setValue(MDL_DIMENSIONS_2D, TVector, id);
        TVector[0] = XX(a);
        TVector[1] = YY(a);
        MD1.setValue(MDL_CRYSTAL_LATTICE_A, TVector, id);
        TVector[0] = XX(b);
        TVector[1] = YY(b);
        MD1.setValue(MDL_CRYSTAL_LATTICE_B, TVector, id);
        TVector[0] = Nshift_dev;
        TVector[1] = Nshift_avg;
        MD1.setValue(MDL_NOISE_PIXEL_LEVEL, TVector, id);
        MD1.setValue(MDL_CRYSTAL_DISAPPEAR_THRE, disappearing_th, id);
        MD1.setValue(MDL_CRYSTAL_ORTHO_PRJ, orthogonal, id);
        MD1.setValue(MDL_CRYSTAL_SHFILE, fn_shift.c_str(), id);
        //MD1.write(fn_crystal, MD_OVERWRITE);
    }
    else
    {
        FILE *fh_param;

        if ((fh_param = fopen(fn_crystal.c_str(), "w")) == NULL)
            REPORT_ERROR(ERR_IO_NOTOPEN,
                         (std::string)"Prog_Project_Parameters::write: There is a problem "
                         "opening the file " + fn_crystal + " for output");

        fprintf(fh_param, "# Crystal dimensions (X, Y)\n");
        fprintf(fh_param, "%d %d\n", crystal_Xdim, crystal_Ydim);
        fprintf(fh_param, "# Crystal a lattice vector (X, Y)\n");
        fprintf(fh_param, "%f %f\n", XX(a), YY(a));
        fprintf(fh_param, "# Crystal b lattice vector (X, Y)\n");
        fprintf(fh_param, "%f %f\n", XX(b), YY(b));

        fprintf(fh_param, "#     Noise (and bias) applied to the magnitude shift\n");
        fprintf(fh_param, "%f ", Nshift_dev);
        if (Nshift_avg != 0)
            fprintf(fh_param, "%f \n", Nshift_avg);
        else
            fprintf(fh_param, "\n");

        fprintf(fh_param, "# Disappearing threshold\n");
        fprintf(fh_param, "%f\n", disappearing_th);

        fprintf(fh_param, "# Orthogonal Projections\n");
        if (orthogonal)
            fprintf(fh_param, "Yes\n");
        else
            fprintf(fh_param, "No\n");

        //   fprintf(fh_param,"# Grid relative size\n");

        fprintf(fh_param, "# File with shifts for each unit cell\n");
        fprintf(fh_param, "%s", fn_shift.c_str());

        fclose(fh_param);
    }
}

/* Project crystal --------------------------------------------------------- */
//#define DEBUG
//#define DEBUG_MORE
void project_crystal(Phantom &phantom, Projection &P,
                     const ParametersProjection &prm,
                     PROJECT_Side_Info &side, const Crystal_Projection_Parameters &prm_crystal,
                     float rot, float tilt, float psi)
{
    SPEED_UP_tempsDouble;
    //Scale crystal output
    //   int aux_crystal_Ydim = ROUND(prm_crystal.crystal_Ydim
    //         *phantom.phantom_scale);
    //   int aux_crystal_Xdim = ROUND(prm_crystal.crystal_Xdim
    //                                       *phantom.phantom_scale);

    // Initialize whole crystal projection
    P().initZeros(prm_crystal.crystal_Ydim, prm_crystal.crystal_Xdim);
    P().setXmippOrigin();

    // Compute lattice vectors in the projection plane
    P.setAngles(rot, tilt, psi);
    Matrix1D<double> proja = P.euler * prm_crystal.a;
    Matrix1D<double> projb = P.euler * prm_crystal.b;

    // Check if orthogonal projections
    // (projXdim,0)'=A*aproj
    // (0,projYdim)'=A*bproj
    Matrix2D<double> Ainv, A, D, Dinv, AuxMat, Daux;
    if (prm_crystal.orthogonal)
    {
        A.resize(2, 2);
        //        A(0, 0) = YY(projb) * XSIZE(P());
        //        A(0, 1) = -XX(projb) * XSIZE(P());
        //        A(1, 0) = -YY(proja) * YSIZE(P());
        //        A(1, 1) = XX(proja) * YSIZE(P());
        A(0, 0) =  YY(projb) * prm.proj_Xdim;
        A(0, 1) = -XX(projb) * prm.proj_Xdim;
        A(1, 0) = -YY(proja) * prm.proj_Ydim;
        A(1, 1) =  XX(proja) * prm.proj_Ydim;
        double nor = 1 / (XX(proja) * YY(projb) - XX(projb) * YY(proja));
        M2x2_BY_CT(A, A, nor);
        A.resize(3, 3);
        A(2, 2) = 1;
        A.inv(Ainv);
        AuxMat.resize(3, 3);
        //Matrix with crystal vectors
        // Delta = (a, b) with a and b crystal column vectors
        AuxMat(0, 0) = XX(prm_crystal.a);
        AuxMat(0, 1) = YY(prm_crystal.a);
        AuxMat(0, 2) = 0.;
        AuxMat(1, 0) = XX(prm_crystal.b);
        AuxMat(1, 1) = YY(prm_crystal.b);
        AuxMat(1, 2) = 0.;
        AuxMat(2, 0) = 0.;
        AuxMat(2, 1) = 0.;
        AuxMat(2, 2) = 1.;
        // Product of Delta(trasposed) and Delta*
        D.resize(3, 3);
        D(0, 0) = XSIZE(P());
        D(0, 1) = 0.;
        D(0, 2) = 0.;
        D(1, 0) = 0.;
        D(1, 1) = YSIZE(P());
        D(1, 2) = 0.;
        D(2, 0) = 0.;
        D(2, 1) = 0.;
        D(2, 2) = 1.;
        Daux = D;
        Daux.inv(D);
        M3x3_BY_M3x3(D, AuxMat, D);
        Dinv.resize(3, 3);
        D.inv(Dinv);
    }
    else
    {
        A.initIdentity(3);
        Ainv.initIdentity(3);
    }

//#define DEBUG
#ifdef DEBUG
    std::cout << "P shape ";
    P().printShape();
    std::cout << std::endl;
    std::cout << "P.euler " << P.euler;
    std::cout << std::endl;
    std::cout << "rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;
    std::cout << "a " << prm_crystal.a.transpose() << std::endl;
    std::cout << "b " << prm_crystal.b.transpose() << std::endl;
    std::cout << "proja " << proja.transpose() << std::endl;
    std::cout << "projb " << projb.transpose() << std::endl;
    std::cout << "proj_Xdim " << prm.proj_Xdim << " proj_Ydim " << prm.proj_Ydim
    << std::endl;
    std::cout << "A\n" << A << std::endl;
    std::cout << "Ainv\n" << Ainv << std::endl;
    std::cout << "D\n" << D << std::endl;
    std::cout << "Dinv\n" << Dinv << std::endl;
#endif

    // Compute aproj and bproj in the deformed projection space
    Matrix1D<double> aprojd = A * proja;
    Matrix1D<double> bprojd = A * projb;
#ifdef DEBUG

    std::cout << "aprojd " << aprojd.transpose() << std::endl;
    std::cout << "bprojd " << bprojd.transpose() << std::endl;
#endif

    // Get rid of all unnecessary components
    Matrix2D<double> A2D = A;
    A2D.resize(2, 2);
    proja.resize(2);
    projb.resize(2);
    aprojd.resize(2);
    bprojd.resize(2);

    // The unit cell projection size as well
    Matrix1D<double> corner1(2), corner2(2);
    // First in the compressed space
    XX(corner1) = FIRST_XMIPP_INDEX(prm.proj_Xdim);
    YY(corner1) = FIRST_XMIPP_INDEX(prm.proj_Ydim);
    XX(corner2) = LAST_XMIPP_INDEX(prm.proj_Xdim);
    YY(corner2) = LAST_XMIPP_INDEX(prm.proj_Ydim);
    // Now deform
#ifdef DEBUG

    std::cout << "corner1 before deformation " << corner1.transpose() << std::endl;
    std::cout << "corner2 before deformation " << corner2.transpose() << std::endl;
#endif

    rectangle_enclosing(corner1, corner2, A2D, corner1, corner2);
#ifdef DEBUG

    std::cout << "corner1 after deformation " << corner1.transpose() << std::endl;
    std::cout << "corner2 after deformation " << corner2.transpose() << std::endl;
#endif

    MultidimArray<double> cell_shiftX, cell_shiftY, cell_shiftZ;
    MultidimArray<int>    cell_inside;
    MultidimArray<double> exp_shifts_matrix_X;
    MultidimArray<double> exp_shifts_matrix_Y;
    MultidimArray<double> exp_shifts_matrix_Z;
    MultidimArray<double> exp_normal_shifts_matrix_X;
    MultidimArray<double> exp_normal_shifts_matrix_Y;
    MultidimArray<double> exp_normal_shifts_matrix_Z;

    fill_cell_positions(P, proja, projb, aprojd, bprojd, corner1, corner2,
                        prm_crystal, cell_shiftX, cell_shiftY, cell_shiftZ, cell_inside,
                        exp_shifts_matrix_X, exp_shifts_matrix_Y, exp_shifts_matrix_Z);

    // Fill a table with all exp shifts
    init_shift_matrix(prm_crystal, cell_inside, exp_shifts_matrix_X,
                      exp_shifts_matrix_Y,
                      exp_shifts_matrix_Z,
                      exp_normal_shifts_matrix_X,
                      exp_normal_shifts_matrix_Y,
                      exp_normal_shifts_matrix_Z,
                      phantom.phantom_scale);

    //std::cout << "prm_crystal : " << prm_crystal. << std::endl;
    //#define DEBUG
#ifdef DEBUG
    std::cout << "prm_crystal ";
    std::cout << std::endl;
    std::cout << "prm_crystal" << prm_crystal.DF_shift << std::endl;
    std::cout << "prm_crystal" << prm_crystal.DF_shift_bool << std::endl;

#endif
#undef DEBUG

    // Prepare matrices to go from uncompressed space to deformed projection
    Matrix2D<double> AE = A * P.euler;   // From uncompressed to deformed
    Matrix2D<double> AEinv = AE.inv(); // From deformed to uncompressed
    // add the shifts to the already compute values
    Matrix1D<double> temp_vect(3);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(exp_shifts_matrix_X)
    {
        //these experimental shift are in phantom
        //coordinates, not into the projection
        //plane, so before adding them we need to project
        temp_vect = AE * vectorR3(exp_shifts_matrix_X(i, j),
                                  exp_shifts_matrix_Y(i, j),
                                  exp_shifts_matrix_Z(i, j));
        //#define DEBUG5
#ifdef DEBUG5

        if (i > 0)
            exp_shifts_matrix_Z(i, j) = 65;
        else
            exp_shifts_matrix_Z(i, j) = 0;
        temp_vect = AE * vectorR3(exp_shifts_matrix_X(i, j),
                                  exp_shifts_matrix_Y(i, j),
                                  exp_shifts_matrix_Z(i, j));
#endif
#undef DEBUG5

        // Add experimental shiftsy
        cell_shiftX(i, j) += XX(temp_vect);
        cell_shiftY(i, j) += YY(temp_vect);
        //entiendo que x e y deban estar en el plano de la proyeccion
        //pero Z!!!!!!!!!!!!!!!!!!!
        cell_shiftZ(i, j) += ZZ(temp_vect);
    }
    //#define DEBUG
#ifdef DEBUG
    std::cout << "Cell inside shape ";
    cell_inside.printShape();
    std::cout << std::endl;
    std::cout << "Cell inside\n" << cell_inside << std::endl;
    std::cout << "Cell shiftX\n" << cell_shiftX << std::endl;
    std::cout << "Cell shiftY\n" << cell_shiftY << std::endl;
    std::cout << "Cell shiftZ\n" << cell_shiftZ << std::endl;
#endif
    #undef DEBUG

    double density_factor = 1.0;
    if (prm_crystal.orthogonal)
    {
        // Remember to compute de density factor
        Matrix1D<double> projection_direction(3);
        (P.euler).getCol(2, projection_direction);
        projection_direction.selfTranspose();
        density_factor = (projection_direction * Dinv).module();
#ifdef DEBUG

        std::cout << "projection_direction" << projection_direction << std::endl;
        std::cout << "projection_direction*A" << projection_direction*A << std::endl;
#endif

    }
#ifdef DEBUG
    std::cout << "X proyectado=" << (AE*vectorR3(1.0, 0.0, 0.0)).transpose() << std::endl;
    std::cout << "Y proyectado=" << (AE*vectorR3(0.0, 1.0, 0.0)).transpose() << std::endl;
    std::cout << "P.euler_shape=" << std::endl;
    std::cout << P.euler << std::endl;
    std::cout << "P.euler="      << P.euler << std::endl;
    std::cout << "AE="           << AE << std::endl;
    std::cout << "AEinv="        << AEinv << std::endl;
    std::cout << "Ainv="            << Ainv << std::endl;
    std::cout << "density_factor="      << density_factor << std::endl;
#endif

    // Project all cells
    FOR_ALL_ELEMENTS_IN_ARRAY2D(cell_inside)
    //      if (cell_inside(i,j) && rnd_unif(0,1)<prm_crystal.disappearing_th) {
    if (cell_inside(i, j))
    {
        // Shift the phantom
        // Remind that displacements are defined in the deformed projection
        // that is why they have to be translated to the Universal
        // coordinate system
        Matrix1D<double> cell_shift(3);
        VECTOR_R3(cell_shift, cell_shiftX(i, j), cell_shiftY(i, j), cell_shiftZ(i, j));
        //SHIFT still pending
        //cell_shift = cell_shift*phantom.phantom_scale;
#ifdef DEBUG

        std::cout << "cell_shift on deformed projection plane "
        << cell_shift.transpose() << std::endl;
#endif

        M3x3_BY_V3x1(cell_shift, AEinv, cell_shift);
#ifdef DEBUG

        std::cout << "cell_shift on real space "
        << cell_shift.transpose() << std::endl;
#endif

        Phantom aux;
        aux = phantom;
        // Now the cell shift is defined in the uncompressed space
        // Any phantom in the projection line will give the same shift
        // in the projection we are interested in the phantom whose center
        // is in the XYplane


        // the phantom is mapped into a surface of different shapes (parabole,cosine, etc)
        Matrix1D<double> normal_vector(3);
        double alpha, beta, gamma;
        Matrix2D<double> angles_matrix, inverse_angles_matrix;
        Matrix2D<double> def_cyl_angles_matrix;
        Matrix2D<double> cyl_angles_matrix;

        // the phantom is rotated only when there exists a shifts file
        if (prm_crystal.DF_shift_bool)
        {
            // for each (h,k) calculate the normal vector and its corresponding rotation matrix
            normal_vector(0) = exp_normal_shifts_matrix_X(i, j);
            normal_vector(1) = exp_normal_shifts_matrix_Y(i, j);
            normal_vector(2) = exp_normal_shifts_matrix_Z(i, j);

            Euler_direction2angles(normal_vector, alpha, beta, gamma);
            gamma = -alpha;
            Euler_angles2matrix(alpha, beta, gamma, angles_matrix);
            inverse_angles_matrix = angles_matrix.inv();

            for (size_t ii = 0; ii < aux.VF.size(); ii++)
                aux.VF[ii]->rotate(angles_matrix);
        }

        cell_shift = cell_shift - ZZ(cell_shift) / ZZ(P.direction) * P.direction;
#ifdef DEBUG

        std::cout << "cell_shift after moving to ground "
        << cell_shift.transpose() << std::endl;
#endif

        aux.shift(XX(cell_shift), YY(cell_shift), ZZ(cell_shift));

        // Project this phantom
        aux.project_to(P, AE, prm_crystal.disappearing_th);
        // Multiply by factor

        P() = P() * density_factor;
#ifdef DEBUG_MORE

        std::cout << "After Projecting ...\n" << aux << std::endl;
        P.write("inter");
        std::cout << "Hit any key\n";
        char c;
        std::cin >> c;
#endif

    }

}
#undef DEBUG
#undef DEBUG_MORE

/* Find crystal limits ----------------------------------------------------- */
#define MIN_MODULE 1e-2
void find_crystal_limits(
    const Matrix1D<double> &proj_corner1, const Matrix1D<double> &proj_corner2,
    const Matrix1D<double> &cell_corner1, const Matrix1D<double> &cell_corner2,
    const Matrix1D<double> &a, const Matrix1D<double> &b,
    int &iamin, int &iamax, int &ibmin, int &ibmax)
{
    if (a.module() < MIN_MODULE || b.module() < MIN_MODULE)
        REPORT_ERROR(ERR_MATRIX_SIZE, "find_crystal_limits: one of the lattice vectors is "
                     "extremely small");

    // Compute area to cover
    double x0 = XX(proj_corner1) + XX(cell_corner1);
    double y0 = YY(proj_corner1) + YY(cell_corner1);
    double xF = XX(proj_corner2) + XX(cell_corner2);
    double yF = YY(proj_corner2) + YY(cell_corner2);

    // Initialize
    SPEED_UP_temps01;
    Matrix2D<double> A(2, 2), Ainv(2, 2);
    A.setCol(0, a);
    A.setCol(1, b);
    M2x2_INV(Ainv, A);
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "A" << A <<std::endl;
    std::cerr << "Ainv" << Ainv <<std::endl;
#endif
#undef DEBUG

    // Now express each corner in the a,b coordinate system
    Matrix1D<double> r(2);
    VECTOR_R2(r, x0, y0);
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "r" << r <<std::endl;
#endif
#undef DEBUG

    M2x2_BY_V2x1(r, Ainv, r);
    iamin = FLOOR(XX(r));
    iamax = CEIL(XX(r));
    ibmin = FLOOR(YY(r));
    ibmax = CEIL(YY(r));
    //#define DEBUG
#ifdef DEBUG

    std::cerr << "r" << r <<std::endl;
    std::cerr << "Ainv" << Ainv <<std::endl;
    std::cerr << "iamin" << iamin <<std::endl;
    std::cerr << "iamax" << iamax <<std::endl;
    std::cerr << "ibmin" << ibmin <<std::endl;
    std::cerr << "ibmax" << ibmax <<std::endl;
#endif
#undef DEBUG

#define CHANGE_COORDS_AND_CHOOSE_CORNERS2D \
    M2x2_BY_V2x1(r,Ainv,r); \
    iamin=XMIPP_MIN(FLOOR(XX(r)),iamin); iamax=XMIPP_MAX(CEIL(XX(r)),iamax); \
    ibmin=XMIPP_MIN(FLOOR(YY(r)),ibmin); ibmax=XMIPP_MAX(CEIL(YY(r)),ibmax);

    VECTOR_R2(r, x0, yF);
    CHANGE_COORDS_AND_CHOOSE_CORNERS2D;
    VECTOR_R2(r, xF, y0);
    CHANGE_COORDS_AND_CHOOSE_CORNERS2D;
    VECTOR_R2(r, xF, yF);
    CHANGE_COORDS_AND_CHOOSE_CORNERS2D;
}

/* Move following spiral --------------------------------------------------- */
/* Given a cell, we study the cells in the inmidiate surrounding, this gives
   a structure of the like

   xxx
   xXx
   xxx

   We will mean by . a non-visited cell, by x a visited one and by o the
   following cell we should visit  given this structure. To identify
   structures we will use the sum of visited cells in rows and columns. The
   following patterns define the spiral movement:

   .o.0 ...0 ...0 .xx2 .xx2 xx.2 xo.1 .o.0 ...0
   .x.1 ox.1 .xx2 .xx2 .xo1 xxo2 xx.2 xx.2 ox.1
   ...0 .x.1 .ox1 .o.0 ...0 ...0 ...0 xx.2 xx.2 etc
   010  020  012  022  021  220  210  220  120

   This function supposes that r is inside the matrix
*/
void move_following_spiral(Matrix1D<double> &r, const MultidimArray<int> &visited)
{
    int r1 = 0, r2 = 0, r3 = 0, c1 = 0, c2 = 0, c3 = 0;

#define x0 STARTINGX(visited)
#define y0 STARTINGY(visited)
#define xF FINISHINGX(visited)
#define yF FINISHINGY(visited)
#define i  (int)YY(r)
#define j  (int)XX(r)
    // Compute row and column sums
    if (i > y0 && j > x0)
    {
        r1 += visited(i - 1, j - 1);
        c1 += visited(i - 1, j - 1);
    }
    if (i > y0)
    {
        r1 += visited(i - 1, j);
        c2 += visited(i - 1, j);
    }
    if (i > y0 && j < xF)
    {
        r1 += visited(i - 1, j + 1);
        c3 += visited(i - 1, j + 1);
    }
    if (j > x0)
    {
        r2 += visited(i  , j - 1);
        c1 += visited(i  , j - 1);
    }
    r2 += visited(i  , j);
    c2 += visited(i  , j);
    if (j < xF)
    {
        r2 += visited(i  , j + 1);
        c3 += visited(i  , j + 1);
    }
    if (i < yF && j > x0)
    {
        r3 += visited(i + 1, j - 1);
        c1 += visited(i + 1, j - 1);
    }
    if (i < yF)
    {
        r3 += visited(i + 1, j);
        c2 += visited(i + 1, j);
    }
    if (i < yF && j < xF)
    {
        r3 += visited(i + 1, j + 1);
        c3 += visited(i + 1, j + 1);
    }

#ifdef DEBUG
    std::cout << r1 << " " << r2 << " " << r3 << " "
    << c1 << " " << c2 << " " << c3 << std::endl;
#endif

    // Decide where to go
    if (r1 == 0 && r2 == 1 && r3 == 0 && c1 == 0 && c2 == 1 && c3 == 0)
    {
        YY(r)--;
    }
    else if (r1 == 0 && r2 == 1 && r3 == 1 && c1 == 0 && c2 == 2 && c3 == 0)
    {
        XX(r)--;
    }
    else if (r1 == 0 && r2 == 2 && r3 == 1 && c1 == 0 && c2 == 1 && c3 == 2)
    {
        YY(r)++;
    }
    else if (r1 == 2 && r2 == 2 && r3 == 0 && c1 == 0 && c2 == 2 && c3 == 2)
    {
        YY(r)++;
    }
    else if (r1 == 2 && r2 == 1 && r3 == 0 && c1 == 0 && c2 == 2 && c3 == 1)
    {
        XX(r)++;
    }
    else if (r1 == 2 && r2 == 2 && r3 == 0 && c1 == 2 && c2 == 2 && c3 == 0)
    {
        XX(r)++;
    }
    else if (r1 == 1 && r2 == 2 && r3 == 0 && c1 == 2 && c2 == 1 && c3 == 0)
    {
        YY(r)--;
    }
    else if (r1 == 0 && r2 == 2 && r3 == 2 && c1 == 2 && c2 == 2 && c3 == 0)
    {
        YY(r)--;
    }
    else if (r1 == 0 && r2 == 1 && r3 == 2 && c1 == 1 && c2 == 2 && c3 == 0)
    {
        XX(r)--;
    }
    else if (r1 == 1 && r2 == 2 && r3 == 2 && c1 == 3 && c2 == 2 && c3 == 0)
    {
        YY(r)--;
    }
    else if (r1 == 0 && r2 == 2 && r3 == 3 && c1 == 1 && c2 == 2 && c3 == 2)
    {
        XX(r)--;
    }
    else if (r1 == 0 && r2 == 2 && r3 == 2 && c1 == 0 && c2 == 2 && c3 == 2)
    {
        XX(r)--;
    }
    else if (r1 == 2 && r2 == 2 && r3 == 1 && c1 == 0 && c2 == 2 && c3 == 3)
    {
        YY(r)++;
    }
    else if (r1 == 3 && r2 == 2 && r3 == 0 && c1 == 2 && c2 == 2 && c3 == 1)
    {
        XX(r)++;
    }
}
#undef i
#undef j
#undef x0
#undef y0
#undef xF
#undef yF

/* Fill cell rotations and shifts ------------------------------------------ */
//#define DEBUG
void fill_cell_positions(Projection &P,
                         Matrix1D<double> &aproj,   Matrix1D<double> &bproj,
                         Matrix1D<double> &aprojd,  Matrix1D<double> &bprojd,
                         Matrix1D<double> &corner1, Matrix1D<double> &corner2,
                         const Crystal_Projection_Parameters &prm_crystal,
                         MultidimArray<double> &cell_shiftX,
                         MultidimArray<double> &cell_shiftY,
                         MultidimArray<double> &cell_shiftZ,
                         MultidimArray<int>    &cell_inside,
                         MultidimArray<double> &exp_shifts_matrix_X,
                         MultidimArray<double> &exp_shifts_matrix_Y,
                         MultidimArray<double> &exp_shifts_matrix_Z)
{

    // Compute crystal limits
    int iamin, iamax, ibmin, ibmax;
    find_crystal_limits(vectorR2(STARTINGX(P()), STARTINGY(P())),
                        vectorR2(FINISHINGX(P()), FINISHINGY(P())),
                        corner1, corner2, aprojd, bprojd, iamin, iamax, ibmin, ibmax);
    //#define DEBUG
#ifdef DEBUG

    P().printShape();
    std::cout << std::endl;
    std::cout << "aprojd=" << aproj.transpose() << std::endl;
    std::cout << "bprojd=" << bproj.transpose() << std::endl;
    std::cout << iamin << " " << iamax << " " << ibmin << " " << ibmax << std::endl;
#endif

    // Compute weight table in the undeformed space
    MultidimArray<double> weight(3, 3);
    STARTINGX(weight) = -1;
    STARTINGY(weight) = -1;
    weight(0, 0) = 0;
    weight(-1, 0) = weight(1, 0) = 1 / aproj.module();
    weight(0, -1) = weight(0, 1) = 1 / bproj.module();
    weight(-1, 1) = weight(1, -1) = 1 / (aproj - bproj).module();
    weight(-1, -1) = weight(1, 1) = 1 / (aproj + bproj).module();

    // Resize all needed matrices
    cell_shiftX.initZeros(ibmax - ibmin + 1, iamax - iamin + 1);
    STARTINGX(cell_shiftX) = iamin;
    STARTINGY(cell_shiftX) = ibmin;
    cell_shiftY.initZeros(cell_shiftX);
    //in this routine cell_shiftZ is set to zero and nothing else
    cell_shiftZ.initZeros(cell_shiftX);
    cell_inside.initZeros(ibmax - ibmin + 1, iamax - iamin + 1);
    STARTINGX(cell_inside) = iamin;
    STARTINGY(cell_inside) = ibmin;

    // Visited has got one cell more than the rest in all directions
    MultidimArray<int> visited;
    int visited_size = XMIPP_MAX(iamax - iamin + 1, ibmax - ibmin + 1) + 2;
    visited.initZeros(visited_size, visited_size);
    STARTINGX(visited) = iamin - (visited_size - (iamax - iamin + 1) + 1) / 2;
    STARTINGY(visited) = ibmin - (visited_size - (ibmax - ibmin + 1) + 1) / 2;
#ifdef DEBUG

    std::cout << "weight=" << weight;
    std::cout << "cell_shiftX shape ";
    cell_shiftX.printShape();
    std::cout << std::endl;
    std::cout << "cell_inside shape ";
    cell_inside.printShape();
    std::cout << std::endl;
    std::cout << "visited shape ";
    visited.printShape();
    std::cout << std::endl;
#endif

    // Perform the crystal deformation starting in the middle of the
    // crystal and going in spirals until the four corners are visited
    // we find one corner
    Matrix1D<double> r(2), ri(2), sh(2);
    r.initZeros();
#define INDEX(r) (int)YY(r),(int)XX(r)

    while (!visited.isCorner(r))
    {
        visited(INDEX(r)) = true;
#ifdef DEBUG

        std::cout << "   Visiting " << r.transpose() << std::endl;
#endif

        // Weight is computed in the undeformed space
        double total_weight = 0;
        double total_shiftX = 0;
        double total_shiftY = 0;
        for (YY(sh) = -1; YY(sh) <= 1; YY(sh)++)
            for (XX(sh) = -1; XX(sh) <= 1; XX(sh)++)
            {
                V2_PLUS_V2(ri, r, sh);
                if (!cell_shiftX.outside(ri))
                {
                    total_weight += weight(INDEX(sh));
                    total_shiftX += weight(INDEX(sh)) * cell_shiftX(INDEX(ri));
                    total_shiftY += weight(INDEX(sh)) * cell_shiftY(INDEX(ri));
                }
            }
        if (total_weight == 0)
            total_weight = 1;
        if (!cell_shiftX.outside(r))
        {
            cell_shiftX(INDEX(r)) = rnd_gaus(prm_crystal.Nshift_avg,
                                             prm_crystal.Nshift_dev) + total_shiftX / total_weight;
            cell_shiftY(INDEX(r)) = rnd_gaus(prm_crystal.Nshift_avg,
                                             prm_crystal.Nshift_dev) + total_shiftY / total_weight;
        }

        // Move to next position
        move_following_spiral(r, visited);
    }

#ifdef DEBUG
    std::cout << "Cell shift X without absolute displacements" << cell_shiftX;
    std::cout << "Cell shift Y without absolute displacements" << cell_shiftY;
#endif

    // The previous shifts are relative to the final position, now
    // express the real final position

    FOR_ALL_ELEMENTS_IN_ARRAY2D(visited)
    {
        if (!cell_shiftX.outside(i, j))
        {
            // Move to final position
            cell_shiftX(i, j) += j * XX(aprojd) + i * XX(bprojd);
            cell_shiftY(i, j) += j * YY(aprojd) + i * YY(bprojd);

            // Check if there is intersection
            Matrix1D<double> auxcorner1(2), auxcorner2(2);
            XX(auxcorner1) = XX(corner1) + cell_shiftX(i, j);
            YY(auxcorner1) = YY(corner1) + cell_shiftY(i, j);
            XX(auxcorner2) = XX(corner2) + cell_shiftX(i, j);
            YY(auxcorner2) = YY(corner2) + cell_shiftY(i, j);

            cell_inside(i, j) =
                (XMIPP_MAX(STARTINGY (P()),YY(auxcorner1))<=
                 XMIPP_MIN(FINISHINGY(P()),YY(auxcorner2))) &&
                (XMIPP_MAX(STARTINGX (P()),XX(auxcorner1))<=
                 XMIPP_MIN(FINISHINGX(P()),XX(auxcorner2)));
            //TEMPORAL FIX FOR PHANTOM AS BIG AS THE WHOLE CRYSTAL
            cell_inside(i, j) = 1;
            //ROBERTO
#ifdef DEBUG

            std::cout << "(i,j)=(" << i << "," << j << ")\n";
            std::cout << "   Projection shape ";
            P().printShape();
            std::cout << std::endl;
            std::cout << "   AuxCorner1 " << auxcorner1.transpose() << std::endl
            << "   Origin     " << cell_shiftX(i, j) << " "
            << cell_shiftY(i, j) << std::endl
            << "   AuxCorner2 " << auxcorner2.transpose() << std::endl;
            std::cout << "   Inside= " << cell_inside(i, j) << std::endl;
#endif

        }
    }
}

/* Fill aux matrix with experimental shifs to add to unit cell
   projection********************************************************/
void init_shift_matrix(const Crystal_Projection_Parameters &prm_crystal,
                       MultidimArray<int>    &cell_inside,
                       MultidimArray<double> &exp_shifts_matrix_X,
                       MultidimArray<double> &exp_shifts_matrix_Y,
                       MultidimArray<double> &exp_shifts_matrix_Z,
                       MultidimArray<double> &exp_normal_shifts_matrix_X,
                       MultidimArray<double> &exp_normal_shifts_matrix_Y,
                       MultidimArray<double> &exp_normal_shifts_matrix_Z,
                       double phantom_scale)
{
    MetaData aux_DF_shift; //crystal_param is cont
    aux_DF_shift = prm_crystal.DF_shift;
    exp_shifts_matrix_X.resize(cell_inside);
    exp_shifts_matrix_X.initZeros();
    exp_shifts_matrix_Y.resize(cell_inside);
    exp_shifts_matrix_Y.initZeros();
    exp_shifts_matrix_Z.resize(cell_inside);
    exp_shifts_matrix_Z.initZeros();

    exp_normal_shifts_matrix_X.resize(cell_inside);
    exp_normal_shifts_matrix_X.initZeros();
    exp_normal_shifts_matrix_Y.resize(cell_inside);
    exp_normal_shifts_matrix_Y.initZeros();
    exp_normal_shifts_matrix_Z.resize(cell_inside);
    exp_normal_shifts_matrix_Z.initZeros();

    //#define DEBUG2
#ifdef DEBUG2

    std::cout << aux_DF_shift;
    std::cout << "exp_shifts_matrix_X shape" << std::endl;
    exp_shifts_matrix_X.printShape();
    std::cout << std::endl;
#endif
#undef DEBUG2
    //fill matrix with docfile data

    FOR_ALL_OBJECTS_IN_METADATA(aux_DF_shift)
    {
        //Check that we are not outside the matrix
        int xcell, ycell;
        aux_DF_shift.getValue(MDL_CRYSTAL_CELLX,xcell,__iter.objId);
        aux_DF_shift.getValue(MDL_CRYSTAL_CELLY,ycell,__iter.objId);
        if (!exp_shifts_matrix_X.outside(xcell,ycell))
        {
            aux_DF_shift.getValue(MDL_SHIFT_X,exp_shifts_matrix_X(ycell, xcell),__iter.objId);
            aux_DF_shift.getValue(MDL_SHIFT_Y,exp_shifts_matrix_Y(ycell, xcell),__iter.objId);
            aux_DF_shift.getValue(MDL_SHIFT_Z,exp_shifts_matrix_Z(ycell, xcell),__iter.objId);

            aux_DF_shift.getValue(MDL_CRYSTAL_SHIFTX,exp_normal_shifts_matrix_X(ycell, xcell),__iter.objId);
            aux_DF_shift.getValue(MDL_CRYSTAL_SHIFTY,exp_normal_shifts_matrix_Y(ycell, xcell),__iter.objId);
            aux_DF_shift.getValue(MDL_CRYSTAL_SHIFTZ,exp_normal_shifts_matrix_Z(ycell, xcell),__iter.objId);
        }
    }
    //#define DEBUG2
#ifdef DEBUG2
    std::cout << "exp_shifts_matrix_X" << exp_shifts_matrix_X;
    std::cout << "exp_shifts_matrix_Y" << exp_shifts_matrix_Y;
    std::cout << "exp_shifts_matrix_Z" << exp_shifts_matrix_Z;
#endif
#undef DEBUG2
}

