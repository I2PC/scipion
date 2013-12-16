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

#include <stdio.h>

#include "symmetries.h"

// Read Symmetry file ======================================================
// crystal symmetry matices from http://cci.lbl.gov/asu_gallery/
int SymList::readSymmetryFile(FileName fn_sym, double accuracy)
{
    int i, j;
    FILE *fpoii;
    char line[80];
    char *auxstr;
    double ang_incr, rot_ang;
    int  fold;
    Matrix2D<double> L(4, 4), R(4, 4);
    Matrix1D<double> axis(3), shift(3);
    int pgGroup, pgOrder;
    std::vector<std::string> fileContent;

    //check if reserved word

    // Open file ---------------------------------------------------------
    if ((fpoii = fopen(fn_sym.c_str(), "r")) == NULL)
    {
        //check if reserved word and return group and order
        if (isSymmetryGroup(fn_sym, pgGroup, pgOrder))
        {
            fillSymmetryClass(fn_sym, pgGroup, pgOrder,fileContent);
        }
        else
            REPORT_ERROR(ERR_IO_NOTOPEN, (std::string)"SymList::read_sym_file:Can't open file: "
                         + " or do not recognize symmetry group" + fn_sym);
    }
    else
    {
        while (fgets(line, 79, fpoii) != NULL)
        {
            if (line[0] == ';' || line[0] == '#' || line[0] == '\0')
                continue;
            fileContent.push_back(line);
        }
        fclose(fpoii);
    }

    //reset space_group
    space_group = 0;
    // Count the number of symmetries ------------------------------------
    true_symNo = 0;
    // count number of axis and mirror planes. It will help to identify
    // the crystallographic symmetry

    int no_axis, no_mirror_planes, no_inversion_points;
    no_axis = no_mirror_planes = no_inversion_points = 0;

    for (size_t n=0; n<fileContent.size(); n++)
    {
        strcpy(line,fileContent[n].c_str());
        auxstr = firstToken(line);
        if (auxstr == NULL)
        {
            std::cout << line;
            std::cout << "Wrong line in symmetry file, the line is skipped\n";
            continue;
        }
        if (strcmp(auxstr, "rot_axis") == 0)
        {
            auxstr = nextToken();
            fold = textToInteger(auxstr);
            true_symNo += (fold - 1);
            no_axis++;
        }
        else if (strcmp(auxstr, "mirror_plane") == 0)
        {
            true_symNo++;
            no_mirror_planes++;
        }
        else if (strcmp(auxstr, "inversion") == 0)
        {
            true_symNo += 1;
            no_inversion_points = 1;
        }
        else if (strcmp(auxstr, "P4212") == 0)
            true_symNo += 7;
        else if (strcmp(auxstr, "P2_122") == 0)
            true_symNo += 3;
        else if (strcmp(auxstr, "P22_12") == 0)
            true_symNo += 3;
    }
    // Ask for memory
    __L.resize(4*true_symNo, 4);
    __R.resize(4*true_symNo, 4);
    __shift.resize(true_symNo, 3);
    __chain_length.resize(true_symNo);
    __chain_length.initConstant(1);

    // Read symmetry parameters
    i = 0;
    shift.initZeros();

    for (size_t n=0; n<fileContent.size(); n++)
    {
        strcpy(line,fileContent[n].c_str());
        auxstr = firstToken(line);
        // Rotational axis ---------------------------------------------------
        if (strcmp(auxstr, "rot_axis") == 0)
        {
            auxstr = nextToken();
            fold = textToInteger(auxstr);
            auxstr = nextToken();
            XX(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            YY(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            ZZ(axis) = textToFloat(auxstr);
            ang_incr = 360. / fold;
            L.initIdentity();
            for (j = 1, rot_ang = ang_incr; j < fold; j++, rot_ang += ang_incr)
            {
                rotation3DMatrix(rot_ang, axis, R);
                setShift(i, shift);
                setMatrices(i++, L, R.transpose());
            }
            __sym_elements++;
            // inversion ------------------------------------------------------
        }
        else if (strcmp(auxstr, "inversion") == 0)
        {
            L.initIdentity();
            L(2, 2) = -1;
            R.initIdentity();
            R(0, 0) = -1.;
            R(1, 1) = -1.;
            R(2, 2) = -1.;
            setShift(i, shift);
            setMatrices(i++, L, R);
            __sym_elements++;
            // mirror plane -------------------------------------------------------------
        }
        else if (strcmp(auxstr, "mirror_plane") == 0)
        {
            auxstr = nextToken();
            XX(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            YY(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            ZZ(axis) = textToFloat(auxstr);
            L.initIdentity();
            L(2, 2) = -1;
            Matrix2D<double> A;
            alignWithZ(axis,A);
            A = A.transpose();
            R = A * L * A.inv();
            setShift(i, shift);
            L.initIdentity();
            setMatrices(i++, L, R);
            __sym_elements++;
            // P4212 -------------------------------------------------------------
        }
        else if (strcmp(auxstr, "P4212") == 0)
        {
            space_group = sym_P42_12;
            accuracy = -1; // Do not compute subgroup
            L.initIdentity();

            // With 0 shift
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = R(1, 1) = -1;
            R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(2, 2) = -1;
            R(0, 1) = R(1, 0) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(2, 2) = R(0, 1) = R(1, 0) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);

            // With 1/2 shift
            VECTOR_R3(shift, 0.5, 0.5, 0);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 1) = -1;
            R(1, 0) = R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(1, 0) = -1;
            R(0, 1) = R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = R(2, 2) = -1;
            R(1, 1) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(1, 1) = R(2, 2) = -1;
            R(0, 0) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);

            __sym_elements++;
        }
        else if (strcmp(auxstr, "P2_122") == 0)
        {
            space_group = sym_P2_122;
            accuracy = -1; // Do not compute subgroup
            L.initIdentity();

            // With 0 shift
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = -1;
            R(1, 1) = -1;
            R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);

            // With 1/2 shift
            VECTOR_R3(shift, 0.5, 0.0, 0.0);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = -1;
            R(1, 1) = 1;
            R(2, 2) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = 1;
            R(1, 1) = -1;
            R(2, 2) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            __sym_elements++;
        }
        else if (strcmp(auxstr, "P22_12") == 0)
        {
            space_group = sym_P22_12;
            accuracy = -1; // Do not compute subgroup
            L.initIdentity();

            // With 0 shift
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = -1;
            R(1, 1) = -1;
            R(2, 2) = 1;
            setShift(i, shift);
            setMatrices(i++, L, R);

            // With 1/2 shift
            VECTOR_R3(shift, 0.0, 0.5, 0.0);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = 1;
            R(1, 1) = -1;
            R(2, 2) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            R.initZeros();
            R(3, 3) = 1;
            R(0, 0) = -1;
            R(1, 1) = 1;
            R(2, 2) = -1;
            setShift(i, shift);
            setMatrices(i++, L, R);
            __sym_elements++;
        }
    }

    if (accuracy > 0)
        computeSubgroup(accuracy);

    //possible crystallographic symmetry
    if (no_axis == 0 && no_mirror_planes == 0 && no_inversion_points == 0 &&
        true_symNo == 7 && space_group == sym_P42_12)
        space_group = sym_P42_12;
    else if (no_axis == 0 && no_mirror_planes == 0 && no_inversion_points == 0 &&
             true_symNo == 3 && space_group == sym_P2_122)
        space_group = sym_P2_122;
    else if (no_axis == 0 && no_mirror_planes == 0 && no_inversion_points == 0 &&
             true_symNo == 3 && space_group == sym_P22_12)
        space_group = sym_P22_12;
    // P4 and P6
    else if (no_axis == 1 && no_mirror_planes == 0 && no_inversion_points == 0 &&
             fabs(R(2, 2) - 1.) < XMIPP_EQUAL_ACCURACY &&
             fabs(R(0, 0) - R(1, 1)) < XMIPP_EQUAL_ACCURACY &&
             fabs(R(0, 1) + R(1, 0)) < XMIPP_EQUAL_ACCURACY)
    {
        switch (true_symNo)
        {
        case(5):
                        space_group = sym_P6;
            break;
        case(3):
                        space_group = sym_P4;
            break;
        default:
            space_group = sym_undefined;
            break;
        }//switch end
    }//end else if (no_axis==1 && no_mirror_planes== 0
    else if (no_axis == 0 && no_inversion_points == 0 && no_mirror_planes == 0)
        space_group = sym_P1;
    else
        space_group = sym_undefined;
    return pgGroup;
}

// Get matrix ==============================================================
void SymList::getMatrices(int i, Matrix2D<double> &L, Matrix2D<double> &R,
                          bool homogeneous)
const
{
    int k, kp, l;
    if (homogeneous)
    {
        L.initZeros(4, 4);
        R.initZeros(4, 4);
        for (k = 4 * i, kp=0; k < 4*i + 4; k++, kp++)
            for (l = 0; l < 4; l++)
            {
                dMij(L,kp, l) = dMij(__L,k, l);
                dMij(R,kp, l) = dMij(__R,k, l);
            }
    }
    else
    {
        L.initZeros(3, 3);
        R.initZeros(3, 3);
        for (k = 4 * i, kp=0; k < 4*i + 3; k++, kp++)
            for (l = 0; l < 3; l++)
            {
                dMij(L,kp, l) = dMij(__L,k, l);
                dMij(R,kp, l) = dMij(__R,k, l);
            }
    }
}

// Set matrix ==============================================================
void SymList::setMatrices(int i, const Matrix2D<double> &L,
                          const Matrix2D<double> &R)
{
    int k, l;
    for (k = 4 * i; k < 4*i + 4; k++)
        for (l = 0; l < 4; l++)
        {
            __L(k, l) = L(k - 4 * i, l);
            __R(k, l) = R(k - 4 * i, l);
        }
}

// Get/Set shift ===========================================================
void SymList::getShift(int i, Matrix1D<double> &shift) const
{
    shift.resize(3);
    XX(shift) = __shift(i, 0);
    YY(shift) = __shift(i, 1);
    ZZ(shift) = __shift(i, 2);
}

void SymList::setShift(int i, const Matrix1D<double> &shift)
{
    if (shift.size() != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "SymList::add_shift: Shift vector is not 3x1");
    __shift(i, 0) = XX(shift);
    __shift(i, 1) = YY(shift);
    __shift(i, 2) = ZZ(shift);
}

void SymList::addShift(const Matrix1D<double> &shift)
{
    if (shift.size() != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "SymList::add_shift: Shift vector is not 3x1");
    int i = MAT_YSIZE(__shift);
    __shift.resize(i + 1, 3);
    setShift(i, shift);
}

// Add matrix ==============================================================
void SymList::addMatrices(const Matrix2D<double> &L, const Matrix2D<double> &R,
                          int chain_length)
{
    if (MAT_XSIZE(L) != 4 || MAT_YSIZE(L) != 4 || MAT_XSIZE(R) != 4 || MAT_YSIZE(R) != 4)
        REPORT_ERROR(ERR_MATRIX_SIZE, "SymList::add_matrix: Transformation matrix is not 4x4");
    if (trueSymsNo() == symsNo())
    {
        __L.resize(MAT_YSIZE(__L) + 4, 4);
        __R.resize(MAT_YSIZE(__R) + 4, 4);
        __chain_length.resize(__chain_length.size() + 1);
    }

    setMatrices(true_symNo, L, R);
    __chain_length(__chain_length.size() - 1) = chain_length;
    true_symNo++;
}

// Compute subgroup ========================================================
bool found_not_tried(const Matrix2D<int> &tried, int &i, int &j,
                     int true_symNo)
{
    i = j = 0;
    size_t n = 0;
    while (n != MAT_YSIZE(tried))
    {
        //       if (tried(i, j) == 0 && !(i >= true_symNo && j >= true_symNo))
        if (dMij(tried,i, j) == 0 && !(i >= true_symNo && j >= true_symNo))
            return true;
        if (i != (int)n)
        {
            // Move downwards
            i++;
        }
        else
        {
            // Move leftwards
            j--;
            if (j == -1)
            {
                n++;
                j = n;
                i = 0;
            }
        }
    }
    return false;
}

//#define DEBUG
void SymList::computeSubgroup(double accuracy)
{
    Matrix2D<double> I(4, 4);
    I.initIdentity();
    Matrix2D<double> L1(4, 4), R1(4, 4), L2(4, 4), R2(4, 4), newL(4, 4), newR(4, 4),identity(4,4);
    Matrix2D<int>    tried(true_symNo, true_symNo);
    Matrix1D<double> shift(3);
    shift.initZeros();
    int i, j;
    int new_chain_length;
    identity.initIdentity();
    while (found_not_tried(tried, i, j, true_symNo))
    {
        tried(i, j) = 1;

        // Form new symmetry matrices
        // if (__chain_length(i)+__chain_length(j)>__sym_elements+2) continue;

        getMatrices(i, L1, R1);
        getMatrices(j, L2, R2);
        newL = L1 * L2;
        newR = R1 * R2;
        new_chain_length = __chain_length(i) + __chain_length(j);

        //if (newL.isIdentity() && newR.isIdentity()) continue;
        //rounding error make newR different from identity
        if (newL.equal(identity, accuracy) &&
            newR.equal(identity, accuracy))
            continue;

        // Try to find it in current ones
        bool found;
        found = false;
        for (int l = 0; l < symsNo(); l++)
        {
            getMatrices(l, L1, R1);
            if (newL.equal(L1, accuracy) && newR.equal(R1, accuracy))
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            //#define DEBUG
#ifdef DEBUG
            static int kjhg=0;
            /* std::cout << "Matrix size " << tried.Xdim() << " "
             << "trying " << i << " " << j << " "
             << "chain length=" << new_chain_length << std::endl;
             std::cout << "Result R Sh\n" << newR << shift;
             */
            //std::cerr << "shift" << __shift <<std::endl;
            std::cerr << "newR "  << kjhg++ << "\n" << newR <<std::endl;
#endif

            addMatrices(newL, newR, new_chain_length);
            addShift(shift);
            tried.resize(MAT_YSIZE(tried) + 1, MAT_XSIZE(tried) + 1);
        }
    }
#ifdef DEBUG
    std::cerr << "__R" << __R <<std::endl;
#endif
#undef DEBUG
}
/** Guess Crystallographic space group.
    Return the group number
    http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen. So
    far it has only been implemented for P1 (1), P2_122 (17) ,
    P22_12,
    P4 (75), P4212 (90) and P6 (168) */

int  SymList::crystallographicSpaceGroup(double mag_a, double mag_b,
        double ang_a2b_deg) const
{

    switch (space_group)
    {
    case sym_undefined:
    case sym_P1:
        return(space_group);
    case sym_P4:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 90) > XMIPP_EQUAL_ACCURACY)
            std::cerr << "\nWARNING: P42 but mag_a != mag_b\n"
            << " or ang_a2b !=90" << std::endl;
        return(space_group);
        break;
    case sym_P2_122:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 90) > XMIPP_EQUAL_ACCURACY)
            std::cerr << "\nWARNING: P2_122 but mag_a != mag_b\n"
            << " or ang_a2b !=90" << std::endl;
        return(space_group);
        break;
    case sym_P22_12:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 90) > XMIPP_EQUAL_ACCURACY)
            std::cerr << "\nWARNING: P22_12 but mag_a != mag_b\n"
            << " or ang_a2b !=90" << std::endl;
        return(space_group);
        break;
    case sym_P42_12:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 90) > XMIPP_EQUAL_ACCURACY)
            std::cerr << "\nWARNING: P42_12 but mag_a != mag_b\n"
            << " or ang_a2b !=90" << std::endl;
        return(space_group);
        break;
    case sym_P6:
        if (fabs((mag_a - mag_b)) > XMIPP_EQUAL_ACCURACY ||
            fabs(ang_a2b_deg - 120.) > XMIPP_EQUAL_ACCURACY)
        {
            std::cerr << "\nWARNING: marked as P6 but mag_a != mag_b\n"
            << "or ang_a2b !=120" << std::endl;
            std::cerr << "\nWARNING: P1 is assumed\n";
            return(sym_P1);
        }
        else
            return(space_group);
        break;
    default:
        std::cerr << "\n Congratulations: you have found a bug in the\n"
        << "routine crystallographic_space_group or\n"
        << "You have called to this rotuine BEFORE reading\n"
        << "the symmetry info" << std::endl;
        exit(0);
        break;
    }//switch(space_group)  end

}//crystallographicSpaceGroup end

// Symmetrize_crystal_vectors==========================================
//A note: Please realize that we are not repeating code here.
//The class SymList deals with symmetries when expressed in
//Cartesian space, that is the basis is orthonormal. Here
//we describe symmetries in the crystallographic way
//that is, the basis and the crystal vectors are the same.
//For same symmetries both representations are almost the same
//but in general they are rather different.

//IMPORTANT: matrix orden should match the one used in "readSymmetryFile"
//if not the wrong angles are assigned to the different matrices

void symmetrizeCrystalVectors(Matrix1D<double> &aint,
                              Matrix1D<double> &bint,
                              Matrix1D<double> &shift,
                              int space_group,
                              int sym_no,
                              const Matrix1D<double> &eprm_aint,
                              const Matrix1D<double> &eprm_bint)
{
    //Notice that we should use R.inv and not R to relate eprm.aint and aint
    shift.initZeros();//I think this init is OK even the vector dim=0
    switch (space_group)
    {
    case(sym_undefined):
                case(sym_P1):
                        XX(aint) =   XX(eprm_aint);
        YY(aint) =                   YY(eprm_aint);
        XX(bint) =   XX(eprm_bint);
        YY(bint) =                   YY(eprm_bint);
        break;
    case(sym_P2):       std::cerr << "\n Group P2 not implemented\n";
        exit(1);
        break;
    case(sym_P2_1):     std::cerr << "\n Group P2_1 not implemented\n";
        exit(1);
        break;
    case(sym_C2):       std::cerr << "\n Group C2 not implemented\n";
        exit(1);
        break;
    case(sym_P222):     std::cerr << "\n Group P222 not implemented\n";
        exit(1);
        break;
    case(sym_P2_122):
                    switch (sym_no)
            {
            case(-1): XX(aint) =   XX(eprm_aint);
                YY(aint) =                   YY(eprm_aint);
                XX(bint) =   XX(eprm_bint);
                YY(bint) =                   YY(eprm_bint);
                break;
            case(0):  XX(aint) = - XX(eprm_aint);
                YY(aint) =                 - YY(eprm_aint);
                XX(bint) = - XX(eprm_bint);
                YY(bint) =                 - YY(eprm_bint);
                break;
            case(1):  XX(aint) = - XX(eprm_aint);
                YY(aint) =                 + YY(eprm_aint);
                XX(bint) = - XX(eprm_bint);
                YY(bint) =                 + YY(eprm_bint);
                VECTOR_R3(shift, 0.5, 0.0, 0.0);
                break;
            case(2):  XX(aint) = + XX(eprm_aint);
                YY(aint) =                 - YY(eprm_aint);
                XX(bint) = + XX(eprm_bint);
                YY(bint) =                 - YY(eprm_bint);
                VECTOR_R3(shift, 0.5, 0.0, 0.0);
                break;
            }//switch P2_122 end
        break;
    case(sym_P22_12):
                    switch (sym_no)
            {
            case(-1): XX(aint) =   XX(eprm_aint);
                YY(aint) =                   YY(eprm_aint);
                XX(bint) =   XX(eprm_bint);
                YY(bint) =                   YY(eprm_bint);
                break;
            case(0):  XX(aint) = - XX(eprm_aint);
                YY(aint) =                 - YY(eprm_aint);
                XX(bint) = - XX(eprm_bint);
                YY(bint) =                 - YY(eprm_bint);
                break;
            case(1):  XX(aint) =   XX(eprm_aint);
                YY(aint) =                 - YY(eprm_aint);
                XX(bint) =   XX(eprm_bint);
                YY(bint) =                 - YY(eprm_bint);
                VECTOR_R3(shift, 0.0, 0.5, 0.0);
                break;
            case(2):  XX(aint) = - XX(eprm_aint);
                YY(aint) =                  YY(eprm_aint);
                XX(bint) = - XX(eprm_bint);
                YY(bint) =                  YY(eprm_bint);
                VECTOR_R3(shift, 0.0, 0.5, 0.0);
                break;
            }//switch P22_12 end
        break;

    case(sym_P22_12_1): std::cerr << "\n Group P22_12_1 not implemented\n";
        exit(1);
        break;
    case(sym_P4):
                    switch (sym_no)
            {
            case(-1): XX(aint) =   XX(eprm_aint);
                YY(aint) =                   YY(eprm_aint);
                XX(bint) =   XX(eprm_bint);
                YY(bint) =                   YY(eprm_bint);
                break;
            case(0):  XX(aint) =                 - YY(eprm_aint);
                YY(aint) =   XX(eprm_aint);
                XX(bint) =                 - YY(eprm_bint);
                YY(bint) =   XX(eprm_bint);
                break;
            case(1):  XX(aint) = - XX(eprm_aint);
                YY(aint) =                 - YY(eprm_aint);
                XX(bint) = - XX(eprm_bint);
                YY(bint) =                 - YY(eprm_bint);
                break;
            case(2):  XX(aint) =                   YY(eprm_aint);
                YY(aint) = - XX(eprm_aint);
                XX(bint) =                   YY(eprm_bint);
                YY(bint) = - XX(eprm_bint);
                break;
            }//switch P4 end
        break;
    case(sym_P422):     REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Group P422 not implemented");
        break;
    case(sym_P42_12):
                    switch (sym_no)
            {
            case(-1): XX(aint) =   XX(eprm_aint);
                YY(aint) =                   YY(eprm_aint);
                XX(bint) =   XX(eprm_bint);
                YY(bint) =                   YY(eprm_bint);
                break;
            case(0):  XX(aint) = - XX(eprm_aint);
                YY(aint) =                 - YY(eprm_aint);
                XX(bint) = - XX(eprm_bint);
                YY(bint) =                 - YY(eprm_bint);
                break;
            case(1):  XX(aint) =                 + YY(eprm_aint);
                YY(aint) = + XX(eprm_aint);
                XX(bint) =                 + YY(eprm_bint);
                YY(bint) = + XX(eprm_bint);
                break;
            case(2):  XX(aint) =                 - YY(eprm_aint);
                YY(aint) = - XX(eprm_aint);
                XX(bint) =                 - YY(eprm_bint);
                YY(bint) = - XX(eprm_bint);
                break;
            case(3):  XX(aint) =                 + YY(eprm_aint);
                YY(aint) = - XX(eprm_aint);
                XX(bint) =                 + YY(eprm_bint);
                YY(bint) = - XX(eprm_bint);
                VECTOR_R3(shift, 0.5, 0.5, 0);
                break;
            case(4):  XX(aint) =          - YY(eprm_aint);
                YY(aint) = + XX(eprm_aint);
                XX(bint) =   - YY(eprm_bint);
                YY(bint) = + XX(eprm_bint);
                VECTOR_R3(shift, 0.5, 0.5, 0);
                break;
            case(5):  XX(aint) = - XX(eprm_aint);
                YY(aint) =          + YY(eprm_aint);
                XX(bint) = - XX(eprm_bint);
                YY(bint) =          + YY(eprm_bint);
                VECTOR_R3(shift, 0.5, 0.5, 0);
                break;
            case(6):  XX(aint) = + XX(eprm_aint);
                YY(aint) =          - YY(eprm_aint);
                XX(bint) = + XX(eprm_bint);
                YY(bint) =          - YY(eprm_bint);
                VECTOR_R3(shift, 0.5, 0.5, 0);
                break;
            default:
                std::cout << "\n Wrong symmetry number "
                "in symmetrize_crystal_vectors, bye" << std::endl;
                exit(1);
                break;


            }//switch P4212 end
        break;
    case(sym_P3):       std::cerr << "\n Group P3 not implemented\n";
        exit(1);
        break;
    case(sym_P312):     std::cerr << "\n Group P312 not implemented\n";
        exit(1);
        break;
    case(sym_P6):
                    switch (sym_no)
            {
            case(-1): XX(aint) =   XX(eprm_aint);
                YY(aint) =                   YY(eprm_aint);
                XX(bint) =   XX(eprm_bint);
                YY(bint) =                   YY(eprm_bint);
                break;
            case(0):  XX(aint) =   XX(eprm_aint) - YY(eprm_aint);
                YY(aint) =   XX(eprm_aint);
                XX(bint) =   XX(eprm_bint) - YY(eprm_bint);
                YY(bint) =   XX(eprm_bint);
                break;
            case(1):  XX(aint) =                 - YY(eprm_aint);
                YY(aint) =   XX(eprm_aint) - YY(eprm_aint);
                XX(bint) =                 - YY(eprm_bint);
                YY(bint) =   XX(eprm_bint) - YY(eprm_bint);
                break;
            case(2):  XX(aint) = - XX(eprm_aint);
                YY(aint) =                 - YY(eprm_aint);
                XX(bint) = - XX(eprm_bint);
                YY(bint) =                 - YY(eprm_bint);
                break;
            case(3):  XX(aint) = - XX(eprm_aint) + YY(eprm_aint);
                YY(aint) = - XX(eprm_aint);
                XX(bint) = - XX(eprm_bint) + YY(eprm_bint);
                YY(bint) = - XX(eprm_bint);
                break;
            case(4):  XX(aint) =                 + YY(eprm_aint);
                YY(aint) = - XX(eprm_aint) + YY(eprm_aint);
                XX(bint) =                 + YY(eprm_bint);
                YY(bint) = - XX(eprm_bint) + YY(eprm_bint);
                break;
            }//switch P6 end
        break;

    case(sym_P622):     std::cerr << "\n Group P622 not implemented\n";
        exit(1);
        break;
    }

}//symmetrizeCrystalVectors end

#define Symmetrize_Vol(X) {\
        for (size_t i=0; i<vol_in.VolumesNo(); i++)\
            X(vol_in(i),vol_in.grid(i),eprm_aint,eprm_bint,mask,i, \
              grid_type);\
    }

// Symmetrize_crystal_volume==========================================
//IMPORTANT: matrix orden should match the one used in "readSymmetryFile"
//if not the wrong angles are assigned to the different matrices
void symmetrizeCrystalVolume(GridVolume &vol_in,
                             const Matrix1D<double> &eprm_aint,
                             const Matrix1D<double> &eprm_bint,
                             int eprm_space_group,
                             const MultidimArray<int> &mask, int grid_type)
{
    //SO FAR ONLY THE GRID CENTERED IN 0,0,0 IS SYMMETRIZED, THE OTHER
    //ONE SINCE REQUIRE INTERPOLATION IS IGNORED
    switch (eprm_space_group)
    {
    case(sym_undefined):
                case(sym_P1):
                        break;
    case(sym_P2):       std::cerr << "\n Group P2 not implemented\n";
        exit(1);
        break;
    case(sym_P2_1):     std::cerr << "\n Group P2_1 not implemented\n";
        exit(1);
        break;
    case(sym_C2):       std::cerr << "\n Group C2 not implemented\n";
        exit(1);
        break;
    case(sym_P222):     std::cerr << "\n Group P222 not implemented\n";
        exit(1);
        break;
    case(sym_P2_122):
                    Symmetrize_Vol(symmetry_P2_122)//already has ;
                    break;
    case(sym_P22_12):
                    Symmetrize_Vol(symmetry_P22_12)//already has ;
                    break;
    case(sym_P22_12_1): std::cerr << "\n Group P22_12_1 not implemented\n";
        exit(1);
        break;
    case(sym_P4):
                    Symmetrize_Vol(symmetry_P4)//already has ;
                    break;
    case(sym_P422):     std::cerr << "\n Group P422 not implemented\n";
        exit(1);
        break;
    case(sym_P42_12):
                    Symmetrize_Vol(symmetry_P42_12)//already has ;
                    break;
    case(sym_P3):       std::cerr << "\n Group P3 not implemented\n";
        exit(1);
        break;
    case(sym_P312):     std::cerr << "\n Group P312 not implemented\n";
        exit(1);
        break;
    case(sym_P6):
                    Symmetrize_Vol(symmetry_P6)//already has ;
                    break;
    case(sym_P622):     std::cerr << "\n Group P622 not implemented\n";
        exit(1);
        break;
    }


}//symmetrizeCrystalVectors end
#define put_inside(j,j_min,j_max,jint)  \
    if( (j) < (j_min) ) { (j) = (j) + (jint);}\
    else if( (j) > (j_max) ) { (j) = (j) - (jint);};

/* Symmetrizes a simple grid with P2_122 symmetry--------------------------*/

void symmetry_P2_122(Image<double> &vol, const SimpleGrid &grid,
                     const Matrix1D<double> &eprm_aint,
                     const Matrix1D<double> &eprm_bint,
                     const MultidimArray<int> &mask, int volume_no,
                     int grid_type)
{

    int ZZ_lowest = (int) ZZ(grid.lowest);
    int YY_lowest = STARTINGY(mask);
    int XX_lowest = STARTINGX(mask);
    int ZZ_highest = (int) ZZ(grid.highest);
    int YY_highest = FINISHINGY(mask);
    int XX_highest = FINISHINGX(mask);

    //if there is an extra slice in the z direction there is no way
    //to calculate the -z slice
    if (ABS(ZZ_lowest) > ABS(ZZ_highest))
        ZZ_lowest = -(ZZ_highest);
    else
        ZZ_highest = ABS(ZZ_lowest);

    while (1)
    {
        if (mask(0, XX_lowest) == 0)
            XX_lowest++;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P2_122, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(0, XX_highest) == 0)
            XX_highest--;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P2_122, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_lowest, 0) == 0)
            YY_lowest++;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P2_122, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_highest, 0) == 0)
            YY_highest--;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P2_122, while(1)" << std::endl;
            exit(0);
        }

    }


    int maxZ, maxY, maxX, minZ, minY, minX;
    int x, y, z;

    int x0, y0, z0;
    int x1, y1, z1;
    int x2, y2, z2;

    int XXaint, YYbint;
    XXaint = (int) XX(eprm_aint);
    YYbint = (int)YY(eprm_bint);

    int XXaint_2;
    XXaint_2 = XXaint / 2;
    int xx, yy, zz;

    if (ABS(XX_lowest) > ABS(XX_highest))
    {
        minX = XX_lowest;
        maxX = 0;
    }
    else
    {
        minX = 0;
        maxX = XX_highest;
    }
    if (ABS(YY_lowest) > ABS(YY_highest))
    {
        minY = YY_lowest;
        maxY = 0;
    }
    else
    {
        minY = 0;
        maxY = YY_highest;
    }
    if (ABS(ZZ_lowest) > ABS(ZZ_highest))
    {
        minZ = ZZ_lowest;
        maxZ = 0;
    }
    else
    {
        minZ = 0;
        maxZ = ZZ_highest;
    }

    //FCC non supported yet
    if (volume_no == 1 && grid_type == FCC)
    {
        std::cerr << "\nSimetries using FCC not implemented\n";
        exit(1);
    }

    for (z = minZ;z <= maxZ;z++)
        for (y = minY;y <= maxY;y++)
            for (x = minX;x <= maxX;x++)
            {
                //sym=-1---------------------------------------------------------
                if (!A2D_ELEM(mask, y, x) || z < ZZ_lowest || z > ZZ_highest)
                    continue;

                //sym=0 ---------------------------------------------------------
                xx = -x;
                yy = -y;
                zz = z;
                //    xx = x; yy=-y; zz=-z;
                if (volume_no == 1)
                {
                    xx--;
                    yy--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }
                x0 = xx;
                y0 = yy;
                z0 = zz;

                //sym=1----------------------------------------------------------
                xx = XXaint_2 - x;
                yy = y;
                zz = -z;
                //    xx = -x; yy= -y+YYbint_2; zz= -z;
                if (volume_no == 1)
                {
                    xx--;
                    zz--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=1 end
                x1 = xx;
                y1 = yy;
                z1 = zz;

                //sym=2----------------------------------------------------------
                xx = XXaint_2 + x;
                yy = -y;
                zz = -z;
                //    xx = -x; yy= y+YYbint_2; zz= -z;
                if (volume_no == 1)
                {
                    yy--;
                    zz--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=2 end
                x2 = xx;
                y2 = yy;
                z2 = zz;

                //only the first simple grid center and the origen is the same
                //point.
                //    if(volume_no==1)
                //    {
                //      switch (grid_type) {
                //  case FCC: std::cerr<< "\nSimetries using FCC not implemented\n";break;
                //  case BCC: x0--;y0--;               z1--;
                //            x2--;y2--;z2--;     y3--;
                //     x4--;          x5--;     z5--;
                //          y6--;z6--;
                //            break;
                //  case CC:  break;
                //      }
                //    }

                VOLVOXEL(vol, z , y , x) = VOLVOXEL(vol, z0, y0, x0) =
                                               VOLVOXEL(vol, z1, y1, x1) = VOLVOXEL(vol, z2, y2, x2) =
                                                                               (VOLVOXEL(vol, z , y , x) + VOLVOXEL(vol, z0, y0, x0) +
                                                                                VOLVOXEL(vol, z1, y1, x1) + VOLVOXEL(vol, z2, y2, x2)) / 4.0;
            }//for end
}//symmetryP2_122 end
/* Symmetrizes a simple grid with P2_122 symmetry--------------------------*/

void symmetry_P22_12(Image<double> &vol, const SimpleGrid &grid,
                     const Matrix1D<double> &eprm_aint,
                     const Matrix1D<double> &eprm_bint,
                     const MultidimArray<int> &mask, int volume_no,
                     int grid_type)
{

    int ZZ_lowest = (int) ZZ(grid.lowest);
    int YY_lowest = STARTINGY(mask);
    int XX_lowest = STARTINGX(mask);
    int ZZ_highest = (int) ZZ(grid.highest);
    int YY_highest = FINISHINGY(mask);
    int XX_highest = FINISHINGX(mask);

    //if there is an extra slice in the z direction there is no way
    //to calculate the -z slice
    if (ABS(ZZ_lowest) > ABS(ZZ_highest))
        ZZ_lowest = -(ZZ_highest);
    else
        ZZ_highest = ABS(ZZ_lowest);

    while (1)
    {
        if (mask(0, XX_lowest) == 0)
            XX_lowest++;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P2_122, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(0, XX_highest) == 0)
            XX_highest--;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P2_122, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_lowest, 0) == 0)
            YY_lowest++;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P2_122, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_highest, 0) == 0)
            YY_highest--;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P2_122, while(1)" << std::endl;
            exit(0);
        }

    }


    int maxZ, maxY, maxX, minZ, minY, minX;
    int x, y, z;

    int x0, y0, z0;
    int x1, y1, z1;
    int x2, y2, z2;

    int XXaint, YYbint;
    XXaint = (int) XX(eprm_aint);
    YYbint = (int)YY(eprm_bint);

    int YYbint_2;
    YYbint_2 = YYbint / 2;
    int xx, yy, zz;

    if (ABS(XX_lowest) > ABS(XX_highest))
    {
        minX = XX_lowest;
        maxX = 0;
    }
    else
    {
        minX = 0;
        maxX = XX_highest;
    }
    if (ABS(YY_lowest) > ABS(YY_highest))
    {
        minY = YY_lowest;
        maxY = 0;
    }
    else
    {
        minY = 0;
        maxY = YY_highest;
    }
    if (ABS(ZZ_lowest) > ABS(ZZ_highest))
    {
        minZ = ZZ_lowest;
        maxZ = 0;
    }
    else
    {
        minZ = 0;
        maxZ = ZZ_highest;
    }

    //FCC non supported yet
    if (volume_no == 1 && grid_type == FCC)
    {
        std::cerr << "\nSimetries using FCC not implemented\n";
        exit(1);
    }

    for (z = minZ;z <= maxZ;z++)
        for (y = minY;y <= maxY;y++)
            for (x = minX;x <= maxX;x++)
            {
                //sym=-1---------------------------------------------------------
                if (!A2D_ELEM(mask, y, x) || z < ZZ_lowest || z > ZZ_highest)
                    continue;

                //sym=0 ---------------------------------------------------------
                xx = -x;
                yy = -y;
                zz = z;
                //    xx = x; yy=-y; zz=-z;
                if (volume_no == 1)
                {
                    xx--;
                    yy--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }
                x0 = xx;
                y0 = yy;
                z0 = zz;

                //sym=1----------------------------------------------------------
                //    xx = XXaint_2-x; yy= y; zz= -z;
                xx = x;
                yy = -y + YYbint_2;
                zz = -z;
                if (volume_no == 1)
                {
                    yy--;
                    zz--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=1 end
                x1 = xx;
                y1 = yy;
                z1 = zz;

                //sym=2----------------------------------------------------------
                //    xx = XXaint_2+x; yy= -y; zz= -z;
                xx = -x;
                yy = y + YYbint_2;
                zz = -z;
                if (volume_no == 1)
                {
                    xx--;
                    zz--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=2 end
                x2 = xx;
                y2 = yy;
                z2 = zz;

                //only the first simple grid center and the origen is the same
                //point.
                //    if(volume_no==1)
                //    {
                //      switch (grid_type) {
                //  case FCC: std::cerr<< "\nSimetries using FCC not implemented\n";break;
                //  case BCC: x0--;y0--;               z1--;
                //            x2--;y2--;z2--;     y3--;
                //     x4--;          x5--;     z5--;
                //          y6--;z6--;
                //            break;
                //  case CC:  break;
                //      }
                //    }

                VOLVOXEL(vol, z , y , x) = VOLVOXEL(vol, z0, y0, x0) =
                                               VOLVOXEL(vol, z1, y1, x1) = VOLVOXEL(vol, z2, y2, x2) =
                                                                               (VOLVOXEL(vol, z , y , x) + VOLVOXEL(vol, z0, y0, x0) +
                                                                                VOLVOXEL(vol, z1, y1, x1) + VOLVOXEL(vol, z2, y2, x2)) / 4.0;
            }//for end
}//symmetryP2_122 end
/* Symmetrizes a simple grid with P4  symmetry --------------------------*/
void symmetry_P4(Image<double> &vol, const SimpleGrid &grid,
                 const Matrix1D<double> &eprm_aint,
                 const Matrix1D<double> &eprm_bint,
                 const MultidimArray<int> &mask, int volume_no, int grid_type)
{
    int ZZ_lowest = (int) ZZ(grid.lowest);
    int YY_lowest = STARTINGY(mask);
    int XX_lowest = STARTINGX(mask);
    int ZZ_highest = (int) ZZ(grid.highest);
    int YY_highest = FINISHINGY(mask);
    int XX_highest = FINISHINGX(mask);

    while (1)
    {
        if (mask(0, XX_lowest) == 0)
            XX_lowest++;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P4, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(0, XX_highest) == 0)
            XX_highest--;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P4, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_lowest, 0) == 0)
            YY_lowest++;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P4, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_highest, 0) == 0)
            YY_highest--;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P4, while(1)" << std::endl;
            exit(0);
        }

    }

    //   int ZZ_lowest =(int) ZZ(grid.lowest);
    //   int YY_lowest =MAX((int) YY(grid.lowest),STARTINGY(mask));
    //   int XX_lowest =MAX((int) XX(grid.lowest),STARTINGX(mask));
    //   int ZZ_highest=(int) ZZ(grid.highest);
    //   int YY_highest=MIN((int) YY(grid.highest),FINISHINGY(mask));
    //   int XX_highest=MIN((int) XX(grid.highest),FINISHINGX(mask));

    int maxZ, maxY, maxX, minZ, minY, minX;
    int x, y, z;

    int x0, y0, z0;
    int x1, y1, z1;
    int x2, y2, z2;

    int XXaint, YYbint;
    XXaint = (int) XX(eprm_aint);
    YYbint = (int)YY(eprm_bint);

    int xx, yy, zz;

    if (ABS(XX_lowest) > ABS(XX_highest))
    {
        minX = XX_lowest;
        maxX = 0;
    }
    else
    {
        minX = 0;
        maxX = XX_highest;
    }
    if (ABS(YY_lowest) > ABS(YY_highest))
    {
        minY = YY_lowest;
        maxY = 0;
    }
    else
    {
        minY = 0;
        maxY = YY_highest;
    }

    minZ = ZZ_lowest;
    maxZ = ZZ_highest;
    //FCC non supported yet
    if (volume_no == 1 && grid_type == FCC)
    {
        std::cerr << "\nSimetries using FCC not implemented\n";
        exit(1);
    }
    for (z = minZ;z <= maxZ;z++)
        for (y = minY;y <= maxY;y++)
            for (x = minX;x <= maxX;x++)
            {
                //sym=-1---------------------------------------------------------
                if (!A2D_ELEM(mask, y, x) || z < ZZ_lowest || z > ZZ_highest)
                    continue;

                //sym=0 ---------------------------------------------------------
                xx = -y;
                yy = x;
                zz = z;
                //only the first simple grid center and the origen is the same
                //point.
                if (volume_no == 1)
                {
                    xx--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }
                x0 = xx;
                y0 = yy;
                z0 = zz;

                //sym=1----------------------------------------------------------
                xx = -x;
                yy = -y;
                zz = z;
                if (volume_no == 1)
                {
                    xx--;
                    yy--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=1 end
                x1 = xx;
                y1 = yy;
                z1 = zz;

                //sym=2----------------------------------------------------------
                xx = y;
                yy = -x;
                zz = z;
                if (volume_no == 1)
                {
                    yy--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=2 end
                x2 = xx;
                y2 = yy;
                z2 = zz;

                VOLVOXEL(vol, z , y , x) = VOLVOXEL(vol, z0, y0, x0) =
                                               VOLVOXEL(vol, z1, y1, x1) = VOLVOXEL(vol, z2, y2, x2) =
                                                                               (VOLVOXEL(vol, z , y , x) + VOLVOXEL(vol, z0, y0, x0) +
                                                                                VOLVOXEL(vol, z1, y1, x1) + VOLVOXEL(vol, z2, y2, x2)) / 4.0;
            }//for end

}
/* Symmetrizes a simple grid with P4212 symmetry--------------------------*/

void symmetry_P42_12(Image<double> &vol, const SimpleGrid &grid,
                     const Matrix1D<double> &eprm_aint,
                     const Matrix1D<double> &eprm_bint,
                     const MultidimArray<int> &mask, int volume_no,
                     int grid_type)
{

    int ZZ_lowest = (int) ZZ(grid.lowest);
    int YY_lowest = STARTINGY(mask);
    int XX_lowest = STARTINGX(mask);
    int ZZ_highest = (int) ZZ(grid.highest);
    int YY_highest = FINISHINGY(mask);
    int XX_highest = FINISHINGX(mask);

    //if there is an extra slice in the z direction there is no way
    //to calculate the -z slice
    if (ABS(ZZ_lowest) > ABS(ZZ_highest))
        ZZ_lowest = -(ZZ_highest);
    else
        ZZ_highest = ABS(ZZ_lowest);

    while (1)
    {
        if (mask(0, XX_lowest) == 0)
            XX_lowest++;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P42_12, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(0, XX_highest) == 0)
            XX_highest--;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P42_12, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_lowest, 0) == 0)
            YY_lowest++;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P42_12, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_highest, 0) == 0)
            YY_highest--;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P42_12, while(1)" << std::endl;
            exit(0);
        }

    }

    //   int ZZ_lowest =(int) ZZ(grid.lowest);
    //   int YY_lowest =MAX((int) YY(grid.lowest),STARTINGY(mask));
    //   int XX_lowest =MAX((int) XX(grid.lowest),STARTINGX(mask));
    //   int ZZ_highest=(int) ZZ(grid.highest);
    //   int YY_highest=MIN((int) YY(grid.highest),FINISHINGY(mask));
    //   int XX_highest=MIN((int) XX(grid.highest),FINISHINGX(mask));

    int maxZ, maxY, maxX, minZ, minY, minX;
    int x, y, z;

    int x0, y0, z0;
    int x1, y1, z1;
    int x2, y2, z2;
    int x3, y3, z3;
    int x4, y4, z4;
    int x5, y5, z5;
    int x6, y6, z6;

    int XXaint, YYbint;
    XXaint = (int) XX(eprm_aint);
    YYbint = (int)YY(eprm_bint);

    int XXaint_2, YYbint_2;
    XXaint_2 = XXaint / 2;
    YYbint_2 = YYbint / 2;
    int xx, yy, zz;

    if (ABS(XX_lowest) > ABS(XX_highest))
    {
        minX = XX_lowest;
        maxX = 0;
    }
    else
    {
        minX = 0;
        maxX = XX_highest;
    }
    if (ABS(YY_lowest) > ABS(YY_highest))
    {
        minY = YY_lowest;
        maxY = 0;
    }
    else
    {
        minY = 0;
        maxY = YY_highest;
    }
    if (ABS(ZZ_lowest) > ABS(ZZ_highest))
    {
        minZ = ZZ_lowest;
        maxZ = 0;
    }
    else
    {
        minZ = 0;
        maxZ = ZZ_highest;
    }

    //FCC non supported yet
    if (volume_no == 1 && grid_type == FCC)
    {
        std::cerr << "\nSimetries using FCC not implemented\n";
        exit(1);
    }

    for (z = minZ;z <= maxZ;z++)
        for (y = minY;y <= maxY;y++)
            for (x = minX;x <= maxX;x++)
            {
                //sym=-1---------------------------------------------------------
                if (!A2D_ELEM(mask, y, x) || z < ZZ_lowest || z > ZZ_highest)
                    continue;

                //sym=0 ---------------------------------------------------------
                xx = -x;
                yy = -y;
                zz = z;
                if (volume_no == 1)
                {
                    xx--;
                    yy--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }
                x0 = xx;
                y0 = yy;
                z0 = zz;

                //sym=1----------------------------------------------------------
                xx = y;
                yy = x;
                zz = -z;//I think z-- is always inside the grid
                //we do not need to check
                if (volume_no == 1)
                {
                    zz--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=1 end
                x1 = xx;
                y1 = yy;
                z1 = zz;

                //sym=2----------------------------------------------------------
                xx = -y;
                yy = -x;
                zz = -z;
                if (volume_no == 1)
                {
                    xx--;
                    yy--;
                    zz--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=2 end
                x2 = xx;
                y2 = yy;
                z2 = zz;

                //sym=3----------------------------------------------------------
                xx = y + XXaint_2;
                yy = -x + YYbint_2;
                zz = z;
                if (volume_no == 1)
                {
                    yy--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=3 end
                x3 = xx;
                y3 = yy;
                z3 = zz;

                //sym=4----------------------------------------------------------
                xx = -y + XXaint_2;
                yy = + x + YYbint_2;
                zz = z;
                if (volume_no == 1)
                {
                    xx--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=4 end
                x4 = xx;
                y4 = yy;
                z4 = zz;
                //sym=5----------------------------------------------------------
                xx = -x + XXaint_2;
                yy = + y + YYbint_2;
                zz = -z;
                if (volume_no == 1)
                {
                    xx--;
                    zz--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=5 end
                x5 = xx;
                y5 = yy;
                z5 = zz;

                //sym=6----------------------------------------------------------
                xx = + x + XXaint_2;
                yy = -y + YYbint_2;
                zz = -z;
                if (volume_no == 1)
                {
                    yy--;
                    zz--;
                }
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=6 end
                x6 = xx;
                y6 = yy;
                z6 = zz;

                //only the first simple grid center and the origen is the same
                //point.
                //    if(volume_no==1)
                //    {
                //      switch (grid_type) {
                //  case FCC: std::cerr<< "\nSimetries using FCC not implemented\n";break;
                //  case BCC: x0--;y0--;               z1--;
                //            x2--;y2--;z2--;     y3--;
                //     x4--;          x5--;     z5--;
                //          y6--;z6--;
                //            break;
                //  case CC:  break;
                //      }
                //    }

                VOLVOXEL(vol, z , y , x) = VOLVOXEL(vol, z0, y0, x0) =
                                               VOLVOXEL(vol, z1, y1, x1) = VOLVOXEL(vol, z2, y2, x2) =
                                                                               VOLVOXEL(vol, z3, y3, x3) = VOLVOXEL(vol, z4, y4, x4) =
                                                                                                               VOLVOXEL(vol, z5, y5, x5) = VOLVOXEL(vol, z6, y6, x6) =
                                                                                                                                               (VOLVOXEL(vol, z , y , x) + VOLVOXEL(vol, z0, y0, x0) +
                                                                                                                                                VOLVOXEL(vol, z1, y1, x1) + VOLVOXEL(vol, z2, y2, x2) +
                                                                                                                                                VOLVOXEL(vol, z3, y3, x3) + VOLVOXEL(vol, z4, y4, x4) +
                                                                                                                                                VOLVOXEL(vol, z5, y5, x5) + VOLVOXEL(vol, z6, y6, x6)) / 8.0;
            }//for end
}//symmetryP42_12 end
/* Symmetrizes a simple grid with P6 symmetry-----------------------------*/
void symmetry_P6(Image<double> &vol, const SimpleGrid &grid,
                 const Matrix1D<double> &eprm_aint,
                 const Matrix1D<double> &eprm_bint,
                 const MultidimArray<int> &mask, int volume_no,
                 int grid_type)
{

    int ZZ_lowest = (int) ZZ(grid.lowest);
    int YY_lowest = STARTINGY(mask);
    int XX_lowest = STARTINGX(mask);
    int ZZ_highest = (int) ZZ(grid.highest);
    int YY_highest = FINISHINGY(mask);
    int XX_highest = FINISHINGX(mask);


    while (1)
    {
        if (mask(0, XX_lowest) == 0)
            XX_lowest++;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P42_12, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(0, XX_highest) == 0)
            XX_highest--;
        else
            break;
        if (XX_lowest == XX_highest)
        {
            std::cerr << "Error in symmetry_P42_12, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_lowest, 0) == 0)
            YY_lowest++;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P42_12, while(1)" << std::endl;
            exit(0);
        }
    }
    while (1)
    {
        if (mask(YY_highest, 0) == 0)
            YY_highest--;
        else
            break;
        if (YY_lowest == YY_highest)
        {
            std::cerr << "Error in symmetry_P42_12, while(1)" << std::endl;
            exit(0);
        }

    }

    //   int ZZ_lowest =(int) ZZ(grid.lowest);
    //   int YY_lowest =MAX((int) YY(grid.lowest),STARTINGY(mask));
    //   int XX_lowest =MAX((int) XX(grid.lowest),STARTINGX(mask));
    //   int ZZ_highest=(int) ZZ(grid.highest);
    //   int YY_highest=MIN((int) YY(grid.highest),FINISHINGY(mask));
    //   int XX_highest=MIN((int) XX(grid.highest),FINISHINGX(mask));

    int maxZ, maxY, maxX, minZ, minY, minX;
    int x, y, z;

    int x0, y0, z0;
    int x1, y1, z1;
    int x2, y2, z2;
    int x3, y3, z3;
    int x4, y4, z4;

    int XXaint, YYbint;
    XXaint = (int) XX(eprm_aint);
    YYbint = (int)YY(eprm_bint);

    int xx, yy, zz;

    if (ABS(XX_lowest) > ABS(XX_highest))
    {
        minX = XX_lowest;
        maxX = 0;
    }
    else
    {
        minX = 0;
        maxX = XX_highest;
    }
    //P6 is tricky. I have decide to apply it to half the volume
    //instead of to 1 sizth. I think the amount of ifs that I save
    //are worth this larger loop

    minY = YY_lowest;
    maxY = YY_highest;

    minZ = ZZ_lowest;
    maxZ = ZZ_highest;

    //FCC non supported yet
    if (volume_no == 1 && grid_type == FCC)
    {
        std::cerr << "\nSimetries using FCC not implemented\n";
        exit(1);
    }

    for (z = minZ;z <= maxZ;z++)
        for (y = minY;y <= maxY;y++)
            for (x = minX;x <= maxX;x++)
            {
                //sym=-1---------------------------------------------------------
                if (!A2D_ELEM(mask, y, x) || z < ZZ_lowest || z > ZZ_highest)
                    continue;

                //sym=0 ---------------------------------------------------------
                xx = x - y;
                yy = x;
                zz = z;
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }
                x0 = xx;
                y0 = yy;
                z0 = zz;

                //sym=1----------------------------------------------------------
                xx = -y;
                yy = x - y;
                zz = z;
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=1 end
                x1 = xx;
                y1 = yy;
                z1 = zz;

                //sym=2----------------------------------------------------------
                xx = -x;
                yy = -y;
                zz = z;
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=2 end
                x2 = xx;
                y2 = yy;
                z2 = zz;

                //sym=3----------------------------------------------------------
                xx = -x + y;
                yy = -x;
                zz = z;
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=3 end
                x3 = xx;
                y3 = yy;
                z3 = zz;

                //sym=4----------------------------------------------------------
                xx = + y;
                yy = -x + y;
                zz = z;
                if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                {
                    put_inside(xx, XX_lowest, XX_highest, XXaint)
                    put_inside(yy, YY_lowest, YY_highest, YYbint)
                    if (!A2D_ELEM(mask, yy, xx) || mask.outside(yy, xx))
                        std::cerr << "ERROR in symmetry_P function"
                        << "after correction spot is still"
                        << "outside mask\a" << std::endl;
                }//sym=4 end
                x4 = xx;
                y4 = yy;
                z4 = zz;

                if (volume_no == 1)
                {
                    switch (grid_type)
                    {
                    case FCC:
                        std::cerr << "\nSimetries using FCC not implemented\n";
                        break;
                        //there is no way to reinforce P6 in the second grid without
                        //interpolation. This is the best we can do.
                    case BCC:
                        x1 = x0 = x;
                        y1 = y0 = y;
                        x2 = x3 = x4 = -x - 1;
                        y2 = y3 = y4 = -y - 1;
                        break;
                    case CC:
                        break;
                    }
                }


                VOLVOXEL(vol, z , y , x) = VOLVOXEL(vol, z0, y0, x0) =
                                               VOLVOXEL(vol, z1, y1, x1) = VOLVOXEL(vol, z2, y2, x2) =
                                                                               VOLVOXEL(vol, z3, y3, x3) = VOLVOXEL(vol, z4, y4, x4) =
                                                                                                               (VOLVOXEL(vol, z , y , x) + VOLVOXEL(vol, z0, y0, x0) +
                                                                                                                VOLVOXEL(vol, z1, y1, x1) + VOLVOXEL(vol, z2, y2, x2) +
                                                                                                                VOLVOXEL(vol, z3, y3, x3) + VOLVOXEL(vol, z4, y4, x4)) / 6.0;

            }//for end

}
#undef wrap_as_Crystal
#undef DEBUG
/** translate string fn_sym to symmetry group, return false
    is translation is not possible. See URL
    http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
     for details  */

bool SymList::isSymmetryGroup(FileName fn_sym, int &pgGroup, int &pgOrder)
{
    char G1,G2,G3='\0',G4;
    char auxChar[3];
    //each case check lenght, check first letter, second, is number
    //Non a point group

    //remove path
    FileName fn_sym_tmp;
    fn_sym_tmp=fn_sym.removeDirectories();
    int mySize=fn_sym_tmp.size();
    bool return_true;
    return_true=false;
    auxChar[2]='\0';
    //size maybe 4 because n maybe a 2 digit number
    if(mySize>4 || mySize<1)
    {
        pgGroup=-1;
        pgOrder=-1;
        return false;
    }
    //get the group character by character
    G1=toupper((fn_sym_tmp.c_str())[0]);
    G2=toupper((fn_sym_tmp.c_str())[1]);
    if (mySize > 2)
    {
        G3=toupper((fn_sym_tmp.c_str())[2]);
        if(mySize > 3)
            G4=toupper((fn_sym.c_str())[3]);
    }
    else
        G4='\0';
    //CN
    if (mySize==2 && G1=='C' && isdigit(G2))
    {
        pgGroup=pg_CN;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    if (mySize==3 && G1=='C' && isdigit(G2) && isdigit(G3))
    {
        pgGroup=pg_CN;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //CI
    else if (mySize==2 && G1=='C' && G2=='I')
    {
        pgGroup=pg_CI;
        pgOrder=-1;
        return_true=true;
    }
    //CS
    else if (mySize==2 && G1=='C' && G2=='S')
    {
        pgGroup=pg_CS;
        pgOrder=-1;
        return_true=true;
    }
    //CNH
    else if (mySize==3 && G1=='C' && isdigit(G2) && G3=='H')
    {
        pgGroup=pg_CNH;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==4 && G1=='C' && isdigit(G2) && isdigit(G3) && G4=='H')
    {
        pgGroup=pg_CNH;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //CNV
    else if (mySize==3 && G1=='C' && isdigit(G2) && G3=='V')
    {
        pgGroup=pg_CNV;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==4 && G1=='C' && isdigit(G2) && isdigit(G3) && G4=='V')
    {
        pgGroup=pg_CNV;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //SN
    else if (mySize==2 && G1=='S' && isdigit(G2) )
    {
        pgGroup=pg_SN;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==3 && G1=='S' && isdigit(G2) && isdigit(G3) )
    {
        pgGroup=pg_SN;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //DN
    else if (mySize==2 && G1=='D' && isdigit(G2) )
    {
        pgGroup=pg_DN;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    if (mySize==3 && G1=='D' && isdigit(G2) && isdigit(G3))
    {
        pgGroup=pg_DN;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //DNV
    else if (mySize==3 && G1=='D' && isdigit(G2) && G3=='V')
    {
        pgGroup=pg_DNV;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==4 && G1=='D' && isdigit(G2) && isdigit(G3) && G4=='V')
    {
        pgGroup=pg_DNV;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //DNH
    else if (mySize==3 && G1=='D' && isdigit(G2) && G3=='H')
    {
        pgGroup=pg_DNH;
        pgOrder=int(G2)-48;
        return_true=true;
    }
    else if (mySize==4 && G1=='D' && isdigit(G2) && isdigit(G3) && G4=='H')
    {
        pgGroup=pg_DNH;
        auxChar[0]=G2;
        auxChar[1]=G3;
        pgOrder=atoi(auxChar);
        return_true=true;
    }
    //T
    else if (mySize==1 && G1=='T')
    {
        pgGroup=pg_T;
        pgOrder=-1;
        return_true=true;
    }
    //TD
    else if (mySize==2 && G1=='T' && G2=='D')
    {
        pgGroup=pg_TD;
        pgOrder=-1;
        return_true=true;
    }
    //TH
    else if (mySize==2 && G1=='T' && G2=='H')
    {
        pgGroup=pg_TH;
        pgOrder=-1;
        return_true=true;
    }
    //O
    else if (mySize==1 && G1=='O')
    {
        pgGroup=pg_O;
        pgOrder=-1;
        return_true=true;
    }
    //OH
    else if (mySize==2 && G1=='O'&& G2=='H')
    {
        pgGroup=pg_OH;
        pgOrder=-1;
        return_true=true;
    }
    //I
    else if (mySize==1 && G1=='I')
    {
        pgGroup=pg_I;
        pgOrder=-1;
        return_true=true;
    }
    //I1
    else if (mySize==2 && G1=='I'&& G2=='1')
    {
        pgGroup=pg_I1;
        pgOrder=-1;
        return_true=true;
    }
    //I2
    else if (mySize==2 && G1=='I'&& G2=='2')
    {
        pgGroup=pg_I2;
        pgOrder=-1;
        return_true=true;
    }
    //I3
    else if (mySize==2 && G1=='I'&& G2=='3')
    {
        pgGroup=pg_I3;
        pgOrder=-1;
        return_true=true;
    }
    //I4
    else if (mySize==2 && G1=='I'&& G2=='4')
    {
        pgGroup=pg_I4;
        pgOrder=-1;
        return_true=true;
    }
    //I5
    else if (mySize==2 && G1=='I'&& G2=='5')
    {
        pgGroup=pg_I5;
        pgOrder=-1;
        return_true=true;
    }
    //IH
    else if (mySize==2 && G1=='I'&& G2=='H')
    {
        pgGroup=pg_IH;
        pgOrder=-1;
        return_true=true;
    }
    //I1H
    else if (mySize==3 && G1=='I'&& G2=='1'&& G3=='H')
    {
        pgGroup=pg_I1H;
        pgOrder=-1;
        return_true=true;
    }
    //I2H
    else if (mySize==3 && G1=='I'&& G2=='2'&& G3=='H')
    {
        pgGroup=pg_I2H;
        pgOrder=-1;
        return_true=true;
    }
    //I3H
    else if (mySize==3 && G1=='I'&& G2=='3'&& G3=='H')
    {
        pgGroup=pg_I3H;
        pgOrder=-1;
        return_true=true;
    }
    //I4H
    else if (mySize==3 && G1=='I'&& G2=='4'&& G3=='H')
    {
        pgGroup=pg_I4H;
        pgOrder=-1;
        return_true=true;
    }
    //I5H
    else if (mySize==3 && G1=='I'&& G2=='5'&& G3=='H')
    {
        pgGroup=pg_I5H;
        pgOrder=-1;
        return_true=true;
    }
    //#define DEBUG7
#ifdef DEBUG7
    std::cerr << "pgGroup" << pgGroup << " pgOrder " << pgOrder << std::endl;
#endif
#undef DEBUG7

    return return_true;
}
void SymList::fillSymmetryClass(const FileName &symmetry, int pgGroup, int pgOrder,
                                std::vector<std::string> &fileContent)
{
    std::ostringstream line1;
    std::ostringstream line2;
    std::ostringstream line3;
    std::ostringstream line4;
    if (pgGroup == pg_CN)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
    }
    else if (pgGroup == pg_CI)
    {
        line1 << "inversion ";
    }
    else if (pgGroup == pg_CS)
    {
        line1 << "mirror_plane 0 0 1";
    }
    else if (pgGroup == pg_CNV)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "mirror_plane 0 1 0";
    }
    else if (pgGroup == pg_CNH)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "mirror_plane 0 0 1";
    }
    else if (pgGroup == pg_SN)
    {
        int order = pgOrder / 2;
        if(2*order != pgOrder)
        {
            std::cerr << "ERROR: order for SN group must be even" << std::endl;
            exit(0);
        }
        line1 << "rot_axis " << order << " 0 0 1";
        line2 << "inversion ";
    }
    else if (pgGroup == pg_DN)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "rot_axis " << "2" << " 1 0 0";
    }
    else if (pgGroup == pg_DNV)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "rot_axis " << "2" << " 1 0 0";
        line3 << "mirror_plane 1 0 0";
    }
    else if (pgGroup == pg_DNH)
    {
        line1 << "rot_axis " << pgOrder << " 0 0 1";
        line2 << "rot_axis " << "2" << " 1 0 0";
        line3 << "mirror_plane 0 0 1";
    }
    else if (pgGroup == pg_T)
    {
        line1 << "rot_axis " << "3" << "  0. 0. 1.";
        line2 << "rot_axis " << "2" << " 0. 0.816496 0.577350";
    }
    else if (pgGroup == pg_TD)
    {
        line1 << "rot_axis " << "3" << "  0. 0. 1.";
        line2 << "rot_axis " << "2" << " 0. 0.816496 0.577350";
        line3 << "mirror_plane 1.4142136 2.4494897 0.0000000";
    }
    else if (pgGroup == pg_TH)
    {
        line1 << "rot_axis " << "3" << "  0. 0. 1.";
        line2 << "rot_axis " << "2" << " 0. -0.816496 -0.577350";
        line3 << "inversion";
    }
    else if (pgGroup == pg_O)
    {
        line1 << "rot_axis " << "3" << "  .5773502  .5773502 .5773502";
        line2 << "rot_axis " << "4" << " 0 0 1";
    }
    else if (pgGroup == pg_OH)
    {
        line1 << "rot_axis " << "3" << "  .5773502  .5773502 .5773502";
        line2 << "rot_axis " << "4" << " 0 0 1";
        line3 << "mirror_plane 0 1 1";
    }
    else if (pgGroup == pg_I || pgGroup == pg_I2)
    {
        line1 << "rot_axis 2  0 0 1";
        line2 << "rot_axis 5  0.525731114  0 0.850650807";
        line3 << "rot_axis 3  0 0.356822076 0.934172364";
    }
    else if (pgGroup == pg_I1)
    {
        line1 << "rot_axis 2  1      0        0";
        line2 << "rot_axis 5 0.85065080702670 0 -0.5257311142635";
        line3 << "rot_axis 3 0.9341723640 0.3568220765 0";
    }
    else if (pgGroup == pg_I3)
    {
        line1 << "rot_axis 2  -0.5257311143 0 0.8506508070";
        line3 << "rot_axis 5  0. 0. 1.";
        line2 << "rot_axis 3  -0.4911234778630044, 0.3568220764705179, 0.7946544753759428";
    }
    else if (pgGroup == pg_I4)
    {
        line1 << "rot_axis 2  0.5257311143 0 0.8506508070";
        line3 << "rot_axis 5  0.8944271932547096 0 0.4472135909903704";
        line2 << "rot_axis 3  0.4911234778630044 0.3568220764705179 0.7946544753759428";
    }
    else if (pgGroup == pg_I5)
    {
        std::cerr << "ERROR: Symmetry pg_I5 not implemented" << std::endl;
        exit(0);
    }
    else if (pgGroup == pg_IH || pgGroup == pg_I2H)
    {
        line1 << "rot_axis 2  0 0 1";
        line2 << "rot_axis 5  0.525731114  0 0.850650807";
        line3 << "rot_axis 3  0 0.356822076 0.934172364";
        line4 << "mirror_plane 1 0 0";
    }
    else if (pgGroup == pg_I1H)
    {
        line1 << "rot_axis 2  1      0        0";
        line2 << "rot_axis 5 0.85065080702670 0 -0.5257311142635";
        line3 << "rot_axis 3 0.9341723640 0.3568220765 0";
        line4 << "mirror_plane 0 0 -1";
    }
    else if (pgGroup == pg_I3H)
    {
        line1 << "rot_axis 2  -0.5257311143 0 0.8506508070";
        line3 << "rot_axis 5  0. 0. 1.";
        line2 << "rot_axis 3  -0.4911234778630044, 0.3568220764705179, 0.7946544753759428";
        line4 << "mirror_plane 0.850650807 0  0.525731114";
    }
    else if (pgGroup == pg_I4H)
    {
        line1 << "rot_axis 2  0.5257311143 0 0.8506508070";
        line3 << "rot_axis 5  0.8944271932547096 0 0.4472135909903704";
        line2 << "rot_axis 3  0.4911234778630044 0.3568220764705179 0.7946544753759428";
        line4 << "mirror_plane 0.850650807 0 -0.525731114";
    }
    else if (pgGroup == pg_I5H)
    {
        std::cerr << "ERROR: Symmetry pg_I5H not implemented" << std::endl;
        exit(0);
    }
    else
    {
        std::cerr << "ERROR: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);
    }
    if (line1.str().size()>0)
        fileContent.push_back(line1.str());
    if (line2.str().size()>0)
        fileContent.push_back(line2.str());
    if (line3.str().size()>0)
        fileContent.push_back(line3.str());
    if (line4.str().size()>0)
        fileContent.push_back(line4.str());
    //#define DEBUG5
#ifdef DEBUG5

    for (int n=0; n<fileContent.size(); n++)
        std::cerr << fileContent[n] << std::endl;
    std::cerr << "fileContent.size()" << fileContent.size() << std::endl;
#endif
    #undef DEBUG5
}
double SymList::nonRedundantProjectionSphere(int pgGroup, int pgOrder)
{
    if (pgGroup == pg_CN)
    {
        return 4.*PI/pgOrder;
    }
    else if (pgGroup == pg_CI)
    {
        return 4.*PI/2.;
    }
    else if (pgGroup == pg_CS)
    {
        return 4.*PI/2.;
    }
    else if (pgGroup == pg_CNV)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_CNH)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_SN)
    {
        return 4.*PI/pgOrder;
    }
    else if (pgGroup == pg_DN)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_DNV)
    {
        return 4.*PI/pgOrder/4;
    }
    else if (pgGroup == pg_DNH)
    {
        return 4.*PI/pgOrder/4;
    }
    else if (pgGroup == pg_T)
    {
        return 4.*PI/12;
    }
    else if (pgGroup == pg_TD)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_TH)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_O)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_OH)
    {
        return 4.*PI/48;
    }
    else if (pgGroup == pg_I || pgGroup == pg_I2)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I1)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I3)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I4)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I5)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_IH || pgGroup == pg_I2H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I1H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I3H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I4H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I5H)
    {
        return 4.*PI/120;
    }
    else
    {
        std::cerr << "ERROR: Symmetry group, order=" << pgGroup
        << " "
        <<  pgOrder
        << "is not known"
        << std::endl;
        exit(0);
    }
}

void SymList::computeDistance(MetaData &md,
                              bool projdir_mode, bool check_mirrors,
                              bool object_rotation)
{
    MDRow row;
    double rot1, tilt1, psi1;
    double rot2, tilt2, psi2;
    double angDistance;
    FOR_ALL_OBJECTS_IN_METADATA(md)
    {
        md.getRow(row,__iter.objId);

        row.getValue(MDL_ANGLE_ROT,rot1);
        row.getValue(MDL_ANGLE_ROT2,rot2);

        row.getValue(MDL_ANGLE_TILT,tilt1);
        row.getValue(MDL_ANGLE_TILT2,tilt2);

        row.getValue(MDL_ANGLE_PSI,psi1);
        row.getValue(MDL_ANGLE_PSI2,psi2);

        angDistance=computeDistance( rot1,  tilt1,  psi1,
                                     rot2,  tilt2,  psi2,
                                     projdir_mode,  check_mirrors,
                                     object_rotation);

        md.setValue(MDL_ANGLE_ROT_DIFF,rot1 - rot2,__iter.objId);
        md.setValue(MDL_ANGLE_TILT_DIFF,tilt1 - tilt2,__iter.objId);
        md.setValue(MDL_ANGLE_PSI_DIFF,psi1 - psi2,__iter.objId);
        md.setValue(MDL_ANGLE_DIFF,angDistance,__iter.objId);
    }

}
double SymList::computeDistance(double rot1, double tilt1, double psi1,
                                double &rot2, double &tilt2, double &psi2,
                                bool projdir_mode, bool check_mirrors,
                                bool object_rotation)
{
    Matrix2D<double> E1, E2;
    Euler_angles2matrix(rot1, tilt1, psi1, E1, false);

    int imax = symsNo() + 1;
    Matrix2D<double>  L(3, 3), R(3, 3);  // A matrix from the list
    double best_ang_dist = 3600;
    double best_rot2=0, best_tilt2=0, best_psi2=0;

    for (int i = 0; i < imax; i++)
    {
        double rot2p, tilt2p, psi2p;
        if (i == 0)
        {
            rot2p = rot2;
            tilt2p = tilt2;
            psi2p = psi2;
        }
        else
        {
            getMatrices(i - 1, L, R, false);
            if (object_rotation)
                Euler_apply_transf(R, L, rot2, tilt2, psi2, rot2p, tilt2p, psi2p);
            else
                Euler_apply_transf(L, R, rot2, tilt2, psi2, rot2p, tilt2p, psi2p);
        }

        double ang_dist = Euler_distanceBetweenAngleSets_fast(E1,rot2p, tilt2p, psi2p,
                          projdir_mode, E2);

        if (ang_dist < best_ang_dist)
        {
            best_rot2 = rot2p;
            best_tilt2 = tilt2p;
            best_psi2 = psi2p;
            best_ang_dist = ang_dist;
        }

        if (check_mirrors)
        {
            Euler_up_down(rot2p, tilt2p, psi2p, rot2p, tilt2p, psi2p);
            double ang_dist_mirror = Euler_distanceBetweenAngleSets_fast(E1,
                                     rot2p, tilt2p, psi2p,projdir_mode, E2);

            if (ang_dist_mirror < best_ang_dist)
            {
                best_rot2 = rot2p;
                best_tilt2 = tilt2p;
                best_psi2 = psi2p;
                best_ang_dist = ang_dist_mirror;
            }

        }
    }
    rot2 = best_rot2;
    tilt2 = best_tilt2;
    psi2 = best_psi2;
    return best_ang_dist;
}

void symmetry_Helical(MultidimArray<double> &Vout, const MultidimArray<double> &Vin, double zHelical, double rotHelical,
		double rot0, MultidimArray<int> *mask)
{
	Vout.initZeros(Vin);
	double izHelical=1.0/zHelical;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(Vin)
	{
		if (mask!=NULL && !A3D_ELEM(*mask,k,i,j))
			continue;
		double rot=atan2((double)i,(double)j)+rot0;
		double rho=sqrt((double)i*i+(double)j*j);
		double l0=ceil((STARTINGZ(Vin)-k)*izHelical);
		double lF=floor((FINISHINGZ(Vin)-k)*izHelical);
		double finalValue=0;
		double L=0;
		for (double l=l0; l<=lF; ++l)
		{
			double kp=k+l*zHelical;
			double rotp=rot+l*rotHelical;
			double ip, jp;
			sincos(rotp,&ip,&jp);
			ip*=rho;
			jp*=rho;
			finalValue+=Vin.interpolatedElement3D(jp,ip,kp,0.0);
			L+=1.0;
		}
		A3D_ELEM(Vout,k,i,j)=finalValue/L;
	}
}
