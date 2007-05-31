/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Roberto Marabini (roberto@mipg.upenn.edu)
 *
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
 *                                      <
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "aph.h"

#include <reconstruction/symmetries.h>
#include <data/args.h>
#include <data/image.h>

#include <fstream>
#include <iomanip>

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
#include <sstream>
#else
#include <strstream.h>
#endif

#define VERBOSE
//#define DEBUG
// APH =====================================================================
void APHFile2D::read(const FileName &fn)
{
    ifstream  fh_aph;
    int       line_no = 1;
    int       hmax = 0, kmax = 0, hmin = 0, kmin = 0;
    string    line;

    // Empties current APH File
    clear();

    // Open file
    fh_aph.open(fn.c_str(), ios::in);
    if (!fh_aph)
        REPORT_ERROR(1601, "aphFile::read: File " + fn + " not found");

    // Read first line and save as title
    fh_aph.peek();
    getline(fh_aph, line);
    astar.resize(2);
    bstar.resize(2);
#if GCC_VERSION < 30300
    istrstream is(line.c_str());
#else
    istringstream is(line.c_str());
#endif
    try
    {
        is >> label >> XX(astar) >> YY(astar) >> XX(bstar) >> YY(bstar)
        >> Xdim >> Ydim >> sampling_rate;
    }
    catch (...)
    {
        REPORT_ERROR(1601, "aphFile::read: Wrong first line in file " + fn);
    }

    // Read each line and keep it in the list of the aphFile object
    fh_aph.peek();
    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (line.length() != 0)
            {
                int h = AtoI(firstToken(line));
                int k = AtoI(nextToken());
                hmax = MAX(hmax, h);
                kmax = MAX(kmax, k);
                hmin = MIN(hmin, h);
                kmin = MIN(kmin, k);
            }
        }
        catch (Xmipp_error)
        {
            cout << "aph File: Line " << line_no << " is skipped due to an error\n";
        }
        line_no++;
    }/* while */
#ifdef DEBUG
    cout << "hmax: " << hmax << " kmax: " << kmax << endl;
    cout << "hmin: " << hmin << " kmin: " << kmin << endl;
#endif

    // Ask for memory
    spots_l.initZeros(kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(spots_l) = hmin;
    STARTINGY(spots_l) = kmin;
    spots_abs.initZeros(kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(spots_abs) = hmin;
    STARTINGY(spots_abs) = kmin;
    spots_arg.initZeros(kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(spots_arg) = hmin;
    STARTINGY(spots_arg) = kmin;
    IQ.initZeros(kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(IQ) = hmin;
    STARTINGY(IQ) = kmin;
    background.initZeros(kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(background) = hmin;
    STARTINGY(background) = kmin;
    CTF.initZeros(background);

    // Read each line (again) and copy values to the matrices
    fh_aph.close();
    fh_aph.clear();
    fh_aph.open(fn.c_str(), ios::in);
    line_no = 1;

    // Read first line and skip it
    getline(fh_aph, line);
    l_is_present = false;
    bool first = true;
    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (line.length() != 0)
            {
                int i = 0;
                int   h         = AtoI(nextToken(line, i));
                int   k         = AtoI(nextToken(line, i));
                double a1, a2, a3, a4, a5, a6;
                a1 = AtoF(nextToken(line, i));
                a2 = AtoF(nextToken(line, i));
                a3 = AtoF(nextToken(line, i));
                a4 = AtoF(nextToken(line, i));
                a5 = AtoF(nextToken(line, i));
                string aux;
                aux = nextToken(line, i);
                if (first)
                {
                    l_is_present = (aux != "");
                    first = false;
                }

                if (l_is_present)
                {
                    a6 = AtoF(aux);
                    spots_l(k, h)    = a1;
                    spots_abs(k, h)  = a2;
                    spots_arg(k, h)  = a3;
                    IQ(k, h)         = (int) a4;
                    background(k, h) = a5;
                    CTF(k, h)        = a6;
                }
                else
                {
                    spots_abs(k, h)  = a1;
                    spots_arg(k, h)  = a2;
                    IQ(k, h)         = (int) a3;
                    background(k, h) = a4;
                    CTF(k, h)        = a5;
                }
//  #define DEBUG
#ifdef DEBUG
                cout << " " << h << " " << k << " " << a1 << " " << a2 << " ";
                cout << a3 << " " << a4 << " " << a5 << " " << a6 << endl;
#endif
//  #undef DEBUG
            }
        }
        catch (Xmipp_error XE)
        {
            cout << XE;
            cout << "aph File: Line " << line_no << " is skipped due to an error\n";
        }
        line_no++;
        fh_aph.peek();
    }/* while */
    // Close file
    fh_aph.close();

    fn_aph = fn;
}/*  APHFile2D::read */

/* ------------------------------------------------------------------------- */
void APHFile2D::write(const FileName &fn) const
{
    ofstream fh;
    fh.open(fn.c_str());
    char aux_char[128];
    if (!fh)
        REPORT_ERROR(1, (string)"APHFile2D::write: Cannot open " +
                     fn + " for output");

    fh << setfill('0') << setw(4) << label << " "
    << XX(astar) << " " << YY(astar) << " "
    << XX(bstar) << " " << YY(bstar) << " "
    << Xdim << " " << Ydim << " " << sampling_rate << endl;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(spots_abs)
    {
        if (spots_abs(i, j) != 0)
        {
            if (l_is_present)
                sprintf(aux_char, " %4d%4d%8.4f%10.1f%7.1f%7d%3d%8.5f%10.1f%7.3f\n",
                        j, i, spots_l(i, j), spots_abs(i, j), spots_arg(i, j),
                        label, IQ(i, j), background(i, j), CTF(i, j));
            else
                sprintf(aux_char, "%4d%4d%8.1f%8.1f%3d%8.1f%8.1f\n",
                        j, i, spots_abs(i, j), spots_arg(i, j),
                        IQ(i, j), background(i, j), CTF(i, j));
            fh << aux_char;
        }
    }
    fh.close();
}

/* ------------------------------------------------------------------------- */
void APHFile2D::clear()
{
    astar.clear();
    bstar.clear();
    Xdim = Ydim = 0;
    sampling_rate = 0;
    spots_abs.clear();
    spots_arg.clear();
    IQ.clear();
    background.clear();
    CTF.clear();
} /*clear*/

/*----transform Xmipp Euler angles into MRC angles *--------------*/
void Euler_to_MRC(double rot, double tilt, double psi,
                  double * mrc_tilt, double * mrc_taxa)
{


    EULER_CLIPPING_RAD(rot, tilt, psi);
    if (tilt == 0) rot = 0;//rot irrelevant if tilt = 0
    else if (rot <= PI / 2.)   *mrc_taxa = PI / 2 - rot;
    else if (rot < PI*3. / 2.) *mrc_taxa = PI * 3. / 2 - rot;
    else if (rot < PI*2.)    *mrc_taxa = PI * 5. / 2 - rot;
    else
    {
        cerr << "\nHORROR: (Euler_to_MRC) Can't find taxa\n)";
        exit(1);
    }
    if ((rot < PI / 2 + 0.1   && rot > PI / 2 - .1) ||
        (rot < PI*3 / 2 + 0.1 && rot > PI*3 / 2 - .1) ||
        (rot < 0.1) ||
        (rot > (PI*2) - .1))
        cerr << "\nWARMING, rot close 0,90,270 or 360 degrees, conversion not reliable\n";


    if ((SGN(tilt) == + 1) && (rot <= PI*3. / 2. && rot > PI / 2.))
    {
        cout << "one\n";
        *mrc_tilt = -tilt;
    }//nrg
    else if ((SGN(tilt) == -1) && (rot > PI*3. / 2.))
    {
        cout << "two\n";
        *mrc_tilt = tilt;
    }//neg
    else if ((SGN(tilt) == -1) && (rot <= PI / 2.))
    {
        cout << "three\n";
        *mrc_tilt = tilt;
    }//neg
    else if ((SGN(tilt) == + 1) && (rot > PI*3. / 2.))
    {
        cout << "four\n";
        *mrc_tilt = tilt;
    }//plus
    else if ((SGN(tilt) == -1) && (rot <= PI*3. / 2. && rot > PI / 2.))
    {
        cout << "five\n";
        *mrc_tilt = -tilt;
    }//plus
    else if ((SGN(tilt) == + 1) && (rot <= PI / 2.))
    {
        cout << "six\n";
        *mrc_tilt = + tilt;
    } //plu
    else
    {
        cerr << "\nHORROR: (Euler_to_MRC) Can't find tilt\n)";
        exit(1);
    }

//cout << "\nDEBUG rot: " << rot<<endl;
//cout << "\nDEBUG tilt: " << tilt<<endl;
//cout << "\nDEBUG psi: " << psi<<endl;
//cout << "\nDEBUG *mrc_tilt: " << *mrc_tilt<<endl;

}

/*----transform Xmipp Euler angles into MRC angles *--------------*/
void MRC_to_Euler(double  mrc_taxa, double  mrc_tilt,
                  double *rot, double *tilt, double *psi)
{
    *psi = 0;
    if (mrc_tilt == 0)
    {
        *tilt = 0.;
        *rot = 0.;
    }//taxa irrelevant if tilt = 0
    else if (mrc_taxa <= PI / 2.)
    {
        *rot = PI / 2 - mrc_taxa;
        *tilt = mrc_tilt;
    }
    else if (mrc_taxa < PI*3. / 2.)
    {
        *rot = PI * 3. / 2 - mrc_taxa;
        *tilt = - mrc_tilt;
    }
    else
    {
        cerr << "\nHORROR: (MRC_to_Euler) taxa bigger than 180\n)";
        exit(1);
    }
    if (*tilt != 0 && ((mrc_taxa < PI / 2 + 0.1   && mrc_taxa > PI / 2 - .1) ||
                       (mrc_taxa < PI + 0.1 && mrc_taxa > PI - .1) ||
                       (mrc_taxa < 0.1)))
        cerr << "\nWARNING, taxa close 0,90 or 180 degrees, conversion not reliable\n";
}

void APHFile2D::copy_reflection(int h, int k, int new_h, int new_k,
                                double new_l, bool conjugate, int sign)
{
    spots_abs(new_k, new_h) = spots_abs(k, h);
    spots_l(new_k, new_h) = new_l;
    double conj = 1;
    if (conjugate) conj = -1;
    if (sign % 2 == 1) spots_arg(new_k, new_h) = conj * spots_arg(k, h) + 180;
    else           spots_arg(new_k, new_h) = conj * spots_arg(k, h);
    IQ(new_k, new_h) = IQ(k, h);
    background(new_k, new_h) = background(k, h);
    CTF(new_k, new_h) = CTF(k, h);
}
#ifdef DISCONTINUATED

/* Generate reflections ---------------------------------------------------- */
void APHFile2D::generate_symmetrical_reflections(int symmetry_group)
{
    // Resize the spot matrices
    int new_kmin = MIN(-FINISHINGY(spots_abs), STARTINGY(spots_abs));
    int new_hmin = MIN(-FINISHINGX(spots_abs), STARTINGX(spots_abs));
    int new_kmax = MAX(-STARTINGY(spots_abs) , FINISHINGY(spots_abs));
    int new_hmax = MAX(-STARTINGX(spots_abs) , FINISHINGX(spots_abs));

    spots_l   .window(new_kmin, new_hmin, new_kmax, new_hmax);
    spots_abs .window(new_kmin, new_hmin, new_kmax, new_hmax);
    spots_arg .window(new_kmin, new_hmin, new_kmax, new_hmax);
    IQ        .window(new_kmin, new_hmin, new_kmax, new_hmax);
    background.window(new_kmin, new_hmin, new_kmax, new_hmax);
    CTF       .window(new_kmin, new_hmin, new_kmax, new_hmax);

    Matrix2D<int> visited;
    visited.initZeros(YSIZE(spots_abs), XSIZE(spots_abs));
    STARTINGY(visited) = STARTINGY(spots_abs);
    STARTINGX(visited) = STARTINGX(spots_abs);

    // Generate symmetrical points
    // Symmetrical points are generated through the multiplication
    //                               [R[0] R[2]  0    0 R[5]]
    //                               [R[1] R[3]  0    0 R[6]]
    // [h' k' l' A' ph']=[h k l A ph][ 0    0   R[4]  0  0  ]
    //                               [ 0    0    0    1  0  ]
    //                               [ 0    0    0    0 R[7]]
    //
    // These R matrices are obtained for each crystallographic group
    // They can be easily obtaine from  origtiltd.for (MRC source code)
    int mode = 1;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(visited)
    if (!visited(i, j))
    {
        visited(i, j) = 1;
        if (spots_abs(i, j) != 0)
        {
            int h = j;
            int k = i;
            switch (symmetry_group)
            {
            case sym_P222_1:
                /* Possible Rs for this group
                          h > 0
                I:[ 1 0 0  1  1 0    0  1] -> ( h, k, l)=(h,k,l)
                R2:[ 1 0 0  1 -1 0    0 -1] -> ( h, k,-l)=conj(h,k,l)
                R3:[ 1 0 0 -1  1 0  180 -1] -> ( h,-k, l)=conj(h,k,l)*(-1)^h
                R3.R2:[ 1 0 0 -1 -1 0 -180  1] -> ( h,-k,-l)=(h,k,l)*(-1)^h
                          h < 0
                R1:[-1 0 0 -1 -1 0    0 -1] -> (-h,-k,-l)=conj(h,k,l)
                R1.R2:[-1 0 0 -1  1 0    0  1] -> (-h,-k, l)=(h,k,l)
                R1.R3:[-1 0 0  1 -1 0 -180  1] -> (-h, k,-l)=(h,k,l)*(-1)^h
                R1.R2.R3:[-1 0 0  1  1 0 -180 -1] -> (-h, k, l)=conj(h,k,l)*(-1)^h

                The special cases of this group are
                Real         Imag
                ----         ----
                (0,2n,l)     (0,2n+1,l)
                (h,k,0)
                (h,0,l)

                I will only use those restrictions valid for all l.
                */
                if ((h == 0 && k % 2 == 0) || (k == 0))
                {
                    double c = cos(DEG2RAD(spots_arg(k, h)));
                    if (c < 0) spots_arg(k, h) = 180;
                    else     spots_arg(k, h) = 0;
                }
                else if (h == 0 && k % 2 == 1)
                {
                    if (spots_l(k, h) != 0)
                    {
                        double s = sin(DEG2RAD(spots_arg(k, h)));
                        if (s < 0) spots_arg(k, h) = -90;
                        else     spots_arg(k, h) = 90;
                    }
                    else if (l_is_present)
                    {
                        spots_arg(k, h) = 0;
                        spots_arg(k, h) = 0;
                    }
                }
// (h,k,0)== real missing
                double l = spots_l(k, h);

                //copy_reflection( h, k,-h,-k,-l,true,0); visited(-k,-h)=1;
                if (mode == 0)
                {
                    //1:copy_reflection( h, k,-h, k, l,true,h); visited( k,-h)=1;
                    //1:copy_reflection(-h, k, h,-k,-l,true,0); visited(-k, h)=1;
                    copy_reflection(h, k, h, -k, -l, false, h);
                    visited(-k, h) = 1;
                    //copy_reflection(h,-k,-h, k, l,true,0); visited( k,-h)=1;
                }
                else
                {
                    copy_reflection(h, k, h, -k, l, false, h);
                    visited(-k, h) = 1;
                    //copy_reflection( h,-k,-h, k,-l,true,0); visited( k,-h)=1;
                    //1:copy_reflection( h,-k,-h, k,-l,true,0); visited( k,-h)=1;
                }
                break;
            }
#ifdef DEBUG
            ImageXmipp save;
            save() = spots_abs;
            save.write("PPPSpots_abs.xmp");
            save() = spots_arg;
            save.write("PPPSpots_arg.xmp");
            type_cast(visited, save());
            save.write("PPPvisited.xmp");
            cout << "Press\n";
            char c;
            cin >> c;
#endif
        }
    }
}
#endif

