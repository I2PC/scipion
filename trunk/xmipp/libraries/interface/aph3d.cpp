/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Roberto Marabini (roberto@mipg.upenn.edu)
 *
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "aph3d.h"

#include <data/args.h>

#include <fstream>

#define VERBOSE
//#define DEBUG
// APH =====================================================================
void APHFile3D::read_from_prepmklcf(const FileName &fn)
{
    ifstream  fh_aph;
    int       line_no = 1;
    int       hmax = 0, kmax = 0, lmax = 0, hmin = 0, kmin = 0, lmin = 0;
    string    line;
    int       h, k, l;

    // Empties current APH File
    clear();
    // Open file
    fh_aph.open(fn.c_str(), ios::in);
    if (!fh_aph)
        REPORT_ERROR(1601, "aphFile::read: File " + fn + " not found");

    // Read first line and skip it
//   fh_aph.peek();
    getline(fh_aph, line);
    try
    {
        getline(fh_aph, line);
        Reduce_phase_accuracy = AtoF(line.c_str() + 47);
        getline(fh_aph, line);
        a = AtoF(line.c_str() + 47);
        getline(fh_aph, line);
        b = AtoF(line.c_str() + 47);
        getline(fh_aph, line);
        gamma = AtoF(line.c_str() + 47);
        getline(fh_aph, line);
        c = AtoF(line.c_str() + 47);
        getline(fh_aph, line);
        resolution = AtoF(line.c_str() + 47);
        getline(fh_aph, line);
        Amp_Scale = AtoF(line.c_str() + 47);
    }
    catch (...)
    {
        REPORT_ERROR(1601, "aph3DFile::read: Wrong first line in file " + fn);
    }

//#define DEBUG_read
#ifdef DEBUG_read
    cout << "Reduce_phase_accuracy: " << Reduce_phase_accuracy << endl;
    cout << "a: " << a << endl;
    cout << "b: " << b << endl;
    cout << "gamma: " << gamma << endl;
    cout << "c: " << c << endl;
    cout << "resolution: " << resolution << endl;
    cout << "Amp_Scale: " << Amp_Scale << endl;
#endif
#undef DEBUG_read

    // look for the begining of the data
    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (string::npos !=
                line.find("   H   K   L      A      P     FOM*100         REJECTS"))
                break;
        }
        catch (Xmipp_error)
        {
            cout << "3Daph File reading error an error\n";
        }
    }/* while */

    // look for maximun and minimum
    line_no = 1;

    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (string::npos != line.find("                                 "))
                continue;
            else if (string::npos != line.find("MKLCF FILE COMPLETED"))
                break;
            h = AtoI(firstToken(line));
            k = AtoI(nextToken());
            l = AtoI(nextToken());
            hmax = MAX(hmax, h);
            kmax = MAX(kmax, k);
            lmax = MAX(lmax, l);
            hmin = MIN(hmin, h);
            kmin = MIN(kmin, k);
            lmin = MIN(lmin, l);
        }
        catch (Xmipp_error)
        {
            cout << "aph File: Line " << line_no << " is skipped due to an error\n";
        }
        line_no++;
    }/* while */
//Space Group
    Space_Group = 0;
    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (string::npos !=
                line.find(" * Space Group ="))
            {
                Space_Group = AtoI(line.c_str() + 16);
                break;
            }
        }
        catch (Xmipp_error)
        {
            cout << "3Daph File reading error an error\n";
        }
    }/* while */

    if (Space_Group == 0)
        REPORT_ERROR(1601, "aphFile::read: File " + fn + " Space Group not found");

#define DEBUG_max
#ifdef DEBUG_max
    cout << "hmax: " << hmax << " kmax: " << kmax << " lmax: " << lmax << endl;
    cout << "hmin: " << hmin << " kmin: " << kmin << " lmin: " << lmin << endl;
#endif
#undef DEBUG_max

    // Ask for memory
    spots_abs.initZeros(lmax - lmin + 1, kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(spots_abs) = hmin;
    STARTINGY(spots_abs) = kmin;
    STARTINGZ(spots_abs) = lmin;
    spots_arg.initZeros(lmax - lmin + 1, kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(spots_arg) = hmin;
    STARTINGY(spots_arg) = kmin;
    STARTINGZ(spots_arg) = lmin;
    FOM.initZeros(lmax - lmin + 1, kmax - kmin + 1, hmax - hmin + 1);
    STARTINGX(FOM) = hmin;
    STARTINGY(FOM) = kmin;
    STARTINGZ(FOM) = lmin;

    // Read each line (again) and copy values to the matrices
    fh_aph.close();
    fh_aph.open(fn.c_str(), ios::in);
    line_no = 1;

    // look for the begining of the data
    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (string::npos !=
                line.find("   H   K   L      A      P     FOM*100         REJECTS"))
                break;
        }
        catch (Xmipp_error)
        {
            cout << "3Daph File reading error an error\n";
        }
    }/* while */

    while (!fh_aph.eof())
    {
        try
        {
            getline(fh_aph, line);
            if (string::npos != line.find("                                 "))
                continue;
            else if (string::npos != line.find("MKLCF FILE COMPLETED"))
                break;
            h     = AtoI(firstToken(line));
            k     = AtoI(nextToken());
            l     = AtoI(nextToken());
            spots_abs(l, k, h)  = AtoF(nextToken());
            spots_arg(l, k, h)  = AtoF(nextToken());
            FOM(l, k, h)       = AtoI(nextToken());
            switch (Space_Group)
            {
            case(1):
                            break;
            case(90)://P4212
                            if (h > k || h < 0 || l < 0)
                    {
                        cerr << "\nHORROR reflection outside the assymetric unit\n"
                        << "(h,k,l)=" << h << " " << k << " " << l << endl;
                        exit(1);
                        break;
                    }
            }//switch end
        }
        catch (Xmipp_error XE)
        {
            cout << XE;
            cout << "aph File: Line " << line_no << " is skipped due to an error\n";
        }
        line_no++;
    }/* while */
    // Close file
    fh_aph.close();

    fn_aph = fn;

}/*  APHFile2D::read */


/* ------------------------------------------------------------------------- */
void APHFile3D::clear()
{
    a = 0.;
    b = 0.;
    gamma = 0.;
    c = 0.;
    resolution = 0.;
    Amp_Scale = 0.;
    spots_abs.clear();
    spots_arg.clear();
    FOM.clear();
    Space_Group = 0;
    FOM.clear();
    fn_aph = "";
} /*clear*/
