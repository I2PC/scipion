/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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

#include "euler.h"

//#define VERBOSE

/* Main routine  ------------------------------------------- */

void ROUT_EULER(const double rot,
                const double tilt,
                const double psi,
                const double d00,
                const double d01,
                const double d10,
                const double d11)
{
    int i;
    Matrix1D<double> w(3);
    Matrix1D<double> new_w(3);
    Matrix1D<double> D_1byw(3);
    matrix2D<double> D(3, 3), D_1(3, 3);
    double module;
    double newrot, newtilt, newpsi;
    Matrix1D<double> fw(3);
    float fnewrot, fnewtilt, fnewpsi;

    D(0, 0) = d00;
    D(0, 1) = d01;
    D(0, 2) = 0.;
    D(1, 0) = d10;
    D(1, 1) = d11;
    D(1, 2) = 0.;
    D(2, 0) = 0.;
    D(2, 1) = 0.;
    D(2, 2) = 1.;
    cout << "D matrix\n" << D;
    D_1 = D.inv();
    cout << "D_1 matrix\n" << D_1;
    Euler_direction(rot, tilt, psi, w);
    cout << "Projection direction\n" << w;

//Euler_direction2angles (w, newrot, newtilt, newpsi);
//  cout << "If everything is OK this should be the original
//  angles" << endl;
//  cout << "rot= "  << newrot
//       << "tilt= " << newtilt
//       << "psi= "  << newpsi << endl;

//fw(0)=(float)w(0);
//fw(1)=(float)w(1);
//fw(2)=(float)w(2);

//Euler_direction( (float)rot, (float)tilt, (float)psi, fw);
//  cout << "Projection direction calculated in float\n" << fw;

//Euler_direction2angles (fw, fnewrot, fnewtilt, fnewpsi);
//  cout << "If everything is OK this should be the original
//  angles calculated in float" << endl;
//  cout << "rot= "  << fnewrot
//       << "tilt= " << fnewtilt
//       << "psi= "  << fnewpsi << endl;


    D_1byw = D_1 * w;
    cout << "Projection direction after streching\n" << D_1byw;

    module = D_1byw.module();
    cout << "module_of_Dbyw= " <<   module << endl;

    new_w = D_1byw / module;
    cout << " D_1byw/module(D_1byw)\n" << new_w;

    Euler_direction2angles(new_w, newrot, newtilt, newpsi);
    if (newtilt == 0.) newrot = rot;
    cout << endl
    << "Old_rot  = " << rot  << " New_rot  = " << newrot
    << endl
    << "Old_tilt = " << tilt << " New_tilt = " << newtilt
    << endl
    << "Old_psi  = " << psi  << " New_psi  = " << newpsi
    << endl;

    Euler_direction(newrot, newtilt, newpsi, w);
    cout << "Projectionn direction made with new angles\n" << w;
    Euler_angles2matrix(rot, tilt, psi, D);
    cout << "Euler Matrix with old angles\n" << D;
    try
    {
        cout << "Inverse Euler Matrix with old angles\n" << D.inv();
    }
    catch (Xmipp_error &XE)
    {
        cout << XE;
        exit(1);
    }

    Euler_angles2matrix(newrot, newtilt, newpsi, D);
    cout << "Euler Matrix with new angles\n" << D;
    try
    {
        cout << "Inverse Euler Matrix with new angles\n" << D.inv();
    }
    catch (Xmipp_error &XE)
    {
        cout << XE;
        exit(1);
    }
}/* ROUT_EULER end */
