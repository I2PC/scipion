/***************************************************************************
 *
 * Authors:     Roberto Marabini
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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


#include <data/args.h>
#include <data/image.h>
#include <data/selfile.h>
#include <data/docfile.h>


/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    DocFile         DFi;
    DocLine Angles;
    Matrix2D<double> Euler(3, 3), Euler_i(3, 3),data(3,3),data1(3,3);
    double rot, tilt,psi;
    double rot1, tilt1,psi1;
    tilt1 = 31.7175;
    rot1=psi1=0.;

    Euler_angles2matrix(rot1, tilt1, psi1, Euler);
    Euler_i=Euler.transpose();
    try
    {

        DFi.read(getParameter(argc, argv, "-i"));

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    DFi.go_beginning();
    DFi.go_first_data_line();
    while (!DFi.eof())
    {
        rot = DFi(0);
        tilt = DFi(1);
        psi = DFi(2);
        rot = realWRAP(rot, -180, 180);
        tilt = realWRAP(tilt, -180, 180);
        psi = realWRAP(psi, -180, 180);
        Euler_angles2matrix(rot, tilt, psi, data);
        data1 = /*Euler * */(data * Euler_i);
        Euler_matrix2angles(data1,rot,tilt,psi);
std::cout << Euler << data << data1;
        DFi.set(0, rot);
        DFi.set(1, tilt);
        DFi.set(2, psi);
        DFi.next_data_line();
    }
DFi.write("kk.doc");
}
