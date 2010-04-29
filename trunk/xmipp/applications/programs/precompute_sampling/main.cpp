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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/fftw.h>
#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/geometry.h>
#include <vector>


/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{

    ImageXmipp test(4,4);
    ImageXmipp rot(4,4);
    Matrix2D<double> I(4,4), AA(3,3);
    Euler_angles2matrix(0,0,-90,AA);

    I.initIdentity();
    test(0,0)=1;
    test(0,1)=2;
    test(0,2)=0;
    test(0,3)=0;

    test(1,0)=3;
    test(1,1)=4;
    test(1,2)=0;
    test(1,3)=0;

    test(2,0)=0;
    test(2,1)=0;
    test(2,2)=0;
    test(2,3)=0;

    test(3,0)=0;
    test(3,1)=0;
    test(3,2)=0;
    test(3,3)=0;

    test().setXmippOrigin();
    //BinaryWedgeMask(test(),-60.,60.,I);
    ////test.write("before.vol");
    std::cerr << "test" <<  test() << std::endl;

    rot=test;
    rot().selfApplyGeometryBSpline(AA,3,IS_NOT_INV,DONT_WRAP,0.);

    std::cerr << "rot real space" <<  rot() << std::endl;
    //rot.write("rotated.vol");
    Matrix2D<std::complex<double> > Faux, Faux2;
    Matrix2D<double> Mr(test()), Mi(test());
    FourierTransform(test(),Faux);
    Mr.setXmippOrigin();
    Mi.setXmippOrigin();
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
    {
        dMij(Mr,i,j) = (dMij(Faux,i,j)).real();
        dMij(Mi,i,j) = (dMij(Faux,i,j)).imag();
    }
    std::cerr << "FourierTransform" <<  Mr*4 << Mi*4 << std::endl;

    //CenterOriginFFT(Faux,true);
    ShiftFFT(Faux, 2, 2);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
    {
        dMij(Mr,i,j) = (dMij(Faux,i,j)).real();
        dMij(Mi,i,j) = (dMij(Faux,i,j)).imag();
    }
    std::cerr << "ShiftFourierTransform" <<  Mr*4 << Mi*4 << std::endl;
    CenterFFT(Faux, true);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
    {
        dMij(Mr,i,j) = (dMij(Faux,i,j)).real();
        dMij(Mi,i,j) = (dMij(Faux,i,j)).imag();
    }
    std::cerr << "CenterShiftFourierTransform" <<  Mr*4 << Mi*4 << std::endl;

    Mr.selfApplyGeometryBSpline(AA,3,IS_NOT_INV,DONT_WRAP,0.);
    Mi.selfApplyGeometryBSpline(AA,3,IS_NOT_INV,DONT_WRAP,0.);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Faux)
    {
        dMij(Faux,i,j) = std::complex<double>(dMij(Mr,i,j), dMij(Mi,i,j));
    }
    std::cerr<<"done rot in Fourier"<< Mr*4 << Mi*4 << std::endl;
    CenterFFT(Faux, false);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
    {
        dMij(Mr,i,j) = (dMij(Faux,i,j)).real();
        dMij(Mi,i,j) = (dMij(Faux,i,j)).imag();
    }
    std::cerr << "Anti_CenterShiftFourierTransform" <<  Mr*4 << Mi*4 << std::endl;
    ShiftFFT(Faux, -2, -2);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(Faux)
    {
        dMij(Mr,i,j) = (dMij(Faux,i,j)).real();
        dMij(Mi,i,j) = (dMij(Faux,i,j)).imag();
    }
    std::cerr << "Anti_ShiftFourierTransform" <<  Mr*4 << Mi*4 << std::endl;


    InverseFourierTransform(Faux,test());
    //test.write("after.vol");
    std::cerr << "rot fourier space" <<  test() << std::endl;
}
