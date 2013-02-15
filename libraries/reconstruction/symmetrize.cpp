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

#include "symmetrize.h"

#include <data/args.h>

/* Read parameters --------------------------------------------------------- */
void ProgSymmetrize::readParams()
{
	XmippMetadataProgram::readParams();
    fn_sym = getParam("--sym");
    do_not_generate_subgroup = checkParam("--no_group");
    wrap = !checkParam("--dont_wrap");
}

/* Usage ------------------------------------------------------------------- */
void ProgSymmetrize::defineParams()
{
    addUsageLine("Symmetrize volumes and images ");
    each_image_produces_an_output=true;
    XmippMetadataProgram::defineParams();
    addParamsLine("    --sym <symmetry>     : For 2D images: a number");
    addParamsLine("                         : For 3D volumes: a symmetry file or point-group description");
    addParamsLine("                         : Valid point-group descriptions are:");
    addParamsLine("                         : C1, Ci, Cs, Cn (from here on n must be an integer number with no more than 2 digits)");
    addParamsLine("                         : Cnv, Cnh, Sn, Dn, Dnv, Dnh, T, Td, Th, O, Oh");
    addParamsLine("                         : I, I1, I2, I3, I4, I5, Ih");
    addParamsLine("                         :+ For a full description of symmetries look at");
    addParamsLine("                         :+ http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry");
    addParamsLine("   [--no_group]          : For 3D volumes: do not generate symmetry subgroup");
    addParamsLine("   [--dont_wrap]         : by default, the image/volume is wrapped");
    addExampleLine("Symmetrize a list of images with 6 fold symmetry",false);
    addExampleLine("   xmipp_symmetrize -i input.sel --sym 6");
    addExampleLine("Symmetrize with i3 symmetry and the volume is not wrapped",false);
    addExampleLine("   xmipp_symmetrize -i input.vol --sym i3 --dont_wrap");
}

/* Show ------------------------------------------------------------------- */
void ProgSymmetrize::show()
{
	if (verbose==0)
		return;
    XmippMetadataProgram::show();
    std::cout
    << "Symmetry: " << fn_sym << std::endl
    << "No group: " << do_not_generate_subgroup << std::endl
    << "Wrap:     " << wrap << std::endl;
}

/* Symmetrize ------------------------------------------------------- */
void symmetrizeVolume(const SymList &SL, const MultidimArray<double> &V_in,
                      MultidimArray<double> &V_out,
                      bool wrap, bool do_outside_avg)
{
    Matrix2D<double> L(4, 4), R(4, 4); // A matrix from the list
    MultidimArray<double> V_aux;
    Matrix1D<double> sh(3);
    double avg = 0.;

    if (do_outside_avg)
    {
        MultidimArray<int> mask;
        mask.resizeNoCopy(V_in);
        mask.setXmippOrigin();
        size_t rad = XMIPP_MIN(XSIZE(V_in), YSIZE(V_in));
        rad = XMIPP_MIN(rad, ZSIZE(V_in));
        BinaryCircularMask(mask, rad / 2, OUTSIDE_MASK);
        double dum;
        computeStats_within_binary_mask(mask, V_in, dum, dum, avg, dum);
    }
    V_out = V_in;

    for (int i = 0; i < SL.symsNo(); i++)
    {
        SL.getMatrices(i, L, R);

        SL.getShift(i, sh);
        R(3, 0) = sh(0) * XSIZE(V_aux);
        R(3, 1) = sh(1) * YSIZE(V_aux);
        R(3, 2) = sh(2) * ZSIZE(V_aux);
        applyGeometry(BSPLINE3, V_aux, V_in, R.transpose(), IS_NOT_INV, wrap, avg);
        arrayByArray(V_out, V_aux, V_out, '+');
    }
    arrayByScalar(V_out, 1.0/(SL.symsNo() + 1.0f), V_out, '*');
}

void symmetrizeImage(int symorder, const MultidimArray<double> &I_in,
                     MultidimArray<double> &I_out,
                     bool wrap, bool do_outside_avg)
{
    double avg = 0.;
    if (do_outside_avg)
    {
        MultidimArray<int> mask;
        mask.resizeNoCopy(I_in);
        mask.setXmippOrigin();
        int rad = XMIPP_MIN(XSIZE(I_in), YSIZE(I_in));
        BinaryCircularMask(mask, rad / 2, OUTSIDE_MASK);
        double dum;
        computeStats_within_binary_mask(mask, I_in, dum, dum, avg, dum);
    }

    I_out = I_in;
    MultidimArray<double> rotatedImg;
    for (int i = 1; i < symorder; i++)
    {
        rotate(BSPLINE3, rotatedImg, I_in, 360.0 / symorder * i,'Z',wrap,avg);
        I_out += rotatedImg;
    }
    I_out *= 1.0/symorder;
}

/* Preprocess ------------------------------------------------------------- */
void ProgSymmetrize::preProcess()
{
    if (!fn_sym.exists() && isdigit(fn_sym[0]))
        symorder=textToInteger(fn_sym);
    else
    {
        double accuracy = (do_not_generate_subgroup) ? -1 : 1e-6;
        SL.readSymmetryFile(fn_sym, accuracy);
        symorder=-1;
    }
}

/* Process image ------------------------------------------------------------- */
void ProgSymmetrize::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    Image<double> Iin;
    Image<double> Iout;
    Iin.readApplyGeo(fnImg, rowIn);
    Iin().setXmippOrigin();
    if (ZSIZE(Iin())==1)
    {
        if (symorder!=-1)
            symmetrizeImage(symorder,Iin(),Iout(),wrap,!wrap);
        else
            REPORT_ERROR(ERR_ARG_MISSING,"The symmetry order is not valid for images");
    }
    else
    {
        if (SL.symsNo()>0)
            symmetrizeVolume(SL,Iin(),Iout(),wrap,!wrap);
        else
            REPORT_ERROR(ERR_ARG_MISSING,"The symmetry description is not valid for volumes");
    }
    Iout.write(fnImgOut);
}
