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
    doMask = false;
    if (checkParam("--mask_in"))
    {
        fn_Maskin = getParam("--mask_in");
        doMask = true;
    }
    helical=(fn_sym=="helical");
    dihedral=(fn_sym=="dihedral");
    helicalDihedral=(fn_sym=="helicalDihedral");
    if (helical || helicalDihedral)
    {
    	Ts=getDoubleParam("--sampling");
        zHelical=getDoubleParam("--helixParams",0)/Ts;
        rotHelical=DEG2RAD(getDoubleParam("--helixParams",1));
        rotPhaseHelical=DEG2RAD(getDoubleParam("--helixParams",2));
    }
    do_not_generate_subgroup = checkParam("--no_group");
    wrap = !checkParam("--dont_wrap");
    sum = checkParam("--sum");
    heightFraction = getDoubleParam("--heightFraction");
    splineOrder = getIntParam("--spline");
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
    addParamsLine("                         : I, I1, I2, I3, I4, I5, Ih, helical, dihedral, helicalDihedral");
    addParamsLine("                         :+ For a full description of symmetries look at");
    addParamsLine("                         :+ http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry");
    addParamsLine("   [--helixParams <z> <rot> <rotPhase=0>]: Helical parameters z(Angstroms), rot(degrees), rotPhase(degrees)");
    addParamsLine("                         :+ V(r,theta,z)=V(r,theta+rotHelix+rotPhase,z+zHelix)");
    addParamsLine("   [--heightFraction <f=0.95>]: Height fraction used for symmetrizing a helix");
    addParamsLine("   [--sampling <T=1>]    : Sampling rate (A/pixel). Used only for helical parameters");
    addParamsLine("   [--no_group]          : For 3D volumes: do not generate symmetry subgroup");
    addParamsLine("   [--dont_wrap]         : by default, the image/volume is wrapped");
    addParamsLine("   [--sum]               : compute the sum of the images/volumes instead of the average. This is useful for symmetrizing pieces");
    addParamsLine("   [--mask_in <fileName>]: symmetrize only in the masked area");
    addParamsLine("   [--spline <order=3>]  : Spline order for the interpolation (valid values are 1 and 3)");
    addExampleLine("Symmetrize a list of images with 6 fold symmetry",false);
    addExampleLine("   xmipp_transform_symmetrize -i input.sel --sym 6");
    addExampleLine("Symmetrize with i3 symmetry and the volume is not wrapped",false);
    addExampleLine("   xmipp_transform_symmetrize -i input.vol --sym i3 --dont_wrap");
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
    << "Wrap:     " << wrap << std::endl
    << "Sum:      " << sum << std::endl
	<< "Spline:   " << splineOrder << std::endl;
    if (doMask)
        std::cout << "mask_in    " << fn_Maskin << std::endl;
    if (helical)
        std::cout << "RotHelical: " << RAD2DEG(rotHelical) << std::endl
        << "RotPhaseHelical: " << RAD2DEG(rotPhaseHelical) << std::endl
        << "ZHelical:   " << zHelical*Ts
		<< "Height fraction: " << heightFraction << std::endl;
}

/* Symmetrize ------------------------------------------------------- */
void symmetrizeVolume(const SymList &SL, const MultidimArray<double> &V_in,
                      MultidimArray<double> &V_out, int spline,
                      bool wrap, bool do_outside_avg, bool sum, bool helical, bool dihedral, bool helicalDihedral,
                      double rotHelical, double rotPhaseHelical, double zHelical, double heightFraction,
                      const MultidimArray<double> * mask)
{
    Matrix2D<double> L(4, 4), R(4, 4); // A matrix from the list
    MultidimArray<double> V_aux;
    Matrix1D<double> sh(3);
    double avg = 0.;

    if (do_outside_avg)
    {
        MultidimArray<int> mask1;
        mask1.resizeNoCopy(V_in);
        mask1.setXmippOrigin();
        size_t rad = XMIPP_MIN(XSIZE(V_in), YSIZE(V_in));
        rad = XMIPP_MIN(rad, ZSIZE(V_in));
        BinaryCircularMask(mask1, rad / 2, OUTSIDE_MASK);
        double dum;
        computeStats_within_binary_mask(mask1, V_in, dum, dum, avg, dum);
    }
    V_out = V_in;

    if (!helical && !dihedral && !helicalDihedral)
    {
    	MultidimArray<double> Bcoeffs;
    	MultidimArray<double> *BcoeffsPtr=NULL;
    	if (spline==3)
    	{
    		produceSplineCoefficients(BSPLINE3, Bcoeffs, V_in);
    		BcoeffsPtr=&Bcoeffs;
    	}
        for (int i = 0; i < SL.symsNo(); i++)
        {
            SL.getMatrices(i, L, R);
            /*FIXME: I do not think this make sense since V_aux is empty. ROB
               SL.getShift(i, sh);
               R(3, 0) = sh(0) * XSIZE(V_aux);
               R(3, 1) = sh(1) * YSIZE(V_aux);
               R(3, 2) = sh(2) * ZSIZE(V_aux);
            */
            applyGeometry(spline, V_aux, V_in, R.transpose(), IS_NOT_INV, wrap, avg, BcoeffsPtr);

            if ( mask==NULL)
                arrayByArray(V_out, V_aux, V_out, '+');
            else               // op1, op2 ,  result
                selfArrayByArrayMask(V_in, V_aux, V_out, '+', mask);
        }

        if (!sum)
            arrayByScalar(V_out, 1.0/(SL.symsNo() + 1.0f), V_out, '*');
    }
    else if (helical)
        symmetry_Helical(V_out,V_in,zHelical,rotHelical,rotPhaseHelical,NULL,false,heightFraction);
    else if (helicalDihedral)
    {
        symmetry_Helical(V_out,V_in,zHelical,rotHelical,rotPhaseHelical,NULL,true,heightFraction);
        MultidimArray<double> Vrotated;
        rotate(spline,Vrotated,V_out,180.0,'X',WRAP);
        V_out+=Vrotated;
        V_out*=0.5;
    }
    else if (dihedral)
    {
    	int zmax=(int)(0.1*ZSIZE(V_in));
        symmetry_Dihedral(V_out,V_in,1,-zmax,zmax,0.5);
    }
}

void symmetrizeImage(int symorder, const MultidimArray<double> &I_in,
                     MultidimArray<double> &I_out, int spline,
                     bool wrap, bool do_outside_avg, bool sum,
                     const MultidimArray<double> * mask)
{
    double avg = 0.;
    if (do_outside_avg)
    {
        MultidimArray<int> mask1;
        mask1.resizeNoCopy(I_in);
        mask1.setXmippOrigin();
        int rad = XMIPP_MIN(XSIZE(I_in), YSIZE(I_in));
        BinaryCircularMask(mask1, rad / 2, OUTSIDE_MASK);
        double dum;
        computeStats_within_binary_mask(mask1, I_in, dum, dum, avg, dum);
    }
    I_out = I_in;
    MultidimArray<double> rotatedImg;
    if ( (mask)!=NULL)
         REPORT_ERROR(ERR_NOT_IMPLEMENTED,"mask symmetrization not implemented for images");
    for (int i = 1; i < symorder; i++)
    {
        rotate(spline, rotatedImg, I_in, 360.0 / symorder * i,'Z',wrap,avg);
        I_out += rotatedImg;
    }
    if (!sum)
        I_out *= 1.0/symorder;
}

/* Preprocess ------------------------------------------------------------- */
void ProgSymmetrize::preProcess()
{
    if (!helical && !dihedral & !helicalDihedral)
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
}

/* Process image ------------------------------------------------------------- */
void ProgSymmetrize::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    Image<double> Iin;
    Image<double> Iout;
    Image<double> mask;
    MultidimArray<double> *mmask=NULL;

    Iin.readApplyGeo(fnImg, rowIn);
    Iin().setXmippOrigin();
    Iout().resize(Iin());

    if (doMask)
    {
    	mask.read(fn_Maskin);
    	mmask=&(mask());
    }

    if (ZSIZE(Iin())==1)
    {
        if (helical)
            REPORT_ERROR(ERR_ARG_INCORRECT,"Helical symmetrization is not meant for images");
        if (symorder!=-1)
            symmetrizeImage(symorder,Iin(),Iout(),splineOrder,wrap,!wrap,sum,mmask);
        else
            REPORT_ERROR(ERR_ARG_MISSING,"The symmetry order is not valid for images");
    }
    else
    {
        if (SL.symsNo()>0 || helical || dihedral || helicalDihedral)
        {
            symmetrizeVolume(SL,Iin(),Iout(),splineOrder,wrap,!wrap,
                             sum,helical,dihedral,helicalDihedral,rotHelical,rotPhaseHelical,zHelical,heightFraction,mmask);
        }
        else
            REPORT_ERROR(ERR_ARG_MISSING,"The symmetry description is not valid for volumes");
    }
    Iout.write(fnImgOut);
}
