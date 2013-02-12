/***************************************************************************
 *
 * Authors:    Sjors Scheres                 (scheres@cnb.csic.es)
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

#include "align2d.h"

#include <data/xmipp_funcs.h>
#include <data/mask.h>
#include <data/filters.h>

// Read arguments ==========================================================
void ProgAlign2d::readParams()
{
    fnSel = getParam("-i");
    fnRoot = getParam("--oroot");
    fnRef = getParam("--ref");
    Niter = getIntParam("--iter");
    dont_mirror = checkParam("--do_not_check_mirrors");
    pspc = checkParam("--pspc");
}

// Show ====================================================================
void ProgAlign2d::show()
{
    if (verbose==0)
        return;
    std::cerr
    << "Input selfile:        " << fnSel       << std::endl
    << "Input reference:      " << fnRef       << std::endl
    << "Output ootname:       " << fnRoot      << std::endl
    << "Number of iterations: " << Niter       << std::endl
    << "Do not check mirrors: " << dont_mirror << std::endl
    << "PSPC:                 " << pspc        << std::endl
    ;
}

// usage ===================================================================
void ProgAlign2d::defineParams()
{
    addUsageLine("Aligns a set of images");
    addParamsLine("  -i <selfile>             : Selfile containing images to be aligned");
    addParamsLine("  --oroot <rootname>       : Output rootname");
    addParamsLine(" [--ref <image=\"\">]      : reference image; if none: pyramidal combination of subset of images");
    addParamsLine(" [--iter <N=5>]            : Number of iterations to perform");
    addParamsLine(" [--do_not_check_mirrors]  : Do not consider mirrors when aligning");
    addParamsLine(" [--pspc]                  : Compute first reference by pspc");
}

// PsPc pyramidal combination of images ========================================
//#define DEBUG
void ProgAlign2d::alignPairs(MetaData &MDin, MetaData &MDout, int level)
{
    // Compute output stack size
    size_t imax=MDin.size()/2;
    size_t remaining=MDin.size()-2*imax;
    size_t Xdim, Ydim, Zdim, Ndim;
    getImageSize(MDin,Xdim,Ydim,Zdim,Ndim);
    if (Zdim!=1 || Ndim!=1)
        REPORT_ERROR(ERR_MATRIX_DIM,"Files in metadata are not 2D images");
    FileName fnOutputStack=fnRoot+"_level_"+integerToString(level,2)+".stk";
    fnOutputStack.deleteFile();
    createEmptyFile(fnOutputStack,Xdim,Ydim,1,imax+remaining,true,WRITE_REPLACE);

    // Align all pairs
    MDIterator mdIter(MDin);
    Image<double> I1, I2;
    Matrix2D<double> M;
    FileName fnOut;
    AlignmentAux aux1;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    std::cerr << "Aligning level " << level << std::endl;
    init_progress_bar(imax);
    for (size_t i=0; i<imax; i++)
    {
        // Read the two input images
        I1.readApplyGeo(MDin,mdIter.objId);
        mdIter.moveNext();
        I2.readApplyGeo(MDin,mdIter.objId);
        mdIter.moveNext();

        // Align images
        MultidimArray<double> &I1m=I1();
        MultidimArray<double> &I2m=I2();
        if (dont_mirror)
            alignImages(I1m,I2m,M);
        else
            alignImagesConsideringMirrors(I1m,I2m,M,aux1,aux2,aux3);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I1m)
            DIRECT_MULTIDIM_ELEM(I1m,n)=0.5*(DIRECT_MULTIDIM_ELEM(I1m,n)+
                                             DIRECT_MULTIDIM_ELEM(I2m,n));
        centerImage(I1m,aux2,aux3);

        // Write to output stack
        fnOut.compose(i+1,fnOutputStack);
        I1.write(fnOut);
        size_t objId=MDout.addObject();
        MDout.setValue(MDL_IMAGE,fnOut,objId);

#ifdef DEBUG
        I1.write("PPPI1.xmp");
        I2.write("PPPI2.xmp");
        std::cout << "M=" << M;
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
#endif
        if (i%100==0)
        	progress_bar(i);
    }
    progress_bar(imax);

    if (remaining)
    {
        I1.readApplyGeo(MDin,mdIter.objId);
        mdIter.moveNext();
        fnOut.compose(imax+1,fnOutputStack);
        I1.write(fnOut);
    }
}
#undef DEBUG

void ProgAlign2d::do_pspc()
{
    int level=0;
    MetaData SFlevel=SF;
    MetaData SFlevel_1;
    while (SFlevel.size()>1)
    {
        alignPairs(SFlevel,SFlevel_1,level);
        if (level>=1)
            unlink((fnRoot+"_level_"+integerToString(level-1,2)+".stk").c_str());
        level++;
        SFlevel=SFlevel_1;
        SFlevel_1.clear();
    }
    size_t objId=SFlevel.firstObject();
    Iref.readApplyGeo(SFlevel,objId);
    unlink((fnRoot+"_level_"+integerToString(level-1,2)+".stk").c_str());
    Iref.write(fnRoot+"_pspc.xmp");
}

void ProgAlign2d::computeMean()
{
    size_t Xdim, Ydim, Zdim, Ndim;
    getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
    if (Zdim!=1 || Ndim!=1)
        REPORT_ERROR(ERR_MATRIX_DIM,"Files in metadata are not 2D images");
    Iref().initZeros(Ydim,Xdim);
    Image<double> I;
    FileName fnImg;
    size_t N=SF.size();
    std::cerr << "Computing average of images ...\n";
    init_progress_bar(N);
    int i=0;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
    	I.readApplyGeo(SF,__iter.objId);
    	Iref()+=I();
    	if ((++i)%100==0)
    		progress_bar(i);
    }
    progress_bar(N);
    Iref()*=1.0/N;
    CorrelationAux aux;
    RotationalCorrelationAux aux2;
    centerImage(Iref(),aux,aux2);
}

// Alignment of all images by iterative refinement  ========================================
//#define DEBUG
void ProgAlign2d::refinement()
{
    Image<double> I;
    size_t N=SF.size();
	double lambda=1.0/(N/2.0);
	double lambdap=1-lambda;
    init_progress_bar(N);
    int i=0;
    MultidimArray<double> &Irefm=Iref();
    Matrix2D<double> M;
    int centerCount=0;
    FileName fnImg;
    AlignmentAux aux1;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
    	SF.getValue(MDL_IMAGE,fnImg,__iter.objId);
    	I.read(fnImg);
    	I().setXmippOrigin();
        MultidimArray<double> Ibackup, Ialigned;
        Ibackup=I();

        // Align images
        MultidimArray<double> &Im=I();
        double corr;
        if (dont_mirror)
            corr=alignImages(Irefm,Im,M);
        else
            corr=alignImagesConsideringMirrors(Irefm,Im,M,aux1,aux2,aux3);
        applyGeometry(LINEAR, Ialigned, Ibackup, M, IS_NOT_INV, WRAP);
        Im=Ialigned;

        // Save alignment
        bool flip;
        double scale, shiftX, shiftY, psi;
        transformationMatrix2Parameters2D(M, flip, scale, shiftX, shiftY, psi);
        SF.setValue(MDL_SHIFT_X,shiftX,__iter.objId);
        SF.setValue(MDL_SHIFT_Y,shiftY,__iter.objId);
        SF.setValue(MDL_ANGLE_PSI,psi,__iter.objId);
        SF.setValue(MDL_FLIP,flip,__iter.objId);
        SF.setValue(MDL_MAXCC,corr,__iter.objId);

        // Update reference
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Irefm)
            DIRECT_MULTIDIM_ELEM(Irefm,n)=lambdap*DIRECT_MULTIDIM_ELEM(Irefm,n)+
                                          lambda*DIRECT_MULTIDIM_ELEM(Im,n);

        // From time to time, recenter the reference
        if ((++centerCount)%10==0)
        {
        	centerImage(Irefm,aux2,aux3);
        	centerCount=0;
        }

        if ((++i)%100==0)
    		progress_bar(i);
#ifdef DEBUG
        Iref.write("PPPref.xmp");
        I.write("PPPexp.xmp");
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
#endif
    }
    progress_bar(N);
    Iref.write(fnRoot+"_ref.xmp");
    SF.write(fnRoot+"_alignment.xmd");
}
#undef DEBUG

// Write out results  ========================================================
void ProgAlign2d::run()
{
    // Read the input metadata
    SF.read(fnSel);

    // Get Reference
    if (fnRef != "")
        Iref.read(fnRef);
    else
    {
    	if (pspc)
    		do_pspc();
    	else
    		computeMean();
    }
    Iref().setXmippOrigin();

    // Circular mask the reference image
    MultidimArray<int> mask;
    mask.resize(Iref());
    mask.setXmippOrigin();
    BinaryCircularMask(mask,XSIZE(Iref())/2);
    apply_binary_mask(mask, Iref(), Iref());

    // Refine
    std::cout << "Refining ...\n";
    for (int n=0; n<Niter; n++)
    {
    	std::cerr << "Refinement iteration " << n << std::endl;
    	refinement();
    }
}
