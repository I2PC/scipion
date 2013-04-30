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

#include <string>
#include <iostream>

#include "args.h"
#include "micrograph.h"
#include "normalize.h"
#include "metadata.h"
#include "xmipp_image_generic.h"

/* Normalizations ---------------------------------------------------------- */
void normalize_OldXmipp(MultidimArray<double> &I)
{
    double mean,std;
    I.computeAvgStdev(mean,std);
    double istd=1.0/std;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I)
    DIRECT_MULTIDIM_ELEM(I,n)=(DIRECT_MULTIDIM_ELEM(I,n)-mean)*istd;
}

void normalize_Near_OldXmipp(MultidimArray<double> &I, const MultidimArray<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    I.computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,
                                    stddevbg);
    I -= avg;
    I /= stddevbg;
}

void normalize_OldXmipp_decomposition(MultidimArray<double> &I, const MultidimArray<int> &bg_mask,
                                      const MultidimArray<double> *mask)
{
    double avgbg, stddevbg, minbg, maxbg;
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,
                                    stddevbg);
    I -= avgbg;
    I /= stddevbg;
    if (mask != NULL)
        I *= *mask;
    normalize_OldXmipp(I);
}

//#define DEBUG
void normalize_tomography(MultidimArray<double> &I, double tilt, double &mui,
                          double &sigmai, bool tiltMask, bool tomography0,
                          double mu0, double sigma0)
{
    // Tilt series are always 2D
    I.checkDimension(2);

    const int L=2;

    // Build a mask using the tilt angle
    I.setXmippOrigin();
    MultidimArray<int> mask;
    mask.initZeros(I);
    int Xdimtilt=(int)XMIPP_MIN(FLOOR(0.5*(XSIZE(I)*cos(DEG2RAD(tilt)))),
                           0.5*(XSIZE(I)-(2*L+1)));
    int N=0;
    for (int i=STARTINGY(I); i<=FINISHINGY(I); i++)
        for (int j=-Xdimtilt; j<=Xdimtilt;j++)
        {
            A2D_ELEM(mask,i,j)=1;
            N++;
        }

    // Estimate the local variance
    MultidimArray<double> localVariance;
    localVariance.initZeros(I);
    double meanVariance=0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    {
        // Center a mask of size 5x5 and estimate the variance within the mask
        double meanPiece=0, variancePiece=0;
        double Npiece=0;
        for (int ii=i-L; ii<=i+L; ii++)
        {
            if (INSIDEY(I,ii))
                for (int jj=j-L; jj<=j+L; jj++)
                {
                    if (INSIDEX(I,jj))
                    {
                        double pixval=A2D_ELEM(I,ii,jj);
                        meanPiece+=pixval;
                        variancePiece+=pixval*pixval;
                        ++Npiece;
                    }
                }
        }
        meanPiece/=Npiece;
        variancePiece=variancePiece/(Npiece-1)-
                      Npiece/(Npiece-1)*meanPiece*meanPiece;
        A2D_ELEM(localVariance,i,j)=variancePiece;
        meanVariance+=variancePiece;
    }
    meanVariance*=1.0/N;

    // Test the hypothesis that the variance in this piece is
    // the same as the variance in the whole image
    double iFu=1/icdf_FSnedecor(4*L*L+4*L,N-1,0.975);
    double iFl=1/icdf_FSnedecor(4*L*L+4*L,N-1,0.025);
    double iMeanVariance=1.0/meanVariance;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(localVariance)
    {
        if (A2D_ELEM(localVariance,i,j)==0)
            A2D_ELEM(mask,i,j)=-2;
        if (A2D_ELEM(mask,i,j)==1)
        {
            double ratio=A2D_ELEM(localVariance,i,j)*iMeanVariance;
            double thl=ratio*iFu;
            double thu=ratio*iFl;
            if (thl>1 || thu<1)
                A2D_ELEM(mask,i,j)=-1;
        }
    }
#ifdef DEBUG
    Image<double> save;
    save()=I;
    save.write("PPP.xmp");
    save()=localVariance;
    save.write("PPPLocalVariance.xmp");
    Image<int> savemask;
    savemask()=mask;
    savemask.write("PPPmask.xmp");
#endif

    // Compute the statistics again in the reduced mask
    double avg, stddev, min, max;
    computeStats_within_binary_mask(mask, I, min, max, avg, stddev);
    double cosTilt=cos(DEG2RAD(tilt));
    double iCosTilt=1.0/cosTilt;
    if (tomography0)
    {
        double iadjustedStddev=1.0/(sigma0*cosTilt);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
        switch (A2D_ELEM(mask,i,j))
        {
        case -2:
            A2D_ELEM(I,i,j)=0;
            break;
        case 0:
            if (tiltMask)
                A2D_ELEM(I,i,j)=0;
            else
                A2D_ELEM(I,i,j)=(A2D_ELEM(I,i,j)*iCosTilt-mu0)*iadjustedStddev;
            break;
        case -1:
        case 1:
            A2D_ELEM(I,i,j)=(A2D_ELEM(I,i,j)*iCosTilt-mu0)*iadjustedStddev;
            break;
        }
    }
    else
    {
        double adjustedStddev=sqrt(meanVariance)*cosTilt;
        adjustedStddev=stddev*cosTilt;
        double iAdjustedStddev=1.0/adjustedStddev;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
        switch (A2D_ELEM(mask,i,j))
        {
        case -2:
            A2D_ELEM(I,i,j)=0;
            break;
        case 0:
            if (tiltMask)
                A2D_ELEM(I,i,j)=0;
            else
                A2D_ELEM(I,i,j)=(A2D_ELEM(I,i,j)-avg)*iAdjustedStddev;
            break;
        case -1:
        case 1:
            A2D_ELEM(I,i,j)=(A2D_ELEM(I,i,j)-avg)*iAdjustedStddev;
            break;
        }
    }

    // Prepare values for returning
    mui=avg;
    sigmai=sqrt(meanVariance);
#ifdef DEBUG

    save()=I;
    save.write("PPPafter.xmp");
    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
#endif
}
#undef DEBUG

void normalize_Michael(MultidimArray<double> &I, const MultidimArray<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    I.computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,
                                    stddevbg);
    if (avgbg > 0)
    {
        I -= avgbg;
        I /= avgbg;
    }
    else
    { // To avoid the contrast inversion
        I -= (avgbg - min);
        I /= (avgbg - min);
    }
}

void normalize_NewXmipp(MultidimArray<double> &I, const MultidimArray<int> &bg_mask)
{
    double avgbg, stddevbg;
    computeAvgStdev_within_binary_mask(bg_mask, I, avgbg, stddevbg);
    double istddevbg=1.0/stddevbg;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I)
    DIRECT_MULTIDIM_ELEM(I,n)=(DIRECT_MULTIDIM_ELEM(I,n)-avgbg)*istddevbg;
}

void normalize_NewXmipp2(MultidimArray<double> &I, const MultidimArray<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    I.computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,
                                    stddevbg);
    double K=1.0/ABS(avg - avgbg);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I)
    DIRECT_MULTIDIM_ELEM(I,n)=(DIRECT_MULTIDIM_ELEM(I,n)-avgbg)*K;
}

void normalize_ramp(MultidimArray<double> &I, MultidimArray<int> &bg_mask)
{
    // Only 2D ramps implemented
    I.checkDimension(2);

    int Npoints=(int)bg_mask.sum();
    if (Npoints<=1)
        return;
    FitPoint *allpoints=new FitPoint[Npoints];

    // Fit a least squares plane through the background pixels
    I.setXmippOrigin();
    bg_mask.setXmippOrigin();
    int idx=0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    {
        if (A2D_ELEM(bg_mask, i, j))
        {
        	FitPoint &p=allpoints[idx++];
            p.x = j;
            p.y = i;
            p.z = A2D_ELEM(I, i, j);
            p.w = 1.;
        }
    }

    double pA, pB, pC;
    least_squares_plane_fit(allpoints, Npoints, pA, pB, pC);
    delete [] allpoints;

    // Substract the plane from the image and compute stddev within mask
    double sum1 = 0;
    double sum2 = 0;
    for (int i=STARTINGY(I); i<=FINISHINGY(I); i++)
    {
        double aux=pB * i + pC;
        for (int j=STARTINGX(I); j<=FINISHINGX(I); j++)
        {
            A2D_ELEM(I, i, j) -= pA * j + aux;
            if (A2D_ELEM(bg_mask, i, j))
            {
                double aux=A2D_ELEM(I,i,j);
                sum1 += aux;
                sum2 += aux*aux;
            }
        }
    }

    double avgbg  = sum1 / Npoints;
    double stddevbg = sqrt(fabs(sum2 / Npoints - avgbg * avgbg) * Npoints / (Npoints - 1));
    if (stddevbg>1e-6)
        I *= 1.0/stddevbg;
}

void normalize_remove_neighbours(MultidimArray<double> &I,
                                 const MultidimArray<int> &bg_mask,
                                 const double &threshold)
{
    double             pA, pB, pC;
    double             avgbg, stddevbg, minbg, maxbg, aux, newstddev;
    double             sum1 = 0.;
    double             sum2 = 0;
    int                N = 0;

    // Only implemented for 2D arrays
    I.checkDimension(2);

    // Fit a least squares plane through the background pixels
    int Npoints=(int)bg_mask.sum();
    FitPoint *allpoints=new FitPoint[Npoints];
    I.setXmippOrigin();

    // Get initial statistics
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,stddevbg);

    // Fit plane through those pixels within +/- threshold*sigma
    int idx=0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    {
        if (A2D_ELEM(bg_mask, i, j))
        {
            if ( fabs(avgbg - A2D_ELEM(I, i, j)) < threshold * stddevbg)
            {
            	FitPoint &p=allpoints[idx++];
                p.x = j;
                p.y = i;
                p.z = A2D_ELEM(I, i, j);
                p.w = 1.;
            }
        }
    }
    least_squares_plane_fit(allpoints, Npoints, pA, pB, pC);
    delete []allpoints;

    // Substract the plane from the image
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    {
        A2D_ELEM(I, i, j) -= pA * j + pB * i + pC;
    }

    // Get std.dev. of the background pixels within +/- threshold*sigma
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    {
        if (A2D_ELEM(bg_mask, i, j))
        {
            if ( ABS(A2D_ELEM(I, i, j)) < threshold * stddevbg)
            {
                N++;
                sum1 +=  (double) A2D_ELEM(I, i, j);
                sum2 += ((double) A2D_ELEM(I, i, j)) *
                        ((double) A2D_ELEM(I, i, j));
            }
        }
    }
    // average and standard deviation
    aux = sum1 / (double) N;
    newstddev = sqrt(fabs(sum2 / N - aux*aux) * N / (N - 1));

    // Replace pixels outside +/- threshold*sigma by samples from
    // a gaussian with avg-plane and newstddev
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
    {
        if (A2D_ELEM(bg_mask, i, j))
        {
            if ( fabs(A2D_ELEM(I, i, j)) > threshold * stddevbg)
            {
                // get local average
                aux = pA * j + pB * i + pC;
                A2D_ELEM(I, i, j)=rnd_gaus(aux, newstddev );
            }
        }
    }

    // Divide the entire image by the new background
    I /= newstddev;
}

void ProgNormalize::defineParams()
{
    each_image_produces_an_output = true;
    allow_apply_geo = true;
    save_metadata_stack = true;
    keep_input_columns = true;
    addUsageLine("Change the range of intensity values of pixels.");
    addUsageLine("In general, most of the methods requires a background to separate "
                 "particles from noise");
    addKeywords("mask,normalization,normalize");
    XmippMetadataProgram::defineParams();
    addParamsLine(" [--method <mth=NewXmipp>]     : Normalizing method.");
    addParamsLine("           where <mth>");
    addParamsLine("           OldXmipp            : I=(I-m(I))/stddev(I)");
    addParamsLine("                               : Avg(I)=0, Stddev(I)=1, does not need background");
    addParamsLine("           Near_OldXmipp       : I=(I-m(I))/stddev(bg)");
    addParamsLine("           NewXmipp            : I=(I-m(bg))/stddev(bg)");
    addParamsLine("                               : Avg(bg)=0, Stddev(I)=1");
    addParamsLine("                               : Positivity constraints can be added in reconstruction");
    addParamsLine("           Tomography          : I=(I-mean(I))/(stddev(I)*cos(tilt))");
    addParamsLine("                               : does not need background, it assumes the tilt series is vertically aligned");
    addParamsLine("                               : Similar to OldXmipp but with an extra division by the cos(tilt)");
    addParamsLine("           Tomography0         : I=(I-mean(I(0 degrees)))/(stddev(I)*cos(tilt))");
    addParamsLine("                               : does not need background");
    addParamsLine("                               : Similar to Tomography but the average at 0 degrees is used for all images");
    addParamsLine("           NewXmipp2           : I=(I-m(bg))/(m(I)-m(bg))");
    addParamsLine("           Michael             : I=(I-m(bg))/stddev(bg)");
    addParamsLine("           None                   : Used for removing only dust");
    addParamsLine("           Random                 : I=aI+b");
    addParamsLine("           Ramp                   : Subtract ramp and then NewXmipp");
    addParamsLine("           Neighbour              : Replace pixels in the background with random noise");
    addParamsLine(" [--invert]                       : Invert contrast.");
    addParamsLine(" [--thr_black_dust <sblack=-3.5>] : Remove black dust particles with sigma threshold sblack.");
    addParamsLine(" [--thr_white_dust <swhite=3.5>]  : Remove white dust particles with sigma threshold swhite.");
    addParamsLine(" [--thr_neigh <value=1.2>]        : Sigma threshold for neighbour removal.");
    addParamsLine(" [--prm <a0> <aF> <b0> <bF>]      : Requires --method Random. I=aI+b.");
    //    addParamsLine("      requires -method Random;");
    addParamsLine(" [--tiltMask]                     : Apply a mask depending on the tilt");
    addParamsLine("                                  : requires --method Tomography or Tomography0");
    addParamsLine(" [--background <mode>] ");
    addParamsLine("        where <mode>");
    addParamsLine("           frame <r>              : Rectangular background of r pixels.");
    addParamsLine("           circle <r>             : Circular background outside radius r.");
    mask_prm.defineParams(this,INT_MASK, "or", "Use an alternative type of background mask.");
    addExampleLine("Normalize using OldXmipp method",false);
    addExampleLine("xmipp_transform_normalize -i images.sel --method OldXmipp",true);
    addExampleLine("Normalize 64x64 images using NewXmipp method",false);
    addExampleLine("xmipp_transform_normalize -i images.sel --method NewXmipp --background circle 29",true);
    addExampleLine("Normalize 64x64 images using NewXmipp method and a crown mask",false);
    addExampleLine("xmipp_transform_normalize -i images.sel --method NewXmipp --mask crown 29 32",true);
    addExampleLine("Normalize a volume to have zero mean and unit variance",false);
    addExampleLine("xmipp_transform_normalize -i volume.vol",true);
}

void ProgNormalize::readParams()
{
    XmippMetadataProgram::readParams();

    // Get normalizing method
    std::string aux;
    aux = getParam("--method");

    if (aux == "OldXmipp")
        method = OLDXMIPP;
    else if (aux == "Near_OldXmipp")
        method = NEAR_OLDXMIPP;
    else if (aux == "NewXmipp")
        method = NEWXMIPP;
    else if (aux == "NewXmipp2")
        method = NEWXMIPP2;
    else if (aux == "Michael")
        method = MICHAEL;
    else if (aux == "Random")
        method = RANDOM;
    else if (aux == "None")
        method = NONE;
    else if (aux == "Ramp")
        method = RAMP;
    else if (aux == "Neighbour")
        method = NEIGHBOUR;
    else if (aux == "Tomography")
        method = TOMOGRAPHY;
    else if (aux == "Tomography0")
        method = TOMOGRAPHY0;
    else
        REPORT_ERROR(ERR_VALUE_INCORRECT, "Normalize: Unknown normalizing method");

    // Invert contrast?
    invert_contrast = checkParam("--invert");

    // Apply a mask depending on the tilt
    tiltMask = checkParam("--tiltMask");

    // Remove dust particles?
    remove_black_dust = checkParam("--thr_black_dust");
    remove_white_dust = checkParam("--thr_white_dust");
    thresh_black_dust = getDoubleParam("--thr_black_dust");
    thresh_white_dust = getDoubleParam("--thr_white_dust");
    thresh_neigh      = getDoubleParam("--thr_neigh");

    // Get background mask
    background_mode = NOBACKGROUND;
    if (method == NEWXMIPP || method == NEWXMIPP2 || method == MICHAEL ||
        method == NEAR_OLDXMIPP || method == RAMP || method == NEIGHBOUR)
    {
        enable_mask = checkParam("--mask");
        if (enable_mask)
        {
            mask_prm.allowed_data_types = INT_MASK;
            mask_prm.readParams(this);
        }
        else
        {
            if (!checkParam("--background"))
                REPORT_ERROR(ERR_ARG_MISSING,
                             "Normalize: --background or --mask parameter required.");

            aux = getParam("--background",0);
            r  =  getIntParam("--background",1);

            if (aux == "frame")
                background_mode = FRAME;
            else if (aux == "circle")
                background_mode = CIRCLE;
            else
                REPORT_ERROR(ERR_VALUE_INCORRECT, "Normalize: Unknown background mode");
        }
    }

    if (method == RANDOM)
    {
        if (!checkParam("--prm"))
            REPORT_ERROR(ERR_ARG_MISSING,
                         "Normalize_parameters: -prm parameter required.");

        a0 = getDoubleParam("--prm", 0);
        aF = getDoubleParam("--prm", 1);
        b0 = getDoubleParam("--prm", 2);
        bF = getDoubleParam("--prm", 3);
    }
}

void ProgNormalize::show()
{
    XmippMetadataProgram::show();

    std::cout << "Normalizing method: ";
    switch (method)
    {
    case OLDXMIPP:
        std::cout << "OldXmipp\n";
        break;
    case NEAR_OLDXMIPP:
        std::cout << "Near_OldXmipp\n";
        break;
    case NEWXMIPP:
        std::cout << "NewXmipp\n";
        break;
    case NEWXMIPP2:
        std::cout << "NewXmipp2\n";
        break;
    case MICHAEL:
        std::cout << "Michael\n";
        break;
    case NONE:
        std::cout << "None\n";
        break;
    case RAMP:
        std::cout << "Ramp\n";
        break;
    case NEIGHBOUR:
        std::cout << "Neighbour\n";
        break;
    case TOMOGRAPHY:
        std::cout << "Tomography\n";
        break;
    case TOMOGRAPHY0:
        std::cout << "Tomography0\n";
        break;
    case RANDOM:
        std::cout << "Random a=[" << a0 << "," << aF << "], " << "b=[" <<
        b0 << "," << bF << "]\n";
        break;
    }

    if (method == NEWXMIPP || method == NEWXMIPP2 ||
        method == NEAR_OLDXMIPP || method == MICHAEL ||
        method == RAMP || method == NEIGHBOUR)
    {
        std::cout << "Background mode: ";
        switch (background_mode)
        {
        case NONE :
            std::cout << "None\n";
            break;
        case FRAME:
            std::cout << "Frame, width=" << r << std::endl;
            std::cout << "Apply transformation to mask: " << apply_geo <<
            std::endl;
            break;
        case CIRCLE:
            std::cout << "Circle, radius=" << r << std::endl;
            std::cout << "Apply transformation to mask: " << apply_geo <<
            std::endl;
            break;
        }
    }

    if (invert_contrast)
        std::cout << "Invert contrast "<< std::endl;

    if (tiltMask)
        std::cout << "Applying a mask depending on tilt "<< std::endl;

    if (remove_black_dust)
        std::cout << "Remove black dust particles, using threshold " <<
        floatToString(thresh_black_dust) << std::endl;

    if (remove_white_dust)
        std::cout << "Remove white dust particles, using threshold " <<
        floatToString(thresh_white_dust) << std::endl;

    if (method == NEWXMIPP && enable_mask)
        mask_prm.show();
}


void ProgNormalize::preProcess()
{
    if (!enable_mask)
    {
        bg_mask.resizeNoCopy(zdimOut, ydimOut, xdimOut);
        bg_mask.setXmippOrigin();

        switch (background_mode)
        {
        case FRAME:
            BinaryFrameMask(bg_mask, xdimOut - 2 * r, ydimOut - 2 * r, zdimOut - 2 * r,
                            OUTSIDE_MASK);
            break;
        case CIRCLE:
            BinaryCircularMask(bg_mask, r, OUTSIDE_MASK);
            break;
        case NOBACKGROUND:
        	break;
        }
    }
    else
    {
        mask_prm.generate_mask(zdimOut, ydimOut, xdimOut);
        bg_mask = mask_prm.imask;
    }
    // backup a copy of the mask for apply_geo mode
    bg_mask_bck = bg_mask;

    //#define DEBUG
#ifdef DEBUG

    Image<int> tt;
    tt()=bg_mask;
    tt.write("PPPmask.xmp");
    std::cerr<<"DEBUG info: written PPPmask.xmp"<<std::endl;
#endif

    // Get the parameters from the 0 degrees
    if (method==TOMOGRAPHY0)
    {
        // Look for the image at 0 degrees
        double bestTilt=1000, tiltTemp;
        FileName bestImage;
        FileName fn_img;
        ImageGeneric Ig;
        MetaData * md = getInputMd();
        FOR_ALL_OBJECTS_IN_METADATA(*md)
        {
            md->getValue(image_label, fn_img, __iter.objId);

            if (fn_img.empty())
                break;

            if (!md->getValue(MDL_ANGLE_TILT,tiltTemp,__iter.objId))
            {
                Ig.readMapped(fn_img);
                tiltTemp = ABS(Ig.tilt());
            }
            if (tiltTemp < bestTilt)
            {
                bestTilt = tiltTemp;
                bestImage = fn_img;
            }
        }
        if (bestImage=="")
            REPORT_ERROR(ERR_VALUE_EMPTY,"Cannot find the image at 0 degrees");

        // Compute the mu0 and sigma0 for this image
        Image<double> I;
        I.read(bestImage);
        normalize_tomography(I(), I.tilt(), mu0, sigma0, tiltMask);
    }
}

void ProgNormalize::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    Image<double> I;
    if (apply_geo)
        I.readApplyGeo(fnImg, rowIn);
    else
        I.read(fnImg);
    I().setXmippOrigin();

    MultidimArray<double> &img=I();

    if (apply_geo)
    {
        Matrix2D<double> A;
        // Applygeo only valid for 2D images for now...
        img.checkDimension(2);

        MultidimArray< double > tmp;
        // get copy of the mask
        tmp.resizeNoCopy(bg_mask_bck);
        typeCast(bg_mask_bck, tmp);

        double outside = DIRECT_A2D_ELEM(tmp, 0, 0);

        // Instead of IS_INV for images use IS_NOT_INV for masks!
        I.getTransformationMatrix(A);
        selfApplyGeometry(BSPLINE3, tmp, A, IS_NOT_INV, DONT_WRAP, outside);

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(bg_mask)
        dAi(bg_mask,n)=(int)round(dAi(tmp,n));
    }

    double a, b;
    if (invert_contrast)
        img *= -1.;

    if (remove_black_dust || remove_white_dust)
    {
        double avg, stddev, min, max, zz;
        img.computeStats(avg, stddev, min, max);

        if ((min - avg) / stddev < thresh_black_dust && remove_black_dust)
        {
            FOR_ALL_ELEMENTS_IN_ARRAY3D(img)
            {
                zz = (A3D_ELEM(img, k, i, j) - avg) / stddev;
                if (zz < thresh_black_dust)
                    A3D_ELEM(img, k, i, j) = rnd_gaus(avg,stddev);
            }
        }

        if ((max - avg) / stddev > thresh_white_dust && remove_white_dust)
        {
            FOR_ALL_ELEMENTS_IN_ARRAY3D(img)
            {
                zz = (A3D_ELEM(img, k, i, j) - avg) / stddev;
                if (zz > thresh_white_dust)
                    A3D_ELEM(img, k, i, j) = rnd_gaus(avg,stddev);
            }
        }
    }

    double mui, sigmai;
    switch (method)
    {
    case OLDXMIPP:
        normalize_OldXmipp(img);
        break;
    case NEAR_OLDXMIPP:
        normalize_Near_OldXmipp(img, bg_mask);
        break;
    case NEWXMIPP:
        normalize_NewXmipp(img, bg_mask);
        break;
    case NEWXMIPP2:
        normalize_NewXmipp2(img, bg_mask);
        break;
    case RAMP:
        normalize_ramp(img, bg_mask);
        break;
    case NEIGHBOUR:
        normalize_remove_neighbours(img, bg_mask, thresh_neigh);
        break;
    case TOMOGRAPHY:
        normalize_tomography(img, I.tilt(), mui, sigmai, tiltMask);
        break;
    case TOMOGRAPHY0:
        normalize_tomography(img, I.tilt(), mui, sigmai, tiltMask,
                             true, mu0, sigma0);
        break;
    case MICHAEL:
        normalize_Michael(img, bg_mask);
        break;
    case RANDOM:
        a = rnd_unif(a0, aF);
        b = rnd_unif(b0, bF);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(img)
        A3D_ELEM(img, k, i, j) = a * A3D_ELEM(img, k, i, j) + b;
        break;
    case NONE:
    	break;
    }

    I.write(fnImgOut);
}
