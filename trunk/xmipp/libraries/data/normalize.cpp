/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include "args.h"
#include "micrograph.h"
#include "mask.h"
#include "normalize.h"
#include "selfile.h"

#include <string>
#include <iostream>

/* Normalizations ---------------------------------------------------------- */
void normalize_OldXmipp(Matrix2D<double> &I)
{
    double avg, stddev, min, max;
    I.computeStats(avg, stddev, min, max);
    I -= avg;
    I /= stddev;
}

void normalize_Near_OldXmipp(Matrix2D<double> &I, const Matrix2D<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    I.computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,
                                     stddevbg);
    I -= avg;
    I /= stddevbg;
}

void normalize_OldXmipp_decomposition(Matrix2D<double> &I, const Matrix2D<int> &bg_mask,
                                     const Matrix2D<double> *mask)
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
void normalize_tomography(Matrix2D<double> &I, double tilt, double &mui,
    double &sigmai, bool tiltMask, bool tomography0, double mu0, double sigma0)
{
    const int L=2;
    double Npiece=(2*L+1)*(2*L+1);

    // Build a mask using the tilt angle
    I.setXmippOrigin();
    Matrix2D<int> mask;
    mask.initZeros(I);
    int Xdimtilt=XMIPP_MIN(FLOOR(0.5*(XSIZE(I)*cos(DEG2RAD(tilt)))),
        0.5*(XSIZE(I)-(2*L+1)));
    double N=0;
    for (int i=STARTINGY(I)+L; i<=FINISHINGY(I)-L; i++)
        for (int j=-Xdimtilt+L; j<=Xdimtilt-L;j++)
        {
            mask(i,j)=1;
            N++;
        }

    // Estimate the local variance
    Matrix2D<double> localVariance;
    localVariance.initZeros(I);
    double meanVariance=0;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
    {
        if (mask(i,j)==0) continue;
        // Center a mask of size 5x5 and estimate the variance within the mask
        double meanPiece=0, variancePiece=0;
        for (int ii=i-L; ii<=i+L; ii++)
            for (int jj=j-L; jj<=j+L; jj++)
            {
                meanPiece+=I(ii,jj);
                variancePiece+=I(ii,jj)*I(ii,jj);
            }
        meanPiece/=Npiece;
        variancePiece=variancePiece/(Npiece-1)-
            Npiece/(Npiece-1)*meanPiece*meanPiece;
        localVariance(i,j)=variancePiece;
        meanVariance+=variancePiece;
    }
    meanVariance/=N;

    // Test the hypothesis that the variance in this piece is
    // the same as the variance in the whole image
    double iFu=1/icdf_FSnedecor(4*L*L+4*L,N-1,0.975);
    double iFl=1/icdf_FSnedecor(4*L*L+4*L,N-1,0.025);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(localVariance)
    {
        if (localVariance(i,j)==0) mask(i,j)=0;
        else
        {
            double ratio=localVariance(i,j)/meanVariance;
            double thl=ratio*iFu;
            double thu=ratio*iFl;
            if (thl>1 || thu<1)
                mask(i,j)=0;
        }
    }
    #ifdef DEBUG
        ImageXmipp save;
        save()=I; save.write("PPP.xmp");
        typeCast(mask,save()); save.write("PPPmask.xmp");
        std::cout << "Press any key\n";
        char c; std::cin >> c;
    #endif

    // Compute the statistics again in the reduced mask
    double avg, stddev, min, max;
    computeStats_within_binary_mask(mask, I, min, max, avg, stddev);
    double cosTilt=cos(DEG2RAD(tilt));
    if (tomography0)
    {
        double adjustedStddev=sigma0*cosTilt;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
            if (!tiltMask || ABS(j)<Xdimtilt)
                I(i,j)=(I(i,j)/cosTilt-mu0)/adjustedStddev;
            else if (tiltMask)
                I(i,j)=0;
    }
    else
    {
        double adjustedStddev=sqrt(meanVariance)*cosTilt;
        adjustedStddev=stddev*cosTilt;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
            if (!tiltMask || ABS(j)<Xdimtilt)
                I(i,j)=(I(i,j)-avg)/adjustedStddev;
            else if (tiltMask)
                I(i,j)=0;
    }
    
    // Prepare values for returning
    mui=avg;
    sigmai=sqrt(meanVariance);
}
#undef DEBUG

void normalize_Michael(Matrix2D<double> &I, const Matrix2D<int> &bg_mask)
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

void normalize_NewXmipp(Matrix2D<double> &I, const Matrix2D<int> &bg_mask)
{
    double avgbg, stddevbg, minbg, maxbg;
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,
                                     stddevbg);
    I -= avgbg;
    I /= stddevbg;
}

void normalize_NewXmipp2(Matrix2D<double> &I, const Matrix2D<int> &bg_mask)
{
    double avg, stddev, min, max;
    double avgbg, stddevbg, minbg, maxbg;
    I.computeStats(avg, stddev, min, max);
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,
                                     stddevbg);
    I -= avgbg;
    I /= ABS(avg - avgbg);
}

void normalize_ramp(Matrix2D<double> &I, const Matrix2D<int> &bg_mask)
{
    fit_point          onepoint;
    std::vector<fit_point>  allpoints;
    double             pA, pB, pC;
    double             avgbg, stddevbg, minbg, maxbg;

    // Fit a least squares plane through the background pixels
    allpoints.clear();
    I.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
    {
        if (MAT_ELEM(bg_mask, i, j))
        {
            onepoint.x = j;
            onepoint.y = i;
            onepoint.z = MAT_ELEM(I, i, j);
            onepoint.w = 1.;
            allpoints.push_back(onepoint);
        }
    }
    least_squares_plane_fit(allpoints, pA, pB, pC);
    // Substract the plane from the image
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
    {
        MAT_ELEM(I, i, j) -= pA * j + pB * i + pC;
    }
    // Divide by the remaining std.dev. in the background region
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,
                                     stddevbg);
    I /= stddevbg;

}

void normalize_remove_neighbours(Matrix2D<double> &I, 
				 const Matrix2D<int> &bg_mask,
                                 const double &threshold)
{
    fit_point          onepoint;
    std::vector<fit_point>  allpoints;
    double             pA, pB, pC;
    double             avgbg, stddevbg, minbg, maxbg, aux, newstddev;
    double             sum1 = 0.;
    double             sum2 = 0;
    int                N = 0;

    // Fit a least squares plane through the background pixels
    allpoints.clear();
    I.setXmippOrigin();
    
    // Get initial statistics
    computeStats_within_binary_mask(bg_mask, I, minbg, maxbg, avgbg,stddevbg);

    // Fit plane through those pixels within +/- threshold*sigma
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
    {
        if (MAT_ELEM(bg_mask, i, j))
        {
	    if ( ABS(avgbg - MAT_ELEM(I, i, j)) < threshold * stddevbg)
	    {
		onepoint.x = j;
		onepoint.y = i;
		onepoint.z = MAT_ELEM(I, i, j);
		onepoint.w = 1.;
		allpoints.push_back(onepoint);
	    }
	}
    }
    least_squares_plane_fit(allpoints, pA, pB, pC);

    // Substract the plane from the image
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
    {
        MAT_ELEM(I, i, j) -= pA * j + pB * i + pC;
    }

    // Get std.dev. of the background pixels within +/- threshold*sigma 
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
    {
        if (MAT_ELEM(bg_mask, i, j))
        {
	    if ( ABS(MAT_ELEM(I, i, j)) < threshold * stddevbg)
	    {
		N++;
		sum1 +=  (double) MAT_ELEM(I, i, j);
		sum2 += ((double) MAT_ELEM(I, i, j)) * 
		    ((double) MAT_ELEM(I, i, j));
	    }
	}
    }
    // average and standard deviation
    aux = sum1 / (double) N;
    newstddev = sqrt(ABS(sum2 / N - aux*aux) * N / (N - 1));

    // Replace pixels outside +/- threshold*sigma by samples from 
    // a gaussian with avg-plane and newstddev
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
    {
        if (MAT_ELEM(bg_mask, i, j))
        {
	    if ( ABS(MAT_ELEM(I, i, j)) > threshold * stddevbg)
	    {
		// get local average
		aux = pA * j + pB * i + pC;
		MAT_ELEM(I, i, j)=rnd_gaus(aux, newstddev );
	    }
	}
    }

    // Divide the entire image by the new background
    I /= newstddev;

}

void Normalize_parameters::read(int argc, char** argv)
{
    Prog_parameters::read(argc, argv);

    // Get normalizing method
    std::string aux;
    aux = getParameter(argc, argv, "-method", "NewXmipp");

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
        REPORT_ERROR(1, "Normalize: Unknown normalizing method");

    // Normalizing a volume
    volume = checkParameter(argc, argv, "-vol");

    // Invert contrast?
    invert_contrast = checkParameter(argc, argv, "-invert");

    // Apply a mask depending on the tilt
    tiltMask = checkParameter(argc, argv, "-tiltMask");

    // Remove dust particles?
    remove_black_dust = checkParameter(argc, argv, "-remove_black_dust");
    remove_white_dust = checkParameter(argc, argv, "-remove_white_dust");
    thresh_black_dust = textToFloat(getParameter(argc, argv, "-thr_black_dust", "-3.5"));
    thresh_white_dust = textToFloat(getParameter(argc, argv, "-thr_white_dust", "3.5"));
    thresh_neigh      = textToFloat(getParameter(argc, argv, "-thr_neigh", "1.2"));

    apply_geo = false;

    // Get background mask
    if (!volume)
    {
        if (method == NEWXMIPP || method == NEWXMIPP2 || method == MICHAEL ||
            method == NEAR_OLDXMIPP || method == RAMP || method == NEIGHBOUR)
        {
            enable_mask = checkParameter(argc, argv, "-mask");
            if (enable_mask)
            {
                mask_prm.allowed_data_types = INT_MASK;
                mask_prm.read(argc, argv);
            }
            else
            {
                enable_mask = false;
                int i = paremeterPosition(argc, argv, "-background");
                if (i + 2 >= argc)
                    REPORT_ERROR(1,
                                 "Normalize: Not enough parameters after -background");

                aux = argv[i + 1];
                r  = textToInteger(argv[i + 2]);

                if (aux == "frame")
                    background_mode = FRAME;
                else if (aux == "circle")
                    background_mode = CIRCLE;
                else
                    REPORT_ERROR(1, "Normalize: Unknown background mode");
            }

            // Default is NOT to apply inverse transformation from image header to the mask
            apply_geo = checkParameter(argc, argv, "-apply_geo");
        }
        else
            background_mode = NONE;

        if (method == RANDOM)
        {
            int i = paremeterPosition(argc, argv, "-prm");
            if (i + 4 >= argc)
                REPORT_ERROR(1,
                             "Normalize_parameters::read: Not enough parameters after -prm");

            a0 = textToFloat(argv[i + 1]);
            aF = textToFloat(argv[i + 2]);
            b0 = textToFloat(argv[i + 3]);
            bF = textToFloat(argv[i + 4]);
        }

        produce_side_info();
    }
}

void Normalize_parameters::produce_side_info()
{
    int Zdim, Ydim, Xdim;
    get_input_size(Zdim, Ydim, Xdim);
    if (!volume)
    {
        if (!enable_mask)
        {
            bg_mask.resize(Ydim, Xdim);
            bg_mask.setXmippOrigin();

            switch (background_mode)
            {
            case FRAME:
                BinaryFrameMask(bg_mask, Xdim - 2 * r, Ydim - 2 * r,
                                OUTSIDE_MASK);
                break;
            case CIRCLE:
                BinaryCircularMask(bg_mask, r, OUTSIDE_MASK);
                break;
            }
        }
        else
        {
            mask_prm.generate_2Dmask(Ydim, Xdim);
            bg_mask = mask_prm.imask2D;
	    // backup a copy of the mask for apply_geo mode
	    bg_mask_bck = bg_mask;
        }
        
        // Get the parameters from the 0 degrees 
        if (method==TOMOGRAPHY0)
        {
            // Look for the image at 0 degrees
            SelFile SF;
            try {
                SF.read(fn_in);
            } catch (Xmipp_error XE)
            {
                REPORT_ERROR(1,(std::string)"There is a problem opening the selfile"+
                    fn_in+". Make sure it is a correct selfile");
            }
            SF.go_first_ACTIVE();
            double bestTilt=1000;
            FileName bestImage;
            while (!SF.eof())
            {
                FileName fn_img=SF.NextImg();
                if (fn_img=="") break;
                ImageXmipp I;
                I.read(fn_img);
                if (ABS(I.tilt())<bestTilt)
                {
                    bestTilt=ABS(I.tilt());
                    bestImage=fn_img;
                }
            }
            if (bestImage=="")
                REPORT_ERROR(1,"Cannot find the image at 0 degrees");
            
            // Compute the mu0 and sigma0 for this image
            ImageXmipp I(bestImage);
            normalize_tomography(I(), I.tilt(), mu0, sigma0, tiltMask);
        }
    }
}

void Normalize_parameters::show()
{
    Prog_parameters::show();
    if (!volume)
    {
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
    }
    else
    {
        std::cout << "Normalizing volumes: " << volume << std::endl;
        if (invert_contrast)
            std::cout << "Invert contrast "<< std::endl;

    }

    if (method == NEWXMIPP && enable_mask)
        mask_prm.show();
}

void Normalize_parameters::usage()
{
    Prog_parameters::usage();

    std::cerr << "NORMALIZATION OF VOLUMES\n"
    << "  [-vol]                    : Activate this mode\n"
    << "  [-invert]                 : Invert contrast \n";

    std::cerr << "NORMALIZATION OF IMAGES\n"
    << "  [-method <mth=NewXmipp>   : Normalizing method. Valid ones are:\n"
    << "                              OldXmipp, Near_OldXmipp, NewXmipp, Tomography, Tomography0\n"
    << "                              NewXmipp2, Michael, None, Random, Ramp, Neighbour\n"
    << "                              Methods NewXmipp, Michael, Near_OldXmipp\n"
    << "                              and Ramp need a background mask:\n"
    << "  [-background frame <r>  | : Rectangular background of r pixels\n"
    << "   -background circle <r> | : Circular background outside radius=r\n"
    << "   -mask <options>]           Use an alternative type of background mask\n"
    << "                               (see xmipp_mask for options) \n"
    << "  [-invert]                 : Invert contrast \n"
    << "  [-remove_black_dust]      : Remove black dust particles \n"
    << "  [-remove_white_dust]      : Remove white dust particles \n"
    << "  [-thr_black_dust=-3.5]    : Sigma threshold for black dust particles \n"
    << "  [-thr_white_dust=3.5]     : Sigma threshold for white dust particles \n"
    << "  [-thr_neigh]              : Sigma threshold for neighbour removal \n"
    << "  [-prm a0 aF b0 bF]        : Only in random mode. y=ax+b\n"
    << "  [-tiltMask]               : Apply a mask depending on the tilt\n";
}

void Normalize_parameters::apply_geo_mask(ImageXmipp& img)
{
    Matrix2D< double > tmp;
    // get copy of the mask 
    tmp.resize(bg_mask_bck);
    typeCast(bg_mask_bck, tmp);

    double outside = DIRECT_MAT_ELEM(tmp, 0, 0);

    // Instead of IS_INV for images use IS_NOT_INV for masks!
    tmp.selfApplyGeometryBSpline(img.get_transformation_matrix(), 3, IS_NOT_INV,
                                DONT_WRAP, outside);

    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(bg_mask)
        dMij(bg_mask,i,j)=ROUND(dMij(tmp,i,j));
}

void Normalize_parameters::apply(ImageXmipp &img)
{
    double a, b;
    if (invert_contrast)
	img() *= -1.;
    
    if (remove_black_dust || remove_white_dust)
    {
        double avg, stddev, min, max, zz;
        img().computeStats(avg, stddev, min, max);

        if ((min - avg) / stddev < thresh_black_dust && remove_black_dust)
        {
            FOR_ALL_ELEMENTS_IN_MATRIX2D(img())
            {
                zz = (img(i, j) - avg) / stddev;
                if (zz < thresh_black_dust)
                    img(i, j) = rnd_gaus(avg,stddev);
            }
        }

        if ((max - avg) / stddev > thresh_white_dust && remove_white_dust)
        {
            FOR_ALL_ELEMENTS_IN_MATRIX2D(img())
            {
                zz = (img(i, j) - avg) / stddev;
                if (zz > thresh_white_dust)
                    img(i, j) = rnd_gaus(avg,stddev);
            }
        }
    }

    double mui, sigmai;
    switch (method)
    {
    case OLDXMIPP:
        normalize_OldXmipp(img());
        break;
    case NEAR_OLDXMIPP:
        normalize_Near_OldXmipp(img(), bg_mask);
        break;
    case NEWXMIPP:
        normalize_NewXmipp(img(), bg_mask);
        break;
    case NEWXMIPP2:
        normalize_NewXmipp2(img(), bg_mask);
        break;
    case RAMP:
        normalize_ramp(img(), bg_mask);
        break;
    case NEIGHBOUR:
        normalize_remove_neighbours(img(), bg_mask, thresh_neigh);
        break;
    case TOMOGRAPHY:
        normalize_tomography(img(), img.tilt(), mui, sigmai, tiltMask);
        break;
    case TOMOGRAPHY0:
        normalize_tomography(img(), img.tilt(), mui, sigmai, tiltMask,
            true, mu0, sigma0);
        break;
    case MICHAEL:
        normalize_Michael(img(), bg_mask);
        break;
    case RANDOM:
        a = rnd_unif(a0, aF);
        b = rnd_unif(b0, bF);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(img())
            img(i, j) = a * img(i, j) + b;
        break;
    }
}
