/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#include "movie_alignment_correlation.h"
#include <data/metadata_extension.h>
#include <data/xmipp_fftw.h>
#include <data/filters.h>

// Read arguments ==========================================================
void ProgMovieAlignmentCorrelation::readParams()
{
    fnMovie = getParam("-i");
    fnOut = getParam("-o");
    fnDark = getParam("--dark");
    fnGain = getParam("--gain");
    maxShift = getDoubleParam("--max_shift");
    Ts = getDoubleParam("--sampling");
    maxFreq = getDoubleParam("--max_freq");
    solverIterations = getIntParam("--solverIterations");
    fnAligned = getParam("--oaligned");
    fnAvg = getParam("--oavg");
    nfirst = getIntParam("--frameRange",0);
    nlast = getIntParam("--frameRange",1);
    xLTcorner= getIntParam("--cropULCorner",0);
    yLTcorner= getIntParam("--cropULCorner",1);
    xDRcorner = getIntParam("--cropDRCorner",0);
    yDRcorner = getIntParam("--cropDRCorner",1);
    show();
}

// Show ====================================================================
void ProgMovieAlignmentCorrelation::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input movie:         " << fnMovie            << std::endl
    << "Output metadata:     " << fnOut              << std::endl
    << "Dark image:          " << fnDark             << std::endl
    << "Gain image:          " << fnGain             << std::endl
    << "Max. Shift:          " << maxShift           << std::endl
    << "Max. Scale:          " << maxFreq            << std::endl
    << "Sampling:            " << Ts                 << std::endl
    << "Solver iterations:   " << solverIterations   << std::endl
    << "Aligned movie:       " << fnAligned          << std::endl
    << "Aligned micrograph:  " << fnAvg              << std::endl
    << "Frame range:         " << nfirst << " " << nlast << std::endl
    << "Crop corners  " << "(" << xLTcorner << ", " << yLTcorner << ") "
    << "(" << xDRcorner << ", " << yDRcorner << ") "
    ;
}

// usage ===================================================================
void ProgMovieAlignmentCorrelation::defineParams()
{
    addUsageLine("Align a set of frames by cross-correlation of the frames");
    addParamsLine("   -i <metadata>               : Metadata with the list of frames to align");
    addParamsLine("  [-o <fn=\"OUT.xmd\">]               : Metadata with the shifts of each frame.");
    addParamsLine("                               : If no filename is given, the input is rewritten");
    addParamsLine("  [--max_shift <s=-1>]         : Maximum shift allowed in pixels");
    addParamsLine("  [--max_freq <s=4>]           : Maximum resolution to align (in Angstroms)");
    addParamsLine("  [--sampling <Ts=1>]          : Sampling rate (A/pixel)");
    addParamsLine("  [--solverIterations <N=2>]   : Number of robust least squares iterations");
    addParamsLine("  [--oaligned <fn=\"\">]       : Give the name of a stack if you want to generate an aligned movie");
    addParamsLine("  [--oavg <fn=\"\">]           : Give the name of a micrograph to generate an aligned micrograph");
    addParamsLine("  [--frameRange <n0=-1> <nF=-1>]  : First and last frame to process, frame numbers start at 0");
    addParamsLine("  [--cropULCorner <x=0> <y=0>]    : crop up left corner (unit=px, index starts at 0)");
    addParamsLine("  [--cropDRCorner <x=-1> <y=-1>]    : crop down right corner (unit=px, index starts at 0), -1 -> no crop");
    addParamsLine("  [--dark <fn=\"\">]           : Dark correction image");
    addParamsLine("  [--gain <fn=\"\">]           : Gain correction image");
    addExampleLine("A typical example",false);
    addExampleLine("xmipp_movie_alignment_correlation -i movie.xmd --oaligned alignedMovie.stk --oavg alignedMicrograph.mrc");
    addSeeAlsoLine("xmipp_movie_optical_alignment_cpu");
}

void computeTotalShift(int iref, int j, const Matrix1D<double> &shiftX, const Matrix1D<double> &shiftY,
                       double &totalShiftX, double &totalShiftY)
{
    totalShiftX=totalShiftY=0;
    if (iref<j)
    {
        for (int jj=j-1; jj>=iref; --jj)
        {
            totalShiftX-=shiftX(jj);
            totalShiftY-=shiftY(jj);
        }
    }
    else if (iref>j)
    {
        for (int jj=j; jj<=iref-1; ++jj)
        {
            totalShiftX+=shiftX(jj);
            totalShiftY+=shiftY(jj);
        }
    }
}

void ProgMovieAlignmentCorrelation::run()
{
    MetaData movie;
    size_t Xdim, Ydim, Zdim, Ndim;

    //if input is an stack create a metadata.
    if (fnMovie.isMetaData())
        movie.read(fnMovie);
    else
    {
        ImageGeneric movieStack;
        movieStack.read(fnMovie,HEADER);
        movieStack.getDimensions(Xdim,Ydim,Zdim,Ndim);
        if (fnMovie.getExtension()=="mrc" and Ndim ==1)
            Ndim = Zdim;
        size_t id;
        FileName fn;
        for (size_t i=0;i<Ndim;i++)
        {
            id = movie.addObject();
            fn.compose(i+FIRST_IMAGE,fnMovie);
            movie.setValue(MDL_IMAGE, fn, id);
        }
    }

    if (nfirst<0)
        nfirst=0;
    if (nlast<0)
        nlast=movie.size();

    // Determine target size of the images
    const double targetOccupancy=0.9; // Set to 1 if you want fmax maps onto 1/(2*newTs)
    newTs=targetOccupancy*maxFreq/2;
    if (newTs<Ts)
        newTs=Ts;
    getImageSize(movie,Xdim, Ydim, Zdim, Ndim);
    if (yDRcorner!=-1)
    {
        Xdim = xDRcorner - xLTcorner +1 ;
        Ydim = yDRcorner - yLTcorner +1 ;
    }

    if (Zdim!=1)
        REPORT_ERROR(ERR_ARG_INCORRECT,"This program is meant to align 2D frames, not 3D");
    double sizeFactor=Ts/newTs;
    newXdim=int(Xdim*sizeFactor);
    newYdim=int(Ydim*sizeFactor);

    // Construct 1D profile of the lowpass filter
    MultidimArray<double> lpf(newXdim/2);
    double iNewXdim=1.0/newXdim;
    double sigma=targetOccupancy/6; // So that from -targetOccupancy to targetOccupancy there is 6 sigma
    double K=-0.5/(sigma*sigma);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(lpf)
    {
        double w=i*iNewXdim;
        A1D_ELEM(lpf,i)=exp(K*(w*w));
    }

    // Compute the Fourier transform of all input images
    if (verbose)
    {
        std::cout << "Computing Fourier transform of frames ..." << std::endl;
        init_progress_bar(movie.size());
    }
    FileName fnFrame;
    Image<double> frame, cropedFrame, reducedFrame, dark, gain;
    int n=0;
    FourierTransformer transformer;
    Matrix1D<double> w(2);
    std::complex<double> zero=0;

    if (fnDark!="")
    {
    	dark.read(fnDark);
        if (yDRcorner!=-1)
            dark().selfWindow(yLTcorner, xLTcorner, yDRcorner, xDRcorner);
    }
    if (fnGain!="")
    {
    	gain.read(fnGain);
        if (yDRcorner!=-1)
            gain().selfWindow(yLTcorner, xLTcorner, yDRcorner, xDRcorner);
        gain()=1.0/gain();
    	double avg=gain().computeAvg();
    	if (isinf(avg) || isnan(avg))
    		REPORT_ERROR(ERR_ARG_INCORRECT,"The input gain image is incorrect, its inverse produces infinite or nan");
    }
    FOR_ALL_OBJECTS_IN_METADATA(movie)
    {

        if (n>=nfirst && n<=nlast)
        {
        	std::cout << fnFrame << std::endl;
            movie.getValue(MDL_IMAGE,fnFrame,__iter.objId);
            if (yDRcorner==-1)
                cropedFrame.read(fnFrame);
            else
            {
                frame.read(fnFrame);
                frame().window(cropedFrame(), yLTcorner, xLTcorner, yDRcorner, xDRcorner);
            }
            if (XSIZE(dark())>0)
            	cropedFrame()-=dark();
            if (XSIZE(gain())>0)
            	cropedFrame()*=gain();

            // Reduce the size of the input frame
            //scaleToSizeFourier(1,newYdim,newXdim,frame(),reducedFrame());
            scaleToSizeFourier(1,newYdim,newXdim,cropedFrame(),reducedFrame());

            // Now do the Fourier transform and filter
            MultidimArray< std::complex<double> > *reducedFrameFourier=new MultidimArray< std::complex<double> >;
            transformer.FourierTransform(reducedFrame(),*reducedFrameFourier,true);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*reducedFrameFourier)
            {
                FFT_IDX2DIGFREQ(i,newYdim,YY(w));
                FFT_IDX2DIGFREQ(j,newXdim,XX(w));
                double wabs=w.module();
                if (wabs>targetOccupancy)
                    A2D_ELEM(*reducedFrameFourier,i,j)=zero;
                else
                    A2D_ELEM(*reducedFrameFourier,i,j)*=lpf.interpolatedElement1D(wabs*newXdim);
            }

            frameFourier.push_back(reducedFrameFourier);
        }
        ++n;
        if (verbose)
            progress_bar(n);
    }
    if (verbose)
        progress_bar(movie.size());

    // Free useless memory
    reducedFrame.clear();
    cropedFrame.clear();
    frame.clear();

    // Now compute all shifts
    size_t N=frameFourier.size();
    Matrix2D<double> A(N*(N-1)/2,N-1);
    Matrix1D<double> bX(N*(N-1)/2), bY(N*(N-1)/2);
    if (verbose)
        std::cout << "Computing shifts between frames ..." << std::endl;
    int idx=0;
    MultidimArray<double> Mcorr;
    Mcorr.resizeNoCopy(newYdim,newXdim);
    Mcorr.setXmippOrigin();
    CorrelationAux aux;
    for (size_t i=0; i<N-1; ++i)
    {
        for (size_t j=i+1; j<N; ++j)
        {
            bestShift(*frameFourier[i],*frameFourier[j],Mcorr,bX(idx),bY(idx),aux,NULL,maxShift);
            if (verbose)
                std::cerr << "Frame " << i+nfirst << " to Frame " << j+nfirst << " -> (" << bX(idx) << "," << bY(idx) << ")\n";
            for (int ij=i; ij<j; ij++)
                A(idx,ij)=1;

            idx++;
        }
        delete frameFourier[i];
    }

    // Finally solve the equation system
    Matrix1D<double> shiftX, shiftY, ex, ey;
    WeightedLeastSquaresHelper helper;
    helper.A=A;
    helper.w.initZeros(VEC_XSIZE(bX));
    helper.w.initConstant(1);

    int it=0;
    double mean, varbX, varbY;
    bX.computeMeanAndStddev(mean,varbX);
    varbX*=varbX;
    bY.computeMeanAndStddev(mean,varbY);
    varbY*=varbY;
    if (verbose)
        std::cout << "\nSolving for the shifts ...\n";
    do
    {
        // Solve the equation system
        helper.b=bX;
        weightedLeastSquares(helper,shiftX);
        helper.b=bY;
        weightedLeastSquares(helper,shiftY);

        // Compute residuals
        ex=bX-A*shiftX;
        ey=bY-A*shiftY;

        // Compute R2
        double mean, vareX, vareY;
        ex.computeMeanAndStddev(mean,vareX);
        vareX*=vareX;
        ey.computeMeanAndStddev(mean,vareY);
        vareY*=vareY;
        double R2x=1-vareX/varbX;
        double R2y=1-vareY/varbY;
        if (verbose)
            std::cout << "Iteration " << it << " R2x=" << R2x << " R2y=" << R2y << std::endl;

        // Identify outliers
        double oldWeightSum=helper.w.sum();
        double stddeveX=sqrt(vareX);
        double stddeveY=sqrt(vareY);
        FOR_ALL_ELEMENTS_IN_MATRIX1D(ex)
        if (fabs(VEC_ELEM(ex,i))>3*stddeveX || fabs(VEC_ELEM(ey,i))>3*stddeveY)
            VEC_ELEM(helper.w,i)=0.0;
        double newWeightSum=helper.w.sum();
        if (newWeightSum==oldWeightSum)
        {
            std::cout << "No outlier found\n";
            break;
        }
        else
            std::cout << "Found " << (int)(oldWeightSum-newWeightSum) << " outliers\n";

        it++;
    }
    while (it<solverIterations);

    // Choose reference image as the minimax of shifts
    double worstShiftEver=1e38;
    int bestIref=-1;
    for (int iref=0; iref<N; ++iref)
    {
        double worstShift=-1;
        for (int j=0; j<N; ++j)
        {
            double totalShiftX, totalShiftY;
            computeTotalShift(iref, j, shiftX, shiftY,totalShiftX, totalShiftY);
            if (fabs(totalShiftX)>worstShift)
                worstShift=fabs(totalShiftX);
            if (fabs(totalShiftX)>worstShift)
                worstShift=fabs(totalShiftX);
        }
        if (worstShift<worstShiftEver)
        {
            worstShiftEver=worstShift;
            bestIref=iref;
        }
    }
    if (verbose)
        std::cout << "Reference frame: " << bestIref+1+nfirst << std::endl;

    // Compute shifts
    int j=0;
    Image<double> shiftedFrame, averageMicrograph;
    Matrix1D<double> shift(2);
    n=0;
    FOR_ALL_OBJECTS_IN_METADATA(movie)
    {
        if (n>=nfirst && n<=nlast)
        {
        	//std::cerr << "DEBUG_ROB: nfirst: " << nfirst << std::endl;
        	//std::cerr << "DEBUG_ROB: nlast: " << nlast << std::endl;
            movie.getValue(MDL_IMAGE,fnFrame,__iter.objId);
            //frame.read(fnFrame);
            if (yDRcorner==-1)
                cropedFrame.read(fnFrame);
            else
            {
                frame.read(fnFrame);
                frame().window(cropedFrame(), yLTcorner, xLTcorner, yDRcorner, xDRcorner);
            }
            if (XSIZE(dark())>0)
            	cropedFrame()-=dark();
            if (XSIZE(gain())>0)
            	cropedFrame()*=gain();

            double totalShiftX, totalShiftY;
            computeTotalShift(bestIref, j, shiftX, shiftY,XX(shift), YY(shift));
            shift/=sizeFactor;
            shift*=-1;
            movie.setValue(MDL_SHIFT_X,XX(shift),__iter.objId);
            movie.setValue(MDL_SHIFT_Y,YY(shift),__iter.objId);

            std::cout << fnFrame << " shiftX=" << XX(shift) << " shiftY=" << YY(shift) << std::endl;
            if (fnAligned!="" || fnAvg!="")
            {
                translate(BSPLINE3,shiftedFrame(),cropedFrame(),shift,WRAP);
                if (fnAligned!="")
                    shiftedFrame.write(fnAligned,j+1,true,WRITE_REPLACE);
                if (fnAvg!="")
                {
                    if (j==0)
                        averageMicrograph()=shiftedFrame();
                    else
                        averageMicrograph()+=shiftedFrame();
                }
            }

            j++;
            movie.setValue(MDL_ENABLED,1,__iter.objId);
        }
        else
        {
            movie.setValue(MDL_ENABLED,-1,__iter.objId);
            movie.setValue(MDL_SHIFT_X,0.0,__iter.objId);
            movie.setValue(MDL_SHIFT_Y,0.0,__iter.objId);
        }
        movie.setValue(MDL_WEIGHT,1.0,__iter.objId);
        n++;
    }
    if (fnAvg!="")
    {
        averageMicrograph()/=N;
        averageMicrograph.write(fnAvg);
    }
    if (fnOut=="")
        std::cerr << "DEBUG_ROB: HORROR: This cannot happend: " << std::endl;
    movie.write(fnOut);
}
