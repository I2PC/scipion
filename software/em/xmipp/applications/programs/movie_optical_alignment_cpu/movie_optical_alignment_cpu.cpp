/***************************************************************************
 * Authors:     Vahid Abrishami (vabrishami@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <vector>
#include <sstream>
#include <fstream>
#include <time.h>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/video/video.hpp"

#ifdef GPU
#include "opencv2/gpu/gpu.hpp"
#endif

#include <data/multidim_array.h>
#include <data/xmipp_image.h>
#include <data/normalize.h>
#include <data/xmipp_fftw.h>

#include <reconstruction/movie_filter_dose.h>

using namespace std;
#ifdef GPU
using namespace cv::gpu;
#endif

#define NOCOMPENSATION   0
#define PRECOMPENSATION  1
#define POSTCOMPENSATION 2

class ProgOpticalAligment: public XmippProgram
{
public:
    FileName fnMovie, fnOut, fnGain, fnDark;
    FileName fnMovieOut, fnMicOut, fnMovieUncOut, fnMicUncOut;
    FileName fnMicInitial; // Filename if writing the initial average micrograph
    MetaData movie;
    int winSize, gpuDevice, nfirst, nlast, numberOfFrames;
    int finalGroupSize;
    bool globalShiftCorr, inMemory;
    double dose0, doseStep, accelerationVoltage, sampling;
    int compensationMode;

    /*****************************/
    /** crop corner **/
    /*****************************/
    /** x left top corner **/
    int xLTcorner;
    /** y left top corner **/
    int yLTcorner;
    /** x right down corner **/
    int xDRcorner;
    /** y right down corner **/
    int yDRcorner;

    // Internal variables
    FileName fnmovieRoot, fnTempStack, fnBaseOut;
    size_t Xdim, Ydim, Zdim, Ndim;
    MultidimArray<float> tempStack;
    Image<float> tempAvg;
	FourierTransformer transformer;
	ProgMovieFilterDose *filterDose;

    void defineParams()
    {
        addUsageLine ("Align movies using optical flow");
        addParamsLine("    -i <inMoviewFnName>           : input movie File Name");
        addParamsLine("   [-o <fn=\"out.xmd\">]           : Metadata with the shifts of each frame.");
        addParamsLine("   [--oavg <fn=\"\">]              : Give the name of a micrograph to generate an aligned micrograph");
        addParamsLine("   [--oavgInitial <fn=\"\">]       : Use this option to save the initial average micrograph ");
        addParamsLine("                                   : before applying any alignment. ");
        addParamsLine("   [--oUnc <fnMic=\"\"> <fnMovie=\"\">] : Give the name of a micrograph and movie to generate an aligned and dose uncompensated micrograph");
        addParamsLine("   [--cropULCorner <x=0> <y=0>]    : crop up left corner (unit=px, index starts at 0)");
        addParamsLine("   [--cropDRCorner <x=-1> <y=-1>]  : crop down right corner (unit=px, index starts at 0), -1 -> no crop");
        addParamsLine("   [--frameRange <n0=-1> <nF=-1>]  : First and last frame to align, frame numbers start at 0");
        addParamsLine("   [--bin <s=-1>]               : Binning factor, it may be any floating number");
        addParamsLine("   [--winSize <int=150>]        : window size for optical flow algorithm");
        addParamsLine("   [--groupSize <int=1>]        : the depth of pyramid for optical flow algorithm");
        addParamsLine("   [--outMovie <fn=\"\">]       : save corrected, dose compensated stack");
        addParamsLine("   [--dark <fn=\"\">]           : Dark correction image");
        addParamsLine("   [--gain <fn=\"\">]           : Gain correction image");
        addParamsLine("   [--inmemory]                 : Do not write a temporary file with the ");
        addParamsLine("   [--doseCorrection <dosePerFrame=0> <Ts=1> <kV=200> <previousDose=0> <mode=\"pre\">] : Set dosePerFrame to 0 if you do not want to correct by the dose");
        addParamsLine("                                : Dose in e/A^2, Sampling rate (Ts) in A. Valid modes are pre and post");
        addParamsLine("                                : Pre compensates before aligning and post after aligning");

#ifdef GPU
        addParamsLine("   [--gpu <int=0>]              : GPU device to be used");
#endif
    }

    void readParams()
    {
        fnMovie = getParam("-i");
        fnOut = getParam("-o");
        fnGain = getParam("--gain");
        fnDark = getParam("--dark");
        fnMicOut = getParam("--oavg");
        fnMicInitial = getParam("--oavgInitial");
        fnMovieOut = getParam("--outMovie");
        if (checkParam("--oUnc"))
        {
            fnMicUncOut = getParam("--oUnc",0);
        	fnMovieUncOut = getParam("--oUnc",1);
        }
        finalGroupSize = getIntParam("--groupSize");
        nfirst = getIntParam("--frameRange",0);
        nlast = getIntParam("--frameRange",1);
        winSize   = getIntParam("--winSize");
        xLTcorner= getIntParam("--cropULCorner",0);
        yLTcorner= getIntParam("--cropULCorner",1);
        xDRcorner = getIntParam("--cropDRCorner",0);
        yDRcorner = getIntParam("--cropDRCorner",1);
        inMemory = checkParam("--inmemory");
        doseStep = getDoubleParam("--doseCorrection",0);
        sampling = getDoubleParam("--doseCorrection",1);
        accelerationVoltage = getDoubleParam("--doseCorrection",2);
        dose0 = getDoubleParam("--doseCorrection",3);
        String aux;
        aux=getParam("--doseCorrection",4);
        if (doseStep==0)
        	compensationMode=NOCOMPENSATION;
        else
        {
        	if (aux=="pre")
        		compensationMode=PRECOMPENSATION;
        	else
        		compensationMode=POSTCOMPENSATION;
        }
#ifdef GPU
        gpuDevice = getIntParam("--gpu");
#endif
    }

    // Save a matrix which is generated by OpenCV
    void saveMat(const FileName& fnOut, const cv::Mat& M)
    {
//    	cv::FileStorage fs(fnOut.c_str(), cv::FileStorage::WRITE);
//    	fs << "flow" << M;
    	FILE *fh=fopen(fnOut.c_str(),"wb");
    	fwrite(&(M.rows),1,sizeof(int),fh);
    	fwrite(&(M.cols),1,sizeof(int),fh);
    	for (int row=0; row<M.rows; ++row)
    		fwrite(&(M.at<float>(row,0)),M.cols,sizeof(float),fh);
    	fclose(fh);
    }

    // Load a matrix which is generated by saveMat
    void readMat(const FileName& fnIn, cv::Mat& M)
    {
//    	cv::FileStorage fs(fnIn.c_str(), cv::FileStorage::READ);
//    	fs["flow"] >> M;
    	FILE *fh=fopen(fnIn.c_str(),"rb");
    	int rows, cols;
    	fread(&rows,1,sizeof(int),fh);
    	fread(&cols,1,sizeof(int),fh);
        M.create(rows, cols, CV_32FC1);
    	for (int row=0; row<rows; ++row)
    		fread(&(M.at<float>(row,0)),cols,sizeof(float),fh);
    	fclose(fh);
    }

    // Converts a XMIPP MultidimArray to OpenCV matrix
    void xmipp2Opencv(const MultidimArray<float> &xmippArray, cv::Mat &opencvMat)
    {
        int h = YSIZE(xmippArray);
        int w = XSIZE(xmippArray);
        opencvMat = cv::Mat::zeros(h, w,CV_32FC1);
    	for (int row=0; row<h; ++row)
    		memcpy(&opencvMat.at<float>(row,0),&DIRECT_A2D_ELEM(xmippArray,row,0),XSIZE(xmippArray)*sizeof(float));
    }

    // Converts an OpenCV matrix to XMIPP MultidimArray
    void opencv2Xmipp(const cv::Mat &opencvMat, MultidimArray<float> &xmippArray)
    {
        int h = opencvMat.rows;
        int w = opencvMat.cols;
        xmippArray.resizeNoCopy(h, w);
    	for (int row=0; row<h; ++row)
    		memcpy(&DIRECT_A2D_ELEM(xmippArray,row,0),&opencvMat.at<float>(row,0),XSIZE(xmippArray)*sizeof(float));
    }

    // Converts an OpenCV float matrix to an OpenCV Uint8 matrix
    void convert2Uint8(cv::Mat opencvDoubleMat, cv::Mat &opencvUintMat)
    {
        double min,max;
        cv::minMaxLoc(opencvDoubleMat, &min, &max);
        opencvDoubleMat.convertTo(opencvUintMat, CV_8U, 255.0/(max - min), -min * 255.0/(max - min));
    }

    /* Apply the window operation if crop option was used */
    void applyWindow(MultidimArray<double> &data)
    {
        if (yDRcorner!=-1)
            data.selfWindow(yLTcorner, xLTcorner, yDRcorner, xDRcorner);
    }

    // Computes the average of a number of frames in movies
    int computeAvg(size_t begin, size_t end, cv::Mat &cvAvgImg)
    {
        Image<float> frame;
        end=std::min(end,movie.size()-1);
        end=std::min(end,(size_t)nlast);
//    	std::cout << "Computing average: " << begin << " " << end << std::endl;

		for (size_t i=begin;i<=end;i++)
		{
			int actualIdx=i-nfirst;
			if (inMemory)
				frame().aliasImageInStack(tempStack,actualIdx);
			else
			{
//                    std::cout << "Reading " << fnTempStack << " " << actualIdx+1 << std::endl;
				frame.read(fnTempStack,DATA,actualIdx+1);
			}

			if (i==begin)
				tempAvg()=frame();
			else
				tempAvg()+=frame();
		}
        tempAvg()/=float(end-begin+1);
        xmipp2Opencv(tempAvg(), cvAvgImg);
        return end-begin+1;
    }

    void evaluateDisplacements(const cv::Mat *flowCurrentGroup, const cv::Mat *flowPreviousGroup, Matrix1D<double> &meanStdDev)
    {
        double sumX=0, sumY=0;
        double absSumX=0, absSumY=0;
        double sqSumX=0, sqSumY=0;
        double valSubtract;
        int h=flowCurrentGroup[0].rows;
        int w=flowCurrentGroup[0].cols;
        for(int i=0;i<h;i++)
            for(int j=0;j<w;j++)
            {
                valSubtract=flowCurrentGroup[0].at<float>(i,j)-flowPreviousGroup[0].at<float>(i,j);
                sumX+=valSubtract;
                absSumX+=abs(valSubtract);
                sqSumX+=valSubtract*valSubtract;
                valSubtract=flowCurrentGroup[1].at<float>(i,j)-flowPreviousGroup[1].at<float>(i,j);
                sumY+=valSubtract;
                absSumY+=abs(valSubtract);
                sqSumY+=valSubtract*valSubtract;
            }
        double n=h*w;
        double avgX=sumX/n;
        double avgY=sumY/n;
        meanStdDev(0)=absSumX/n;
        meanStdDev(1)=sqrt(std::max(sqSumX/n-avgX*avgX,0.0));
        meanStdDev(2)=absSumY/n;
        meanStdDev(3)=sqrt(std::max(sqSumY/n-avgY*avgY,0.0));
    }

    void removeFlows(int level)
    {
    	int g=0;
    	FileName fnAuxX=fnmovieRoot+fnBaseOut.removeLastExtension()+formatString("flowx_%d_%d.flow",level,g);
    	FileName fnAuxY=fnmovieRoot+fnBaseOut.removeLastExtension()+formatString("flowy_%d_%d.flow",level,g++);
    	while (fileExists(fnAuxX))
    	{
    		deleteFile(fnAuxX);
    		deleteFile(fnAuxY);
    		fnAuxX=fnmovieRoot+fnBaseOut.removeLastExtension()+formatString("flowx_%d_%d.flow",level,g);
    		fnAuxY=fnmovieRoot+fnBaseOut.removeLastExtension()+formatString("flowy_%d_%d.flow",level,g++);
    	}
    }

    void produceSideInfo()
    {
        fnmovieRoot=fnMovie.getDir();
        fnBaseOut=fnOut.removeDirectories();

        // If input is an stack create a metadata.
        if (fnMovie.isMetaData())
        {
            movie.read(fnMovie);
            getImageSize(movie, Xdim, Ydim, Zdim, Ndim);
        }
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

        // Get corners
        if (yDRcorner!=-1)
        {
            Xdim = xDRcorner - xLTcorner +1 ;
            Ydim = yDRcorner - yLTcorner +1 ;
        }
        if (Zdim!=1)
            REPORT_ERROR(ERR_ARG_INCORRECT,"This program is meant to align 2D frames, not 3D");

        // Dark and gain images
        Image<double> dark, gain;
        if (fnDark!="")
        {
            dark.read(fnDark);
            applyWindow(dark());
        }

        if (fnGain!="")
        {
            gain.read(fnGain);
            applyWindow(gain());
            gain() = 1.0/gain();
            double avg = gain().computeAvg();
            if (isinf(avg) || isnan(avg))
                REPORT_ERROR(ERR_ARG_INCORRECT,
                             "The input gain image is incorrect, its inverse produces infinite or nan");
        }

        // Count the number of frames to process
        numberOfFrames=Ndim;
        if (nfirst<0)
            nfirst=0;
        if (nlast<0)
            nlast=numberOfFrames-1;
        numberOfFrames=nlast-nfirst+1;

        // Initialize the stack for the output movie
        if (!fnMovieOut.isEmpty())
            createEmptyFile(fnMovieOut, Xdim, Ydim, 1, numberOfFrames, true, WRITE_REPLACE);
        if (!fnMovieUncOut.isEmpty())
            createEmptyFile(fnMovieUncOut, Xdim, Ydim, 1, numberOfFrames, true, WRITE_REPLACE);

        // Prepare stack
        if (inMemory)
        	tempStack.resizeNoCopy(numberOfFrames,1,Ydim,Xdim);
        else
        {
            fnTempStack=fnmovieRoot+"tmpMovie.stk";
            createEmptyFile(fnTempStack,Xdim,Ydim,1,numberOfFrames,true);
        }

        int currentFrameInIdx=0,currentFrameOutIdx=0;
        FileName fnFrame;
        Image<double> frameImage, translatedImage;
    	MultidimArray<float> Ifloat;
        Matrix1D<double> shift(2);
    	filterDose=new ProgMovieFilterDose(accelerationVoltage);
    	filterDose->pixel_size=sampling;
    	MultidimArray< std::complex<double> > FFTI;
        FOR_ALL_OBJECTS_IN_METADATA(movie)
        {
        	if (currentFrameInIdx>=nfirst && currentFrameInIdx<=nlast)
        	{
                movie.getValue(MDL_IMAGE, fnFrame, __iter.objId);
                frameImage.read(fnFrame);
                applyWindow(frameImage());

                if (XSIZE(dark())>0)
                	frameImage()-=dark();
                if (XSIZE(gain())>0)
                	frameImage()*=gain();

                if (movie.containsLabel(MDL_SHIFT_X))
                {
                	movie.getValue(MDL_SHIFT_X, XX(shift), __iter.objId);
                	movie.getValue(MDL_SHIFT_Y, YY(shift), __iter.objId);
                    translate(LINEAR, translatedImage(), frameImage(), shift, WRAP);
                    frameImage()=translatedImage();
                }

                if (fnMovieUncOut!="")
                	frameImage.write(fnMovieUncOut, currentFrameOutIdx+1, true, WRITE_REPLACE);

                if (compensationMode==PRECOMPENSATION)
                {
					transformer.FourierTransform(frameImage(), FFTI, false);
					filterDose->applyDoseFilterToImage(YSIZE(frameImage()), XSIZE(frameImage()), FFTI,
													  dose0+currentFrameInIdx*doseStep, dose0+(currentFrameInIdx+1)*doseStep);
					transformer.inverseFourierTransform();
                }
                if (inMemory)
                {
                	typeCast(frameImage(),Ifloat);
                	memcpy(&DIRECT_NZYX_ELEM(tempStack,currentFrameOutIdx,0,0,0),&Ifloat(0,0),MULTIDIM_SIZE(Ifloat)*sizeof(float));
                }
                else
                	frameImage.write(fnTempStack, currentFrameOutIdx+1, true, WRITE_REPLACE);
                currentFrameOutIdx++;
        	}
        	currentFrameInIdx++;
        }
    }

    void run()
    {
    	produceSideInfo();

        Matrix1D<double> meanStdev;
    	Image<float> undeformedGroupAverage, uncompensatedMic;

#ifdef GPU
        // Matrices required in GPU part
        GpuMat d_flowx, d_flowy, d_currentReference8, d_currentGroupAverage8;
#else
        cv::Mat flow;
#endif

        // Matrices required by Opencv
        cv::Mat cvCurrentReference, cvNewReference, cvCurrentGroupAverage, cvUndeformedGroupAverage, cvCurrentGroupAverage8, cvCurrentReference8;
        cv::Mat flowxCurrentGroup, flowyCurrentGroup, flowxPreviousGroup, flowyPreviousGroup;
        cv::Mat flowCurrentGroup[]={flowxCurrentGroup, flowyCurrentGroup};
        cv::Mat flowPreviousGroup[]={flowxPreviousGroup, flowyPreviousGroup};

        meanStdev.initZeros(4);
#ifdef GPU
        // Object for optical flow
        FarnebackOpticalFlow d_calc;
        setDevice(gpuDevice);

        // Initialize the parameters for optical flow structure
        d_calc.numLevels=6;
        d_calc.pyrScale=0.5;
        d_calc.fastPyramids=true;
        d_calc.winSize=winSize;
        d_calc.numIters=1;
        d_calc.polyN=5;
        d_calc.polySigma=1.1;
        d_calc.flags=0;
#endif

        computeAvg(nfirst, nlast, cvCurrentReference);
        if (!fnMicInitial.isEmpty())
            tempAvg.write(fnMicInitial);

        cout << "Frames " << nfirst << " to " << nlast << " under processing ..." << std::endl;

        int numberOfGroups=2;
        int levelNum=int(ceil(log(double(numberOfFrames)/finalGroupSize)/log(2.0))), levelCounter=0;
        MetaData MDout; // To save plot information
        while (levelCounter<levelNum)
        {
            convert2Uint8(cvCurrentReference,cvCurrentReference8);
    #ifdef GPU
    		d_currentReference8.upload(cvCurrentReference8);
    #endif
    		int currentGroupSize=int(ceil((float)numberOfFrames/numberOfGroups));

            bool lastLevel = levelCounter==levelNum-1;
            // avgStep to hold the sum of aligned frames of each group at each step
            cvNewReference=cv::Mat::zeros(Ydim, Xdim, CV_32FC1);

            cout << "Level " << levelCounter << "/" << levelNum-1
                 << " of the pyramid is under processing" << std::endl;

            // Compute time for each level
            clock_t tStart = clock();

            for (int currentGroup=0; currentGroup<numberOfGroups; currentGroup++)
            {
				int NimgsInAvg=computeAvg(currentGroup*currentGroupSize+nfirst, (currentGroup+1)*currentGroupSize+nfirst-1, cvCurrentGroupAverage);
                convert2Uint8(cvCurrentGroupAverage,cvCurrentGroupAverage8);

                if (numberOfGroups>2)
                {
                    FileName flowXFileName=fnmovieRoot+fnBaseOut.removeLastExtension()+formatString("flowx_%d_%d.flow",levelCounter-1,currentGroup/2);
                    readMat(flowXFileName.c_str(), flowCurrentGroup[0]);
//                    std::cout << "Reading flow " << flowXFileName << std::endl;

                    FileName flowYFileName=fnmovieRoot+fnBaseOut.removeLastExtension()+formatString("flowy_%d_%d.flow",levelCounter-1,currentGroup/2); // ***
                    readMat(flowYFileName.c_str(), flowCurrentGroup[1]);
                }

#ifdef GPU
                d_currentGroupAverage8.upload(cvCurrentGroupAverage8);

                if (numberOfGroups>2)
                {
                    d_flowx.upload(flowCurrentGroup[0]);
                    d_flowy.upload(flowCurrentGroup[1]);

                    d_calc.flags=cv::OPTFLOW_USE_INITIAL_FLOW;
                }
                d_calc(d_currentReference8, d_currentGroupAverage8, d_flowx, d_flowy);

                d_flowx.download(flowCurrentGroup[0]);
                d_flowy.download(flowCurrentGroup[1]);
                d_currentReference8.release();
                d_currentGroupAverage8.release();
                d_flowx.release();
                d_flowy.release();
#else
                int ofFlags=0;
                // Check if we should use the flows from the previous steps
                if (numberOfGroups==2)
                {
                    flowxCurrentGroup=cv::Mat::zeros(Ydim, Xdim, CV_32FC1);
                    flowyCurrentGroup=cv::Mat::zeros(Ydim, Xdim, CV_32FC1);
                    flow=cv::Mat::zeros(Ydim, Xdim, CV_32FC1);
                }
                else
                    ofFlags=cv::OPTFLOW_USE_INITIAL_FLOW;

                merge(flowCurrentGroup,2,flow);
                calcOpticalFlowFarneback(cvCurrentReference8, cvCurrentGroupAverage8, flow, 0.5, 6, winSize, 1, 5, 1.1, ofFlags);
                split(flow, flowCurrentGroup);
#endif

                // Save the flows if we are in the last step
                if (lastLevel)
                {
                    if (currentGroup > 0)
                    {
                    	evaluateDisplacements(flowCurrentGroup,flowPreviousGroup,meanStdev);
                        size_t id=MDout.addObject();
                        MDout.setValue(MDL_OPTICALFLOW_MEANX, double(meanStdev(0)), id);
                        MDout.setValue(MDL_OPTICALFLOW_MEANY, double(meanStdev(2)), id);
                        MDout.setValue(MDL_OPTICALFLOW_STDX, double(meanStdev(1)), id);
                        MDout.setValue(MDL_OPTICALFLOW_STDY, double(meanStdev(3)), id);
                    }
                    flowCurrentGroup[0].copyTo(flowPreviousGroup[0]);
                    flowCurrentGroup[1].copyTo(flowPreviousGroup[1]);
                }
                else
                {
					FileName flowXFileName=fnmovieRoot+fnBaseOut.removeLastExtension()+formatString("flowx_%d_%d.flow",levelCounter,currentGroup);
					FileName flowYFileName=fnmovieRoot+fnBaseOut.removeLastExtension()+formatString("flowy_%d_%d.flow",levelCounter,currentGroup);
					saveMat(flowXFileName.c_str(), flowCurrentGroup[0]);
					saveMat(flowYFileName.c_str(), flowCurrentGroup[1]);
//                    std::cout << "Saving flow " << flowXFileName << std::endl;
                }

                // Update cvNewReference
                for (float row = 0; row < flowCurrentGroup[0].rows; row++ )
                    for (float col = 0; col < flowCurrentGroup[0].cols; col++ )
                    {
						flowCurrentGroup[0].at<float>(row,col) += col;
						flowCurrentGroup[1].at<float>(row,col) += row;
                    }
                cv::remap(cvCurrentGroupAverage, cvUndeformedGroupAverage, flowCurrentGroup[0], flowCurrentGroup[1], cv::INTER_CUBIC);
                if (lastLevel)
                {
                    if (compensationMode==POSTCOMPENSATION)
                    {
                    	MultidimArray<float> If;
                    	MultidimArray<double> Id;
                    	MultidimArray<std::complex<double> > fId;
						opencv2Xmipp(cvUndeformedGroupAverage,If);
						typeCast(If,Id);
    					transformer.FourierTransform(Id, fId, false);
    					filterDose->applyDoseFilterToImage(YSIZE(Id), XSIZE(Id), fId,
    													  dose0+currentGroup*doseStep, dose0+(currentGroup+1)*doseStep);
    					transformer.inverseFourierTransform();
    					typeCast(Id,If);
    					xmipp2Opencv(If,cvUndeformedGroupAverage);
                    }
                	if (!fnMovieOut.isEmpty())
					{
						opencv2Xmipp(cvUndeformedGroupAverage,undeformedGroupAverage());
	                	// std::cout << "Writing frame " << fnMovieOut << " " << currentGroup+1 << std::endl;
						undeformedGroupAverage.write(fnMovieOut, currentGroup+1, true, WRITE_REPLACE);
					}
                	if (!fnMovieUncOut.isEmpty())
                	{
                        Image<float> frame;
                        cv::Mat cvFrame, cvUndeformedFrame;
                        frame.read(formatString("%d@%s",currentGroup+1,fnMovieUncOut.c_str()));
                        xmipp2Opencv(frame(), cvFrame);
                        cv::remap(cvFrame, cvUndeformedFrame, flowCurrentGroup[0], flowCurrentGroup[1], cv::INTER_CUBIC);
						opencv2Xmipp(cvUndeformedFrame,frame());
						frame.write(fnMovieUncOut, currentGroup+1, true, WRITE_REPLACE);

						if (currentGroup==0)
							uncompensatedMic()=frame();
						else
							uncompensatedMic()+=frame();
                	}
                }
                cvNewReference+=cvUndeformedGroupAverage;
            }

            cvCurrentReference=cvNewReference;
            cvCurrentReference*=1.0/numberOfGroups;
            printf("Processing time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
            numberOfGroups=std::min(2*numberOfGroups,numberOfFrames);
            levelCounter++;
            if (levelCounter>1)
            	removeFlows(levelCounter-2);
        }
        MDout.write(fnOut, MD_APPEND);
    	removeFlows(levelNum-1);
    	removeFlows(levelNum);

        if (!fnMicOut.isEmpty())
        {
            opencv2Xmipp(cvCurrentReference, tempAvg());
            tempAvg.write(fnMicOut);
        }

        if (!fnMicUncOut.isEmpty())
        	uncompensatedMic.write(fnMicUncOut);

        if (!inMemory)
        	deleteFile(fnTempStack);
    }
};

RUN_XMIPP_PROGRAM(ProgOpticalAligment)
