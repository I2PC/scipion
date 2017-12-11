/***************************************************************************
 * Authors:     Mohsen Kazemi (mkazemi@cnb.csic.es)
 *				Javier Vargas (javier.vargasbalbuena@mcgill.ca)
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

#include "../mpi_volume_homogenizer/volume_homogenizer.h"
#define numParams 3
//# define DEBUG

ProgVolumeHomogenizer::ProgVolumeHomogenizer()
{
	rank = 0;
	Nprocessors = 1;
}

void ProgVolumeHomogenizer::readParams()
{
    fnVol         	= getParam("-i");
    fnRef         	= getParam("-ref");
    fnSetOfImgIn  	= getParam("-img");
    fnSetOfImgOut 	= getParam("-o");
    winSize       	= getIntParam("--winSize");
    cutFreq       	= getDoubleParam("--cutFreq");
    pyr_scale	  	= getDoubleParam("--pyr_scale");
    levels        	= getIntParam("--levels");
    iterations      = getIntParam("--iterations");
    addToInput 		= checkParam("--auto");
}

void ProgVolumeHomogenizer::defineParams()
{
    addUsageLine("Correcting the set of input images of the input volume with ");
    addUsageLine("respect to the reference map and using optical flow algorithm");
    addParamsLine("  -i <selfile>            		: Selfile with input volume. This is the volume that we aimed to deform its related particles.");
    addParamsLine("                           		: Then we will use its deformed (corrected) images to merge with reference map images.");
    addParamsLine("  -ref <selfile>           		: Selfile with reference volume. This is the volume that we use as the reference in OF algorithm.");
    addParamsLine("                           		: This is the volume that we will be deformed");
    addParamsLine("  -img <selfile>           		: Selfile with input aligned images. These images are related to the input map");
    addParamsLine(" [-o <rootname=\"\">]      		: Output fileName");
    addParamsLine("                           		: Output contains the list of corrected (deformed) images and their angles with respect to the reference map");
    addParamsLine("                           		: NOTE: you do not need to define any extension for the output root name. ");
    addParamsLine("                           		: This program will automatically create one .stk and one .xmd");
    addParamsLine(" [--cutFreq <cutFreq=0.5>] 		: cut-off frequency to use for low-pass filtering of input and reference volumes. This is digital frequency");
    addParamsLine(" [--pyr_scale <pyr_scale=0.5>] 	: parameter, specifying the image scale (<1) to build pyramids for each image.");
    addParamsLine("									: pyr_scale=0.5 means a classical pyramid, where each next layer is twice smaller than the previous one.");
    addParamsLine(" [--levels <levels=2>] 			: number of pyramid layers including the initial image; levels=1 means");
    addParamsLine("							  		: that no extra layers are created and only the original images are used.");
    addParamsLine(" [--winSize <winSize=50>] 		: averaging window size; larger values increase the algorithm robustness to image noise and give more chances");
    addParamsLine("							  		: for fast motion detection, but yield more blurred motion field.");
    addParamsLine(" [--iterations <iterations=10>] 	: number of iterations the algorithm does at each pyramid level.");
    addParamsLine(" [--auto] 						: Automatic determination of pyr_scale, levels and iterations");

    addExampleLine("xmipp_volume_homogenizer  -i volume1.vol -ref volume2.vol -img particlesOfVolume1.xmd -o outputParticles");
}

// Converts a XMIPP MultidimArray to OpenCV matrix
void ProgVolumeHomogenizer::xmipp2Opencv(const MultidimArray<double> &xmippArray, cv::Mat &opencvMat)
{
    int h = YSIZE(xmippArray);
    int w = XSIZE(xmippArray);
    opencvMat = cv::Mat::zeros(h, w,CV_32FC1);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    opencvMat.at<float>(i,j) = (float) DIRECT_A2D_ELEM(xmippArray,i,j);
}

// Converts an OpenCV float matrix to an OpenCV Uint8 matrix
void ProgVolumeHomogenizer::convert2Uint8(cv::Mat opencvDoubleMat, cv::Mat &opencvUintMat)
{
    double min,max;
    cv::minMaxLoc(opencvDoubleMat, &min, &max);
    opencvDoubleMat.convertTo(opencvUintMat, CV_8U, 255.0/(max - min), -min * 255.0/(max - min));
}

// Converts an OpenCV matrix to XMIPP MultidimArray
void ProgVolumeHomogenizer::opencv2Xmipp(const cv::Mat &opencvMat, MultidimArray<double> &xmippArray)
{
    int h = opencvMat.rows;
    int w = opencvMat.cols;
    xmippArray.initZeros(h, w);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    DIRECT_A2D_ELEM(xmippArray,i,j) = opencvMat.at<float>(i,j);
}

void ProgVolumeHomogenizer::parameterList(size_t imgSize, int maxParamsTrial, float* paramList)
{
	srand (time(NULL));
	int offset = 0;
	for (int i=0; i < maxParamsTrial; i++)
	{
		if (i==0)
		{
			paramList[offset] = pyr_scale;
			paramList[1+offset] = levels;
			paramList[2+offset] = iterations;
		}
		else
		{
			paramList[offset] = (float(0.5*rand()) / RAND_MAX + 0.5); // fraction of image per pyramid 0.5, 2, winSize, 15
			paramList[1+offset] = (rand() % 10+1); // number of pyramids (int)
			paramList[2+offset] = (rand() % 15 + 1); // number of iterations  (int)
		}

		offset +=  numParams;

	}

#ifdef DEBUG
	std::cout << "\n" << std::endl;
	for (int i=0; i < maxParamsTrial*numParams; i++)
	{
		std::cout << paramList[i] << std::endl;
	}
#endif

}


void ProgVolumeHomogenizer::run()
{
	FileName fn_proj, fnIn;
	Image<double> inV, refV;
	Image<double> imgIn, imgCorr;
	MetaData setOfImgIn, setOfImgOut;
	double rot, tilt, psi, xShift, yShift;
	bool flip;
	Projection projIn, projRef;
	cv::Mat ProjIn, ProjRef, ProjIn8, ProjRef8, ImgIn;
	cv::Mat flow, ImgCorr;
	cv::Mat planes[]={flow, flow};
	size_t objId;
	int count = 0;
	int projIdx = FIRST_IMAGE;
    FileName stackName = fnSetOfImgOut.removeAllExtensions() + ".stk";
    FileName mdName = fnSetOfImgOut.removeAllExtensions() + ".xmd";
    MDRow rowInput;

	inV.read(fnVol);
	inV().setXmippOrigin();

	refV.read(fnRef);
	refV().setXmippOrigin();

	setOfImgIn.read(fnSetOfImgIn);
	setOfImgIn.removeDisabled();

	fn_proj = stackName;
	createEmptyFile(fn_proj, XSIZE(inV()), YSIZE(inV()), 1, setOfImgIn.size(), true, WRITE_OVERWRITE);

	//filtering the input and reference volumes
    FourierFilter Filter;
    Filter.FilterBand=LOWPASS;
    Filter.FilterShape=RAISED_COSINE;
    Filter.w1=cutFreq;
    Filter.generateMask(inV());
    Filter.do_generate_3dmask=true;
	Filter.applyMaskSpace(inV());
	Filter.applyMaskSpace(refV());

	if (addToInput) //Do automatic OF parameter search
	{
		int maxImagsCheck = 5;
		int maxParamsTrial = 100;
		cv::Mat corr;
		double bestCorr = -1e6;
		float* paramList = new float[numParams*maxParamsTrial];

		std::cout << "Estimating optimal parameters for OF alignment" << std::endl;
		init_progress_bar(maxImagsCheck);

		size_t Xdim = 0;
		size_t Ydim = 0;
		size_t Zdim = 0;
		size_t Ndim = 0;

		ApplyGeoParams p;
		p.only_apply_shifts = true;
		imgIn.readApplyGeo(setOfImgIn, 1, p);
		imgIn.getDimensions(Xdim,Ydim,Zdim,Ndim);

		parameterList(Xdim, maxParamsTrial, paramList);

		for (size_t i = 0; i < maxImagsCheck; i++)
		{
			imgIn.readApplyGeo(setOfImgIn, i+1, p);

			setOfImgIn.getRow(rowInput,i+1);
			rowInput.getValue(MDL_ANGLE_ROT, rot);
			rowInput.getValue(MDL_ANGLE_TILT, tilt);
			rowInput.getValue(MDL_ANGLE_PSI, psi);
			if (rowInput.containsLabel(MDL_FLIP))
				rowInput.getValue(MDL_FLIP, flip);

			//Reprojection from input and reference volumes to calculate optical flow (OF)
			projectVolume(inV(), projIn, YSIZE(inV()), XSIZE(inV()), rot, tilt, psi);
			projectVolume(refV(), projRef, YSIZE(refV()), XSIZE(refV()), rot, tilt, psi);

			if (flip)
			{
				projIn.mirrorX();
				projRef.mirrorX();
			}

			//Preparing input data to use for OF algorithm
			xmipp2Opencv(projIn(), ImgIn); //We wamt to compare the projIn() with the projRef() imgIn here is only a holder to compare
			xmipp2Opencv(projIn(), ProjIn);
			xmipp2Opencv(projRef(), ProjRef);

			convert2Uint8(ProjIn,ProjIn8);
			convert2Uint8(ProjRef,ProjRef8);

			int offset = 0;
			for (size_t j = 0; j < maxParamsTrial; ++j)
			{
				cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, paramList[0+offset], paramList[1+offset],paramList[2+offset], paramList[3+j], 7, 1.5, 0);
				cv::split(flow, planes);

				for( int row = 0; row < planes[0].rows; row++ )
				    for( int col = 0; col < planes[0].cols; col++ )
				    {
				        planes[0].at<float>(row,col) += (float)col;
				        planes[1].at<float>(row,col) += (float)row;
				    }

				//applying the flow on the input image to use as the corrected one
				cv::remap(ImgIn, ImgCorr, planes[0], planes[1], cv::INTER_CUBIC);
				cv::matchTemplate(ImgCorr, ProjRef, corr, cv::TM_CCOEFF_NORMED);

				if( corr.at<float>(0,0) > bestCorr)
				{
					bestCorr = corr.at<float>(0,0);

					pyr_scale = paramList[0+offset];
					levels = paramList[1+offset];
					iterations = paramList[2+offset];
				}

				offset += numParams;

			}

			progress_bar(i+1);

		}


			std::cout << "Current parameters for OF: " << std::endl;
			std::cout << "fitness value : " << bestCorr << std::endl;
			std::cout << "Scale : " << pyr_scale << " " <<  " levels : "<< levels <<" WinSize : " << winSize << " Iterations : " << iterations << std::endl;

	}

	//calculating progress time
	size_t maxNImg = setOfImgIn.size();

	if (rank == 0)
		init_progress_bar(maxNImg);

	for (size_t i = 0; i < maxNImg; i++)
	{
		if ((i+1) % Nprocessors == rank)
		{
			ApplyGeoParams p;
			p.only_apply_shifts = true;
			imgIn.readApplyGeo(setOfImgIn, i+1, p);

			setOfImgIn.getRow(rowInput,i+1);
			rowInput.getValue(MDL_ANGLE_ROT, rot);
			rowInput.getValue(MDL_ANGLE_TILT, tilt);
			rowInput.getValue(MDL_ANGLE_PSI, psi);
			if (rowInput.containsLabel(MDL_FLIP))
				rowInput.getValue(MDL_FLIP, flip);

			//Reprojection from input and reference volumes to calculate optical flow (OF)
			projectVolume(inV(), projIn, YSIZE(inV()), XSIZE(inV()), rot, tilt, psi);
			projectVolume(refV(), projRef, YSIZE(refV()), XSIZE(refV()), rot, tilt, psi);

			if (flip)
			{
				projIn.mirrorX();
				projRef.mirrorX();
			}

			//preparing input data to use for OF algorithm
			xmipp2Opencv(imgIn(), ImgIn);
			xmipp2Opencv(projIn(), ProjIn);
			xmipp2Opencv(projRef(), ProjRef);
			convert2Uint8(ProjIn,ProjIn8);
			convert2Uint8(ProjRef,ProjRef8);

			/*OF algorithm/*
			Parameters:
			----------
			prev    	– first 8-bit single-channel input image.
			next		– second input image of the same size and the same type as prev.
			flow 		– computed flow image that has the same size as prev and type CV_32FC2.
			pyr_scale 	– parameter, specifying the image scale (<1) to build pyramids for each image; pyr_scale=0.5 means a classical pyramid, where each next layer is twice smaller than the previous one.
			levels 		– number of pyramid layers including the initial image; levels=1 means that no extra layers are created and only the original images are used.
			winsize 	– averaging window size; larger values increase the algorithm robustness to image noise and give more chances for fast motion detection, but yield more blurred motion field.
			iterations 	– number of iterations the algorithm does at each pyramid level.
			poly_n 		– size of the pixel neighborhood used to find polynomial expansion in each pixel; larger values mean that the image will be approximated with smoother surfaces, yielding more robust algorithm and more blurred motion field, typically poly_n =5 or 7.
			poly_sigma 	– standard deviation of the Gaussian that is used to smooth derivatives used as a basis for the polynomial expansion; for poly_n=5, you can set poly_sigma=1.1, for poly_n=7, a good value would be poly_sigma=1.5.
			flags 		–
			 */

			//cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, 0.9, 5, winSize, 15, 7, 1.5, 0);
			//Moshen: cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, 0.5, 6, winSize, 3, 5, 1.1, 0);
			//JV: cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, 0.5, 2, winSize, 15, 7, 1.5, 0);


			cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, pyr_scale, levels, winSize, iterations, 7, 1.5, 0);
			cv::split(flow, planes);

			for( int row = 0; row < planes[0].rows; row++ )
			    for( int col = 0; col < planes[0].cols; col++ )
			    {
			        planes[0].at<float>(row,col) += (float)col;
			        planes[1].at<float>(row,col) += (float)row;
			    }

			//applying the flow on the input image to use as the corrected one
			cv::remap(ImgIn, ImgCorr, planes[0], planes[1], cv::INTER_CUBIC);

			//preparing output data to use for xmipp
			opencv2Xmipp(ImgCorr, imgCorr());

#ifdef DEBUG
			if (i==5)
			{
				imgIn.write("test_raw.spi"); // Raw image to deform Not filtered or modifed
				projIn.write("test_projIn.spi"); //map projection of conformation to deform (Noise free and filter)
				projRef.write("test_projRef.spi"); //map projection of conformation used as reference (deform to this conformation) (Noise free and filter)

				cv::remap(ProjIn, ImgCorr, planes[0], planes[1], cv::INTER_CUBIC); //apply OF to the map projetion to see something
				opencv2Xmipp(ImgCorr, imgCorr()); //transform to Xmipp
				imgCorr.write("test_ProjIn_Deformed.spi"); //  Image deformed

				Image<double> imgFlowX, imgFlowY;
				opencv2Xmipp(planes[0], imgFlowX());
				opencv2Xmipp(planes[1], imgFlowY());
				imgFlowX.write("test_FlowX.spi");
				imgFlowY.write("test_FlowY.spi");
			}
#endif

			//filling output metaData
			fn_proj.compose(i+1, stackName);

			rowInput.setValue(MDL_IMAGE, fn_proj);
			rowInput.setValue(MDL_SHIFT_X, 0.0);
			rowInput.setValue(MDL_SHIFT_Y, 0.0);
			rowInput.setValue(MDL_SHIFT_Z, 0.0);

			imgCorr.write(fn_proj, i+1,true, WRITE_OVERWRITE);

			mdPartialParticles.addRow(rowInput);
			rowInput.clear();
		}

		if (rank==0)
			progress_bar(i+1);
	}

	synchronize();
	gatherResults();

	if (rank == 0)
		mdPartialParticles.write(mdName);

}


