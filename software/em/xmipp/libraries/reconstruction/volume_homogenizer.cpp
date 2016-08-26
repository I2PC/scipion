/***************************************************************************
 * Authors:     Mohsen Kazemi (mkazemi@cnb.csic.es)
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

#include "volume_homogenizer.h"


ProgVolumeHomogenizer::ProgVolumeHomogenizer()
{
	rank = 0;
	Nprocessors = 1;
}

void ProgVolumeHomogenizer::readParams()
{
    fnVol         = getParam("-i");
    fnRef         = getParam("-ref");
    fnSetOfImgIn  = getParam("-img");
    fnSetOfImgOut = getParam("-o");
    winSize       = getIntParam("--winSize");
    cutFreq       = getDoubleParam("--cutFreq");
    //for test
    //fnTestOut1    = getParam("-o1");
    //fnTestOut2    = getParam("-o2");
    //fnTestOut3    = getParam("-o3");
}

void ProgVolumeHomogenizer::defineParams()
{
    addUsageLine("Corrects the set of input images of the input volume with respect to the reference map and using optical flow algorithm");
    addParamsLine("  -i <selfile>             : Selfile with input volume. This is the volume with better resolution. then we will use its reformed images to join with reference map images.");
    addParamsLine("  -ref <selfile>           : Selfile with reference volume. This is the volume that we aimed to improve its resolution");
    addParamsLine("  -img <selfile>           : Selfile with input aligned images. These images are related to the input map");
    addParamsLine(" [-o <rootname=\"\">]      : Output rootname");
    addParamsLine("                           : rootname contains the list of corrected images and their angles with respect to the reference map");
    addParamsLine("                           : NOTE: you do not need to define any extension for the output root name. Program will automatically create one .stk and one .xmd");
    addParamsLine(" [--winSize <winSize=20>]  : window size for optical flow algorithm");
    addParamsLine(" [--cutFreq <cutFreq=0.5>] : cut-off frequency to use for low-pass filtering of input and reference volumes. This is digital frequency");
    //for test
    //addParamsLine(" [-o1<rootname=\"\">]      : Output rootname for test and save applied OF on reprojections of input");
    //addParamsLine(" [-o2 <rootname=\"\">]     : Output test .stk input reprojection rootname");
    //addParamsLine(" [-o3 <rootname=\"\">]     : Output test .stk reference reprojection rootname");
    //example
    /////////////////////////////////////////////////////addExampleLine("xmipp_transform_filter  -i volume.vol -o volumeFiltered.vol -f band_pass 0.1 0.3");
    //////////////tanzime khotoot va tanzimate fasele dar barname
}

// Converts a XMIPP MultidimArray to OpenCV matrix
void ProgVolumeHomogenizer::xmipp2Opencv(const MultidimArray<double> &xmippArray, cv::Mat &opencvMat)
{
    int h = YSIZE(xmippArray);
    int w = XSIZE(xmippArray);
    opencvMat = cv::Mat::zeros(h, w,CV_32FC1);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(xmippArray)
    opencvMat.at<float>(i,j) = DIRECT_A2D_ELEM(xmippArray,i,j);
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

void ProgVolumeHomogenizer::run()
{
	FileName fn_proj, fnIn;
	Image<double> inV, refV;
	Image<double> imgIn;
	MetaData setOfImgIn, setOfImgOut;
	double rot, tilt, psi, xShift, yShift;
	bool flip;
	Projection projIn, projRef, projCorr;
	cv::Mat ProjIn, ProjRef, ProjIn8, ProjRef8, ImgIn;
	cv::Mat flow, ProjCorr;
	cv::Mat planes[]={flow, flow};
	size_t objId;
	int count = 0;
	int projIdx = FIRST_IMAGE;
    FileName stackName = fnSetOfImgOut.removeAllExtensions() + ".stk";
    FileName mdName = fnSetOfImgOut.removeAllExtensions() + ".xmd";
    FourierFilter filter;

	inV.read(fnVol);
	inV().setXmippOrigin();

	refV.read(fnRef);
	refV().setXmippOrigin();

	setOfImgIn.read(fnSetOfImgIn);
	setOfImgIn.removeDisabled();

	fn_proj = stackName;
	createEmptyFile(fn_proj, XSIZE(inV()), YSIZE(inV()), 1, setOfImgIn.size(), true, WRITE_OVERWRITE);

	//for test
	//cv::Mat ProjCorr_test;
	//Projection projCorr_test;
	//FileName fn_proj_test;
	//MetaData setOfImgOut_test;
	//FileName stackName_test = fnTestOut1.removeAllExtensions() + ".stk";
	//FileName mdName_test = fnTestOut1.removeAllExtensions() + ".xmd";
	//fn_proj_test = stackName_test;
	//createEmptyFile(fn_proj_test, XSIZE(inV()), YSIZE(inV()), 1, setOfImgIn.size(), true, WRITE_OVERWRITE);


	//filtering the input and reference volumes
	filter.FilterShape = RAISED_COSINE;
	filter.FilterBand = LOWPASS;
	filter.w1 = cutFreq;
	filter.applyMaskSpace(inV());
	filter.applyMaskSpace(refV());

	//for test
	//inV.write("FinV.vol");
	//refV.write("FrefV.vol");


	//MetaData mdSetOfImgInSort;
	//size_t sz = setOfImgIn.size(), maxNImg;
	//mdSetOfImgInSort.sort(setOfImgIn,MDL_IMAGE_IDX,true,-1,0);
	//mdSetOfImgInSort.getValue(MDL_IMAGE_IDX,maxNImg,sz);
	size_t maxNImg = setOfImgIn.size();
	std::cerr << maxNImg<<"@@@@@@@@@22"<<std::endl;

	//calculating progress time
	//std::cerr<<"calculating OF and remapping for each new reprojection......\n";
	if (rank == 0)
		init_progress_bar(maxNImg);


	//size_t i = 0;
	//FOR_ALL_OBJECTS_IN_METADATA (setOfImgIn)
	for (size_t i = 1; i <= maxNImg; i++)
	{
		if ((i+1) % Nprocessors == rank)
		{

			/* if I use FOR_ALL_Object... I have to uncomment this part
			//get the coordinate of each image for reprojection purpose
			//********WARNING: be careful about applying shift, for experimental data it needs to be checked
			ApplyGeoParams p;
			p.only_apply_shifts = true;
			imgIn.readApplyGeo(setOfImgIn, __iter.objId, p);
			//setOfImgIn.getValue(MDL_IMAGE, fnIn, __iter.objId);
			//imgIn.read(fnIn);

			setOfImgIn.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
			setOfImgIn.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
			setOfImgIn.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
			//if (setOfImgIn.containsLabel(MDL_SHIFT_X))
			//	setOfImgIn.getValue(MDL_SHIFT_X, xShift, __iter.objId);
			//if (setOfImgIn.containsLabel(MDL_SHIFT_Y))
			//	setOfImgIn.getValue(MDL_SHIFT_Y, yShift, __iter.objId);
			if (setOfImgIn.containsLabel(MDL_FLIP))
				setOfImgIn.getValue(MDL_FLIP, flip, __iter.objId);

			//reprojection from input and reference volumes to calculate optical flow (OF)
			projectVolume(inV(), projIn, YSIZE(inV()), XSIZE(inV()), rot, tilt, psi);
			projectVolume(refV(), projRef, YSIZE(refV()), XSIZE(refV()), rot, tilt, psi);
			if (flip)
			{
				projIn.mirrorX();
				projRef.mirrorX();
			}

			//for test
			//projIn.write(fnTestOut2, ALL_IMAGES, true, WRITE_APPEND);
			//projRef.write(fnTestOut3, ALL_IMAGES, true, WRITE_APPEND);

			//preparing input data to use for OF algorithm
			xmipp2Opencv(imgIn(), ImgIn);
			xmipp2Opencv(projIn(), ProjIn);
			xmipp2Opencv(projRef(), ProjRef);
			convert2Uint8(ProjIn,ProjIn8);
			convert2Uint8(ProjRef,ProjRef8);

			//OF algorithm
			cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, 0.5, 6, winSize, 3, 5, 1.1, 0);
			cv::split(flow, planes);

			//for test
			//if (count == 0)
			//{
			//	std::cout<<"count: "<<count<<std::endl;
			//	std::cout<<"planes[0]"<<planes[0]<<"\n"<<std::endl;
			//	std::cout<<"planes[1]"<<planes[1]<<"\n"<<std::endl;
			//}

			for( int row = 0; row < planes[0].rows; row++ )
			    for( int col = 0; col < planes[0].cols; col++ )
			    {
			        planes[0].at<float>(row,col) += (float)col;
			        planes[1].at<float>(row,col) += (float)row;
			    }

			//applying the flow on the input image to use as the corrected one
			cv::remap(ImgIn, ProjCorr, planes[0], planes[1], cv::INTER_CUBIC);

			//for test
			//cv::remap(ProjIn, ProjCorr_test, planes[0], planes[1], cv::INTER_CUBIC);
			//opencv2Xmipp(ProjCorr_test, projCorr_test());
			//fn_proj_test.compose(projIdx, stackName_test);
			//objId = setOfImgOut_test.addObject();
			//setOfImgOut_test.setValue(MDL_IMAGE, fn_proj_test, objId);
			//setOfImgOut_test.setValue(MDL_ANGLE_ROT, rot, objId);
			//setOfImgOut_test.setValue(MDL_ANGLE_TILT, tilt, objId);
			//setOfImgOut_test.setValue(MDL_ANGLE_PSI, psi, objId);
			//setOfImgOut_test.setValue(MDL_ENABLED, 1, objId);
			//projCorr_test.write(fn_proj_test, projIdx, true, WRITE_OVERWRITE);


			//preparing output data to use for xmipp
			opencv2Xmipp(ProjCorr, projCorr());
			*/
			//get the coordinate of each image for reprojection purpose
			//********WARNING: be careful about applying shift, for experimental data it needs to be checked
			ApplyGeoParams p;
			p.only_apply_shifts = true;
			imgIn.readApplyGeo(setOfImgIn, i, p);
			//setOfImgIn.getValue(MDL_IMAGE, fnIn, __iter.objId);
			//imgIn.read(fnIn);

			setOfImgIn.getValue(MDL_ANGLE_ROT, rot, i);
			setOfImgIn.getValue(MDL_ANGLE_TILT, tilt, i);
			setOfImgIn.getValue(MDL_ANGLE_PSI, psi, i);
			//if (setOfImgIn.containsLabel(MDL_SHIFT_X))
			//	setOfImgIn.getValue(MDL_SHIFT_X, xShift, __iter.objId);
			//if (setOfImgIn.containsLabel(MDL_SHIFT_Y))
			//	setOfImgIn.getValue(MDL_SHIFT_Y, yShift, __iter.objId);
			if (setOfImgIn.containsLabel(MDL_FLIP))
				setOfImgIn.getValue(MDL_FLIP, flip, i);

			//reprojection from input and reference volumes to calculate optical flow (OF)
			projectVolume(inV(), projIn, YSIZE(inV()), XSIZE(inV()), rot, tilt, psi);
			projectVolume(refV(), projRef, YSIZE(refV()), XSIZE(refV()), rot, tilt, psi);
			if (flip)
			{
				projIn.mirrorX();
				projRef.mirrorX();
			}

			//for test
			//projIn.write(fnTestOut2, ALL_IMAGES, true, WRITE_APPEND);
			//projRef.write(fnTestOut3, ALL_IMAGES, true, WRITE_APPEND);

			//preparing input data to use for OF algorithm
			xmipp2Opencv(imgIn(), ImgIn);
			xmipp2Opencv(projIn(), ProjIn);
			xmipp2Opencv(projRef(), ProjRef);
			convert2Uint8(ProjIn,ProjIn8);
			convert2Uint8(ProjRef,ProjRef8);

			//OF algorithm
			cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, 0.5, 6, winSize, 3, 5, 1.1, 0);
			cv::split(flow, planes);

			//for test
			//if (count == 0)
			//{
			//	std::cout<<"count: "<<count<<std::endl;
			//	std::cout<<"planes[0]"<<planes[0]<<"\n"<<std::endl;
			//	std::cout<<"planes[1]"<<planes[1]<<"\n"<<std::endl;
			//}

			for( int row = 0; row < planes[0].rows; row++ )
			    for( int col = 0; col < planes[0].cols; col++ )
			    {
			        planes[0].at<float>(row,col) += (float)col;
			        planes[1].at<float>(row,col) += (float)row;
			    }

			//applying the flow on the input image to use as the corrected one
			cv::remap(ImgIn, ProjCorr, planes[0], planes[1], cv::INTER_CUBIC);

			//for test
			//cv::remap(ProjIn, ProjCorr_test, planes[0], planes[1], cv::INTER_CUBIC);
			//opencv2Xmipp(ProjCorr_test, projCorr_test());
			//fn_proj_test.compose(projIdx, stackName_test);
			//objId = setOfImgOut_test.addObject();
			//setOfImgOut_test.setValue(MDL_IMAGE, fn_proj_test, objId);
			//setOfImgOut_test.setValue(MDL_ANGLE_ROT, rot, objId);
			//setOfImgOut_test.setValue(MDL_ANGLE_TILT, tilt, objId);
			//setOfImgOut_test.setValue(MDL_ANGLE_PSI, psi, objId);
			//setOfImgOut_test.setValue(MDL_ENABLED, 1, objId);
			//projCorr_test.write(fn_proj_test, projIdx, true, WRITE_OVERWRITE);


			//preparing output data to use for xmipp
			opencv2Xmipp(ProjCorr, projCorr());


			/* change to mdPartialParticle
			//preparing output metadata included corrected images and angles
			fn_proj.compose(projIdx, stackName);
			objId = setOfImgOut.addObject();
			setOfImgOut.setValue(MDL_IMAGE, fn_proj, objId);
			setOfImgOut.setValue(MDL_ANGLE_ROT, rot, objId);
			setOfImgOut.setValue(MDL_ANGLE_TILT, tilt, objId);
			setOfImgOut.setValue(MDL_ANGLE_PSI, psi, objId);
			setOfImgOut.setValue(MDL_ENABLED, 1, objId);
			*/
			fn_proj.compose(projIdx, stackName);
			objId = mdPartialParticles.addObject();
			mdPartialParticles.setValue(MDL_IMAGE, fn_proj, objId);
			mdPartialParticles.setValue(MDL_ANGLE_ROT, rot, objId);
			mdPartialParticles.setValue(MDL_ANGLE_TILT, tilt, objId);
			mdPartialParticles.setValue(MDL_ANGLE_PSI, psi, objId);
			mdPartialParticles.setValue(MDL_ENABLED, 1, objId);


			projCorr.write(fn_proj, projIdx, true, WRITE_OVERWRITE);
			projIdx++;
		}
		//if (i <= maxNImg)
		//	i++;

		//calculating progress time
		//count++;
		//if (count%100 == 0)
		//	progress_bar(count);

		if (rank==0)
			progress_bar(i+1);
	}

	synchronize();
	gatherResults();


	if (rank == 0)
		mdPartialParticles.write(mdName);

	//setOfImgOut.write(mdName);

	//for test
	//setOfImgOut_test.write(mdName_test);

	//for test I have disabled this line...finally I think I have to uncomment it
	//progress_bar(maxNImg);
}


