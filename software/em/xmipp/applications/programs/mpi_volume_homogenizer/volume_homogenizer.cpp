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

#include "../mpi_volume_homogenizer/volume_homogenizer.h"


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
}

void ProgVolumeHomogenizer::defineParams()
{
    addUsageLine("Correcting the set of input images of the input volume with ");
    addUsageLine("respect to the reference map and using optical flow algorithm");
    addParamsLine("  -i <selfile>             : Selfile with input volume. This is the volume that we aimed to deform its related particles.");
    addParamsLine("                           : Then we will use its deformed (corrected) images to merge with reference map images.");
    addParamsLine("  -ref <selfile>           : Selfile with reference volume. This is the volume that we use as the reference in OF algorithm.");
    addParamsLine("                           : This is the volume that we will be deformed");
    addParamsLine("  -img <selfile>           : Selfile with input aligned images. These images are related to the input map");
    addParamsLine(" [-o <rootname=\"\">]      : Output fileName");
    addParamsLine("                           : Output contains the list of corrected (deformed) images and their angles with respect to the reference map");
    addParamsLine("                           : NOTE: you do not need to define any extension for the output root name. ");
    addParamsLine("                           : This program will automatically create one .stk and one .xmd");
    addParamsLine(" [--winSize <winSize=50>]  : window size for optical flow algorithm");
    addParamsLine(" [--cutFreq <cutFreq=0.5>] : cut-off frequency to use for low-pass filtering of input and reference volumes. This is digital frequency");
    addExampleLine("xmipp_volume_homogenizer  -i volume1.vol -ref volume2.vol -img particlesOfVolume1.xmd -o outputParticles");
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
	filter.FilterShape = RAISED_COSINE;
	filter.FilterBand = LOWPASS;
	filter.w1 = cutFreq;
	filter.applyMaskSpace(inV());
	filter.applyMaskSpace(refV());

	//calculating progress time
	size_t maxNImg = setOfImgIn.size();
	if (rank == 0)
		init_progress_bar(maxNImg);

	for (size_t i = 0; i < maxNImg; i++)
	{
		if ((i+1) % Nprocessors == rank)
		{
			//get the coordinate of each image for reprojection purpose
			//********WARNING: be careful about applying shift, for experimental data it needs to be checked
			ApplyGeoParams p;
			p.only_apply_shifts = true;
			imgIn.readApplyGeo(setOfImgIn, i+1, p);

			setOfImgIn.getRow(rowInput,i+1);
			rowInput.getValue(MDL_ANGLE_ROT, rot);
			rowInput.getValue(MDL_ANGLE_TILT, tilt);
			rowInput.getValue(MDL_ANGLE_PSI, psi);
			if (rowInput.containsLabel(MDL_FLIP))
				rowInput.getValue(MDL_FLIP, flip);

			//reprojection from input and reference volumes to calculate optical flow (OF)
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

			//OF algorithm
			cv::calcOpticalFlowFarneback(ProjRef8, ProjIn8, flow, 0.5, 6, winSize, 3, 5, 1.1, 0);
			cv::split(flow, planes);

			for( int row = 0; row < planes[0].rows; row++ )
			    for( int col = 0; col < planes[0].cols; col++ )
			    {
			        planes[0].at<float>(row,col) += (float)col;
			        planes[1].at<float>(row,col) += (float)row;
			    }

			//applying the flow on the input image to use as the corrected one
			cv::remap(ImgIn, ProjCorr, planes[0], planes[1], cv::INTER_CUBIC);

			//preparing output data to use for xmipp
			opencv2Xmipp(ProjCorr, projCorr());

			//filling output metaData
			fn_proj.compose(i+1, stackName);
			rowInput.setValue(MDL_IMAGE, fn_proj);
			projCorr.write(fn_proj, i+1,true, WRITE_OVERWRITE);

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


