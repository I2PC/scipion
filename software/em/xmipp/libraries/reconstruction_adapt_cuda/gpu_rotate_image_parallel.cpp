/***************************************************************************
 *
 * Authors:    Amaya Jimenez      ajimenez@cnb.csic.es (2002)
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

#include "gpu_rotate_image_parallel.h"



void transformationMatrixCudaParallel(const MDRow &imageGeo, Matrix2D<double> &A, bool isVol, int splineDegree){

    double psi = 0, rot = 0., tilt = 0., shiftX = 0., shiftY = 0., shiftZ=0., scale=1.;
    bool flip;
    int dim = A.Xdim() - 1;

    imageGeo.getValue(MDL_ANGLE_PSI, psi);
    imageGeo.getValue(MDL_SHIFT_X, shiftX);
    imageGeo.getValue(MDL_SHIFT_Y, shiftY);
    imageGeo.getValue(MDL_SCALE, scale);
    imageGeo.getValue(MDL_FLIP, flip);

    psi = realWRAP(psi, 0., 360.);

    //AJ to make shifts in the same direction as transform_geometry
    if (splineDegree<BSPLINE3){
    	shiftX=-shiftX;
    	shiftY=-shiftY;
    }

    if(isVol){
    	imageGeo.getValue(MDL_ANGLE_ROT, rot);
    	imageGeo.getValue(MDL_ANGLE_TILT, tilt);
    	imageGeo.getValue(MDL_SHIFT_Z, shiftZ);
    	//AJ to make shifts in the same direction as transform_geometry
    	if(splineDegree<BSPLINE3){
    		shiftZ=-shiftZ;
    	}
    	//AJ the - sign in the angles is to make rotations in the same direction as transform_geometry
    	Euler_angles2matrix(-rot, -tilt, -psi, A, true);
    	dMij(A, 2, dim) = shiftZ;
    }else{
    	//AJ the - sign in the angle is to make rotations in the same direction as transform_geometry
    	rotation2DMatrix(-psi, A, true);
    }
    dMij(A, 0, dim) = shiftX;
    dMij(A, 1, dim) = shiftY;

    if (scale != 1.){
    	if (scale==0.) // Protection against badly formed metadatas
    		scale=1.0;
        if (dim == 2){
            M3x3_BY_CT(A, A, (double)1/scale);
        }
        else if (dim == 3){
        	M4x4_BY_CT(A, A, (double)1/scale);
        }
        dMij(A, dim, dim) = 1.;
    }

    if (flip){
    	dMij(A, 0, 0) *= -1.;
    	dMij(A, 0, 1) *= -1.;
        if (dim == 3)
        	dMij(A, 0, 2) *= -1.;
    }

}

//CUDA functions
float** cuda_rotate_image_parallel(float **image, size_t *Xdim,
		size_t *Ydim, size_t *Zdim, size_t mdInSize, double** transform_matrix, int interp, int wrap, int first_call){

	float **rotated_image;

	if (interp == BSPLINE3){
		for(int i=0; i<mdInSize; i++){
			first_call++;
			cuda_rotate_image_bspline(image[i], rotated_image[i], Xdim[i], Ydim[i], Zdim[i], transform_matrix[i], wrap, first_call);
		}
	}
	else
	{
	    if (interp == NEAREST)
	    	interp = CUDA_NEAREST;
	    else if (interp == LINEAR)
	    	interp = CUDA_LINEAR;
	    else
	    	REPORT_ERROR(ERR_PARAM_INCORRECT,"Non supported interpolation type");

		rotated_image = cuda_rotate_image_linear_nearest_parallel(image, Xdim, Ydim, Zdim, transform_matrix, interp, wrap, first_call);

	    /*for(int i=0; i<mdInSize; i++){
	    	first_call++;
	    	cuda_rotate_image_linear_nearest(image[i], rotated_image[i], Xdim[i], Ydim[i], Zdim[i], transform_matrix[i], interp, wrap, first_call);
	    }*/

	}

	return rotated_image;



}


// Read arguments ==========================================================
void ProgGpuRotateImageParallel::readParams()
{

    fn_sel = getParam("-i");
    fn_out = getParam("-o");

    String degree = getParam("--interp");
    if (degree == "spline")
    	splineDegree = BSPLINE3;
    else if (degree == "linear")
        splineDegree = LINEAR;
    else if (degree == "nearest")
        splineDegree =  NEAREST;
    else
     	REPORT_ERROR(ERR_PARAM_INCORRECT,"Non supported interpolation type");

    String mode = getParam("--wrap");
    if (mode == "border")
     	wrapMode = CUDA_BORDER;
    else if (mode == "clamp")
        wrapMode = CUDA_CLAMP;
    else if (mode == "wrap")
        wrapMode =  CUDA_WRAP;
    else
        REPORT_ERROR(ERR_PARAM_INCORRECT,"Non supported wrap mode");

    first_call=0;


}

// Show ====================================================================

void ProgGpuRotateImageParallel::show()
{
    std::cout
	<< "Input:          " << fn_sel    << std::endl
	<< "Output:          " << fn_out    << std::endl
    ;
}

// usage ===================================================================
void ProgGpuRotateImageParallel::defineParams()
{

	addParamsLine(" -i <input_file>      : Input file.");
	addParamsLine(" -o  <output_file>   : Output file.");
    addUsageLine("Computes the rotation of an image with CUDA in GPU");

    addParamsLine("[--interp <interpolation_type=spline>] : Interpolation type to be used. ");
    addParamsLine("      where <interpolation_type>");
    addParamsLine("        spline          : Use spline interpolation");
    addParamsLine("        linear          : Use bilinear/trilinear interpolation");
    addParamsLine("        nearest         : Use nearest neighbor interpolation");

    addParamsLine("[--wrap <wrap_type=border>] : Wrapping type to be used");
    addParamsLine("      where <wrap_type>");
    addParamsLine("        border          : Padding with zero values");
    addParamsLine("        clamp           : Padding with the border pixel values");
    addParamsLine("        wrap            : Padding with the image repeated, only valid for non-spline interpolation ");
}

/*
void ProgGpuRotateImageParallel::preProcess()
{

	//If zdimOut greater than 1, is a volume
    isVol = false;
    dim = isVol ? 3 : 2;

    //if(wrapMode==CUDA_WRAP &&  isVol)
    		//REPORT_ERROR(ERR_PARAM_INCORRECT,"Non supported wrap mode with volumes");

    if(wrapMode==CUDA_WRAP &&  splineDegree == BSPLINE3)
        	REPORT_ERROR(ERR_PARAM_INCORRECT,"Non supported wrap mode spline interpolation");

    A.initIdentity(dim+1);


    if(use_md!=USE_MD){
    	if(use_md==USE_BOTH){
    		isMd=true;
    	}else if(use_md==USE_COMMAND){
    		isMd=false;
    	}
    	MDRow rowGeo;
    	if (checkParam("--shift")){
    		rowGeo.setValue(MDL_SHIFT_X, getDoubleParam("--shift", 0));
    		rowGeo.setValue(MDL_SHIFT_Y, getDoubleParam("--shift", 1));
    		if(isVol){
    			rowGeo.setValue(MDL_SHIFT_Z, getDoubleParam("--shift", 2));
    		}
    	}

    	if (checkParam("--scale")){
    		rowGeo.setValue(MDL_SCALE, getDoubleParam("--scale"));
    	}

    	if (checkParam("--rotate")){
    		rowGeo.setValue(MDL_ANGLE_PSI, getDoubleParam("--rotate"));
    	}


    	//Here A has the information about shift and rotation if they are supplied by command line
    	transformationMatrixCudaParallel(rowGeo, A, isVol, splineDegree);

    	if (flip){
    		MAT_ELEM(A, 0, 0) *= -1.;
    		MAT_ELEM(A, 0, 1) *= -1.;
    		if (dim == 3)
    			MAT_ELEM(A, 0, 2) *= -1.;
    	}

    }//end if(use_md!=USE_MD)
    else if(use_md==USE_MD){
    	isMd=true;
    }
}
*/


//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgGpuRotateImageParallel::run()
{

	Image<float> Iref, Iout;
	size_t *Xdim, *Ydim, *Zdim, *Ndim;
	float **original_image_gpu;
	double **transform_matrix;

	MultidimArray<float> rotated_image;
	float **rotated_image_gpu;

	//If zdimOut greater than 1, is a volume
    isVol = false;
    dim = isVol ? 3 : 2;

	//Read input metadataFile
	FileName fnImg, fnImgOut;
	MDRow rowIn;
	size_t objIndex = 0;

	SF.read(fn_sel,NULL);

	size_t mdInSize = SF.size();
	original_image_gpu = new float *[(int)mdInSize];
	//rotated_image_gpu = new float *[(int)mdInSize];
	transform_matrix = new double *[(int)mdInSize];
	Xdim = new size_t [(int)mdInSize];
	Ydim = new size_t [(int)mdInSize];
	Zdim = new size_t [(int)mdInSize];
	Ndim = new size_t [(int)mdInSize];

	FOR_ALL_OBJECTS_IN_METADATA(SF){

		++objIndex; //increment for composing starting at 1
		SF.getRow(rowIn, objIndex);
		rowIn.getValue(MDL_IMAGE, fnImg);
		std::cerr << objIndex << ". Input image: " << fnImg << std::endl;
		Iref.read(fnImg);
		Iref.getDimensions(Xdim[objIndex-1], Ydim[objIndex-1], Zdim[objIndex-1], Ndim[objIndex-1]);
		original_image_gpu[objIndex-1] = new float[Xdim[objIndex-1]*Ydim[objIndex-1]];
		//rotated_image_gpu[objIndex-1] = new float[Xdim[objIndex-1]*Ydim[objIndex-1]];
		memcpy(original_image_gpu[objIndex-1], MULTIDIM_ARRAY(Iref()), Xdim[objIndex-1]*Ydim[objIndex-1]*sizeof(float));

		R.initIdentity(dim+1);
		transformationMatrixCudaParallel(rowIn, R, isVol, splineDegree);
		transform_matrix[objIndex-1] = new double[(dim+1)*(dim+1)];
		memcpy(transform_matrix[objIndex-1], MATRIX2D_ARRAY(R), (dim+1)*(dim+1)*sizeof(double));

		/*rotated_image.coreAllocate(Ndim[objIndex-1], Zdim[objIndex-1], Ydim[objIndex-1], Xdim[objIndex-1]);
		rotated_image_gpu = MULTIDIM_ARRAY(rotated_image);
		memcpy(rotated_image_gpu, original_image_gpu[objIndex-1], Xdim[objIndex-1]*Ydim[objIndex-1]*sizeof(float));

		Iout() = rotated_image;
		fnImgOut.compose("test", (int)objIndex, "jpg");
		Iout.write(fnImgOut);

		rotated_image.coreDeallocate();*/

	}

	rotated_image_gpu = cuda_rotate_image_parallel(original_image_gpu, Xdim, Ydim, Zdim, mdInSize, transform_matrix, splineDegree, wrapMode, first_call);

	objIndex = 0;
	FOR_ALL_OBJECTS_IN_METADATA(SF){

		++objIndex; //increment for composing starting at 1
		rotated_image.coreAllocate(Ndim[objIndex-1], Zdim[objIndex-1], Ydim[objIndex-1], Xdim[objIndex-1]);
		memcpy(MULTIDIM_ARRAY(rotated_image), rotated_image_gpu[objIndex-1], Xdim[objIndex-1]*Ydim[objIndex-1]*sizeof(float));

		Iout() = rotated_image;
		fnImgOut.compose("test", (int)objIndex, "jpg");
		Iout.write(fnImgOut);

		rotated_image.coreDeallocate();

	}

	for( int i=0 ; i < objIndex; i++){
		delete[] original_image_gpu[i];
		delete[] transform_matrix[i];
	}
	delete[] original_image_gpu;
	delete[] transform_matrix;
	delete[] Xdim;
	delete[] Ydim;
	delete[] Zdim;
	delete[] Ndim;

	freeRotatedImage(rotated_image_gpu);

}

