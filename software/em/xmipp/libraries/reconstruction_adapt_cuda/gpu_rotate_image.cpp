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

#include "gpu_rotate_image.h"

#include <reconstruction_cuda/cuda_gpu_rotate_image.h>
#include <reconstruction_cuda/cuda_gpu_rotate_image_v2.h>
#include <data/xmipp_image.h>
#include "data/xmipp_program.h"
#include <data/args.h>
#include "data/matrix2d.h"


//#define PI ((float)3.141592653589793)


// Read arguments ==========================================================
void ProgGpuRotateImage::readParams()
{
	fnRef = getParam("-i");
	fnOut = getParam("-o");
	ang = getParam("-ang");
	//angX = getParam("-angX");
	//angY = getParam("-angY");
	//angZ = getParam("-angZ");
	interp = getParam("-interp");

}

// Show ====================================================================
void ProgGpuRotateImage::show()
{
    std::cout
	<< "Input:          " << fnRef    << std::endl
	<< "Output:          " << fnOut    << std::endl
    << "Rotation angle for images: " << ang       << std::endl
	//<< "Rotation angle for volumes in X axis: " << angX       << std::endl
	//<< "Rotation angle for volumes in Y axis: " << angY       << std::endl
	//<< "Rotation angle for volumes in Z axis: " << angZ       << std::endl
	<< "The user can choose the interpolation method: 0 - Point, 1 - Linear, 2 - Cubic " << interp    << std::endl
    ;
}

// usage ===================================================================
void ProgGpuRotateImage::defineParams()
{
    addUsageLine("Computes the rotation of an image or a volume with CUDA in GPU");
    addParamsLine("   -i <input_file>        : Input image or volume");
    addParamsLine("   -o <output_file>        : Output image or volume");
    addParamsLine("   -ang <angle_im>     : For images: Rotation angle in degrees");
    //addParamsLine("   -angX <angle_volx>    : For volumes: Rotation angle in X axis in degrees");
    //addParamsLine("   -angY <angle_voly>    : For volumes: Rotation angle in Y axis in degrees");
    //addParamsLine("   -angZ <angle_volz>    : For volumes: Rotation angle in Z axis in degrees");
    addParamsLine("   -interp <interpolation>        : Interpolation method: 0 - Point, 1 - Linear, 2 - Cubic");

}

//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgGpuRotateImage::run()
{

	double angIm;
	int interpol;
	//int angVolX, angVolY, angVolZ;
	angIm = ang.getNumber();
	interpol = interp.getNumber();


	float radIm = PI*(float)angIm/180.0;
    std::cout << "Inside run with deg " << angIm << " and rad " << radIm << std::endl;

    Image<float> Iref, Iout;
    Iref.read(fnRef);
    size_t Xdim, Ydim, Zdim, Ndim;
    //Iref.getSize();
    Iref.getDimensions(Xdim, Ydim, Zdim, Ndim);
    std::cout << "Xdim = " << Xdim << " Ydim = " << Ydim << " Zdim = " << Zdim << " Ndim = " << Ndim << std::endl;
    //if (Zdim>1 || Ndim>1){
    if (Ndim>1){
    	REPORT_ERROR(ERR_MATRIX_DIM,"Problem with image dimensions");
    }

    double alfa = radIm;
    double beta, chi;
    if(Zdim==1){
    	beta=0;
    	chi=0;
    }else{
    	beta=radIm;
    	chi=radIm;
    }
    Matrix2D<double> rot_matrix(3,3), rotZ(3,3), rotY(3,3), rotX(3,3);

    rotZ(0,0)=cos(alfa);
    rotZ(0,1)=-1*sin(alfa);
    rotZ(1,0)=sin(alfa);
    rotZ(1,1)=cos(alfa);
    rotZ(2,2)=1;

    rotY(0,0)=cos(beta);
    rotY(0,2)=sin(beta);
    rotY(1,1)=1;
    rotY(2,0)=-1*sin(beta);
    rotY(2,2)=cos(beta);

    rotX(0,0)=1;
    rotX(1,1)=cos(chi);
    rotX(1,2)=-1*sin(beta);
    rotX(2,1)=sin(beta);
    rotX(2,2)=cos(chi);

    rot_matrix=rotZ*(rotY*rotX);
    std::cout << "rotZ " << rotZ << std::endl;
    std::cout << "rotY " << rotY << std::endl;
    std::cout << "rotX " << rotX << std::endl;
    std::cout << "rot_matrix " << rot_matrix << std::endl;
    double rot_vector[9];
    int aux=0;
    for(int i=0; i<3; i++){
    	for(int j=0; j<3; j++){
    		rot_vector[aux]=rot_matrix(i,j);
    		aux++;
    	}
    }
    for (int i=0; i++; i<9){
    	std::cout << "i " << i << " rot_vector " << rot_vector << std::endl;
    }


    MultidimArray<float> &original_image=Iref();

    MultidimArray<float> rotated_image(Zdim,Ydim,Xdim);

    float *original_image_gpu = MULTIDIM_ARRAY(original_image);
    float *rotated_image_gpu = MULTIDIM_ARRAY(rotated_image);

    if (interpol<2){
    	cuda_rotate_image(original_image_gpu, rotated_image_gpu, Xdim, Ydim, Zdim, rot_vector, interpol);
    }else if(interpol==2){
    	cuda_rotate_image_v2(original_image_gpu, rotated_image_gpu, Xdim, Ydim, Zdim, rot_vector);
    }

    Iout() = rotated_image;
    Iout.write(fnOut);


}

