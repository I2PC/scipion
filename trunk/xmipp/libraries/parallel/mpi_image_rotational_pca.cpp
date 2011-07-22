/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2002)
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

// Translated from MATLAB code by Yoel Shkolnisky

#include "mpi_image_rotational_pca.h"
#include <data/metadata_extension.h>
#include <mpi.h>

// Empty constructor =======================================================
ProgImageRotationalPCA::ProgImageRotationalPCA(int argc, char **argv)
{
	node=new MpiNode(argc,argv);
    if (!node->isMaster())
    	verbose=0;
}

// MPI destructor
ProgImageRotationalPCA::~ProgImageRotationalPCA()
{
	delete node;
}

// Read arguments ==========================================================
void ProgImageRotationalPCA::readParams()
{
	fnIn = getParam("-i");
	fnRoot = getParam("--oroot");
	Neigen = getIntParam("--eigenvectors");
	Nits = getIntParam("--iterations");
    max_shift_change = getDoubleParam("--max_shift_change");
    psi_step = getDoubleParam("--psi_step");
    shift_step = getDoubleParam("--shift_step");
}

// Show ====================================================================
void ProgImageRotationalPCA::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input:               " << fnIn << std::endl
    << "Output root:         " << fnRoot << std::endl
    << "Number eigenvectors: " << Neigen << std::endl
    << "Number iterations:   " << Nits << std::endl
    << "Psi step:            " << psi_step << std::endl
    << "Max shift change:    " << max_shift_change << " step: " << shift_step << std::endl
    ;
}

// usage ===================================================================
void ProgImageRotationalPCA::defineParams()
{
    addUsageLine("Makes a rotational invariant representation of the image collection");
    addParamsLine("    -i <selfile>               : Selfile with experimental images");
    addParamsLine("   --oroot <rootname>          : Rootname for output");
    addParamsLine("  [--eigenvectors <N=200>]     : Number of eigenvectors");
    addParamsLine("  [--iterations <N=2>]         : Number of iterations");
    addParamsLine("  [--max_shift_change <r=0>]   : Maximum change allowed in shift");
    addParamsLine("  [--psi_step <ang=1>]         : Step in psi in degrees");
    addParamsLine("  [--shift_step <r=1>]         : Step in shift in pixels");
    addExampleLine("Typical use:",false);
    addExampleLine("xmipp_mpi_image_rotational_pca -i images.stk --oroot images_eigen");
}

// Produce side info =====================================================
void ProgImageRotationalPCA::produceSideInfo()
{
	MD.read(fnIn);
	Nimg=MD.size();
	int Ydim, Zdim;
	size_t Ndim;
	ImgSize(MD, Xdim, Ydim, Zdim, Ndim);
	Nangles=floor(360.0/psi_step);
	Nshifts=(2*max_shift_change+1)/shift_step;
	Nshifts*=Nshifts;

	Hblock.resizeNoCopy(Nangles*Nshifts,Neigen+2);
	W.resizeNoCopy(Xdim*Xdim, Neigen+2);
	Wnode.resizeNoCopy(Xdim*Xdim,Neigen+2);
	Wtranspose.resizeNoCopy(Neigen+2,Xdim*Ydim);

	if (node->isMaster())
	{
		// F (#images*#shifts*#angles) x (#eigenvectors+2)*(its+1)
		createEmptyFileWithGivenLength(fnRoot+"_matrixF.raw",Nimg*Nshifts*Nangles*(Neigen+2)*(Nits+1)*sizeof(double));

		// Initialize with random numbers between -1 and 1
		FOR_ALL_ELEMENTS_IN_MATRIX2D(W)
		MAT_ELEM(W,i,j)=rnd_unif(-1.0,1.0);

		// Send to workers
		MPI_Bcast(&MAT_ELEM(W,0,0),MAT_XSIZE(W)*MAT_YSIZE(W),MPI_DOUBLE,0,MPI_COMM_WORLD);
	}
	else
		// Receive W
		MPI_Bcast(&MAT_ELEM(W,0,0),MAT_XSIZE(W)*MAT_YSIZE(W),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

// Apply T ================================================================
void ProgImageRotationalPCA::applyT()
{
	W.initZeros();
	Wnode.initZeros();

	FileName fnImg;
	size_t idx=0;
	FOR_ALL_OBJECTS_IN_METADATA(MD)
	{
		// Process image only if it is in my rank
		if ((++idx)%(node->size)!=node->rank)
			continue;

		// Read image
		MD.getValue(MDL_IMAGE,fnImg,__iter.objId);
		I.read(fnImg);
		const MultidimArray<double> &mI=I();

		// Locate the corresponding index in Matrix H
		// and copy a block in memory to speed up calculations
		size_t Hidx=(idx-1)*Nangles*Nshifts;
		memcpy(&MAT_ELEM(Hblock,0,0),&MAT_ELEM(H,Hidx,0),MAT_XSIZE(Hblock)*MAT_YSIZE(Hblock)*sizeof(double));

		// For each rotation and shift
		int block_idx=0;
		for (double psi=0; psi<360; psi+=psi_step)
		{
			rotation2DMatrix(psi,A,true);
			for (double y=-max_shift_change; y<max_shift_change; y+=shift_step)
			{
				MAT_ELEM(A,1,2)=y;
				for (double x=-max_shift_change; x<max_shift_change; x+=shift_step, ++block_idx)
				{
					MAT_ELEM(A,0,2)=x;

					// Rotate and shift image
					applyGeometry(1,Iaux,mI,A,IS_INV,true);

					// Update Wnode
					FOR_ALL_ELEMENTS_IN_MATRIX2D(Wnode)
					MAT_ELEM(Wnode,i,j)+=DIRECT_MULTIDIM_ELEM(Iaux,i)*MAT_ELEM(Hblock,block_idx,j);
				}
			}
		}
	}
    MPI_Allreduce(MATRIX2D_ARRAY(Wnode), MATRIX2D_ARRAY(W), MAT_XSIZE(W)*MAT_YSIZE(W),
    		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

// Apply T ================================================================
void ProgImageRotationalPCA::applyTt()
{
	// Compute W transpose to accelerate memory access
	FOR_ALL_ELEMENTS_IN_MATRIX2D(Wtranspose)
    MAT_ELEM(Wtranspose,i,j) = MAT_ELEM(W,j,i);

	FileName fnImg;
	size_t idx=0;
	FOR_ALL_OBJECTS_IN_METADATA(MD)
	{
		// Process image only if it is in my rank
		if ((++idx)%(node->size)!=node->rank)
			continue;

		// Read image
		MD.getValue(MDL_IMAGE,fnImg,__iter.objId);
		I.read(fnImg);
		const MultidimArray<double> &mI=I();

		// For each rotation and shift
		int block_idx=0;
		for (double psi=0; psi<360; psi+=psi_step)
		{
			rotation2DMatrix(psi,A,true);
			for (double y=-max_shift_change; y<max_shift_change; y+=shift_step)
			{
				MAT_ELEM(A,1,2)=y;
				for (double x=-max_shift_change; x<max_shift_change; x+=shift_step, ++block_idx)
				{
					MAT_ELEM(A,0,2)=x;

					// Rotate and shift image
					applyGeometry(1,Iaux,mI,A,IS_INV,true);

					// Update Hblock
					for (int j=0; j<MAT_XSIZE(Hblock); j++)
					{
						double dotproduct=0;
						const double *ptrIaux=MULTIDIM_ARRAY(Iaux);
						const double *ptrWtranspose=&MAT_ELEM(Wtranspose,j,0);
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iaux)
						dotproduct+=(*ptrIaux++)*(*ptrWtranspose++);
						MAT_ELEM(Hblock,block_idx,j)=dotproduct;
					}
				}
			}
		}

		// Locate the corresponding index in Matrix H
		// and copy block to disk
		size_t Hidx=(idx-1)*Nangles*Nshifts;
		memcpy(&MAT_ELEM(H,Hidx,0),&MAT_ELEM(Hblock,0,0),MAT_XSIZE(Hblock)*MAT_YSIZE(Hblock)*sizeof(double));// 72720
	}
}

// Run ====================================================================
void ProgImageRotationalPCA::run()
{
    show();
    produceSideInfo();

    // Compute matrix F:
	// Set H pointing to the first block of F
	H.mapToFile(fnRoot+"_matrixF.raw",Nimg*Nshifts*Nangles,Neigen+2);
	applyTt(); // H=Tt(W)
	for (int it=0; it<Nits; it++)
	{
		applyT(); // W=T(H)
		H.mapToFile(fnRoot+"_matrixF.raw",MAT_YSIZE(H),MAT_XSIZE(H),
				    (it+1)*MAT_YSIZE(H)*MAT_XSIZE(H)*sizeof(double)); // Move H to next block
		applyTt(); // H=Tt(W)
	}
}
