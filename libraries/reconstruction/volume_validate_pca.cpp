/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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


#include "volume_validate_pca.h"

// Define params
void ProgVolumeValidationPCA::defineParams()
{
    //usage
    addUsageLine("Validate an obtained volume from a set of class averages");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input classes");
    addParamsLine("   -o <md_file>                : Metadata file with output information");
    addParamsLine("  [ --vol <file=\"\">]         : Input volume");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--numVols <N=5>]            : Number of intermediate volumes to generate");
    addParamsLine("  [--numClasses <N=8>]         : Number of classes to generate the intermediate volumes");
}


// Read arguments ==========================================================
void ProgVolumeValidationPCA::readParams()
{
    fnClasses = getParam("-i");
    fnOut = getParam("-o");
    fnSym = getParam("--sym");
    NVols = getIntParam("--numVols");
    NClasses = getIntParam("--numClasses");
}

// Show ====================================================================
void ProgVolumeValidationPCA::show()
{
    if (verbose > 0)
    {
        std::cout << "Input classes metadata           : "  << fnClasses  << std::endl;
        std::cout << "Output metadata                  : "  << fnOut      << std::endl;
        std::cout << "Number of intermediate volumes to generate : "  << NVols      << std::endl;
        std::cout << "Number of classes to be used               : "  << NClasses   << std::endl;
        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
    }
}

void ProgVolumeValidationPCA::produceSideinfo()
{
    mdClasses.read(fnClasses);
    mdClasses.removeDisabled();
    getImageSize(mdClasses,xdim, ydim, zdim, ndim);
}

void ProgVolumeValidationPCA::reconstructCurrent()
{
    String args=formatString("-i %s -o %s --operate random_subset %d --mode overwrite",fnClasses.c_str(),fnAngles.c_str(),NClasses);
    String cmd=(String)"xmipp_metadata_utilities "+args;
    system(cmd.c_str());

    reconstruct();
}

void ProgVolumeValidationPCA::reconstruct()
{
    //number of threads in the reconstruction
    int Nthr = 4;

    String args=formatString("-i %s -o %s --sym %s --weight --thr %d",fnAngles.c_str(),fnVol.c_str(),fnSym.c_str(),Nthr);
    String cmd=(String)"xmipp_reconstruct_fourier "+args;
    system(cmd.c_str());

    args=formatString("-i %s --mask circular %d -v 0",fnVol.c_str(),-xdim/2);
    cmd=(String)"xmipp_transform_mask "+args;
    system(cmd.c_str());

    args=formatString("-i %s --select below 0 --substitute value 0 -v 0",fnVol.c_str());
    cmd=(String)"xmipp_transform_threshold "+args;
    system(cmd.c_str());
}

void ProgVolumeValidationPCA::run()
{
    show();
    produceSideinfo();
    int index=0;

    std::stringstream ss;

    for( int index = 0; index< NVols; index++)
    {
        ss << index;
        fnAngles=fnClasses.removeAllExtensions()+ss.str();
        fnAngles+=".xmd";
        ss.str(std::string());

        ss << index;
        fnVol=fnClasses.removeAllExtensions()+ss.str();
        fnVol+=".vol";

        reconstructCurrent();
        ss.str(std::string());
    }

    Matrix2D<double> X; // Input data
    Matrix2D<double> X2; // Input data
    Image<double> img;
	Matrix1D<double> img1D;

    fnAngles=fnClasses;
    fnVol=fnClasses.removeAllExtensions();
    fnVol+=".vol";
    reconstruct();
	Matrix2D<double> C2;
	img.read(fnVol);
	img().getAliasAsRowVector(img1D);
	X2.resizeNoCopy(1,img().getSize());
	X2.setRow(index,img1D);
    matrixOperation_AAt(X2,C2);
    std::cout << "C2 : " << std::sqrt(C2(0,0)) << std::endl;

    for(index = 0; index< NVols; index++)
    {
		ss << index;
		fnVol=fnClasses.removeAllExtensions()+ss.str();
		fnVol+=".vol";
		img.read(fnVol);
		img().resize(img().getSize());

		if (index==0)
			X.resizeNoCopy(NVols,img().getSize());

		img().getAliasAsRowVector(img1D);
		X.setRow(index,img1D);
		ss.str(std::string());
    }

	subtractColumnMeans(X);
	Matrix2D<double> C, M;
	matrixOperation_AAt(X,C);
	Matrix1D<double> lambda;

	Matrix2D<double> U, V;
	Matrix1D<double> D;
	C.svd(U,D,V);

    std::cout << "Parameter C2: " << ((C2(0,0))/D.sum2())  << std::endl;
}
