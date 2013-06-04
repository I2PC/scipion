/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

#include "matrix_dimred.h"
#include <data/matrix2d.h>

ProgDimRed::ProgDimRed()
{
	algorithm=NULL;
}

void ProgDimRed::readParams()
{
    fnIn = getParam("-i");
    fnOut = getParam("-o");
    dimRefMethod = getParam("-m");
    outputDim  = getIntParam("--dout");
    dimEstMethod = getParam("--dout",1);

    if (dimRefMethod=="LTSA" || dimRefMethod=="LLTSA" || dimRefMethod=="LPP" || dimRefMethod=="LE" || dimRefMethod=="HLLE")
    	kNN=getIntParam("-m",1);
    if (dimRefMethod=="DM" || dimRefMethod=="kPCA")
    	sigma=getDoubleParam("-m",1);
    if (dimRefMethod=="LPP" || dimRefMethod=="LE")
    	sigma=getDoubleParam("-m",2);
    if (dimRefMethod=="DM")
    	t=getDoubleParam("-m",2);
    if (dimRefMethod=="pPCA")
    	Niter=getIntParam("-m",1);
}

// Show ====================================================================
void ProgDimRed::show()
{
    if (verbose>0)
        std::cout
        << "Input metadata file:    " << fnIn          << std::endl
        << "Output metadata:        " << fnOut         << std::endl
        << "Dim Red Method:         " << dimRefMethod  << std::endl
        << "Dimension out:          " << outputDim     << std::endl
        ;
    if (dimRefMethod=="LTSA" || dimRefMethod=="LLTSA" || dimRefMethod=="LPP" || dimRefMethod=="LE" || dimRefMethod=="HLLE")
    	std::cout << "k=" << kNN << std::endl;
    if (dimRefMethod=="DM" || dimRefMethod=="kPCA" || dimRefMethod=="LPP" || dimRefMethod=="LE")
    	std::cout << "sigma=" << sigma << std::endl;
    if (dimRefMethod=="DM")
    	std::cout << "t=" << t << std::endl;
    if (dimRefMethod=="pPCA")
    	std::cout << "Niter=" << Niter << std::endl;
}

// usage ===================================================================
void ProgDimRed::defineParams()
{
    processDefaultComment("-i","-i <input>");
    processDefaultComment("-o","[-o <output=\"\">]");
    addParamsLine("  [-m <dimRefMethod=PCA>]: Dimensionality Reduction method selected");
    addParamsLine("      where <dimRefMethod>");
    addParamsLine("             PCA            : Principal Component Analysis");
    addParamsLine("             LTSA <k=12>    : Local Tangent Space Alignment, k=number of nearest neighbours");
    addParamsLine("             DM <s=1> <t=1> : Diffusion map, t=Markov random walk, s=kernel sigma");
    addParamsLine("             LLTSA <k=12>   : Linear Local Tangent Space Alignment, k=number of nearest neighbours");
    addParamsLine("             LPP <k=12> <s=1> : Linearity Preserving Projection, k=number of nearest neighbours, s=kernel sigma");
    addParamsLine("             kPCA <s=1>     : Kernel PCA, s=kernel sigma");
    addParamsLine("             pPCA <n=200>   : Probabilistic PCA, n=number of iterations");
    addParamsLine("             LE <k=7> <s=1> : Laplacian Eigenmap, k=number of nearest neighbours, s=kernel sigma");
    addParamsLine("             HLLE <k=12>    : Hessian Locally Linear Embedding, k=number of nearest neighbours");
    addParamsLine("  [--dout <d=2> <method=CorrDim>] : Output dimension. Set to -1 for automatic estimation with a specific method");
    addParamsLine("       where <method>");
    addParamsLine("                  CorrDim: Correlation dimension");
    addParamsLine("                  MLE: Maximum Likelihood Estimate");
}

// Produce Side info  ====================================================================
void ProgDimRed::produceSideInfo()
{
    if (dimRefMethod=="PCA")
    {
    	algorithm=&algorithmPCA;
    } else if (dimRefMethod=="LTSA")
    {
    	algorithm=&algorithmLTSA;
    	algorithmLTSA.setSpecificParameters(kNN);
    } else if (dimRefMethod=="DM")
    {
    	algorithm=&algorithmDiffusionMaps;
    	algorithmDiffusionMaps.setSpecificParameters(t,sigma);
    } else if (dimRefMethod=="LLTSA")
    {
    	algorithm=&algorithmLLTSASCG;
    	algorithmLLTSASCG.setSpecificParameters(kNN);
    } else if (dimRefMethod=="LPP")
    {
    	algorithm=&algorithmLPP;
    	algorithmLPP.setSpecificParameters(kNN,sigma);
    } else if (dimRefMethod=="kPCA")
    {
    	algorithm=&algorithmKernelPCA;
    	algorithmKernelPCA.setSpecificParameters(sigma);
    } else if (dimRefMethod=="pPCA")
    {
    	algorithm=&algorithmProbabilisticPCA;
    	algorithmProbabilisticPCA.setSpecificParameters(Niter);
    } else if (dimRefMethod=="LE")
    {
    	algorithm=&algorithmLaplacianEigenmap;
    	algorithmLaplacianEigenmap.setSpecificParameters(sigma,kNN);
    } else if (dimRefMethod=="HLLE")
    {
    	algorithm=&algorithmHessianLLE;
    	algorithmHessianLLE.setSpecificParameters(kNN);
    }

    algorithm->setOutputDimensionality(outputDim);
}

// Estimate dimension
void ProgDimRed::estimateDimension()
{
	outputDim=intrinsicDimensionality(X, dimEstMethod, false, algorithm->distance);
	std::cout << "Estimated dimensionality: " << outputDim << std::endl;
}

void ProgMatrixDimRed::readParams()
{
	ProgDimRed::readParams();
    inputDim  = getIntParam("--din");
    Nsamples  = getIntParam("--samples");
}

void ProgMatrixDimRed::show()
{
    if (verbose>0)
    {
    	ProgDimRed::show();
        std::cerr
        << "Dimension in:           " << inputDim      << std::endl
        << "Number of samples:      " << Nsamples      << std::endl
        ;
    }
}

void ProgMatrixDimRed::defineParams()
{
	addUsageLine("This program takes an input matrix, whose rows are individual observations and ");
	addUsageLine("projects each sample onto a lower dimensional space using the selected method");
	setDefaultComment("-i","Input matrix with data. Each observation is a row.");
	setDefaultComment("-o","Output matrix with projected data");
	ProgDimRed::defineParams();
    addParamsLine("   --din <d>             : Input dimension");
    addParamsLine("   --samples <N>         : Number of observations in the input matrix");
    addExampleLine("xmipp_matrix_dimred -i matrixIn.txt -o matrixOut.txt --din 30 --dout 2 --samples 1000");
}

// Produce side info  ======================================================
void ProgMatrixDimRed::produceSideInfo()
{
	ProgDimRed::produceSideInfo();
	X.resizeNoCopy(Nsamples,inputDim);
	X.read(fnIn);
	algorithm->setInputData(X);
}

// Run  ====================================================================
void ProgMatrixDimRed::run()
{
    show();
    produceSideInfo();
    if (outputDim<0)
    	estimateDimension();
    algorithm->reduceDimensionality();
    algorithm->getReducedData().write(fnOut);
}
