/***************************************************************************
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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


#include "program_image_residuals.h"
#include <data/filters.h>

void ProgImageResiduals::defineParams()
{
    each_image_produces_an_output = true;
    addUsageLine("Analyze image residuals");
    XmippMetadataProgram::defineParams();
    addExampleLine("xmipp_image_residuals -i residuals.stk -o autocorrelations.stk --save_metadata_stack autocorrelations.xmd");
}

void ProgImageResiduals::readParams()
{
    XmippMetadataProgram::readParams();
}

void ProgImageResiduals::preProcess()
{
	i=0;
	resmean.initZeros(mdInSize);
	resvar=resmean;
}

void ProgImageResiduals::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    Image<double> img;
    img.read(fnImg);
    covarianceMatrix(img(), R);
    IR()=R;
    IR.write(fnImgOut);

    img().computeAvgStdev(A1D_ELEM(resmean,i),A1D_ELEM(resvar,i));
    i++;
}

// See formula (25) of Cherian 2013
void updateRavg(MetaData &mdR, Matrix2D<double> &Ravg)
{
	FileName fnR;
	Matrix2D<double> R, Rinv, newRavg;
	Image<double> IR;
	newRavg.initZeros(Ravg);
	FOR_ALL_OBJECTS_IN_METADATA(mdR)
	{
		mdR.getValue(MDL_IMAGE,fnR,__iter.objId);
		IR.read(fnR);
		IR().copy(R);

		R+=Ravg;
		R*=0.5;
		R.inv(Rinv);
		newRavg+=Rinv;
	}
	newRavg*=1.0/mdR.size();
	newRavg.inv(Ravg);
}

double computeCovarianceMatrixDivergence(const Matrix2D<double> &C1, const Matrix2D<double> &C2)
{
	Matrix1D<double> D;
	Matrix2D<double> P;
	Matrix2D<double> C;
	C=C1+C2;
	C*=0.5;
	firstEigs(C, MAT_XSIZE(C), D, P, false);
	double retval=0;
	FOR_ALL_ELEMENTS_IN_MATRIX1D(D)
	{
		double l=VEC_ELEM(D,i);
		if (l>1e-14)
			retval+=log(l);
	}

	C=C1*C2;
	firstEigs(C, MAT_XSIZE(C), D, P, false);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(D)
	{
		double l=VEC_ELEM(D,i);
		if (l>1e-14)
			retval-=0.5*log(l);
	}
	return retval;
}

void ProgImageResiduals::postProcess()
{
	MetaData mdR(fn_out);

	// Ravg
	Matrix2D<double> Ravg;
	Ravg.initIdentity(MAT_XSIZE(R));
	std::cerr << "Calculating covariance centroid ..." << std::endl;
	init_progress_bar(10);
	for (int i=0; i<10; i++)
	{
		progress_bar(i);
		updateRavg(mdR,Ravg);
	}
	progress_bar(10);

	// Calculate the zscore of the mean and stddev
	double mean, stddev;
	resmean.computeAvgStdev(mean,stddev);
	resmean-=mean;
	resmean/=stddev;
	resvar.computeAvgStdev(mean,stddev);
	resvar-=mean;
	resvar/=stddev;

	Image<double> IR;
	FileName fnR;
	Matrix2D<double> R;
	std::cerr << "Calculating covariance divergence ..." << std::endl;
	init_progress_bar(mdR.size());
	size_t n=0;
	FOR_ALL_OBJECTS_IN_METADATA(mdR)
	{
		mdR.getValue(MDL_IMAGE,fnR,__iter.objId);
		IR.read(fnR);
		IR().copy(R);

		double d=computeCovarianceMatrixDivergence(Ravg,R);
		mdR.setValue(MDL_ZSCORE_RESMEAN,fabs(A1D_ELEM(resmean,n)),__iter.objId);
		mdR.setValue(MDL_ZSCORE_RESVAR,fabs(A1D_ELEM(resvar,n)),__iter.objId);
		mdR.setValue(MDL_ZSCORE_RESCOV,d,__iter.objId);
		n++;
		if (n%100==0)
			progress_bar(n);
	}
	progress_bar(mdR.size());
	mdR.renameColumn(MDL_IMAGE,MDL_IMAGE2);
	mdR.write(fn_out.replaceExtension("xmd"));
}
