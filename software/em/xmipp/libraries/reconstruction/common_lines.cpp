/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "common_lines.h"

#include <data/args.h>
#include <data/mask.h>
#include <data/metadata_extension.h>
#include <data/filters.h>
#include <reconstruction/fourier_filter.h>
#include <reconstruction/radon.h>
#include <fstream>
#include <iomanip>

/* Common line ------------------------------------------------------------- */
CommonLine::CommonLine()
{
    angi = angj = 0;
    distanceij = -1;
}

/* Read parameters --------------------------------------------------------- */
void ProgCommonLine::readParams() {
	fn_sel = getParam("-i");
	fn_out = getParam("-o");
	outputStyle = getParam("--ostyle");
	lpf = getDoubleParam("--lpf");
	hpf = getDoubleParam("--hpf");
	stepAng = getDoubleParam("--stepAng");
	mem = getDoubleParam("--mem");
	Nthr = getIntParam("--thr");
	scaleDistance = checkParam("--scaleDistance");
	outlierFraction = getDoubleParam("--outlierFraction");
	Nmpi = 1;
}

/* Usage ------------------------------------------------------------------- */
void ProgCommonLine::defineParams() {
	addUsageLine(
			"For every pair of images in the metadata, find the common line");
	addParamsLine("    -i <file_in>         : input selfile");
	addParamsLine("   [-o <file_out=\"commonlines.txt\">]  : Output filename");
	addParamsLine("   [--ostyle <style=matrix>]  : Output style");
	addParamsLine("          where <style>");
	addParamsLine("                matrix: Text file with the indexes of the corresponding common lines");
	addParamsLine("                line_entries: Text file with a line for each matrix entry");
	addParamsLine("   [--lpf <f=0.01>]      : low pass frequency (0<lpf<0.5)");
	addParamsLine("   [--hpf <f=0.35>]      : high pass frequency (lpf<hpf<0.5)");
	addParamsLine("                         : Set both frequencies to -1 for no filter");
	addParamsLine("   [--stepAng <s=1>]     : angular step");
	addParamsLine("   [--mem <m=1>]         : float number with the memory available in Gb");
	addParamsLine("   [--thr <thr=1>]       : number of threads for each process");
	addParamsLine("   [--scaleDistance]     : scale the output distance to 16-bit integers");
	addParamsLine("   [--outlierFraction <f=0.25>] : Fraction of expected outliers in the detection of common lines");
}

/* Side info --------------------------------------------------------------- */
void ProgCommonLine::produceSideInfo()
{
    SF.read(fn_sel);
    Nimg = SF.size();

	// Compute the number of images in each block
	size_t Ydim, Zdim, Ndim;
	getImageSize(SF, Xdim, Ydim, Zdim, Ndim);
	Nblock = FLOOR(sqrt(mem*pow(2.0,30.0)/(2*Ydim*(360/stepAng)*sizeof(double))));
	Nblock = XMIPP_MIN(Nblock,CEIL(((float)Nimg)/Nmpi));

	// Ask for memory for the common line matrix
	CommonLine dummy;
	for (int i = 0; i < Nimg * Nimg; i++)
		CLmatrix.push_back(dummy);
}

/* Show -------------------------------------------------------------------- */
void ProgCommonLine::show() {
	std::cout << "File in:       " << fn_sel << std::endl << "File out:      "
			<< fn_out << std::endl << "Lowpass:       " << lpf << std::endl
			<< "Highpass:      " << hpf << std::endl << "StepAng:       "
			<< stepAng << std::endl << "Memory(Gb):    " << mem << std::endl
			<< "Block size:    " << Nblock << std::endl << "N.Threads:     "
			<< Nthr << std::endl;
}

/* Get and prepare block --------------------------------------------------- */
struct ThreadPrepareImages {
	int myThreadID;
	ProgCommonLine * parent;
	MetaData *SFi;
	std::vector<MultidimArray<std::complex<double> > > *blockRTFs;
	std::vector<MultidimArray<double> > *blockRTs;
};

void * threadPrepareImages(void * args)
{
    ThreadPrepareImages * master = (ThreadPrepareImages *) args;
    ProgCommonLine * parent = master->parent;
    MetaData SFi = *(master->SFi);
    size_t Ydim, Xdim, Zdim, Ndim;
    getImageSize(SFi, Xdim, Ydim, Zdim, Ndim);

    MultidimArray<int> mask;
    mask.resize(Ydim, Xdim);
    mask.setXmippOrigin();
    BinaryCircularMask(mask, Xdim / 2, OUTSIDE_MASK);
    int NInsideMask = (int)(XSIZE(mask) * YSIZE(mask) - mask.sum());

    FourierFilter Filter;
    Filter.w1 = -1;
    if (parent->lpf > 0 && parent->hpf > 0)
    {
        Filter.FilterBand = BANDPASS;
        Filter.w1 = parent->lpf;
        Filter.w2 = parent->hpf;
    }
    else if (parent->lpf > 0)
    {
        Filter.FilterBand = LOWPASS;
        Filter.w1 = parent->lpf;
    }
    else if (parent->hpf > 0)
    {
        Filter.FilterBand = HIGHPASS;
        Filter.w1 = parent->hpf;
    }
    Filter.raised_w = Filter.w1 / 3;

	int ii = 0;
	bool first = true;
	Image<double> I;
	FileName fnImg;
	MultidimArray<double> RT, linei(Xdim);
	FourierTransformer transformer;
	transformer.setReal(linei);
	MultidimArray<std::complex<double> >&mlineiFourier=transformer.fFourier;
	MultidimArray<std::complex<double> > RTFourier;
	FOR_ALL_OBJECTS_IN_METADATA(SFi)
	{
		if ((ii + 1) % parent->Nthr == master->myThreadID) {
			I.readApplyGeo(SFi, __iter.objId);
			I().setXmippOrigin();
			MultidimArray<double> &mI = I();

            // Bandpass filter images
            if (Filter.w1 > 0)
            {
                if (first)
                {
                    Filter.generateMask(mI);
                    first = false;
                }
                Filter.applyMaskSpace(mI);
            }

			// Cut the image outside the largest circle
			// And compute the DC value inside the mask
			double meanInside = 0;
			FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
				if (A2D_ELEM(mask,i,j))A2D_ELEM(mI,i,j)=0;
				else
				meanInside+=A2D_ELEM(mI,i,j);
			meanInside /= NInsideMask;

            // Substract the mean inside the circle
            FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
            if (!A2D_ELEM(mask,i,j))
                A2D_ELEM(mI,i,j)-=meanInside;

			// Compute the Radon transform
			Radon_Transform(I(), parent->stepAng, RT);

			// Normalize each line in the Radon Transform so that the
			// multiplication of any two lines is actually their correlation index
			RTFourier.resize(YSIZE(RT), XSIZE(mlineiFourier));
			for (size_t i = 0; i < YSIZE(RT); i++) {
				memcpy(&(DIRECT_A1D_ELEM(linei,0)),&DIRECT_A2D_ELEM(RT,i,0),XSIZE(linei)*sizeof(double));
				linei.statisticsAdjust(0,1);
				transformer.FourierTransform();
			    transformer.fReal->initZeros();
			    transformer.inverseFourierTransform();
				memcpy(&DIRECT_A2D_ELEM(RTFourier,i,0),&(DIRECT_A1D_ELEM(mlineiFourier,0)),XSIZE(mlineiFourier)*sizeof(std::complex<double>));

				memcpy(&DIRECT_A2D_ELEM(RT,i,0),&(DIRECT_A1D_ELEM(linei,0)),XSIZE(linei)*sizeof(double));
			}

			(*(master->blockRTFs))[ii] = RTFourier;
			(*(master->blockRTs))[ii] = RT;
		}
		ii++;
	}
	return NULL;
}

void ProgCommonLine::getAndPrepareBlock(int i,
		std::vector<MultidimArray<std::complex<double> > > &blockRTFs,
		std::vector<MultidimArray<double> > &blockRTs) {
	// Get the selfile
	MetaData SFi;
	SFi.selectPart(SF, i * Nblock, Nblock);

	// Ask for space for all the block images
	int jmax = SFi.size();
	MultidimArray<std::complex<double> > dummyF;
	MultidimArray<double> dummy;
	for (int j = 0; j < jmax; j++)
	{
		blockRTFs.push_back(dummyF);
		blockRTs.push_back(dummy);
	}

	// Read and preprocess the images
	pthread_t * th_ids = new pthread_t[Nthr];
	ThreadPrepareImages * th_args = new ThreadPrepareImages[Nthr];
	for (int nt = 0; nt < Nthr; nt++) {
		// Passing parameters to each thread
		th_args[nt].parent = this;
		th_args[nt].myThreadID = nt;
		th_args[nt].SFi = &SFi;
		th_args[nt].blockRTFs = &blockRTFs;
		th_args[nt].blockRTs = &blockRTs;
		pthread_create((th_ids + nt), NULL, threadPrepareImages,
				(void *) (th_args + nt));
	}

	// Waiting for threads to finish
	for (int nt = 0; nt < Nthr; nt++)
		pthread_join(*(th_ids + nt), NULL);

    // Threads structures are not needed any more
    delete []th_ids;
    delete []th_args;
}

/* Process block ----------------------------------------------------------- */
void commonLineTwoImages(
		std::vector<MultidimArray<std::complex<double> > > &RTFsi,
		std::vector<MultidimArray<double> > &RTsi,
		int idxi,
		std::vector<MultidimArray<std::complex<double> > >&RTFsj,
		std::vector<MultidimArray<double> > &RTsj,
		int idxj,
		ProgCommonLine *parent, CommonLine &result, FourierTransformer &transformer) {
	MultidimArray<std::complex<double> > &RTFi = RTFsi[idxi];
	MultidimArray<std::complex<double> > &RTFj = RTFsj[idxj];
	MultidimArray<double> &RTi = RTsi[idxi];
	MultidimArray<double> &RTj = RTsj[idxj];

	result.distanceij = 1e60;
	result.angj = -1;
	MultidimArray<std::complex<double> > lineFi, lineFj;
	MultidimArray<double> linei, linej;
	MultidimArray<double> correlationFunction;
	int jmax;
	for (size_t ii = 0; ii < YSIZE(RTFi) / 2 + 1; ii++) {
		lineFi.aliasRow(RTFi,ii);
		linei.aliasRow(RTi,ii);

		for (size_t jj=0; jj<YSIZE(RTFj); jj++)
		{
			lineFj.aliasRow(RTFj,jj);
			linej.aliasRow(RTj,jj);

			// Compute distance between the two lines
			fast_correlation_vector(lineFi, lineFj, correlationFunction, transformer);
			correlationFunction.maxIndex(jmax);
			double distance=1-A1D_ELEM(correlationFunction,jmax);

			// Check if this is the best match
			if (distance<result.distanceij)
			{
				result.distanceij=distance;
				result.angi=ii;
				result.angj=jj;
				result.jmax=jmax;
			}
		}
	}
	result.angi *= parent->stepAng;
	result.angj *= parent->stepAng;
}

struct ThreadCompareImages {
	int myThreadID;
	ProgCommonLine * parent;
	int i;
	int j;
	std::vector<MultidimArray<std::complex<double> > > *RTFsi;
	std::vector<MultidimArray<std::complex<double> > > *RTFsj;
	std::vector<MultidimArray<double> > *RTsi;
	std::vector<MultidimArray<double> > *RTsj;
}
;

void * threadCompareImages(void * args)
{
    ThreadCompareImages * master = (ThreadCompareImages *) args;
    ProgCommonLine * parent = master->parent;

	int blockIsize = master->RTFsi->size();
	int blockJsize = master->RTFsj->size();
	FourierTransformer transformer;
	MultidimArray<double> linei;
	linei.resize(parent->Xdim);
	transformer.setReal(linei);
	for (int i = 0; i < blockIsize; i++) {
		long int ii = parent->Nblock * master->i + i;
		for (int j = 0; j < blockJsize; j++) {
			// Check if this two images have to be compared
			long int jj = parent->Nblock * master->j + j;
			if (ii >= jj)
				continue;
			if ((ii * blockJsize + jj + 1) % parent->Nthr != master->myThreadID)
				continue;

			// Effectively compare the two images
			long int idx_ij = ii * parent->Nimg + jj;
			commonLineTwoImages(
					*(master->RTFsi), *(master->RTsi), i,
					*(master->RTFsj), *(master->RTsj), j,
					parent, parent->CLmatrix[idx_ij], transformer);

			// Compute the symmetric element
			long int idx_ji = jj * parent->Nimg + ii;
			parent->CLmatrix[idx_ji].distanceij = parent->CLmatrix[idx_ij].distanceij;
			parent->CLmatrix[idx_ji].angi = parent->CLmatrix[idx_ij].angj;
			parent->CLmatrix[idx_ji].angj = parent->CLmatrix[idx_ij].angi;
			parent->CLmatrix[idx_ji].jmax = -parent->CLmatrix[idx_ij].jmax;
		}
	}
	return NULL;
}

void ProgCommonLine::processBlock(int i, int j)
{
    if (i > j)
        return;

	// Preprocess each one of the selfiles
	std::vector<MultidimArray<std::complex<double> > > RTFsi, RTFsj;
	std::vector<MultidimArray<double> > RTsi, RTsj;
	getAndPrepareBlock(i, RTFsi, RTsi);
	if (i != j)
		getAndPrepareBlock(j, RTFsj, RTsj);

	// Compare all versus all
	// Read and preprocess the images
	pthread_t * th_ids = new pthread_t[Nthr];
	ThreadCompareImages * th_args = new ThreadCompareImages[Nthr];
	for (int nt = 0; nt < Nthr; nt++) {
		// Passing parameters to each thread
		th_args[nt].parent = this;
		th_args[nt].myThreadID = nt;
		th_args[nt].i = i;
		th_args[nt].j = j;
		th_args[nt].RTFsi = &RTFsi;
		th_args[nt].RTsi = &RTsi;
		if (i != j)
		{
			th_args[nt].RTFsj = &RTFsj;
			th_args[nt].RTsj = &RTsj;
		}
		else
		{
			th_args[nt].RTFsj = &RTFsi;
			th_args[nt].RTsj = &RTsi;
		}
		pthread_create((th_ids + nt), NULL, threadCompareImages,
				(void *) (th_args + nt));
	}

    // Waiting for threads to finish
    for (int nt = 0; nt < Nthr; nt++)
        pthread_join(*(th_ids + nt), NULL);

    // Threads structures are not needed any more
    delete (th_ids);
    delete (th_args);
}

/* Qualify common lines ---------------------------------------------------- */
void ProgCommonLine::qualifyCommonLines()
{
    std::vector<double> distance;
    size_t imax=CLmatrix.size();
    distance.reserve(imax);
    double dmax=imax;
    for (size_t i=0; i<imax; ++i)
    	distance.push_back(CLmatrix[i].distanceij);

    std::sort(distance.begin(),distance.end());

    std::vector<double>::iterator dbegin=distance.begin();
    std::vector<double>::iterator low;
    for (size_t i=0; i<imax; ++i)
    {
    	low=std::lower_bound(dbegin, distance.end(), CLmatrix[i].distanceij);
    	CLmatrix[i].percentile=1.0-(double)(low-dbegin)/dmax;
    }
}

/* Write results ----------------------------------------------------------- */
void ProgCommonLine::writeResults()
{
    // Look for the minimum and maximum of the common line matrix
    double minVal = 2, maxVal = -2;
    if (scaleDistance)
    {
        for (int i = 0; i < Nimg; i++)
            for (int j = 0; j < Nimg; j++)
            {
                double val = CLmatrix[i * Nimg + j].distanceij;
                if (val > 0)
                {
                    minVal = XMIPP_MIN(minVal,val);
                    maxVal = XMIPP_MAX(maxVal,val);
                }
            }
    }

	// Write the common line matrix
	if (outputStyle == "matrix") {
		Matrix2D<int> CL(Nimg,Nimg);
		CL.initConstant(-1);
		for (int j = 1; j < Nimg; j++)
			for (int i = 0; i < j; i++) {
				int ii = i * Nimg + j;
				MAT_ELEM(CL,i,j)= (int)round(CLmatrix[ii].angi/stepAng);
				MAT_ELEM(CL,j,i)= (int)round(CLmatrix[ii].angj/stepAng);
			}
		CL.write(fn_out);
	} else {
		std::ofstream fh_out;
		fh_out.open(fn_out.c_str());
		if (!fh_out)
			REPORT_ERROR(ERR_IO_NOWRITE, fn_out);
		for (int j = 1; j < Nimg; j++)
			for (int i = 0; i < j; i++) {
				int ii = i * Nimg + j;
				if (CLmatrix[ii].distanceij > 0) {
					fh_out << j << " " << i << " ";
					if (scaleDistance)
						fh_out << round(65535*(CLmatrix[ii].distanceij-minVal)/
								(maxVal-minVal)) << " ";
					else
						fh_out << CLmatrix[ii].distanceij << " ";
					fh_out << round(CLmatrix[ii].angi/stepAng) << " "
						   << round(CLmatrix[ii].angj/stepAng) << " "
						   << CLmatrix[ii].jmax << " "
						   << CLmatrix[ii].percentile << std::endl;
				}
			}
		fh_out.close();
	}

	// Write the aligned images
    int idx=0;
	FOR_ALL_OBJECTS_IN_METADATA(SF)
	{
		SF.setValue(MDL_SHIFT_X,-shift(idx++)/2,__iter.objId); // *** FIXME: COSS Why /2?
		SF.setValue(MDL_SHIFT_Y,-shift(idx++)/2,__iter.objId);
	}
	SF.write(fn_out.insertBeforeExtension("_aligned_images"));
}

/* Main program ------------------------------------------------------------ */
void ProgCommonLine::solveForShifts()
{
	WeightedLeastSquaresHelper solver;
    int nmax=Nimg*(Nimg-1)/2;
	solver.A.resize(nmax,Nimg*2);
	solver.w.resize(nmax);
	solver.b.resize(nmax);
	int n=0;
	for (int j = 1; j < Nimg; j++)
		for (int i = 0; i < j; i++, n++) {
			int ii = i * Nimg + j;
			if (CLmatrix[ii].distanceij > 0) {
				const CommonLine& cl=CLmatrix[ii];
				double sini, cosi, sinj, cosj;
				sincos(DEG2RAD(cl.angi),&sini,&cosi);
				sincos(DEG2RAD(cl.angj),&sinj,&cosj);
				int idxi=2*i;
				int idxj=2*j;
				MAT_ELEM(solver.A,n,idxi)=cosi;
				MAT_ELEM(solver.A,n,idxi+1)=sini;
				MAT_ELEM(solver.A,n,idxj)=-cosj;
				MAT_ELEM(solver.A,n,idxj+1)=-sinj;
				VEC_ELEM(solver.b,n)=cl.jmax;
				VEC_ELEM(solver.w,n)=cl.percentile;
			}
		}

	const double maxShiftError=1;
	ransacWeightedLeastSquares(solver,shift,maxShiftError,10000,outlierFraction,Nthr);
}

/* Main program ------------------------------------------------------------ */
void ProgCommonLine::run()
{
	produceSideInfo();

    // Process all blocks
    int maxBlock = CEIL(((float)Nimg)/Nblock);
    if (rank == 0)
        init_progress_bar(maxBlock * maxBlock);
    for (int i = 0; i < maxBlock; i++)
        for (int j = 0; j < maxBlock; j++)
        {
            int numBlock = i * maxBlock + j;
            if ((numBlock + 1) % Nmpi == rank)
                processBlock(i, j);
            if (rank == 0)
                progress_bar(i * maxBlock + j);
        }
    if (rank == 0)
    {
        progress_bar(maxBlock * maxBlock);
        qualifyCommonLines();
        solveForShifts();
        writeResults();
    }
}


/***************** Porting code from Yoel Shkolnisky **********************************
 *
 */

void
randomQuaternions(int k, DMatrix &quaternions)
{
    //  function q=qrand(K)
    //  %
    //  % q=qrand(K)
    //  %
    //  % Generate K random uniformly distributed quaternions.
    //  % Each quaternions is a four-elements column vector. Returns a matrix of
    //  % size 4xK.
    //  %
    //  % The 3-sphere S^3 in R^4 is a double cover of the rotation group SO(3),
    //  % SO(3) = RP^3.
    //  % We identify unit norm quaternions a^2+b^2+c^2+d^2=1 with group elements.
    //  % The antipodal points (-a,-b,-c,-d) and (a,b,c,d) are identified as the
    //  % same group elements, so we take a>=0.
    //  %
    quaternions.initGaussian(4, k);
    //quaternions.resize(4, k);
    //quaternions.initRandom(0., 1.);
    saveMatrix("random_numbers.txt", quaternions);

    DVector q(4);

    for (int j = 0; j < k; ++j)
    {
        //Get the j-column that is the quaternion
        quaternions.getCol(j, q);
        q.selfNormalize();
        if (XX(q) < 0)
            q *= -1;
        quaternions.setCol(j, q);
    }
    saveMatrix("random_quaternions.txt", quaternions);
}

void saveMatrix(const char *fn, DMatrix &matrix)
{
    std::fstream fs;
    fs.open(fn, std::fstream::trunc | std::fstream::out);
    fs.precision(16);

    for (size_t j = 0; j < MAT_YSIZE(matrix); ++j)
    {
        for (size_t i = 0; i < MAT_XSIZE(matrix); ++i)
        {
            fs << std::scientific << dMij(matrix, j, i) << "\t";
            //std::cerr << dAij(matrix, j, i) << " ";
        }
        fs << std::endl;
        //std::cerr << std::endl;
    }
    fs.close();
    //std::cerr << "DEBUG_JM: leaving saveMatrix" <<std::endl;
}

void quaternionToMatrix(const DVector &q, DMatrix &rotMatrix)
{
    //   function rot_matrix = q_to_rot(q)
    //   %
    //   % Convert a quaternion into a rotation matrix.
    //   %
#define q0 dMi(q, 0)
#define q1 dMi(q, 1)
#define q2 dMi(q, 2)
#define q3 dMi(q, 3)

    rotMatrix.resize(3, 3);

    rotMatrix(0, 0) = POW2(q0) + POW2(q1) - POW2(q2) - POW2(q3);
    rotMatrix(0, 1) = 2 * q1 * q2 - 2 * q0 * q3;
    rotMatrix(0, 2) = 2 * q0 * q2 + 2 * q1 * q3;

    rotMatrix(1, 0) = 2 * q1 * q2 + 2 * q0 * q3;
    rotMatrix(1, 1) = POW2(q0) - POW2(q1) + POW2(q2) - POW2(q3);
    rotMatrix(1, 2) = -2 * q0 * q1 + 2 * q2 * q3;

    rotMatrix(2, 0) = -2 * q0 * q2 + 2 * q1 * q3;
    rotMatrix(2, 1) = 2 * q0 * q1 + 2 * q2 * q3;
    rotMatrix(2, 2) = POW2(q0) - POW2(q1) - POW2(q2) + POW2(q3);

    //saveMatrix("rot_matrix.txt", rotMatrix);
}



#define AXIS(var) DVector var(3); var.initZeros()


// function [idx1,idx2,Z3]=commonline_q(q,k1,k2,n_theta)
// %
// % Given a Kx4 array of quaternions, find the common lines between
// % projections k1 and k2. Each projection has n_theta Fourier rays.
// %
// % idx1 is the index of the common line in projection k1. idx2 is the index
// % of the common line in projection k2. Z3 is the direction vector of the
// % common line.
// %
// % Yoel Shkolnisky, October 2008.

void quaternionCommonLines(const DMatrix &quaternions, CommonLineInfo &clInfo)
// clInfo serve as input param: k1, k2 and nRays
// and as output: idx1, idx2 and vector
{
    DVector qa1(4);
    DVector qa2(4); //quaternions k1 and k2
    DMatrix  R1, R2; //rotation matrixes

    for (int i = 0; i < 4; ++i)
    {
        qa1(i) = quaternions(i, clInfo.getImage(0));
        qa2(i) = quaternions(i, clInfo.getImage(1));
    }
    quaternionToMatrix(qa1, R1);
    quaternionToMatrix(qa2, R2);
    //Should be the inversve, but transpose is fine here
    //   R1.selfInverse();
    //   R2.selfInverse();
    R1 = R1.transpose();
    R2 = R2.transpose();

    // % Rotated coordinate system is the columns of the rotation matrix.
    // % We use explicit multiplication to show which column corresponds to x,
    // % which to y, and which to z. See commonline_euler for the difference.
    // X1=R1*([1 0 0].');
    // Y1=R1*([0 1 0].');
    // Z1=R1*([0 0 1].');
    //
    // X2=R2*([1 0 0].');
    // Y2=R2*([0 1 0].');
    // Z2=R2*([0 0 1].');
    //
    // Z3=[Z1(2)*Z2(3)-Z1(3)*Z2(2);...
    //     Z1(3)*Z2(1)-Z1(1)*Z2(3);...
    //     Z1(1)*Z2(2)-Z1(2)*Z2(1)];
    DVector X1, Y1, Z1, X2, Y2, Z2;
    DVector &Z3 = clInfo.vector;
    R1.getCol(0, X1);
    R1.getCol(1, Y1);
    R1.getCol(2, Z1);

    R2.getCol(0, X2);
    R2.getCol(1, Y2);
    R2.getCol(2, Z2);

    XX(Z3) = YY(Z1) * ZZ(Z2) - ZZ(Z1) * YY(Z2);
    YY(Z3) = ZZ(Z1) * XX(Z2) - XX(Z1) * ZZ(Z2);
    ZZ(Z3) = XX(Z1) * YY(Z2) - YY(Z1) * XX(Z2);

    double normZ3 = Z3.module();

    if (normZ3 < 1.0e-8)
        REPORT_ERROR(ERR_NUMERICAL, "GCAR:normTooSmall','Images have same orientation");

    Z3 /= normZ3;
    DVector Z3_t(Z3);
    Z3_t.selfTranspose();

    //   % Compute coordinates of the common-line in each local coordinate system.
    //   XY1=[X1 Y1];
    //   XY2=[X2 Y2];
    //   c1=(Z3.')*XY1;
    //   c2=(Z3.')*XY2;
    DMatrix XY1(3, 2), XY2(3, 2);
    XY1.setCol(0, X1);
    XY1.setCol(1, Y1);
    XY2.setCol(0, X2);
    XY2.setCol(1, Y2);
    DVector c1 = Z3_t * XY1;
    DVector c2 = Z3_t * XY2;
    // % Verify that the common-line is indeed common to both planes. The
    // % following warning should never happen! Just to make sure nothing went
    // % terribly wrong.
    // ev1=XY1*c1(:)-Z3;
    // ev2=XY2*c2(:)-Z3;
    c1.selfTranspose();
    c2.selfTranspose();
    DVector ev1 = XY1 * c1 - Z3;
    DVector ev2 = XY2 * c2 - Z3;

    double normEv1 = ev1.module() / normZ3;
    double normEv2 = ev2.module() / normZ3;

    if (normEv1 > 1.0e-12 || normEv2 > 1.0e-12)
        REPORT_ERROR(ERR_NUMERICAL,
                     formatString("GCAR:largeErrors: Common line is not common. Error1 = %f, Error2 = %f", normEv1, normEv2));
    //   % Compute angle of the common line at each projection's coordinate system
    //   theta1=atan2(c1(2),c1(1));
    //   theta2=atan2(c2(2),c2(1));
    //
    //   PI=4*atan(1.0);
    //   theta1=theta1+PI; % Shift from [-pi,pi] to [0,2*pi].
    //   theta2=theta2+PI;
    //std::cerr << "DEBUG_JM: c1: " << c1 <<std::endl;
    //std::cerr << "DEBUG_JM: c2: " << c2 << std::endl;
    double theta1 = atan2(YY(c1), XX(c1)) + PI;
    double theta2 = atan2(YY(c2), XX(c2)) + PI;
    //
    clInfo.setIndex(0, theta1);
    clInfo.setIndex(1, theta2);
}//function quaternionCommonLines


// function [clmatrix,clcorr]=clmatrix_cheat_q(q,n_theta)
// %
// % Build common lines matrix using the true quaternions corresponding to the
// % projections orientations. Each projection has n_theta rays.
// % clcorr is set to 1.0e-8.
// %
// % Yoel Shkolnisky, October 2008.
//
// N=size(q,2);
//
// clmatrix=zeros(N);   % common lines matrix
// clcorr=zeros(N);     % correlation coefficient for ach common line
//
// for k1=1:N-1
//     for k2=k1+1:N
//         [idx1,idx2]=commonline_q(q,k1,k2,n_theta);
//         clmatrix(k1,k2)=idx1+1;
//         clmatrix(k2,k1)=idx2+1;
//         clcorr(k1,k2)=1.0e-8;
//     end
// end
void commonlineMatrixCheat(const DMatrix &quaternions, size_t nRays,
                           DMatrix &clMatrix, DMatrix &clCorr)
{
    int n = quaternions.Xdim();

    clMatrix.initConstant(n, n, 0);
    clCorr.initZeros(n, n);

    CommonLineInfo clInfo(nRays);

    for (int i = 0; i < n - 1; ++i)
        for (int j = i + 1; j < n; ++j)
        {
            clInfo.setImages(i, j);
            quaternionCommonLines(quaternions, clInfo);
            dMij(clMatrix, i, j) = clInfo.getIndex(0);
            dMij(clMatrix, j, i) = clInfo.getIndex(1);
            dMij(clCorr, i, j) = 1.0e-8;
        }

}//function commonlineMatrixCheat



void anglesRotationMatrix(const DMatrix &clMatrix, size_t nRays, int clI, int clJ,
                          const DVector &Q1, const DVector &Q2,
                          DMatrix &R)
{

    double alpha1 = TWOPI * clI / nRays;
    double alpha2 = TWOPI * clJ / nRays;
    DVector N1(3);
    vectorProduct(Q1, Q2, N1);
    N1.selfNormalize();

    DMatrix U(3, 3);

    sincos(alpha1, &dMij(U, 1, 0), &dMij(U, 0, 0));
    sincos(alpha2, &dMij(U, 1, 1), &dMij(U, 0, 1));
    dMij(U, 2, 2) = 1.;

    if (U.det3x3() < 0)
        dMij(U, 2, 2) = -1.; //% Make sure we have a right-handed system)

    DMatrix Q(3, 3);
    Q.setCol(0, Q1);
    Q.setCol(1, Q2);
    Q.setCol(2, N1);

    U.selfInverse();
    DMatrix T = Q * U;
    T.svd(U, N1, Q);
    R = U * Q.transpose();
}//function

#define  EPS 1.0e-13 // % Should be 1.0e-13 after fixing XXX below.
#define  MAX_COND 1000 // % Largest allowed condition number for the system of equations
#define cl(i, j) dMij(clMatrix, k##i, k##j)

/** Negative output means error
 * -101 Triangle too small
 */
int tripletRotationMatrix(const DMatrix &clMatrix, size_t nRays,
                          int k1, int k2, int k3, DMatrix &R)
{
    DMatrix G(3, 3), Q;
    // % Find Q12, Q13, Q23
    double factor = TWOPI / nRays;
    double a = cos(factor * (cl(3,2)-cl(3,1)));
    double b = cos(factor * (cl(2,3)-cl(2,1)));
    double c = cos(factor * (cl(1,3)-cl(1,2)));

    if (1 + 2 * a * b * c - (POW2(a) + POW2(b) + POW2(c))<1.0e-5)
        return SMALL_TRIANGLE;

    G.initIdentity();
    // Compute G matrix according to (4.2)
    dMij(G, 0, 1) = a;
    dMij(G, 0, 2) = b;
    dMij(G, 1, 0) = a;
    dMij(G, 1, 2) = c;
    dMij(G, 2, 0) = b;
    dMij(G, 2, 1) = c;

    cholesky(G, Q);

    //Our cholesky gives us the lower triangular matrix
    //so we need to transpose to have the same as in Yoel code
    //that's why we take now rows instead of columns
    DVector Q12, Q13, Q23;
    Q.getRow(0, Q23);
    Q.getRow(1, Q13);
    Q.getRow(2, Q12);
    Q12.setCol();
    Q13.setCol();
    Q23.setCol();

    DMatrix R1, R2;
    anglesRotationMatrix(clMatrix, nRays, (int)cl(1, 2), (int)cl(1, 3), Q12, Q13, R1);
    anglesRotationMatrix(clMatrix, nRays, (int)cl(2, 1), (int)cl(2, 3), Q12, Q23, R2);
    // Compute rotation matrix according to (4.6)
    R = R1.transpose() * R2;

    return 0;
}//function tripleRotationMatrix

/** Helper function to set the 2x2 upper left corner of a rotation
 * matrix into a bigger matrix indexed by images
 */
void putRotationMatrix(const DMatrix &R, int k1, int k2, DMatrix &syncMatrix)
{
#define SET_POINT(x, y) dMij(syncMatrix, x+2*k1, y+2*k2) = dMij(R, x, y)
    SET_POINT(0, 0);
    SET_POINT(0, 1);
    SET_POINT(1, 0);
    SET_POINT(1, 1);
}

//%
//% Construct the CryoEM synchronization matrix, given a common lines matrix
//% clmatrix, that was constructed using angular resolution of L radial lines
//% per image.
//%
//% refq (optional) are the quaternions used to computed the common lines
//% matrix.
//%
//% Yoel Shkolnisky, August 2010.
void computeSyncMatrix(const DMatrix &clMatrix, size_t nRays, DMatrix &sMatrix, DMatrix * pQuaternions)
{
    int K = clMatrix.Xdim();
    DMatrix J, R, S(3,3);
    J.initIdentity(3);
    dMij(J, 2, 2) = -1;
    int kk;
    int result = 0;

    sMatrix.initIdentity(2*K);

    for (int k1 = 0; k1 < K - 1; ++k1)
    {
        for (int k2 = k1 + 1; k2 < K; ++k2)
        {
            kk = 0; //% Each triplet (k1,k2,k3) defines some rotation matrix from
            //                % image k1 to image k2. kk counts the number of such
            //                % rotations, that is, the number of thrid images k3 that give
            //                % rise to a rotation from image k1 to image k2. Not all
            //                % images k3 give rise to such a rotation, since, for example
            //                % the resuling triangle can be too small.

            //For now average all rotation matrixes, maybe later on
            // we will need to store all of them for a better estimation
            // than just average
            S.initZeros();

            for (int k3 = 0; k3 < K; ++k3)
                if (k3 != k1 && k3 != k2)
                {
                    result = tripletRotationMatrix(clMatrix, nRays, k1, k2, k3, R);
                    if (result >= 0)
                    {
                      S += R;
                      ++kk;
                    }
                    std::cerr << "DEBUG_JM: triplet: " << formatString("%d %d, %d", k1, k2, k3) << std::endl;
                    std::cerr << "DEBUG_JM: R: " << R << std::endl;
                }
            if (kk > 0)
              S /= kk;
            putRotationMatrix(S, k1, k2, sMatrix);
            putRotationMatrix(S.transpose(), k2, k1, sMatrix);
        }
    }
}//function computeSyncMatrix

void rotationsFromSyncMatrix(const DMatrix &sMatrix, DMatrix * pQuaternions)
{
    int K = sMatrix.Xdim() / 2;
    int Kx3 = 3 * K;

    DMatrix U, V, V1(3, K), V2(3, K);
    DVector W;
    IVector indexes;
    sMatrix.svd(U, W, V);
    indexes.resizeNoCopy(W);
    indexes.enumerate();

    //std::cerr << "DEBUG_JM: W: " << W << std::endl;

#define SORTED_INDEX(i) dMi(indexes, i)
#define SORTED_ELEM(M, i) dMi(M, SORTED_INDEX(i))

    //Order by insertion sort
    int iAux;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(W)
        for (int j = i; j > 0 && SORTED_ELEM(W, j) > SORTED_ELEM(W, j-1); --j)
          VEC_SWAP(indexes, j, j-1, iAux);

//    std::cerr << "DEBUG_JM: VEC_XSIZE(indexes): " << VEC_XSIZE(indexes) << std::endl;
//    std::cerr << "DEBUG_JM: 10 bigger eigenvalues:" <<std::endl;
//    std::cerr << "DEBUG_JM: indexes: " << indexes << std::endl;
//    for (int i = 0; i < VEC_XSIZE(indexes); ++i)
//      std::cerr << SORTED_ELEM(W, i) << std::endl;

    for (int i = 0; i < K; ++i)
    {
      iAux = 2 * i;
      for (int j = 0; j < 3; ++j)
      {
        dMij(V1, j, i) = dMij(V, iAux, SORTED_INDEX(j));
        dMij(V2, j, i) = dMij(V, iAux+1, SORTED_INDEX(j));
      }
    }

    std::cerr << "DEBUG_JM: V1: " << V1 << std::endl;
    std::cerr << "DEBUG_JM: V2: " << V2 << std::endl;

    //% 3*K equations in 9 variables (3 x 3 matrix entries).
    DMatrix equations(Kx3, 9);
    int index;

    for (int k = 0; k < K; ++k)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          index = 3 * i + j;
          iAux = 3 * k;
          dMij(equations, iAux, index) = dMij(V1, i, k) * dMij(V1, j, k);
          dMij(equations, iAux+1, index) = dMij(V2, i, k) * dMij(V2, j, k);
          dMij(equations, iAux+2, index) = dMij(V1, i, k) * dMij(V2, j, k);
        }
//    std::cerr << "DEBUG_JM: equations: " << equations << std::endl;
    PseudoInverseHelper pih;

    DMatrix &trunc_equations = pih.A;
    trunc_equations.resizeNoCopy(Kx3, 6);
    //Use vector W to select columns from equations
    int cols[6] = {0, 1, 2, 4, 5, 8};
    for (int i = 0; i < 6; ++i)
    {
      equations.getCol(cols[i], W);
      trunc_equations.setCol(i, W);
    }

    DVector &b = pih.b;
    b.initConstant(Kx3, 1.);
    for (int i = 2; i < Kx3; i+=3)
      dMi(b, i) = 0.;

    //% Find the least squares approximation
    DVector ATA_vec;
    DMatrix ATA(3, 3);

    //C
    solveLinearSystem(pih, ATA_vec);

//    std::cerr << "DEBUG_JM: trunc_equations: " << trunc_equations << std::endl;
//    saveMatrix("truncated.txt", trunc_equations);
//    std::cerr << "DEBUG_JM: b: " << b << std::endl;
//    std::cerr << "DEBUG_JM: ATA_vec: " << ATA_vec << std::endl;

    dMij(ATA, 0, 0) = dMi(ATA_vec, 0);
    dMij(ATA, 0, 1) = dMi(ATA_vec, 1);
    dMij(ATA, 0, 2) = dMi(ATA_vec, 2);
    dMij(ATA, 1, 0) = dMi(ATA_vec, 1);
    dMij(ATA, 1, 1) = dMi(ATA_vec, 3);
    dMij(ATA, 1, 2) = dMi(ATA_vec, 4);
    dMij(ATA, 2, 0) = dMi(ATA_vec, 2);
    dMij(ATA, 2, 1) = dMi(ATA_vec, 4);
    dMij(ATA, 2, 2) = dMi(ATA_vec, 5);

    std::cerr << "DEBUG_JM: ATA: " << ATA << std::endl;

    DMatrix A, R1, R2;
    cholesky(ATA, A);
    A = A.transpose();

    std::cerr << "DEBUG_JM: A: " << A << std::endl;

    R1 = A * V1;
    R2 = A * V2;
    //DMatrix R3(R1);
    DVector v1, v2, v3(3);
    std::vector<DMatrix> rotations(K);

    MetaData MD("images.xmd");
    MDIterator it(MD);
    for (int i = 0; i < K; ++i)
    {
      DMatrix &R = rotations[i];
      R.resizeNoCopy(3, 3);
      R1.getCol(i, v1);
      R2.getCol(i, v2);
      vectorProduct(v1, v2, v3);
      //R3.setCol(i, v3);
      R.setCol(0, v1);
      R.setCol(1, v2);
      R.setCol(2, v3);
      //% Enforce R to be a rotation (in case the error is large)
      R.svd(U, W, V);
      R = U * V.transpose();

      double rot, tilt, psi;
      //
      Euler_matrix2angles(R.transpose(), rot, tilt, psi);
      MD.setValue(MDL_ANGLE_ROT,rot,it.objId);
      MD.setValue(MDL_ANGLE_TILT,tilt,it.objId);
      MD.setValue(MDL_ANGLE_PSI,psi,it.objId);
      it.moveNext();

      std::cerr << "DEBUG_JM: R" << i << " : " << R << std::endl;
    }
    MD.write("imagesReconstruction.xmd");
}//function rotationsFromSyncMatrix

