/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2009)
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

#include "class_averages.h"
#include <data/docfile.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/polar.h>

/* VQProjection basics ---------------------------------------------------- */
void VQProjection::updateProjection(const Matrix2D<double> &I,
    double corrCode, int idx)
{
    Pupdate+=I;
    nextListImg.push_back(idx);
    if (!classicalMultiref)
        nextClassCorr.push_back(corrCode);
}

int VQProjection::sendMPI(int d_rank)
{
    int size;
    
    // Projection
    // Matrix2D<double> P;
    MPI_Send( &(P.xdim), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(P.ydim), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( P.data, P.yxdim, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
    
    // Update for next iteration
    // Matrix2D<double> Pupdate;
    MPI_Send( &(Pupdate.xdim), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(Pupdate.ydim), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( Pupdate.data, Pupdate.yxdim, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
     	   
    MPI_Send( &(rotationalCorr.xdim), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &rotationalCorr.data,rotationalCorr.xdim , MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );

    size = currentListImg.size(); 
    MPI_Send( &size, 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &currentListImg[0], size, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    
    size = nextListImg.size(); 
    MPI_Send( &size, 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &nextListImg[0], size, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
     
    // Correlations of the next class members
    // std::vector<double> nextClassCorr;
    size = nextClassCorr.size();
    MPI_Send( &size, 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &nextClassCorr[0], size, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
    
    // Correlations of the next non-class members
    // std::vector<double> nextNonClassCorr;
    size = nextNonClassCorr.size();
    MPI_Send( &size, 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD ); 
    MPI_Send( &nextNonClassCorr[0], size, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
       
    // Histogram of the correlations of the current class members
    MPI_Send( &(histClass.xdim), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(histClass.hmin), 1, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(histClass.hmax), 1, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(histClass.step_size), 1, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(histClass.no_samples), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( histClass.data, histClass.xdim, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );

    // Histogram of the correlations of the current non-class members
    MPI_Send( &(histNonClass.xdim), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(histNonClass.hmin), 1, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(histNonClass.hmax), 1, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(histNonClass.step_size), 1, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &(histNonClass.no_samples), 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( histNonClass.data, histNonClass.xdim, MPI_DOUBLE, d_rank, 0, MPI_COMM_WORLD );

    // List of neighbour indexes
    // std::vector<int> neighboursIdx;
    size = neighboursIdx.size();
    MPI_Send( &size, 1, MPI_INT, d_rank, 0, MPI_COMM_WORLD );
    MPI_Send( &neighboursIdx[0], neighboursIdx.size(), MPI_INT, d_rank, 0, MPI_COMM_WORLD );
}
    
int VQProjection::receiveMPI(int s_rank)
{    
    // Projection
    // Matrix2D<double> P;
    int sizex,sizey,size;
    MPI_Recv( &sizex, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL ); 
    MPI_Recv( &sizey, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );  
    P.resize( sizey, sizex );
    P.setXmippOrigin();
    MPI_Recv( P.data, sizey*sizex, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );
    
    // Update for next iteration
    // Matrix2D<double> Pupdate;
    MPI_Recv( &sizex, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );  
    MPI_Recv( &sizey, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );  
    Pupdate.resize( sizey, sizex );
    Pupdate.setXmippOrigin();
    MPI_Recv( Pupdate.data, sizey*sizex, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );

    // Rotational correlation for best_rotation
    // Matrix1D<double> rotationalCorr;
    MPI_Recv( &size, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    rotationalCorr.resize( size );
    rotationalCorr.setXmippOrigin();
    MPI_Recv( rotationalCorr.data, size, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );

    // local_transformer depends on rotationalCorr
    normalizedPolarFourierTransform(P,polarFourierP,false,XSIZE(P)/5,XSIZE(P)/2,plans,1);
    local_transformer.setReal(rotationalCorr);

    // List of images assigned
    // std::vector<int> currentListImg;
    MPI_Recv( &size, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    currentListImg.resize( size );
    MPI_Recv( &currentListImg[0], size, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );

    // List of images assigned
    // std::vector<int> nextListImg;
    MPI_Recv( &size, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    nextListImg.resize( size );
    MPI_Recv( &nextListImg[0], size, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    
    // Correlations of the next class members
    // std::vector<double> nextClassCorr;
    MPI_Recv( &size, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    nextClassCorr.resize( size );
    MPI_Recv( &nextClassCorr[0], size, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );
    
    // Correlations of the next non-class members
    // std::vector<double> nextNonClassCorr;
    MPI_Recv( &size, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL ); 
    nextNonClassCorr.resize( size );
    MPI_Recv( &nextNonClassCorr[0], size, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );
     
    double d_value;
    int i_value;
    
    // Histogram of the correlations of the current class members
    MPI_Recv( &size, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    MPI_Recv( &d_value, 1, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL);
    histClass.hmin = d_value;
    MPI_Recv( &d_value, 1, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );
    histClass.hmax = d_value;
    MPI_Recv( &d_value, 1, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );
    histClass.step_size = d_value;
    MPI_Recv( &i_value, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    histClass.no_samples = i_value;
    histClass.resize( size );
    MPI_Recv( histClass.data, size, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );

    // Histogram of the correlations of the current non-class members
    MPI_Recv( &size, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    MPI_Recv( &d_value, 1, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );
    histNonClass.hmin = d_value;
    MPI_Recv( &d_value, 1, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );
    histNonClass.hmax = d_value;
    MPI_Recv( &d_value, 1, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );
    histNonClass.step_size = d_value;
    MPI_Recv( &i_value, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    histNonClass.no_samples = i_value;
    histNonClass.resize( size );
    MPI_Recv( histNonClass.data, size, MPI_DOUBLE, s_rank, 0, MPI_COMM_WORLD, NULL );

    // List of neighbour indexes
    // std::vector<int> neighboursIdx;
    MPI_Recv( &size, 1, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
    neighboursIdx.resize( size );
    MPI_Recv( &neighboursIdx[0], size, MPI_INT, s_rank, 0, MPI_COMM_WORLD, NULL );
}

//#define DEBUG
void VQProjection::transferUpdate()
{
    if (nextListImg.size()>0)
    {
        double iNq=1.0/nextListImg.size();
        Pupdate*=iNq;
        Pupdate.statisticsAdjust(0,1);
        P=Pupdate;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(P)
            if (!(*mask)(i,j)) P(i,j)=0;
        Pupdate.initZeros(P);
        computeTransforms();
        std::sort(nextListImg.begin(),nextListImg.end() );
	currentListImg=nextListImg;
        nextListImg.clear();
    
        if (!classicalMultiref)
        {
            Matrix1D<double> classCorr, nonClassCorr;
            classCorr.initZeros(nextClassCorr.size());
            nonClassCorr.initZeros(nextNonClassCorr.size());
            FOR_ALL_ELEMENTS_IN_MATRIX1D(classCorr)
                classCorr(i)=nextClassCorr[i];
            FOR_ALL_ELEMENTS_IN_MATRIX1D(nonClassCorr)
                nonClassCorr(i)=nextNonClassCorr[i];
            nextClassCorr.clear();
            nextNonClassCorr.clear();
            double minC, maxC; classCorr.computeDoubleMinMax(minC,maxC);
            double minN, maxN; nonClassCorr.computeDoubleMinMax(minN,maxN);
            double c0=XMIPP_MIN(minC,minN);
            double cF=XMIPP_MAX(maxC,maxN);
            compute_hist(classCorr,histClass,c0,cF,200);
            compute_hist(nonClassCorr,histNonClass,c0,cF,200);
            histClass+=1; // Laplace correction
            histNonClass+=1;
            histClass/=histClass.sum();
            histNonClass/=histNonClass.sum();
        }

        #ifdef DEBUG
            histClass.write("PPPclass.txt");
            histNonClass.write("PPPnonClass.txt");
            std::cout << "Histograms have been written. Press any key\n";
            char c; std::cin >> c;
        #endif
    }
    else
    {
        currentListImg.clear();
        P.initZeros();
    }
}
#undef DEBUG

void VQProjection::computeTransforms()
{
    // Make sure the image is centered
    centerImage(P);

    // Compute the polar Fourier transform of the full image
    normalizedPolarFourierTransform(P,polarFourierP,false,
        XSIZE(P)/5,XSIZE(P)/2,plans,1);
    rotationalCorr.resize(2*polarFourierP.getSampleNoOuterRing()-1);
    local_transformer.setReal(rotationalCorr);
}

/* Show -------------------------------------------------------------------- */
void VQProjection::show() const
{
    ImageXmipp save;
    save()=P;
    save.write("PPPcode.xmp");
    for (int i=0; i<10; i++)
        std::cout << nextListImg[i] << " ";
    std::cout << std::endl;
    std::cout << "Images have been saved. Press any key\n";
    char c; std::cin >> c;
}

double fastCorrentropy(const Matrix2D<double> &x, const Matrix2D<double> &y,
    double sigma, const GaussianInterpolator &G, const Matrix2D<int> &mask)
{
    double retvalxy=0, retvalxx=0, retvalyy=0;
    double isigma=1.0/sigma;
    double minx, maxx, avgx, stddevx;
    double miny, maxy, avgy, stddevy;
    computeStats_within_binary_mask(mask, x, minx, maxx, avgx, stddevx);
    computeStats_within_binary_mask(mask, y, miny, maxy, avgy, stddevy);
    int maskSum=0;
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(x)
    {
        if (DIRECT_MAT_ELEM(mask,i,j))
        {
            retvalxy+=G.getValue(isigma*(DIRECT_MAT_ELEM(x,i,j)-DIRECT_MAT_ELEM(y,i,j)));
            maskSum++;
        }
    }
    return (retvalxy/maskSum);
}

//#define DEBUG
void VQProjection::fitBasic(Matrix2D<double> &I,
    double sigma, double &corrCode)
{
    Matrix2D<double> ARS, ASR, R;
    ARS.initIdentity(3);
    ASR.initIdentity(3);
    Matrix2D<double> IauxSR=I, IauxRS=I;

    // Align the image with the node
    for (int i=0; i<2; i++)
    {
        double shiftX, shiftY;
	
        // Shift then rotate	
        best_shift(P,IauxSR,shiftX,shiftY);
        ASR(0,2)+=shiftX;
        ASR(1,2)+=shiftY;
        applyGeometry(IauxSR,ASR,I,IS_NOT_INV,WRAP);
        
        Polar< std::complex<double> > polarFourierI;
	normalizedPolarFourierTransform(
		IauxSR,
		polarFourierI,
		true,
                XSIZE(P)/5,
		XSIZE(P)/2,
                plans,
                1);
        
        double bestRot = best_rotation(polarFourierP,polarFourierI,
            local_transformer);
	R=rotation2DMatrix(-bestRot);
        ASR=R*ASR;
        applyGeometry(IauxSR,ASR,I,IS_NOT_INV,WRAP);

        // Rotate then shift
	normalizedPolarFourierTransform(
		IauxRS,
		polarFourierI,
		true,
                XSIZE(P)/5,
		XSIZE(P)/2,
                plans,
                1);
        bestRot = best_rotation(polarFourierP,polarFourierI,
            local_transformer);
	R=rotation2DMatrix(-bestRot);
        ARS=R*ARS;
        applyGeometry(IauxRS,ARS,I,IS_NOT_INV,WRAP);

        best_shift(P,IauxRS,shiftX,shiftY);
        ARS(0,2)+=shiftX;
        ARS(1,2)+=shiftY;
        applyGeometry(IauxRS,ARS,I,IS_NOT_INV,WRAP);
    }
    
    //applyGeometryBSpline(IauxRS,ARS,I,3,IS_NOT_INV,WRAP);
    //applyGeometryBSpline(IauxSR,ASR,I,3,IS_NOT_INV,WRAP);

    applyGeometry(IauxRS,ARS,I,IS_NOT_INV,WRAP);
    applyGeometry(IauxSR,ASR,I,IS_NOT_INV,WRAP);

    // Compute the correntropy
    double corrCodeRS, corrCodeSR;
    if (useCorrelation)
    {
        corrCodeRS=correlation(P,IauxRS,mask);
        corrCodeSR=correlation(P,IauxSR,mask);
    }
    else
    {
        double sigmaPI=sigma;
        if (useFixedCorrentropy) sigmaPI*=sqrt(1+1.0/currentListImg.size());
        corrCodeRS=fastCorrentropy(P,IauxRS,sigmaPI,
            *gaussianInterpolator,*mask);
        corrCodeSR=fastCorrentropy(P,IauxSR,sigmaPI,
            *gaussianInterpolator,*mask);
    }
    if (corrCodeRS>corrCodeSR)
    {
        I=IauxRS;
        corrCode=corrCodeRS;
    }
    else
    {
        I=IauxSR;
        corrCode=corrCodeSR;
    }

    #ifdef DEBUG
        ImageXmipp save;
        save()=P; save.write("PPPI1.xmp");
        save()=I;   save.write("PPPI2.xmp");
        save()=P-I; save.write("PPPdiff.xmp");
        std::cout << "sigma=" << sigma << " corr=" << corrCode
                  << ". Press" << std::endl;
        char c; std::cin >> c;
    #endif
}
#undef DEBUG

void VQProjection::fit(Matrix2D<double> &I,
    double sigma, bool noMirror, double &corrCode, double &likelihood)
{
    // Try this image
    Matrix2D<double> Iaux=I;
    double corrCodeAux;
    fitBasic(Iaux,sigma,corrCodeAux);

    Matrix2D<double> bestImg=Iaux;
    double bestCorrCode=corrCodeAux;

    // Try its mirror
    if (!noMirror)
    {
        Iaux=I;
        Iaux.selfReverseX();  
	Iaux.setXmippOrigin();

        fitBasic(Iaux,sigma,corrCodeAux);

        if (corrCodeAux>bestCorrCode)
        {
            bestImg=Iaux;
            bestCorrCode=corrCodeAux;
        }
    }
    I=bestImg;
    corrCode=bestCorrCode;

    likelihood=0;
    if (!classicalMultiref)
    {
        // Find the likelihood
        int idx;
        histClass.val2index(corrCode,idx);
        if (idx<0) return;
        double likelihoodClass=0;
        for (int i=0; i<=idx; i++) likelihoodClass+=histClass(i);
        histNonClass.val2index(corrCode,idx);
        if (idx<0) return;
        double likelihoodNonClass=0;
        for (int i=0; i<=idx; i++) likelihoodNonClass+=histNonClass(i);
        likelihood=likelihoodClass*likelihoodNonClass;
    }
}

/* Look for K neighbours in a list ----------------------------------------- */
//#define DEBUG
void VQProjection::lookForNeighbours(const std::vector<VQProjection *> listP,
    double sigma, bool noMirror, int K)
{
    int Q=listP.size();
    neighboursIdx.clear();
    if (K==Q)
    {
        for (int q=0; q<Q; q++)
            neighboursIdx.push_back(q);
    }
    else
    {
        Matrix1D<double> distanceCode;
        distanceCode.initZeros(Q);
        double likelihood;
        for (int q=0; q<Q; q++)
        {
            if (listP[q]==this)
                distanceCode(q)=1;
            else
            {
                Matrix2D<double> I=listP[q]->P;
		fit(I,sigma,noMirror,distanceCode(q),likelihood);
            }
        }

        Matrix1D<int> idx=distanceCode.indexSort();
        for (int k=0; k<K; k++)
            neighboursIdx.push_back(idx(Q-k-1)-1);
    }
    #ifdef DEBUG
        ImageXmipp save;
        save()=P; save.write("PPPthis.xmp");
        for (int k=0; k<K; k++)
        {
            save()=listP[neighboursIdx[k]]->P;
			save.write((std::string)"PPPneigh"+integerToString(k,1));
        }
        std::cout << "Neighbours saved. Press any key\n";
        char c; std::cin >> c;
    #endif
}
#undef DEBUG

/* VQ initialization ------------------------------------------------ */
//#define DEBUG
void VQ::initialize(SelFile &_SF, int _Niter, int _Nneighbours,
    double _PminSize, std::vector< Matrix2D<double> > _codes0, int _Ncodes0,
    bool _noMirror, bool _verbose, bool _corrSplit, bool _useCorrelation,
    bool _useFixedCorrentropy, bool _classicalMultiref, bool _alignImages,
    bool _fast, int rank)
{
    // Only "parent" worker prints
    if( rank == 0 )
    	std::cout << "Initializing ...\n";
    
    gaussianInterpolator.initialize(6,60000,false);
    SF=&_SF;
    Niter=_Niter;
    Nneighbours=_Nneighbours;
    PminSize=_PminSize;
    noMirror=_noMirror;
    verbose=_verbose;
    corrSplit=_corrSplit;
    useCorrelation=_useCorrelation;
    useFixedCorrentropy=_useFixedCorrentropy;
    classicalMultiref=_classicalMultiref;
    alignImages=_alignImages;
    fast=_fast;
    int Ydim, Xdim;
    int mpi_size;
    
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
    SF->ImgSize(Ydim,Xdim);
    
    // Prepare mask for evaluating the noise outside
    mask.resize(Ydim,Xdim);
    mask.setXmippOrigin();
    BinaryCircularMask(mask,Xdim/2, INNER_MASK);
	
    // Start with _Ncodes0 codevectors
    for (int q=0; q<_Ncodes0; q++)
    {
    	P.push_back(new VQProjection());
        P[q]->gaussianInterpolator=&gaussianInterpolator;
        P[q]->plans=NULL;
        P[q]->mask=&mask;
        P[q]->useCorrelation=useCorrelation;
        P[q]->useFixedCorrentropy=useFixedCorrentropy;
        P[q]->classicalMultiref=classicalMultiref;
	if (_codes0.size()!=0)
        {
            P[q]->Pupdate=_codes0[q];
            P[q]->nextListImg.push_back(-1);
        } else            
            P[q]->Pupdate.initZeros(Ydim,Xdim);
        P[q]->Pupdate.setXmippOrigin();
        if (_codes0.size()!=0)
            P[q]->transferUpdate();
    }	
	
    // Estimate sigma and if no previous classes have been given,
    // assign randomly
    sigma = 0;
    int N=SF->ImgNo();
    std::vector<int> initNodeAssign( N, -1 );
    
    if( rank == 0 )
    	init_progress_bar(N);
    
    for( int n=0; n < N ; n++ )
    {
    	ImageXmipp I;
	FileName auxFN = SF->NextImg();
	SFv.push_back( auxFN );
			
    	if( n % mpi_size == rank )
	{
	    I.read( auxFN );
	    I().setXmippOrigin();
	    I().statisticsAdjust(0,1);

	    // Measure the variance of the signal outside a circular mask
            double min_val, max_val, avg, stddev;
            computeStats_within_binary_mask(mask,I(),
                min_val, max_val, avg, stddev);
            sigma+=stddev;

            // Put it randomly in one of the classes
	    if (_codes0.size()==0)
            {
            	//int q=ROUND(rnd_unif(0,_Ncodes0-1));
              	int q = n%_Ncodes0;

	        initNodeAssign[n] = q;
            	P[q]->updateProjection(I(),0,n);
            }

            if (n%100==0 && rank==0) 
		progress_bar(n);
	}
    }
    if( rank == 0 )
    	progress_bar(N);

    // Put sigma in common
    double sigmaAux=0;
    MPI_Allreduce( &sigma, &sigmaAux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
    sigma=sigmaAux;	
    sigma/=N;
    sigma*=sqrt(2.0);

    // Really make the assignment if some initial images were given
    if (_codes0.size()!=0)
    {
        // Really make the assignment once sigma is known
        if( rank == 0 )
        {
    	    std::cout << "Making assignment ...\n";
            init_progress_bar(N);
        }
        for( int n=0; n < N ; n++ )
        {
    	    if( n % mpi_size == rank )
	    {
    	        ImageXmipp I;
		I.read( SFv[n] );
		I().setXmippOrigin();
		I().statisticsAdjust(0,1);

                double bestCorr=-2;
                int q;
                for (int qp=0; qp<_Ncodes0; qp++)
                {
   	            Matrix2D<double> Iaux;
	            double corrCode, likelihood;
                    Iaux=I();
                    P[qp]->fit(Iaux, sigma, noMirror, corrCode, likelihood);
                    if (corrCode>bestCorr)
                    {
                        bestCorr=corrCode;
                        q=qp;
                    }
                }

	        initNodeAssign[n] = q;
            	P[q]->updateProjection(I(),0,n);
	    }
            if (n%100==0 && rank==0) 
		progress_bar(n);
        }
        if( rank == 0 )
    	    progress_bar(N);
    }

    // Put assignment in common
    std::vector<int> auxNodeAssign( N, -1);
    MPI_Allreduce( &(initNodeAssign[0]), &(auxNodeAssign[0]), N, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    initNodeAssign = auxNodeAssign;

    // Update for next iteration
    // Matrix2D<double> Pupdate;
    Matrix2D<double> PupdateAux;
    PupdateAux.resize( P[0]->Pupdate.ydim, P[0]->Pupdate.xdim );
    PupdateAux.setXmippOrigin( );
        
    for (int q=0; q<_Ncodes0; q++)
    {
    	PupdateAux.initZeros();
        MPI_Allreduce( P[q]->Pupdate.data, PupdateAux.data, P[q]->Pupdate.yxdim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        P[q]->Pupdate = PupdateAux;

    	// Share nextClassCorr and nextNonClassCorr
    	int oldSizeClassCorr =  P[q]->nextClassCorr.size();
    	int oldSizeNextListImg = P[q]->nextListImg.size();

    	std::vector<double> oldListClassCorr = P[q]->nextClassCorr;
    	std::vector<int> oldListNextListImg = P[q]->nextListImg;
	    
    	for( int ranks = 0 ; ranks < mpi_size ; ranks ++ )
    	{
    		if( ranks == rank )
		{
		    MPI_Bcast( &oldSizeClassCorr, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListClassCorr[0]), oldSizeClassCorr, MPI_DOUBLE, rank, MPI_COMM_WORLD );    
	
		    MPI_Bcast( &oldSizeNextListImg, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListNextListImg[0]),oldSizeNextListImg, MPI_INT, rank, MPI_COMM_WORLD );
		}
		else
		{
	    	    int size;
		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
		    std::vector<double> aux(size, 0);
		    MPI_Bcast( &(aux[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );
				
		    for( int element = 0; element < size; element ++ )
		    {
	                P[q]->nextClassCorr.push_back(aux[element]);
		    }
	    
    	    	    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
	    	    std::vector<int> aux2(size, 0);
	    	    MPI_Bcast( &(aux2[0]), size, MPI_INT, ranks, MPI_COMM_WORLD );
			
	    	    for( int element = 0; element < size; element ++ )
	    	    {
                        P[q]->nextListImg.push_back(aux2[element]);
	    	    }
	        }
    	 }
    }
			
    for (int q=0; q<_Ncodes0; q++)
    	P[q]->transferUpdate();
		
    // Now compute the histogram of corr values
    if (!classicalMultiref) {
        if( rank == 0 )
    	    init_progress_bar(N);

        for (int n=0; n<N; n++)
        {
    	    ImageXmipp I;
    	    if( n % mpi_size == rank )
	    {
	        I.read( SFv[n] );
	        I().setXmippOrigin();
	        I().statisticsAdjust(0,1);

   	        Matrix2D<double> Iaux;

	        int q=initNodeAssign[n];

	        double corrCode, likelihood;
                Iaux=I();
                P[q]->fit(Iaux, sigma, noMirror, corrCode, likelihood);
                P[q]->updateProjection(Iaux,corrCode,n);
                updateNonCode(I(),q);

	        for (int qp=0; qp<_Ncodes0; qp++)
	        {
        	    if (q==qp) continue;
        	    Matrix2D<double> Iaux=I();
		    double corrNonCode, nonCodeLikelihood;
		    P[qp]->fit(Iaux,sigma,noMirror,corrNonCode,nonCodeLikelihood);
		    P[qp]->updateNonProjection(corrNonCode);
	        }
                if( rank == 0 )
		    if (n%100==0) progress_bar(n);
    	    }
        }

        if( rank == 0 )
    	    progress_bar(N);
    
        for (int q=0; q<_Ncodes0; q++)
        {
    	    PupdateAux.initZeros();
            MPI_Allreduce( P[q]->Pupdate.data, PupdateAux.data, P[q]->Pupdate.yxdim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            P[q]->Pupdate = PupdateAux;

    	    // Share nextClassCorr and nextNonClassCorr
    	    int oldSizeClassCorr =  P[q]->nextClassCorr.size();
    	    int oldSizeNextListImg = P[q]->nextListImg.size();
	    int oldSizeNextNonClassCorr = P[q]->nextNonClassCorr.size();

    	    std::vector<double> oldListClassCorr = P[q]->nextClassCorr;
    	    std::vector<int> oldListNextListImg = P[q]->nextListImg;
	    std::vector<double> oldListNextNonClassCorr = P[q]->nextNonClassCorr;    

    	    for( int ranks = 0 ; ranks < mpi_size ; ranks ++ )
    	    {
    	        if( ranks == rank )
	        {
		    MPI_Bcast( &oldSizeClassCorr, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListClassCorr[0]), oldSizeClassCorr, MPI_DOUBLE, rank, MPI_COMM_WORLD );    

		    MPI_Bcast( &oldSizeNextListImg, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListNextListImg[0]),oldSizeNextListImg, MPI_INT, rank, MPI_COMM_WORLD );

		    MPI_Bcast( &oldSizeNextNonClassCorr, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListNextNonClassCorr[0]), oldSizeNextNonClassCorr, MPI_DOUBLE, rank, MPI_COMM_WORLD );   
	        }
	        else
	        {
	    	    int size;
		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
		    std::vector<double> aux(size, 0);
		    MPI_Bcast( &(aux[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );

		    for( int element = 0; element < size; element ++ )
		    {
	                P[q]->nextClassCorr.push_back(aux[element]);
		    }

    	    	    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
	    	    std::vector<int> aux2(size, 0);
	    	    MPI_Bcast( &(aux2[0]), size, MPI_INT, ranks, MPI_COMM_WORLD );

	    	    for( int element = 0; element < size; element ++ )
	    	    {
                        P[q]->nextListImg.push_back(aux2[element]);
	    	    }

		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
	    	    std::vector<double> aux3(size, 0);
	    	    MPI_Bcast( &(aux3[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );

	    	    for( int element = 0; element < size; element ++ )
	    	    {
                        P[q]->nextNonClassCorr.push_back(aux3[element]);
	    	    }
	        }
    	     }
        }

        for (int q=0; q<_Ncodes0; q++)
    	    P[q]->transferUpdate();
    }
}
#undef DEBUG

/* VQ write --------------------------------------------------------- */
//#define DEBUG
void VQ::write(const FileName &fnRoot, bool final) const
{
    int Q=P.size();
    int Nimg=SFv.size();
    SelFile SFout;
    DocFile DFout;
    Matrix1D<double> aux, Nq, dq;
    aux.resize(2);
    Nq.resize(Q);
    dq.resize(Q);
    DFout.append_comment("Number of imgs assigned, Fraction of imgs assigned");
    for (int q=0; q<Q; q++)
    {
        FileName fnOut;
        fnOut.compose(fnRoot,q,"xmp");
        ImageXmipp I;
        #ifdef DEBUG
            std::cout << "About to write node " << q << " address=" << &P[q] << std::endl;
        #endif
        I()=P[q]->P;
        I.write(fnOut);
        SFout.insert(fnOut);
        DFout.append_comment(fnOut);
        Nq(q)=aux(0)=P[q]->currentListImg.size();
        aux(1)=aux(0)/Nimg;
        DFout.append_data_line(aux);
    }
    SFout.write(fnRoot+".sel");
    DFout.write(fnRoot+".doc");
    
    // Make the selfiles of each class
    for (int q=0; q<Q; q++)
    {
        SelFile SFq;
        SFq.clear();
        int imax=P[q]->currentListImg.size();
        for (int i=0; i<imax; i++)
        {
            if (alignImages && final)
            {
                ImageXmipp I;
                I.read(SFv[P[q]->currentListImg[i]]);
                I().setXmippOrigin();
                
                double corrCode, likelihood;
                P[q]->fit(I(), sigma, noMirror, corrCode, likelihood);
                I.write();
            }
            SFq.insert(SFv[P[q]->currentListImg[i]]);
        }
        SFq.write(fnRoot+integerToString(q,6)+".sel");
    }
}
#undef DEBUG

void VQ::lookNode(Matrix2D<double> &I, int idx, int oldnode,
    int &newnode, double &corrCode, double &likelihood)
{
    Matrix2D<double> Ibackup;
    Ibackup=I;
    
    int Q=P.size();
    double bestCorrCode=-1, bestLikelihood=-1;
    int bestq=-1;
    Matrix2D<double> bestImg;
    Matrix1D<double> corrCodeList;
    corrCodeList.initZeros(Q);
    int Nimg=SFv.size();
    for (int q=0; q<Q; q++)
    {
        // Check if q is neighbour of the oldnode
        bool proceed=false;
        bool neighbour=false;
        if (oldnode!=-1)
        {
            int imax=P[oldnode]->neighboursIdx.size();
            for (int i=0; i<imax; i++)
            {
                if (P[oldnode]->neighboursIdx[i]==q)
                {
                    proceed=true;
                    neighbour=true;
                    break;
                }
            }
            if (!proceed && fast)
            {
                double threshold=3.0*P[oldnode]->currentListImg.size();
                threshold=XMIPP_MAX(threshold,1000);
                threshold=(double)(XMIPP_MIN(threshold,Nimg))/Nimg;
                proceed=(rnd_unif(0,1)<threshold);
            } else
                proceed=true;
        }
        else
        {
            proceed=true;
            neighbour=true;
        }

        if (proceed) {
            // Try this image
            Matrix2D<double> Iaux=I;
            double likelihood;
	    P[q]->fit(Iaux,sigma,noMirror,corrCodeList(q),likelihood);
            if ((!classicalMultiref &&
                    ((likelihood>ABS(bestLikelihood)) ||
                    (ABS(likelihood)>ABS(bestLikelihood) && bestLikelihood<0)))
                || (classicalMultiref && corrCodeList(q)>bestCorrCode)
                || bestCorrCode<0)
            {
                if (neighbour)
                {
                    bestq=q;
                    bestImg=Iaux;
                    bestCorrCode=corrCodeList(q);
                    bestLikelihood=likelihood;
                }
            }
        }
    }

    I=bestImg;
    newnode=bestq;
    corrCode=bestCorrCode;
    likelihood=bestLikelihood;

    // Assign it to the new node and remove it from the rest
    // of nodes if it was among the best
    P[newnode]->updateProjection(I,corrCode,idx);
    if (!classicalMultiref)
        for (int q=0; q<Q; q++)
    	    if (q!=newnode && corrCodeList(q)>0)
	        P[q]->updateNonProjection(corrCodeList(q));
}

void VQ::transferUpdates()
{
    int Q=P.size();
    for (int q=0; q<Q; q++)
        P[q]->transferUpdate();
}

void VQ::updateNonCode(Matrix2D<double> &I, int newnode)
{
    int Q=P.size();
    
    for (int q=0; q<Q; q++)
    {
        if (q==newnode) continue;
        
        // Try this image
        Matrix2D<double> Iaux=I;
        double corrNonCode, nonCodeLikelihood;
        P[q]->fit(Iaux,sigma,noMirror,corrNonCode,nonCodeLikelihood);
        P[q]->updateNonProjection(corrNonCode);
    }
}

/* Run VQ ------------------------------------------------------------------ */
//#define DEBUG
void VQ::run(const FileName &fnOut, int level, int rank)
{
    int N=SF->ImgNo();
    int Q=P.size();
        
    if( rank == 0 )
    	std::cout << "Quantizing with " << Q << " codes...\n";
    
    Matrix1D<int> oldAssignment;
    oldAssignment.resize(N);
    oldAssignment.initConstant(-1);
        
    Matrix1D<int> newAssignment;
    newAssignment.resize(N);
    
    Matrix1D<int> auxAssignment;
    auxAssignment.resize(N);
      
    Matrix2D<double> aux_matrix;   
    int n=0;
    
    bool goOn=true;
    
    int mpi_size;
    
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
	 
    while (goOn)
    {
        if( rank == 0 )
        {
	    std::cout << "Iteration " << n << " ...\n";
	    std::cerr << "Iteration " << n << " ...\n";
            if (verbose)
        	for (int q=0; q<Q; q++)
        	    std::cout << "Before q=" << q
                              << " Nq=" << P[q]->currentListImg.size()
                              << std::endl;
            init_progress_bar(N);
        }
        
	SF->go_first_ACTIVE();
        
	int i=0;
        
	newAssignment.initConstant(-1);
        

        int K=XMIPP_MIN(Nneighbours+1,Q);
        if (K==0) K=Q;
	for (int q=0; q<Q; q++)
        {
            P[q]->nextListImg.clear();
            P[q]->lookForNeighbours(P, sigma, noMirror, K);
        }
        
	ImageXmipp I;
        int node;
	double corrCode, likelihood;
	double corrCodeSum=0;
        int progressStep=XMIPP_MAX(1,N/60);
	while (!SF->eof()) {
            if( i % mpi_size == rank )    
	    {
	        I.read(SF->NextImg());
                I().setXmippOrigin();
                I().statisticsAdjust(0,1);	    
	    
                lookNode(I(),i,oldAssignment(i),node,corrCode,likelihood);
        	newAssignment(i)=node;
                corrCodeSum+=corrCode;
            
	    	if( rank == 0 && (i%progressStep==0)) progress_bar(i);
	    }
            else
	    	SF->NextImg();
	    
	    i++;
        }
	
	
	// All nodes need to know the "global" value for corrSum in order to proceed correctly
	double corrCodeAux=0;
	
	MPI_Allreduce( &corrCodeSum, &corrCodeAux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	corrCodeSum = corrCodeAux;
	
	double aux_double;
	
	for (int q=0; q<Q; q++)
        {
	    aux_matrix.initZeros( YSIZE(I()), XSIZE(I()) );
	
	    MPI_Allreduce( P[q]->Pupdate.data, aux_matrix.data, P[q]->Pupdate.yxdim, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD );
	    P[q]->Pupdate = aux_matrix;
	    P[q]->Pupdate.setXmippOrigin();
	    
	    // Share nextClassCorr and nextNonClassCorr
	    int oldSizeClassCorr =  P[q]->nextClassCorr.size();
	    int oldSizeNonClassCorr = P[q]->nextNonClassCorr.size();
	    int oldSizeNextListImg = P[q]->nextListImg.size();
	    
	    std::vector<double> oldListClassCorr = P[q]->nextClassCorr;
	    std::vector<double> oldListNonClassCorr = P[q]->nextNonClassCorr;
	    std::vector<int> oldListNextListImg = P[q]->nextListImg;
	    
	    for( int ranks = 0 ; ranks < mpi_size ; ranks ++ )
	    {
	    	if( ranks == rank ) 
	     	{
			MPI_Bcast( &oldSizeClassCorr, 1, MPI_INT, rank, MPI_COMM_WORLD );
			MPI_Bcast( &(oldListClassCorr[0]), oldSizeClassCorr, MPI_DOUBLE, rank, MPI_COMM_WORLD );    
            	
			MPI_Bcast( &oldSizeNonClassCorr, 1, MPI_INT, rank, MPI_COMM_WORLD );
			MPI_Bcast( &(oldListNonClassCorr[0]), oldSizeNonClassCorr, MPI_DOUBLE, rank, MPI_COMM_WORLD );    
			
			MPI_Bcast( &oldSizeNextListImg, 1, MPI_INT, rank, MPI_COMM_WORLD );
			MPI_Bcast( &(oldListNextListImg[0]),oldSizeNextListImg, MPI_INT, rank, MPI_COMM_WORLD );
		}
		else
		{
			int size;
			MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
			std::vector<double> aux(size, 0);
			MPI_Bcast( &(aux[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );
			
			for( int element = 0; element < size; element ++ )
			{
                                P[q]->nextClassCorr.push_back(aux[element]);
			}
			
			MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
			std::vector<double> aux2(size, 0);
			MPI_Bcast( &(aux2[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );
			
			for( int element = 0; element < size; element ++ )
			{
                                P[q]->nextNonClassCorr.push_back(aux2[element]);
			}
			
			MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
			std::vector<int> auxint(size, 0);
			MPI_Bcast( &(auxint[0]), size, MPI_INT, ranks, MPI_COMM_WORLD );
			
			for( int i = 0 ; i < size ; i++ )
			{
				P[q]->nextListImg.push_back( auxint[i] );
			}
		}
	    }
	    
	
	}
	
	if( rank == 0 )
	{
            progress_bar(N);
            std::cout << "\nAverage correlation with input vectors="
                      << corrCodeSum/N << std::endl;
	}

	MPI_Allreduce( newAssignment.data, auxAssignment.data, N , MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	
	newAssignment = auxAssignment;
	
	int Nchanges=0;
        if (n>0)
        {
            FOR_ALL_ELEMENTS_IN_MATRIX1D(newAssignment)
                if (newAssignment(i)!=oldAssignment(i))
                    Nchanges++;
            if( rank == 0 )
	    	std::cout << "Number of assignment changes=" << Nchanges
                      << std::endl;
        }
        	
	oldAssignment=newAssignment;

        transferUpdates();
	
	if (rank == 0 && verbose)
            for (int q=0; q<Q; q++)
            {
        	std::cout << "After q=" << q << " Nq="
                          << P[q]->currentListImg.size()
                          << std::endl;
            }
	
        // Check if there are empty nodes
        bool smallNodes;
        do
        {
            smallNodes=false;
            int largestNode=-1, sizeLargestNode=-1,
                smallNode=-1, sizeSmallestNode=N+1;
            for (int q=0; q<Q; q++)
            {
		if (P[q]->currentListImg.size()<sizeSmallestNode)
                {
                    smallNode=q;
                    sizeSmallestNode=P[q]->currentListImg.size();
                }
                if ((int)(P[q]->currentListImg.size())>sizeLargestNode)
                {
                    sizeLargestNode=P[q]->currentListImg.size();
                    largestNode=q;
                }
            }
            if (sizeSmallestNode<PminSize*N/Q*0.01 )
            {
                if (rank==0 && verbose)
                    std::cout << "Splitting node " << largestNode
                              << " by overwriting " << smallNode << std::endl;
                smallNodes=true;
                
                // Clear the old assignment of the images in the small node
                for (int ii=0; ii<P[smallNode]->currentListImg.size(); ii++)
                    oldAssignment(P[smallNode]->currentListImg[ii])=-1;
                std::vector<int> toReassign;
                for (int ii=0; ii<P[largestNode]->currentListImg.size(); ii++)
                {
                    newAssignment(P[largestNode]->currentListImg[ii])=-1;
                    toReassign.push_back(
                        P[largestNode]->currentListImg[ii]);
                }
                
                // Now split the largest node
                #ifdef DEBUG
                    if (rank==0)
                    {
                        std::cout << "Largest node address: " << P[largestNode] << std::endl;
                        std::cout << "Small node address: " << P[smallNode] << std::endl;
                    }
                #endif
                VQProjection *node1=new VQProjection;
                VQProjection *node2=new VQProjection;
                node1->useCorrelation=useCorrelation;
                node2->useCorrelation=useCorrelation;
                node1->useFixedCorrentropy=useFixedCorrentropy;
                node2->useFixedCorrentropy=useFixedCorrentropy;
                node1->classicalMultiref=classicalMultiref;
                node2->classicalMultiref=classicalMultiref;
                #ifdef DEBUG
                    if (rank==0) {
                        std::cout << "Node1 address: " << node1 << std::endl;
                        std::cout << "Node2 address: " << node2 << std::endl;
                    }
                #endif
		std::vector<int> finalAssignment;
                splitNode(P[largestNode],node1,node2,rank, finalAssignment);
                #ifdef DEBUG
                    if (rank==0) {
                        std::cout << "Deleting Largest node address: " << P[largestNode] << std::endl;
                        std::cout << "Deleting Small node address: " << P[smallNode] << std::endl;
                    }
                #endif
                delete P[largestNode];
                delete P[smallNode];
                P[largestNode]=node1;
                P[smallNode]=node2;
                
                for (int i=0; i<toReassign.size(); i++)
                {
		    newAssignment(toReassign[i]) = -1;

		    for( int j=0; j<finalAssignment.size(); j+=2)
		    {
		    	if( finalAssignment[j] == toReassign[i])
			{
			     if (finalAssignment[j+1]==1)
			         newAssignment(toReassign[i])=largestNode;
			     else
			         newAssignment(toReassign[i])=smallNode;
			     break;
			}
		    }
                    
                    // If not assigned yet, reassign it to any of the two
                    // outcoming nodes
                    if (newAssignment(toReassign[i])==-1)
                    {
                        ImageXmipp I;
                        I.read(SFv[toReassign[i]]);
                        I().setXmippOrigin();
                    
                	Matrix2D<double> Iaux1;
                	Iaux1=I();
                	double corrCode1, likelihood1;
                	P[largestNode]->fit(Iaux1, sigma, noMirror,
                            corrCode1, likelihood1);

                	Matrix2D<double> Iaux2;
                	Iaux2=I();
                	double corrCode2, likelihood2;
                	P[smallNode]->fit(Iaux2, sigma, noMirror,
                            corrCode2, likelihood2);

                	if ((!classicalMultiref &&
                                (likelihood1>likelihood2 && likelihood1>0 ||
                                likelihood1<0 && likelihood2<0 &&
                                   ABS(likelihood1)>ABS(likelihood2)))
                             || (classicalMultiref && corrCode1>corrCode2)
                            )
                	{
                            P[largestNode]->updateProjection(Iaux1,
                                corrCode1,toReassign[i]);
                            newAssignment(toReassign[i])=largestNode;
                            if (!classicalMultiref)
                                P[smallNode]->updateNonProjection(corrCode2);
                	}
                	else
                	{
                            P[smallNode]->updateProjection(Iaux2,corrCode2,
                                toReassign[i]);
                            newAssignment(toReassign[i])=smallNode;
                            if (!classicalMultiref)
                                P[largestNode]->updateNonProjection(corrCode1);
                	}
                    }
		}
            }
        } while (smallNodes);
        
        currentAssignment=newAssignment;
        if (rank==0)
	    write(fnOut+"_level_"+integerToString(level,2)+"_");

        if (n>0 && Nchanges<0.005*N && Q>1 || n>=(Niter-1)) 
		goOn=false;
        n++;
    }
    
    std::sort(P.begin(),P.end(),SDescendingClusterSort());
}

/* Clean ------------------------------------------------------------------- */
int VQ::cleanEmptyNodes()
{
    int retval=0;
    std::vector<VQProjection *>::iterator ptr=P.begin();
    while (ptr!=P.end())
        if ((*ptr)->currentListImg.size()==0)
        {
            ptr=P.erase(ptr);
            retval++;
        }
        else
            ptr++;
    return retval;
}

/* Split ------------------------------------------------------------------- */
//#define DEBUG
void VQ::splitNode(VQProjection *node,
    VQProjection *&node1, VQProjection *&node2, int rank,
    std::vector<int> &finalAssignment) const
{
    bool finish=true;
    std::vector<VQProjection *> toDelete; 
    Matrix1D<int> newAssignment;

    int mpi_size;
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
    
    int Ninitial=node->currentListImg.size();
	
    do
    {
        finish=true;
        
        node1->gaussianInterpolator=&gaussianInterpolator;
        node1->plans=NULL;
        node1->mask=&mask;
        node1->Pupdate.initZeros(YSIZE(node->Pupdate),
                                 XSIZE(node->Pupdate));
        node1->Pupdate.setXmippOrigin();
        node1->neighboursIdx=node->neighboursIdx;
        node1->P=node->P;
        node1->classicalMultiref=false;

        node2->gaussianInterpolator=&gaussianInterpolator;
        node2->plans=NULL;
        node2->mask=&mask;
        node2->Pupdate.initZeros(YSIZE(node->Pupdate),
                                 XSIZE(node->Pupdate));
        node2->Pupdate.setXmippOrigin();
        node2->neighboursIdx=node->neighboursIdx;
        node2->P=node->P;
        node2->classicalMultiref=false;

        int imax=node->currentListImg.size();

	Matrix1D<int> oldAssignment;
        oldAssignment.resize(imax);
        oldAssignment.initConstant(-1);
        newAssignment.resize(imax);
        Matrix1D<int> auxAssignment;
        auxAssignment.resize(imax);
        Matrix1D<double> corrList;
        corrList.initZeros(imax);
	
	for (int it=0; it<Niter; it++)
        {
            node1->nextListImg.clear();
            node2->nextListImg.clear();
            
	    if( rank == 0 )
	    {	
	    	std::cerr << "Split iteration " << it << std::endl;
            	init_progress_bar(imax);
            }
	    
	    double corrThreshold;
            if (it==1 && (corrSplit || imax<0.5*Ninitial))
            {
                histogram1D hist;
                compute_hist(corrList,hist,100);
                corrThreshold=hist.percentil(50);
            }
	    
	    newAssignment.initZeros();
	    
            for (int i=0; i<imax; i++)
            {
	    	if( i % mpi_size == rank )
		{
                    // Read image
                    ImageXmipp I;
	            I.read(SFv[node->currentListImg[i]]);
                    I().setXmippOrigin();
                    I().statisticsAdjust(0,1);

                    if (it==0 && corrSplit)
                    {
                	double corrCode, likelihood;
                	node->fit(I(), sigma, noMirror, corrCode, likelihood);
                	corrList(i)=corrCode;	    
	            }
                    else if (it==0 && !corrSplit)
                    {
                	double corrCode, likelihood;
                	node->fit(I(), sigma, noMirror, corrCode, likelihood);
                	if( rnd_unif(0,1) < 0.5 ) 
			{
                            node1->updateProjection(I(),corrCode,
                        	node->currentListImg[i]);
                            newAssignment(i)=1;
                            node2->updateNonProjection(corrCode);
                	}
                	else
                	{
                            node2->updateProjection(I(),corrCode,
                        	node->currentListImg[i]);
                            newAssignment(i)=2;
                            node1->updateNonProjection(corrCode);
                	}
                    }
                    else if (it==1 && (corrSplit || imax<0.5*Ninitial))
                    {
                	double corrCode, likelihood;
                	node->fit(I(), sigma, noMirror, corrCode, likelihood);
                	if( corrCode < corrThreshold ) 
			{
                            node1->updateProjection(I(),corrCode,
                        	node->currentListImg[i]);
                            newAssignment(i)=1;
                            node2->updateNonProjection(corrCode);
                	}
                	else
                	{
                            node2->updateProjection(I(),corrCode,
                        	node->currentListImg[i]);
                            newAssignment(i)=2;
                            node1->updateNonProjection(corrCode);
                	}
                    }
                    else
                    {
                	Matrix2D<double> Iaux1;
                	Iaux1=I();
                	double corrCode1, likelihood1;
                	node1->fit(Iaux1, sigma, noMirror, corrCode1, likelihood1);

                	Matrix2D<double> Iaux2;
                	Iaux2=I();
                	double corrCode2, likelihood2;
                	node2->fit(Iaux2, sigma, noMirror, corrCode2, likelihood2);

                        #ifdef DEBUG
                            if (rank==0 && false) std::cout << "Previously Assigned to " << newAssignment(i) << std::endl;
                        #endif

                	if (likelihood1>likelihood2 && likelihood1>0 ||
                                likelihood1<0 && likelihood2<0 &&
                                   ABS(likelihood1)>ABS(likelihood2))
                	{
                            node1->updateProjection(Iaux1,corrCode1,node->currentListImg[i]);
                            newAssignment(i)=1;
                            node2->updateNonProjection(corrCode2);
                	}
                	else
                	{
                            node2->updateProjection(Iaux2,corrCode2,node->currentListImg[i]);
                            newAssignment(i)=2;
                            node1->updateNonProjection(corrCode1);
                	}

                        #ifdef DEBUG
                            if (rank==0 && false)
                            {
                                std::cout << "Lik1=" << likelihood1 << " Lik2="
                                          << likelihood2 << std::endl;
                                std::cout << "Cor1=" << corrCode1 << " Cor2="
                                          << corrCode2 << std::endl;
                                std::cout << "Assigned to " << newAssignment(i) << std::endl;
                                ImageXmipp save;
                                save()=I(); save.write("PPPI.xmp");
                                save()=Iaux1; save.write("PPPIaux1.xmp");
                                save()=Iaux2; save.write("PPPIaux2.xmp");
                                std::cout << "Press any key\n";
                                char c; std::cin >> c;
                            }
                        #endif
                    }
		    if (rank==0 && i%25==0) progress_bar(i);
		}
            }
	    
	    // Share among other mpi workers
	    MPI_Allreduce( newAssignment.data, auxAssignment.data, imax , MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	
	    if (it==0 && corrSplit)
            {
             	Matrix1D<double> corrListAux;
        	corrListAux.initZeros(imax);
	  	
		MPI_Allreduce( corrList.data, corrListAux.data, imax, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            	corrList = corrListAux;
	    }
	    
	    newAssignment = auxAssignment;
	    Matrix2D<double> PupdateAux;
    	    PupdateAux.resize( node1->Pupdate.ydim, node1->Pupdate.xdim );
    	    PupdateAux.setXmippOrigin( );
   	    
	    PupdateAux.initZeros();
	    MPI_Allreduce( node1->Pupdate.data, PupdateAux.data, node1->Pupdate.yxdim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            node1->Pupdate = PupdateAux;
	    
	    PupdateAux.initZeros();
	    MPI_Allreduce( node2->Pupdate.data, PupdateAux.data, node2->Pupdate.yxdim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            node2->Pupdate = PupdateAux;

	    // Share nextClassCorr and nextNonClassCorr
	    int oldSizeNonClassCorr1 = node1->nextNonClassCorr.size();
	    int oldSizeNonClassCorr2 = node2->nextNonClassCorr.size();
	    int oldSizeClassCorr1 =  node1->nextClassCorr.size();
	    int oldSizeClassCorr2 =  node2->nextClassCorr.size();
	    int oldSizeNextListImg1 = node1->nextListImg.size();
	    int oldSizeNextListImg2 = node2->nextListImg.size();
	    
	    std::vector<double> oldListNonClassCorr1 = node1->nextNonClassCorr;
	    std::vector<double> oldListNonClassCorr2 = node2->nextNonClassCorr;
	    std::vector<double> oldListClassCorr1 = node1->nextClassCorr;
	    std::vector<double> oldListClassCorr2 = node2->nextClassCorr;
	    std::vector<int> oldListNextListImg1 = node1->nextListImg;
	    std::vector<int> oldListNextListImg2 = node2->nextListImg;
	    
	    for( int ranks = 0 ; ranks < mpi_size ; ranks ++ )
	    {
	    	if( ranks == rank ) 
	     	{
		    MPI_Bcast( &oldSizeNonClassCorr1, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListNonClassCorr1[0]), oldSizeNonClassCorr1, MPI_DOUBLE, rank, MPI_COMM_WORLD );    

		    MPI_Bcast( &oldSizeNonClassCorr2, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListNonClassCorr2[0]), oldSizeNonClassCorr2, MPI_DOUBLE, rank, MPI_COMM_WORLD );    

		    MPI_Bcast( &oldSizeClassCorr1, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListClassCorr1[0]), oldSizeClassCorr1, MPI_DOUBLE, rank, MPI_COMM_WORLD );    

		    MPI_Bcast( &oldSizeClassCorr2, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListClassCorr2[0]), oldSizeClassCorr2, MPI_DOUBLE, rank, MPI_COMM_WORLD );    

		    MPI_Bcast( &oldSizeNextListImg1, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListNextListImg1[0]), oldSizeNextListImg1, MPI_INT, rank, MPI_COMM_WORLD );    

		    MPI_Bcast( &oldSizeNextListImg2, 1, MPI_INT, rank, MPI_COMM_WORLD );
		    MPI_Bcast( &(oldListNextListImg2[0]), oldSizeNextListImg2, MPI_INT, rank, MPI_COMM_WORLD );    
		}
		else
		{
		    int size, element;
		    std::vector<double> daux;
		    std::vector<int> iaux;

		    // Receive and gather nextNonClassCorr values
		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
		    daux.resize(size, 0);
		    MPI_Bcast( &(daux[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );

		    for( element = 0; element < size; element ++ )
                        node1->nextNonClassCorr.push_back(daux[element]);

		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
		    daux.resize(size, 0);
		    MPI_Bcast( &(daux[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );

		    for( element = 0; element < size; element ++ )
                        node2->nextNonClassCorr.push_back(daux[element]);

		    // Receive and gather nextClassCorr values
		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
		    daux.resize(size, 0);
		    MPI_Bcast( &(daux[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );

		    for( element = 0; element < size; element ++ )
                        node1->nextClassCorr.push_back(daux[element]);

		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
		    daux.resize(size, 0);
		    MPI_Bcast( &(daux[0]), size, MPI_DOUBLE, ranks, MPI_COMM_WORLD );

		    for( element = 0; element < size; element ++ )
                        node2->nextClassCorr.push_back(daux[element]);

		    // Receive and gather nextListImg values
		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
		    iaux.resize(size, 0);
		    MPI_Bcast( &(iaux[0]), size, MPI_INT, ranks, MPI_COMM_WORLD );

		    for( element = 0; element < size; element ++ )
                        node1->nextListImg.push_back(iaux[element]);

		    MPI_Bcast( &size, 1, MPI_INT ,ranks, MPI_COMM_WORLD );
		    iaux.resize(size, 0);
		    MPI_Bcast( &(iaux[0]), size, MPI_INT, ranks, MPI_COMM_WORLD );

		    for( element = 0; element < size; element ++ )
                        node2->nextListImg.push_back(iaux[element]);
		}
	    }
	   
	    if( rank == 0 )
	    	progress_bar(imax);
	    
            if (imax<0.5*Ninitial && it==1) break;
            
	    int Nchanges;
            if (it>=1 || !corrSplit)
            {
                node1->transferUpdate();
                node2->transferUpdate();

                if( rank == 0 && verbose)
                {
	            std::cout
                        << "  Split iteration " << it << std::endl
                        << "  node1 Nq=" << node1->currentListImg.size()
                        << std::endl
                        << "  node2 Nq=" << node2->currentListImg.size()
                        << std::endl;
                    ImageXmipp save;
                    save()=node1->P; save.write("PPPnode1.xmp");
                    save()=node2->P; save.write("PPPnode2.xmp");
                }

                if (it>=2 || (!corrSplit && it>=1))
                {
                    Nchanges=0;
                    FOR_ALL_ELEMENTS_IN_MATRIX1D(newAssignment)
                        if (newAssignment(i)!=oldAssignment(i))
                            Nchanges++;
                    if( rank == 0 && verbose)
	    	        std::cout << "Number of assignment split changes=" << Nchanges << std::endl;
                }

                // Check if one of the nodes is too small
                if (node1->currentListImg.size()<PminSize*0.01*imax/2 ||
                    node2->currentListImg.size()<PminSize*0.01*imax/2)
                    if (it>1 || (!corrSplit && it>0)) break;

                oldAssignment=newAssignment;
            }
	    if (Nchanges<0.005*imax && it>1) break;
        }

        if (node1->currentListImg.size()<(PminSize*0.01*imax/2))
        {
	    if( rank == 0 && verbose)
            	std::cout << "Removing node1, it's too small "
                          << node1->currentListImg.size() << " "
                          << PminSize*0.01*imax/2 << "...\n";
            if (node1!=node) delete node1;
            node1=new VQProjection();
            node1->useCorrelation=useCorrelation;
            node1->useFixedCorrentropy=useFixedCorrentropy;
            node1->classicalMultiref=false;
            toDelete.push_back(node2);
            node=node2;
            node2=new VQProjection();
            node2->useCorrelation=useCorrelation;
            node2->useFixedCorrentropy=useFixedCorrentropy;
            node2->classicalMultiref=false;
            finish=false;
        }
        else if (node2->currentListImg.size()<(PminSize*0.01*imax/2))
        {
	    if( rank == 0 && verbose)
            	std::cout << "Removing node2, it's too small "
                          << node2->currentListImg.size() << " "
                          << PminSize*0.01*imax/2 << "...\n";
            if (node2!=node) delete node2;
            node2=new VQProjection();
            node2->useCorrelation=useCorrelation;
            node2->useFixedCorrentropy=useFixedCorrentropy;
            node2->classicalMultiref=false;
            toDelete.push_back(node1);
            node=node1;
            node1=new VQProjection();
            node1->useCorrelation=useCorrelation;
            node1->useFixedCorrentropy=useFixedCorrentropy;
            node1->classicalMultiref=false;
            finish=false;
        }
    } while (!finish);
    for (int i=0; i<toDelete.size(); i++)
        if (toDelete[i]!=node) delete toDelete[i];
   
    for (int i=0; i<node->currentListImg.size(); i++)
    {
        finalAssignment.push_back( node->currentListImg[i] ); 
	finalAssignment.push_back( newAssignment(i) );
    }    
    node1->classicalMultiref=classicalMultiref;
    node2->classicalMultiref=classicalMultiref;
}
#undef DEBUG

void VQ::splitFirstNode(int rank) {
    std::sort(P.begin(),P.end(),SDescendingClusterSort());
    int Q=P.size();
    P.push_back(new VQProjection());
    P.push_back(new VQProjection());
    P[Q]->useCorrelation=useCorrelation;
    P[Q]->useFixedCorrentropy=useFixedCorrentropy;
    P[Q]->classicalMultiref=classicalMultiref;
    P[Q+1]->useCorrelation=useCorrelation;
    P[Q+1]->useFixedCorrentropy=useFixedCorrentropy;
    P[Q+1]->classicalMultiref=classicalMultiref;
    std::vector<int> finalAssignment;
    splitNode(P[0],P[Q],P[Q+1],rank, finalAssignment);
    delete P[0];
    P[0]=NULL;
    P.erase(P.begin());    
}

/* VQPrm I/O --------------------------------------------------------------- */
void Prog_VQ_prm::read(int argc, char** argv)
{
    fnSel=getParameter(argc,argv,"-i");
    fnOut=getParameter(argc,argv,"-o","class");
    fnCodes0=getParameter(argc,argv,"-codesSel0","");
    Niter=textToInteger(getParameter(argc,argv,"-iter","20"));
    Nneighbours=textToInteger(getParameter(argc,argv,"-neigh","4"));
    Ncodes0=textToInteger(getParameter(argc,argv,"-codes0","2"));
    Ncodes=textToInteger(getParameter(argc,argv,"-codes","16"));
    PminSize=textToFloat(getParameter(argc,argv,"-minsize","20"));
    noMirror=checkParameter(argc,argv,"-no_mirror");
    verbose=checkParameter(argc,argv,"-verbose");
    corrSplit=checkParameter(argc,argv,"-corr_split");
    fast=checkParameter(argc,argv,"-fast");
    useCorrelation=checkParameter(argc,argv,"-useCorrelation");
    useFixedCorrentropy=checkParameter(argc,argv,"-useFixedCorrentropy");
    classicalMultiref=checkParameter(argc,argv,"-classicalMultiref");
    alignImages=checkParameter(argc,argv,"-alignImages");
}

void Prog_VQ_prm::show() const
{
    std::cout << "Input images:            " << fnSel               << std::endl
              << "Output images:           " << fnOut               << std::endl
              << "Iterations:              " << Niter               << std::endl
              << "CodesSel0:               " << fnCodes0            << std::endl
              << "Codes0:                  " << Ncodes0             << std::endl
              << "Codes:                   " << Ncodes              << std::endl
              << "Neighbours:              " << Nneighbours         << std::endl
              << "Minimum node size:       " << PminSize            << std::endl
              << "No mirror:               " << noMirror            << std::endl
              << "Verbose:                 " << verbose             << std::endl
              << "Corr Split:              " << corrSplit           << std::endl
              << "Fast:                    " << fast                << std::endl
              << "Use Correlation:         " << useCorrelation      << std::endl
              << "Use Fixed Correntropy:   " << useFixedCorrentropy << std::endl
              << "Classical Multiref:      " << classicalMultiref   << std::endl
              << "Align images:            " << alignImages         << std::endl;
}
    
void Prog_VQ_prm::usage() const
{
    std::cerr << "Usage:\n"
              << "    -i <selfile>         : Selfile with the input images\n"
              << "   [-o <root=class>]     : Output rootname, by default, class\n"
              << "   [-iter <N=20>]        : Number of iterations\n"
              << "   [-codes0 <N=2]        : Initial number of code vectors\n"
              << "   [-codesSel0 <selfile>]: Selfile with initial code vectors\n"
              << "   [-codes <N=16>]       : Final number of code vectors\n"
              << "   [-neigh <N=4>]        : Number of neighbour code vectors\n"
              << "                         : Set -1 for all\n"
              << "   [-minsize <N=20>]     : Percentage minimum node size\n"
              << "   [-no_mirror]          : Do not check mirrors\n"
              << "   [-verbose]            : Verbose\n"
              << "   [-corr_split]         : Correlation split\n"
              << "   [-fast]               : Fast calculations, suboptimal\n"
              << "   [-useCorrelation]     : Instead of correntropy\n"
              << "   [-useFixedCorrentropy]: Instead of correntropy\n"
              << "   [-classicalMultiref]  : Instead of enhanced clustering\n"
              << "   [-alignImages]        : Align the original images at the end\n"
    ;
}
    
void Prog_VQ_prm::produce_side_info(int rank)
{
    SF.read(fnSel);
    std::vector< Matrix2D<double> > codes0;
    if (fnCodes0!="")
    {
        SelFile SFCodes(fnCodes0);
        while (!SFCodes.eof())
        {
            ImageXmipp I;
            I.read(SFCodes.NextImg());
            I().setXmippOrigin();
            codes0.push_back(I());
        }
        Ncodes0=codes0.size();
    }
    vq.initialize(SF,Niter,Nneighbours,PminSize,
        codes0,Ncodes0,noMirror,verbose,corrSplit,useCorrelation,
        useFixedCorrentropy,classicalMultiref,alignImages,fast,rank);
}

void Prog_VQ_prm::run(int rank)
{
    int level=0;
    vq.run(fnOut,level,rank);
    
    int Q=vq.P.size();
    
    while (Q<Ncodes)
    {
        if( rank == 0 )
		std::cout << "Spliting nodes ...\n";
        
        int Nclean=vq.cleanEmptyNodes();
	int Nsplits=XMIPP_MIN(Q,Ncodes-Q)+Nclean;
		
	for (int i=0; i<Nsplits; i++)
            vq.splitFirstNode(rank);
        
	Q=vq.P.size();	
        level++;
        vq.run(fnOut,level,rank);
    }
    if (rank==0)
    {
	std::sort(vq.P.begin(),vq.P.end(),SDescendingClusterSort());
    	vq.write(fnOut,true);
    }
}

/* Main -------------------------------------------------------------------- */
int main(int argc, char** argv)
{
    Prog_VQ_prm prm;
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    try {
        prm.read(argc,argv);
    } catch (Xmipp_error XE)
    {
    	if (rank==0)
	{
            std::cout << XE << std::endl;
            prm.usage();
	}
	
	MPI_Finalize();
        return 1;
    }

    try {
        if (rank==0) prm.show();
	prm.produce_side_info(rank);
        prm.run(rank);
    } catch (Xmipp_error XE)
    {
        std::cout << XE << std::endl;
	MPI_Finalize();
        return 1;
    }
    MPI_Finalize();
    return 0;
}
