#include "npe.h"
#include <data/matrix2d.h>
#include <data/matrix1d.h>


void NPE::setSpecificParameters(int k){
	this->k=k;
}

/*
 * Reduce dimensionality method based on NPE algorithm
 */
void NPE::reduceDimensionality()
{

	X->write("dimred/X.txt");
	if (MAT_YSIZE(*X) <= MAT_XSIZE(*X))
			REPORT_ERROR(ERR_MATRIX_DIM, "Number of samples should be higher than number of dimensions.");

	int n = MAT_YSIZE(*X);

	//Find nearest neighbours
	Matrix2D<double> D;
	Matrix2D<int> idx;
	kNearestNeighbours(*X,k,idx,D,distance,false);
	D.write("dimred/D.txt");

	Matrix2D<double> W(k,n), Xi, C, M;
	Matrix1D<double> wi;

	PseudoInverseHelper h;

	h.b.resizeNoCopy(k);
	h.b.initConstant(1);

	double tol=1e-5;
	for (int ip=0; ip<n; ++ip)
	{

		extractNearestNeighbours(*X, idx, ip, Xi);

		// Move the origin to observation ip
		FOR_ALL_ELEMENTS_IN_MATRIX2D(Xi)
			MAT_ELEM(Xi,i,j)-=MAT_ELEM(*X,ip,j);

		matrixOperation_AAt(Xi,C);

		double K=C.trace()*tol;
		for (int ii=0; ii<MAT_XSIZE(C); ++ii)
			MAT_ELEM(C,ii,ii)+=K;

		h.A=C;

		solveLinearSystem(h,wi);
		wi/=wi.sum();

		W.setCol(ip,wi);
	}


	W.write("dimred/W.txt");
	std::cout << "Finished loop" << std::endl;
	idx.write("dimred/idx.txt");
	Xi.write("dimred/Xi.txt");

	//Find the sparse cost matrix
	M.initIdentity(1000);
	Matrix1D<int> neighboursi;


	for(int i=0;i<n; ++i)
	{
		// Get wi for this observation
		W.getCol(i,wi);

		// Get the neighbours of this observation
		idx.getRow(i,neighboursi);

		for(size_t p1=0;p1<VEC_XSIZE(neighboursi); p1++)
		{
			int j1=VEC_ELEM(neighboursi,p1); // j is the index of the neighbour
			double w1=VEC_ELEM(wi,p1);
			MAT_ELEM(M,i,j1)-=w1;
			MAT_ELEM(M,j1,i)-=w1;
			for (size_t p2=0; p2<VEC_XSIZE(neighboursi); p2++)
			{
				int j2=VEC_ELEM(neighboursi,p2);
				MAT_ELEM(M,j1,j2)+=w1*VEC_ELEM(wi,p2);
			}

		}
	}
	M.write("dimred/M.txt");
	std::cout << "Finished sparse cost matrix" << std::endl;

	//Check symmetry
	Matrix2D<double> DP, WP;
	M.write("dimred/MFinal.txt");

	X->transpose();
	matrixOperation_XtAX_symmetric(*X,M,WP);
	std::cout << "WP-> " << MAT_YSIZE(WP) << "x" << MAT_XSIZE(WP) << std::endl;
	matrixOperation_AtA(*X, DP);

	std::cout << "Symmetry checked" << std::endl;
	DP.write("dimred/DP.txt");
	WP.write("dimred/WP.txt");


	//Solve eigenvector problem
	Matrix2D<double> Peigvec, eigvector;
	Matrix1D<double> Deigval;
	generalizedEigs(WP,DP,Deigval,Peigvec);

	//Sort eigenvalues

	Matrix1D<int> idx2;
	Deigval.indexSort(idx2);

	eigvector.resizeNoCopy(MAT_YSIZE(Peigvec),outputDim);
	for(size_t j =0;j<outputDim;++j){
		int idxj=VEC_ELEM(idx2,j)-1;
		for(size_t i=0;i<MAT_YSIZE(Peigvec);++i)
			MAT_ELEM(eigvector, i, j)=MAT_ELEM(Peigvec,i,idxj);
	}

	eigvector.operator *=(-1);
	//Compute results
	Y=*X*eigvector;

	Y.write("dimred/Y.txt");
	std::cout << "Reduction completed" << std::endl;

}
