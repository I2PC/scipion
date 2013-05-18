/***************************************************************************
 *
 * Authors:    Clara M. Mejias Perez      clrmejias7@gmail.com (2013)
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

#include "chartingmanifold.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <data/xmipp_image.h>

    /*Para copiar la matriz3D
     * Image<double> save;
    save()=CZZ;
    save.write("dimred/CZZ.vol");*/

//Column Max
	void columnMax(Matrix2D<double> &A, Matrix1D<double> &MAXA){
		MAXA.resizeNoCopy(MAT_XSIZE(A));
        for (int j=0; j<MAT_XSIZE(A); ++j){
       		double max=-1e38;
       		for(int i=0; i<MAT_YSIZE(A); i++){
       			if (MAT_ELEM(A,i,j)>max)
       				max=MAT_ELEM(A,i,j);
       		}
       		VEC_ELEM(MAXA,j)=max;
       	}
	}

//Sumatory of the columns
	//NOTA IMPORTANTE: antes no me daba problemas pero ahora salta error. He creado una funcion (colSum) en matrix2d.h basandome en rowSum y funciona bien
	void columnSum(const Matrix2D<double> &A, Matrix1D<double> &S){
		S.resizeNoCopy(MAT_XSIZE(A));
		double sum=.0;
		for (int i=0; i<MAT_XSIZE(A); i++){
			for (int j=0; j<MAT_YSIZE(A); j++)
				sum+=MAT_ELEM(A,i,j);
			VEC_ELEM(S,i)=sum;
		}
	}

//Compute Matrix2D^n
	void matrixPow (const Matrix2D<double> &A, Matrix2D<double> &A2, int n){
		A2.resizeNoCopy(MAT_YSIZE(A),MAT_XSIZE(A));
		FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
				MAT_ELEM(A2,i,j)=pow(MAT_ELEM(A,i,j),n);
	}

//Mean of the elements in a Matrix1D
double meanVector (const Matrix1D<double> &A){
	double vm;
	FOR_ALL_ELEMENTS_IN_MATRIX1D(A)
		vm=vm+VEC_ELEM(A,i);
	vm=vm/VEC_XSIZE(A);
	return (vm);
}

//Array Multriply --> bsxfun(@times) Aij*Vj
void arrayMultiply(const Matrix1D<double> &A, const Matrix2D<double> &B, Matrix2D<double> &ANS){
	ANS.resizeNoCopy(MAT_YSIZE(B),MAT_XSIZE(B));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(B){
		MAT_ELEM(ANS,i,j)=MAT_ELEM(B,i,j)*VEC_ELEM(A,j);
	}
}

//Array Multriply --> bsxfun(@times) Akij*V1ij
void arrayMultiply3D(const MultidimArray<double> &A, const MultidimArray<double> &B, MultidimArray<double> &ANS){
	Matrix2D<double>AUXB;
	AUXB.resizeNoCopy(YSIZE(B), XSIZE(B));
	ANS.resizeNoCopy(ZSIZE(A),YSIZE(A),XSIZE(A));
	FOR_ALL_ELEMENTS_IN_ARRAY3D(B)
		MAT_ELEM(AUXB,i,j)=A3D_ELEM(B,k,i,j);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(ANS)
		A3D_ELEM(ANS,k,i,j)=A3D_ELEM(A,k,i,j)*MAT_ELEM(AUXB,i,j);
}

//Sumatory of Rows
void rowSum(const Matrix2D<double> &A, Matrix1D<double> &S){
	S.resizeNoCopy(MAT_YSIZE(A));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(A){
		VEC_ELEM(S,i)+=MAT_ELEM(A,i,j);
	}
}

//Array division
void arrayDivision(const Matrix2D<double> &Q, const Matrix1D<double> P, Matrix2D<double> &D){
	FOR_ALL_ELEMENTS_IN_MATRIX2D(Q)
		MAT_ELEM(D,i,j)=MAT_ELEM(Q,i,j)/VEC_ELEM(P,j);
}

//bsxfun(@times,M,N) --> Multiplicar mij*nij
void elementsMultiply (const Matrix2D<double>&X, const Matrix2D<double>&N, Matrix2D<double> &ANS){
	ANS.resize(MAT_YSIZE(X),MAT_XSIZE(X));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(X)
		MAT_ELEM(ANS,i,j)=MAT_ELEM(X,i,j)*MAT_ELEM(N,i,j);
}

//Para 3D bsxfun(@times,M,N) --> Multiplicar mkij*nkij
void elementsMultiply3D (const MultidimArray<double>&X, const MultidimArray<double>&N, MultidimArray<double> &ANS){
	ANS.resizeNoCopy(ZSIZE(X),YSIZE(X),XSIZE(X));
	FOR_ALL_ELEMENTS_IN_ARRAY3D(ANS)
		A3D_ELEM(ANS,k,i,j)=A3D_ELEM(X,k,i,j)*A3D_ELEM(N,k,i,j);
}

//Covariance Matrix
void covariance (const Matrix2D<double> &B, Matrix2D<double> &COV){
	Matrix2D<double> C, Ct, M;
	Matrix1D<double>R;
	int y=MAT_YSIZE(B);
	int x=MAT_XSIZE(B);
	B.computeRowMeans(R);
	M.resizeNoCopy(y,x);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
		MAT_ELEM(M,i,j)=VEC_ELEM(R,i);
	C=B-M;
	Ct=C.transpose();
	COV=(1./(x-1))*C*Ct;
}

//Standard deviation of the Matrix (ROWS)
void stdRows (const Matrix2D<double> &B, Matrix1D<double>&ANS){
	double sum,desv;
	Matrix1D<double> M;
	int y=MAT_YSIZE(B);
	int x=MAT_XSIZE(B);
	ANS.resizeNoCopy(y);
	B.computeRowMeans(M);
	for (int i=0; i<y; ++i){
		sum=0;
		for(int j=0; j<x ;++j){
			sum=sum+pow(MAT_ELEM(B,i,j)-VEC_ELEM(M,j),2);
		}
		sum=sum/(x-1);
		desv=sqrt(sum);
		VEC_ELEM(ANS,i)=desv;
	}
}

//bsxfun(@plus,VEC,M2D) M2D(i,j)+VEC(i) (For all j);
void vecPlusMatrix(const Matrix2D<double> &M, const Matrix1D<double> &V, Matrix2D<double> &ANS){
	ANS.resizeNoCopy(MAT_YSIZE(M),MAT_XSIZE(M));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(M)
		MAT_ELEM(ANS,i,j)=MAT_ELEM(M,i,j)+VEC_ELEM(V,i);
}

void claraWrite(const Matrix2D<double> &X, const FileName &fn)
{
	std::ofstream fhOut;
	fhOut.open(fn.c_str());
	for (int i=0; i<MAT_YSIZE(X); ++i)
	{
		for (int j=0; j<MAT_XSIZE(X); ++j)
			fhOut << MAT_ELEM(X,i,j) << " ";
		fhOut << std::endl;
	}
	fhOut.close();
}

void claraWrite1D (const Matrix1D<double> &V, const FileName &fn){
	std::ofstream fhOut;
	fhOut.open(fn.c_str());
	for (int i=0; i<VEC_XSIZE(V); ++i)
	{
		fhOut << VEC_ELEM(V,i) << " ";
		fhOut << std::endl;
	}
	fhOut.close();
}

void vectorInverse(const Matrix1D<double> &A, Matrix1D<double> &I){
	I.resizeNoCopy(VEC_XSIZE(A));
	FOR_ALL_ELEMENTS_IN_MATRIX1D(I){
		VEC_ELEM(I,i)=1/VEC_ELEM(A,i);
	}
}

// mmpca runs EM algorithm and computes local factor analyzers
void mppca (const Matrix2D<double> &X, int no_dims, int no_analyzers, double tol, int maxiter, double minstd,
            MultidimArray<double> &LX, Matrix2D<double> &MX, Matrix1D<double> &PX)
{

    double pi=3.141592653589;
	double m,lik,sr;
	Matrix1D<double> cc,MM,PII,RS,AUXR,AUXR1,DCC,AUXMAX,MAXR,SUMR,LOGRS,r,AUXT,MCC,srz,srx,RNDM,SAUXPX;
	Matrix2D<double> SS,CC,RND1,MUL,X2,Xt,AUX1,R,L,LTPI,LTPIL,ILTPIL,ID,AUXD,CCI,COVZ,DELTA,MEANZ,AUXR2,AUX3,DELTA2,E;
    Matrix2D<double>z,srxz,m1,m2,rz,A,AUXM2,EL1, AUXDM,IAUXDM,LXij,AUXPX;
	MultidimArray<double> RNDM2,CZZ, ZZ;

	//NOTA: AL ENTRAR EN MMPCA CAMBIA LAS DIMENSIONES DE X(DxN) a X(NxD)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	A=X.transpose();

    // Initialize some  variables
    double D=MAT_YSIZE(A);       /* row size */
    double N=MAT_XSIZE(A);       /* column size */
    double epsilon=1e-9;
    double minvar=pow(minstd,2);

    // Randomly initialize factor analyzers
    A.computeRowMeans(MM);
    covariance(A,SS);

    RNDM2.resize(no_analyzers,D,no_dims);
    RNDM2.initRandom(0,1);
    LX=minstd*RNDM2;

    cholesky(SS,CC);
    CC=CC.transpose();
    RND1.initGaussian(D, no_analyzers);
    claraWrite(RND1,"dimred/RND1.txt");
    MUL=CC*RND1;
    vecPlusMatrix(MUL,MM,MX);
    claraWrite(MX,"dimred/MX.txt");
    AUXT.resizeNoCopy(D);
    AUXT.initConstant(1);
    CC.getDiagonal(MCC);
    m=meanVector(MCC);
    PX=2*m*AUXT;

    /*catch
    stdRows(A,cc);
    arrayMultiply(cc,RND1,AUX1);
    vecPlusMatrix(AUX1,MM,MX);
    PX.resizeNoCopy(D);
    meancc=meanVector(cc);
    for (int i=0;i<D;++i){
       	PX(i)=2*meancc;

    }
*/
    MM.initZeros();
    SS.initZeros();
    CC.initZeros();
    matrixPow(A,X2,2);

    //Initialize some variables
    double cnst=-(D/2)*log(2*pi);
    lik=-1e300;

    R.initZeros(no_analyzers,N);
    CZZ.initZeros(no_analyzers,no_dims,no_dims);
    ZZ.initZeros(no_analyzers,no_dims,N);
    E.resizeNoCopy(no_analyzers,N);
    L.resizeNoCopy(D,no_dims);
    ID.initIdentity(no_dims);
	AUXR2.resizeNoCopy(no_dims,N);
    AUXD.resizeNoCopy(D,N);
    Matrix1D<double>MXk;
    MXk.resizeNoCopy(D);
    AUXR2.resizeNoCopy(no_dims,N);

  	for (int i=0; i<maxiter; i++)
    {
        //E step
    	vectorInverse(PX,PII);
        for (int k=0; k< no_analyzers; ++k){
        	FOR_ALL_ELEMENTS_IN_MATRIX2D(L)
        		MAT_ELEM(L,i,j)=A3D_ELEM(LX, k, i, j);
        	arrayMultiply(PII,L,LTPI);
          	LTPI=LTPI.transpose();
          	LTPIL=LTPI*L;
        	ILTPIL=ID+LTPIL;
          	cholesky(ILTPIL,CC);
          	CC=CC.transpose();
        	CC.inv(CCI);
        	COVZ=CCI*CCI.transpose();
            FOR_ALL_ELEMENTS_IN_MATRIX1D(MXk)
            	VEC_ELEM(MXk,i)=MAT_ELEM(MX,i,k);
        	FOR_ALL_ELEMENTS_IN_MATRIX2D(AUXD)
        		MAT_ELEM(AUXD,i,j)=VEC_ELEM(MXk,i);
        	DELTA=A-AUXD;
        	MEANZ=((ID-(LTPIL*COVZ))*LTPI)*DELTA;
           	FOR_ALL_ELEMENTS_IN_MATRIX2D(COVZ)
        		A3D_ELEM(CZZ, k, i, j) = MAT_ELEM(COVZ,i,j);
        	FOR_ALL_ELEMENTS_IN_MATRIX2D(MEANZ)
        	    A3D_ELEM(ZZ, k, i, j) = MAT_ELEM(MEANZ,i,j);
        	AUX3=ILTPIL*MEANZ;
        	FOR_ALL_ELEMENTS_IN_MATRIX2D(AUXR2)
        		MAT_ELEM(AUXR2,i,j)=MAT_ELEM(MEANZ,i,j)*MAT_ELEM(AUX3,i,j);
           	AUXR2.colSum(AUXR1);
           	CC.getDiagonal(DCC);
           	matrixPow(DELTA,DELTA2,2);
        	AUXR=-0.5*(PII.transpose()*DELTA2-AUXR1)-log(DCC.sum());
           	FOR_ALL_ELEMENTS_IN_MATRIX1D(AUXR)
        		R(k,i)=AUXR(i);
        }

        // Compute responsibilities of datapoints to the clusters
        FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
        	MAT_ELEM(R,i,j)+=(cnst+.5*log(PI));
        MAXR.resizeNoCopy(MAT_XSIZE(R));
        columnMax(R,MAXR);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
        	MAT_ELEM(E,i,j)=MAT_ELEM(R,i,j)-VEC_ELEM(MAXR,j);
        R.initZeros(R);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
        	MAT_ELEM(R,i,j)=exp(MAT_ELEM(E,i,j));
        R.colSum(RS);
        arrayDivision(R,RS,R);

        // Update likelihood of estimation
        double oldlik=lik;
        LOGRS.initZeros(N);
        RS.initZeros(RS);
        MAXR.initZeros(MAXR);
        columnMax(R,MAXR);
        R.colSum(RS);
        FOR_ALL_ELEMENTS_IN_MATRIX1D(LOGRS)
        	VEC_ELEM(LOGRS,i)=log(VEC_ELEM(RS,i));
        SUMR=MAXR+LOGRS;
        lik=0;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(SUMR)
        	lik+=VEC_ELEM(SUMR,i);

        if (fabs(oldlik - lik)<tol)
            break;

        // M step
    	m2.resizeNoCopy(no_dims+1,no_dims+1);
        AUXM2.resizeNoCopy(no_dims,no_dims);
        IAUXDM.resizeNoCopy(no_dims,no_dims);
        PX.initZeros();
        z.resizeNoCopy(no_dims,N);
        m1.resizeNoCopy(N,no_dims+1);
        r.initZeros(N);
        MXk.initZeros(D);
        LXij.resizeNoCopy(D,no_dims);
        AUXDM.resizeNoCopy(m2);
        for(int k=0 ;k<no_analyzers; k++){
        	FOR_ALL_ELEMENTS_IN_MATRIX1D(r)
				VEC_ELEM(r,i)=MAT_ELEM(R,k,i);
           	FOR_ALL_ELEMENTS_IN_MATRIX2D(z)
             	MAT_ELEM(z,i,j)=A3D_ELEM(ZZ, k, i, j);
        	arrayMultiply(r,z,rz);
        	sr=r.sum();
          	rowSum(rz,srz);
        	srxz=A*rz.transpose();
        	r.isCol();
        	srx=A*r;
        	/*Aunque compile no funciona de esta forma, por lo que lo he vuelto a escribir con la otra
        	 * for (int i=0; i<D; i++){
        		MAT_ELEM(m1,i,N-1)=VEC_ELEM(srx,i);
        		memcpy(&MAT_ELEM(m1,i,0),&MAT_ELEM(srxz,i,0),MAT_XSIZE(srxz)*sizeof(double));
        	}*/
        	for (int i=0; i<no_dims+1; i++){
        		for (int j=0; j<no_dims; j++){
        			MAT_ELEM(m1,i,j)=MAT_ELEM(srxz,i,j);
					MAT_ELEM(m1,i,no_dims-1)=VEC_ELEM(srx,i);
        		}
        	}
        	FOR_ALL_ELEMENTS_IN_MATRIX2D(AUXM2)
        	    MAT_ELEM(AUXM2,i,j)=A3D_ELEM(CZZ, k, i, j);
        	EL1=sr*AUXM2+z*rz.transpose();
           	/*Mismo problema que con m1
           	 *
           	 for (int i=0; i<no_dims+1; i++){
           		MAT_ELEM(m2,i,no_dims)=VEC_ELEM(srz,i);
        		memcpy(&MAT_ELEM(m2,i,0),&MAT_ELEM(EL1,i,0),no_dims*sizeof(double));
           	}
           	memcpy(&MAT_ELEM(m2,no_dims,0),&VEC_ELEM(srz,0),no_dims*sizeof(double));*/

        	for (int j=0; j<no_dims+1; j++){
        		for (int i=0; i<no_dims; i++){
        			MAT_ELEM(m2,i,j)=MAT_ELEM(EL1,i,j);
        			MAT_ELEM(m2,i,no_dims)=VEC_ELEM(srz,i);
        		}
        			MAT_ELEM(m2,no_dims,j)=VEC_ELEM(srz,j);
        	}
        	MAT_ELEM(m2,no_dims,no_dims)=sr;
           	FOR_ALL_ELEMENTS_IN_MATRIX2D(m2)
           		MAT_ELEM(AUXDM,i,j)=MAT_ELEM(m2,i,j)+epsilon*rnd_gaus(0,1);
           	AUXDM.inv(IAUXDM);
           	m1=m1*IAUXDM;
           	{FOR_ALL_ELEMENTS_IN_ARRAY3D(LX)
           		A3D_ELEM(LX, k, i, j)=MAT_ELEM(m1,i,j);}
           	FOR_ALL_ELEMENTS_IN_MATRIX2D(MX)
           		MAT_ELEM(MX,i,j)=MAT_ELEM(m1,i,no_dims+1);
           	FOR_ALL_ELEMENTS_IN_MATRIX2D(LXij)
           		MAT_ELEM(LXij,i,j)=A3D_ELEM(LX,k,i,j);
           	AUXPX.resizeNoCopy(LXij);
           	FOR_ALL_ELEMENTS_IN_MATRIX2D(AUXPX)
           		MAT_ELEM(AUXPX,i,j)=MAT_ELEM(LXij,i,j)*MAT_ELEM(srxz,i,j);
           	rowSum(AUXPX,SAUXPX);
           	FOR_ALL_ELEMENTS_IN_MATRIX1D(MXk)
           		VEC_ELEM(MXk,i)=MAT_ELEM(MX,i,k);
        	PX = PX + X2 * r - SAUXPX - MXk* srx;
        }

        FOR_ALL_ELEMENTS_IN_MATRIX1D(PX){
        	if (minvar>VEC_ELEM(PX,i)/N)
        		VEC_ELEM(PX,i)=minvar;
        	else
        		VEC_ELEM(PX,i)=VEC_ELEM(PX,i)/N;
        }
      FOR_ALL_ELEMENTS_IN_MATRIX1D(PX)
      	  VEC_ELEM(PX,i)=meanVector(PX);
    }

}

//Inter MoFA using information from EM algorithm in MPPCA
void intermfa(const Matrix2D<double> &X, MultidimArray<double> &LX, Matrix2D<double> &MX, Matrix1D<double> &PX, Matrix2D<double> &R,
		      MultidimArray<double> &Z, int no_analyzers, int no_dims)
{

	Matrix1D<double> pii,AUXR1,DCC,AUXR,LOGpii,MAXR,RS;
	Matrix2D<double> l,ltpi,ltpil,ID,iltpil,CC,CCI,COVZ, MEANZ,DELTA,AUX3,AUXR2,DELTA2,A,AUXMX,RMAXR,ER,RRS;

	A=X.transpose();

	pii=1/PX;
	//Initialize some parameters
	int D=MAT_YSIZE(A);		//row size
    int N=MAT_XSIZE(A);   	//column size
    double cnst=-D/2*log(2*PI);
    LX.resize(no_analyzers,D,no_dims);
    R.initZeros(no_analyzers,N);
    Z.initZeros(no_analyzers,no_dims,N);
    ltpi.resizeNoCopy(no_dims,D);
    l.resizeNoCopy(D,no_dims);
    Matrix2D<double>repl;
    repl.resizeNoCopy(VEC_XSIZE(pii),no_dims);
	AUXR2.resizeNoCopy(no_dims,N);

    //Estimate the cluster centers based on input (= E-step)
    for(int k=0; k<no_analyzers; ++k){
    	FOR_ALL_ELEMENTS_IN_MATRIX2D(l)
    	    MAT_ELEM(l,i,j)=A3D_ELEM(LX, k, i, j);
    	FOR_ALL_ELEMENTS_IN_MATRIX2D(repl)
    		MAT_ELEM(repl,i,j)=VEC_ELEM(pii,i);
    	FOR_ALL_ELEMENTS_IN_MATRIX2D(ltpi)
    		MAT_ELEM(ltpi,i,j)=MAT_ELEM(repl,i,j)*MAT_ELEM(l,i,j);
    	ltpil=ltpi * l;
    	ID.initIdentity(no_dims);
    	iltpil=ID+ltpil;
    	cholesky(iltpil,CC);
    	CC.inv(CCI);
    	COVZ=CCI*CCI.transpose();
    	AUXMX.resizeNoCopy(D,N);
    	FOR_ALL_ELEMENTS_IN_MATRIX2D(AUXMX){
    		for (int n=k; n<N; n++)
    			MAT_ELEM(AUXMX,i,j)=MAT_ELEM(MX,i,n);
    	}
    	DELTA= A-AUXMX;
    	MEANZ=((ID-ltpil*COVZ)*ltpi)*DELTA;
    	FOR_ALL_ELEMENTS_IN_MATRIX2D(MEANZ)
    	    A3D_ELEM(Z, k, i, j) = MAT_ELEM(MEANZ,i,j);
    	AUX3=iltpil*MEANZ;
    	FOR_ALL_ELEMENTS_IN_MATRIX2D(MEANZ)
      		MAT_ELEM(AUXR2,i,j)=MAT_ELEM(MEANZ,i,j)*MAT_ELEM(AUX3,i,j);
    	AUXR2.colSum(AUXR1);
    	CC.getDiagonal(DCC);
    	matrixPow(DELTA,DELTA2,2);
    	AUXR=-0.5*(pii.transpose()*DELTA2-AUXR1)-log(DCC.sum());
    	FOR_ALL_ELEMENTS_IN_MATRIX1D(AUXR)
    	    MAT_ELEM(R,k,i)=VEC_ELEM(AUXR,i);
        }

    //Compute responsibilities of clusters to points
    LOGpii.resizeNoCopy(VEC_XSIZE(pii));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(LOGpii)
    	VEC_ELEM(LOGpii,i)=log(VEC_ELEM(pii,i));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
       MAT_ELEM(R,i,j)=MAT_ELEM(R,i,j)+cnst+.5*VEC_ELEM(LOGpii,i);
    RMAXR.resizeNoCopy(MAT_YSIZE(R),MAT_XSIZE(R));
    ER.resizeNoCopy(MAT_YSIZE(R),MAT_XSIZE(R));
    columnMax(R,MAXR);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(RMAXR)
    	MAT_ELEM(RMAXR,i,j)=VEC_ELEM(MAXR,j);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(ER)
    	MAT_ELEM(ER,i,j)= MAT_ELEM(R,i,j)-MAT_ELEM(RMAXR,i,j);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
    	MAT_ELEM(R,i,j)=exp(MAT_ELEM(ER,i,j));
    RRS.resizeNoCopy(MAT_YSIZE(R),MAT_XSIZE(R));
    R.colSum(RS);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(RRS)
    	MAT_ELEM(RRS,i,j)=VEC_ELEM(RS,j);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(R)
    	MAT_ELEM(R,i,j)=MAT_ELEM(R,i,j)/MAT_ELEM(RRS,i,j);
}

void ChartingManifold::setSpecificParameters(int max_iterations, int no_analyzers)
{
    this->max_iterations=max_iterations;
    this->no_analyzers=no_analyzers;
}

// m=Z(:,i,j)*Z(:,i,j)^t
void multiplyZijZijt(const MultidimArray<double>& Z, size_t i, size_t j, Matrix2D<double> &m)
{
    m.resizeNoCopy(ZSIZE(Z),ZSIZE(Z));
    for (size_t k1=0; k1<ZSIZE(Z); ++k1)
        for (size_t k2=0; k2<ZSIZE(Z); ++k2)
            MAT_ELEM(m,k1,k2)=DIRECT_A3D_ELEM(Z,k1,i,j)*DIRECT_A3D_ELEM(Z,k2,i,j);
}

void ChartingManifold::reduceDimensionality()
{
    Matrix1D<double> XrowMeans, PX,AUXU, lambdaD;
    Matrix2D<double> MX,D,Ds,R,TRA,U, Ut,lambda, V, Vv, mappedX;
    MultidimArray<double> LX,Z,ZP,Znew,R3,AAA;
    //Initialize some parameters
    double tol=1e-10;
    double min_std=1e-3;
    double kf=no_analyzers*(outputDim+1);

    claraWrite(*X,"dimred/X.txt");

    //Construct MFA model on data
    mppca(*X,outputDim,no_analyzers,tol,max_iterations,min_std,LX,MX,PX);
    intermfa(*X,LX,MX,PX,R,Z,no_analyzers, outputDim);

    //Adds last entry = 1 in posterior mean to handle means of factor analyzers
    Znew.resizeNoCopy(ZSIZE(Z),YSIZE(Z)+1,XSIZE(Z));
    FOR_ALL_ELEMENTS_IN_ARRAY3D(Z)
    	A3D_ELEM(Znew,k,i,j)=A3D_ELEM(Z,k,i,j);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(Znew)
    	A3D_ELEM(Znew, k,outputDim, j) = 1;
    Znew.write("dimred/Zn.txt");
    ZP.resizeNoCopy(XSIZE(Znew),YSIZE(Znew),ZSIZE(Znew));
    FOR_ALL_ELEMENTS_IN_ARRAY3D(ZP)
		A3D_ELEM(ZP, k,i,j) =A3D_ELEM(Znew,j,i,k);

    // Construct blockdiagonal matrix D
    D.initZeros((outputDim+1)*no_analyzers,(outputDim+1)*no_analyzers);
    for (int i=0;i<no_analyzers;++i)
    {
        Ds.initZeros(outputDim+1,outputDim+1);
        for (size_t j=1; j<MAT_XSIZE(*X); ++j)
        {
            multiplyZijZijt(Z,i,j,TRA);
            double Rij=MAT_ELEM(R,i,j);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(Ds)
            MAT_ELEM(Ds,i,j)+=Rij*MAT_ELEM(TRA,i,j);
        }
    }

    // Construct responsibility weighted local representation matrix U
    R3.setNdim(MAT_YSIZE(*X));
    R3.setYdim(no_analyzers);
    R3.setXdim(1);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(R3)
    	A3D_ELEM(R3,k,i,j)=MAT_ELEM(R,i,j);
    U.resizeNoCopy(kf,MAT_YSIZE(*X));
    AUXU.resizeNoCopy(kf*MAT_YSIZE(*X));
    arrayMultiply3D(Z,R3,AAA);
    //NOTA: Desde aqui tarda mucho tiempo
    for (int n=0; n<VEC_XSIZE(AUXU); n++){
    	for (int k=0; k<NSIZE(AAA); ++k){
    		for(int j=0; j<XSIZE(AAA); ++j){
    			for (int i=0; i<YSIZE(AAA); ++i)
    				VEC_ELEM(AUXU,n)=A3D_ELEM(AAA,k,i,j);
    		}
    	}
    }
    for (int n=0; n<VEC_XSIZE(AUXU); n++){
    	for (int j=0; j<MAT_YSIZE(*X); j++){
    		for (int i=0; i<kf; i++)
    			MAT_ELEM(U,i,j)=VEC_ELEM(AUXU,n);
        }
    }

    Ut=U.transpose();
    /*Problema: NECESITO QUE V SEA 2D y en matrix2d.cpp es 1D!!!
    generalizedEigs(D-(Ut*U),Ut*U,V,lambda);
*/
    lambda.getDiagonal(lambdaD);
    Vv.resizeNoCopy(MAT_YSIZE(V),VEC_XSIZE(lambdaD)-2);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(V){
    	for (size_t k=0; k<VEC_XSIZE(lambdaD)-2; k++)
    		MAT_ELEM(V,i,j)=MAT_ELEM(V,i,k);
    }
    mappedX=U*V;
}
