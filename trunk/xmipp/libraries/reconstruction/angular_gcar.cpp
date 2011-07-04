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

#include "angular_gcar.h"

// Read arguments ==========================================================
void ProgAngularGCAR::readParams()
{
}

// Show ====================================================================
void ProgAngularGCAR::show()
{
    if (!verbose)
        return;
}

// usage ===================================================================
void ProgAngularGCAR::defineParams()
{
	addParamsLine("-i <file>   : Input file");
}

// Produce side information ================================================
void ProgAngularGCAR::run()
{
    Matrix2D<double> q;

    //Matrix2D<double> rotM;
    //Matrix1D<double> qK1;
	//Matrix1D<double> qK2;

    //int idx1;
    //int idx2;

    //Matrix2D<int> clmatrix;
	//Matrix2D<int> clcorr;

    Matrix2D<double> PHI;
	//Matrix2D<double> PHI2;

    //Matrix2D<double> PHIRef;

    int K = 100;
    int L = 360;

    //qrand(q,K);
    //std::ofstream fout("/home/xmippuser/Desktop/salida_qrand.txt");
    //fout << q;
    //fout.close();

    //q.getCol(0,qK1);
	//q.getCol(0,qK2);
	//std::cout << "qk1 :" << qK1 << std::endl;
	//qToRot(qK1,rotM);
	//std::cout << rotM << std::endl;

    //prm.commonLineQ(q,0,1,L,idx1,idx2);
    //std::cout << "idx1 = " << idx1 << ", idx2 = " << idx2 << std::endl;

    //clmatrixCheatQ(q,L,clmatrix,clcorr);
    //fout.open("/home/xmippuser/Desktop/salida_clmatrixCheatQ.txt");
    //fout << clmatrix;
    //fout.close();
    //PHI.initZeros(K*L,3);
    //PHI.read("/home/xmippuser/Desktop/PHI.txt");
    //std::cout << PHI;
    //cryoCosmetify(PHI,PHI2,K,L);
    //std::ofstream fout("/home/xmippuser/Desktop/salida_PHI2.txt");
    //fout << PHI2;
	//fout.close();

    int nEigs = 20;
    SparseMatrix2D W;
    W.loadMatrix("/home/xmippuser/Desktop/cosa.txt");
    cryoOrientationsSdp(W,nEigs,PHI);

	//q.initZeros(4,100);
	//q.read("/home/xmippuser/Desktop/Q_Q2S2.txt");
    //Q2S2(q,PHIRef,L);
    //fout.open("/home/xmippuser/Desktop/salida_PHIRef.txt");
    //fout<<PHIRef;
    //fout.close();
}

void ProgAngularGCAR::qrand(Matrix2D<double>& q, int K)
{
	Matrix1D<double> l2Norm;
	l2Norm.initZeros(K);
	q.initZeros(4,K);

	//generar numeros aleatorios segun una normal
	//FOR_ALL_ELLEMENTS_IN_MATRIX2d(q)
	//	q(i,j) = rnd_gaus(double a, double b);
	q.read("/home/xmippuser/Desktop/matrix_q.txt");

	FOR_ALL_ELEMENTS_IN_MATRIX1D(l2Norm)
	{
		l2Norm(i) = sqrt(pow(q(0,i),2) + pow(q(1,i),2) + pow(q(2,i),2) + pow(q(3,i),2));
		q(0,i) = q(0,i)/l2Norm(i);
		q(1,i) = q(1,i)/l2Norm(i);
		q(2,i) = q(2,i)/l2Norm(i);
		q(3,i) = q(3,i)/l2Norm(i);
	}

	for(int k=0; k<K; k++)
	{
		if(q(0,k)<0)
		{
			q(0,k) = -q(0,k);
			q(1,k) = -q(1,k);
			q(2,k) = -q(2,k);
			q(3,k) = -q(3,k);
		}
	}
}

/// qToRot
void ProgAngularGCAR::qToRot(const Matrix1D<double>& q, Matrix2D<double> &rotMatrix)
{
	rotMatrix.resizeNoCopy(3,3);

	double q0 = VEC_ELEM(q,0);
	double cuaQ0 = pow(q0,2);
	double q1 = VEC_ELEM(q,1);
	double cuaQ1 = pow(q1,2);
	double q2 = VEC_ELEM(q,2);
	double cuaQ2 = pow(q2,2);
	double q3 = VEC_ELEM(q,3);
	double cuaQ3 = pow(q3,2);

	MAT_ELEM(rotMatrix,0,0) = cuaQ0 + cuaQ1 - cuaQ2 - cuaQ3;
	MAT_ELEM(rotMatrix,0,1) = 2*q1*q2 - 2*q0*q3;
	MAT_ELEM(rotMatrix,0,2) = 2*q0*q2 + 2*q1*q3;

	MAT_ELEM(rotMatrix,1,0) = 2*q1*q2 + 2*q0*q3;
	MAT_ELEM(rotMatrix,1,1) = cuaQ0 - cuaQ1 + cuaQ2 - cuaQ3;
	MAT_ELEM(rotMatrix,1,2) = -2*q0*q1 + 2*q2*q3;

	MAT_ELEM(rotMatrix,2,0) = -2*q0*q2 + 2*q1*q3;
	MAT_ELEM(rotMatrix,2,1) = 2*q0*q1 + 2*q2*q3;
	MAT_ELEM(rotMatrix,2,2) = cuaQ0 - cuaQ1 - cuaQ2 + cuaQ3;
}

/// commonLineQ
void ProgAngularGCAR::commonLineQ(const Matrix2D<double>& q, int k1, int k2, int nTheta, int& idx1, int& idx2)
{
	Matrix2D<double> R1;
	Matrix2D<double> R2;
	Matrix2D<double> invR1;
	Matrix2D<double> invR2;
	Matrix1D<double> qK1;
	Matrix1D<double> qK2;
	Matrix1D<double> X1;
	Matrix1D<double> Y1;
	Matrix1D<double> Z1;
	Matrix1D<double> X2;
	Matrix1D<double> Y2;
	Matrix1D<double> Z2;
	Matrix1D<double> Z3;
	Matrix2D<double> XY1;
	Matrix2D<double> XY2;
	Matrix1D<double> c1;
	Matrix1D<double> c2;
	Matrix1D<double> ev1;
	Matrix1D<double> ev2;
	double idx1d, idx2d;

	q.getCol(k1,qK1);
	q.getCol(k2,qK2);

	qToRot(qK1,R1);
	qToRot(qK2,R2);

	R1.inv(invR1);
	R2.inv(invR2);

	//std::cout << "R1 " << R1 << std::endl;
	//std::cout << "R2 " << R2 << std::endl;

	invR1.getCol(0,X1);
	invR1.getCol(1,Y1);
	invR1.getCol(2,Z1);
	invR2.getCol(0,X2);
	invR2.getCol(1,Y2);
	invR2.getCol(2,Z2);

	Z3.resizeNoCopy(3);
	VEC_ELEM(Z3,0) = VEC_ELEM(Z1,1)*VEC_ELEM(Z2,2) - VEC_ELEM(Z1,2)*VEC_ELEM(Z2,1);
	VEC_ELEM(Z3,1) = VEC_ELEM(Z1,2)*VEC_ELEM(Z2,0) - VEC_ELEM(Z1,0)*VEC_ELEM(Z2,2);
	VEC_ELEM(Z3,2) = VEC_ELEM(Z1,0)*VEC_ELEM(Z2,1) - VEC_ELEM(Z1,1)*VEC_ELEM(Z2,0);

	//std::cout << "X1 " << X1 << std::endl;
	//std::cout << "Y1 " << Y1 << std::endl;
	//std::cout << "Z1 " << Z1 << std::endl;
	//std::cout << "X2 " << X2 << std::endl;
	//std::cout << "Y2 " << Y2 << std::endl;
	//std::cout << "Z2 " << Z2 << std::endl;
	//std::cout << "Z3" << Z3 << std::endl;
	//getchar();

	double normaZ3 = Z3.module();
	if(normaZ3<1.0e-8)
		std::cerr << "GCAR:normTooSmall -> Images have same orientation" << std::endl;

	Z3 = Z3/normaZ3;

	XY1.resizeNoCopy(3,2);
	XY2.resizeNoCopy(3,2);
	XY1.setCol(0,X1);
	XY1.setCol(1,Y1);
	XY2.setCol(0,X2);
	XY2.setCol(1,Y2);

	c1 = Z3.transpose()*XY1;
	c1.setCol();
	c2 = Z3.transpose()*XY2;
	c2.setCol();

	//std::cout << "c1" << c1 << std::endl;
	//std::cout << "c2" << c2 << std::endl;
	//std::cout << "Z3" << Z3 << std::endl;
	//std::cout << "xy1" << XY1 << std::endl;
	//std::cout << "xy2" << XY2 << std::endl;

	ev1 = XY1*c1-Z3.transpose();
	ev2 = XY2*c2-Z3.transpose();

	//std::cout << "ev1" << XY1*c1 << std::endl;
	//std::cout << "ev2" << XY1*c1 << std::endl;

	if(ev1.module()/normaZ3>1.0e-12 || ev2.module()/normaZ3>1.0e-12)
		std::cerr << "GCAR:largeErrors -> Common line is not common. Error1 = "
				  << ev1.module()/normaZ3 << ", Error2 = " << ev2.module()/normaZ3 << std::endl;

	double theta1 = atan2(VEC_ELEM(c1,1),VEC_ELEM(c1,0)) + PI;
	double theta2 = atan2(VEC_ELEM(c2,1),VEC_ELEM(c2,0)) + PI;
	//std::cout << "theta1 " << theta1 << std::endl;
	//std::cout << "theta2 " << theta2 << std::endl;
	idx1d = theta1/(2*PI)*nTheta;
	idx2d = theta2/(2*PI)*nTheta;
	//std::cout << "idx1 " << idx1d << std::endl;
	//std::cout << "idx2 " << idx2d << std::endl;
	idx1 = ((int) round(idx1d))%nTheta;
	idx2 = ((int) round(idx2d))%nTheta;

}

void ProgAngularGCAR::clmatrixCheatQ(const Matrix2D<double>& q, int nTheta, Matrix2D<int>& clmatrix, Matrix2D<int>&  clcorr)
{
	int N;
	int idx1;
	int idx2;

	N = q.Xdim();

	clmatrix.initZeros(N,N);
	clcorr.initZeros(N,N);

	for(int k1=0; k1<N-1; k1++)
		for(int k2=k1+1; k2<N; k2++)
		{
			commonLineQ(q,k1,k2,nTheta,idx1,idx2);
			MAT_ELEM(clmatrix,k1,k2) = idx1+1;
			MAT_ELEM(clmatrix,k2,k1) = idx2+1;
			MAT_ELEM(clcorr,k1,k2) = 1.0e-8;
		}

}

void ProgAngularGCAR::cryoCosmetify(const Matrix2D<double>& PHI, Matrix2D<double>& PHI2, int K, int L)
{
	Matrix2D<double> circ;
	Matrix1D<double> equiAngles;
	Matrix1D<double> col;
	Matrix1D<double> col2;
	Matrix1D<double> col3;
	Matrix2D<double> U;
	Matrix2D<double> V;
	Matrix1D<double> D;
	Matrix1D<double> diag;
	Matrix1D<double> cosTheta;
	Matrix1D<double> sinTheta;
	Matrix1D<double> theta;
	Matrix2D<double> PHICopia;
	double avgCosTheta;
	double avgSinTheta;
	double avgTheta;
	double nr;
	int nPositive;
	int orientation;
	int minPos, maxPos;

	circ.initZeros(L,3);

	equiAngles.initZeros(L);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(equiAngles)
	{
		VEC_ELEM(equiAngles,i) = 2*PI*i/L;
	}
	equiAngles.setCol();

	PHI2.resizeNoCopy(K*L,3);
	PHI2.initZeros(K*L,3);


	for(int k=1; k<=K; k++)
	{
		PHICopia = PHI;
		//std::cout << "dimx " << PHI.Xdim() << "dimy" << PHI.Ydim() << std::endl;
		//std::cout << "dimx " << PHICopia.Xdim() << "dimy" << PHICopia.Ydim() << std::endl;
		PHICopia.submatrix((k-1)*L,0,(k*L)-1,PHI.Xdim()-1);

		PHICopia.getCol(0,col);
		circ.setCol(0,col);


		PHICopia.getCol(1,col);
		circ.setCol(1,col);

		PHICopia.getCol(2,col);
		circ.setCol(2,col);

		//std::cout << "circ " << circ << std::endl;

		svdcmp(circ,U,D,V);
		//std::cout<< "U "<<U<<std::endl;
		//std::cout<< "D "<<D<<std::endl;
		//std::cout<< "V "<<V<<std::endl;

		D.minIndex(minPos);
		D.maxIndex(maxPos);

		V.getCol(minPos,col);
		V.getCol(3-minPos-maxPos,col2);
		V.getCol(maxPos,col3);

		V.setCol(0,col3);
		V.setCol(1,col2);
		V.setCol(2,col);

		//std::cout<< "V "<<V<<std::endl;

		V.getCol(0,col);
		cosTheta = circ*col;
		V.getCol(1,col);
		sinTheta = circ*col;
		//std::cout<<"sin"<<sinTheta<<std::endl;
		//std::cout<<"cos"<<cosTheta<<std::endl;

		theta.initZeros(cosTheta.vdim);
		FOR_ALL_ELEMENTS_IN_MATRIX1D(cosTheta)
		{
			nr = sqrt(pow(VEC_ELEM(cosTheta,i),2)+pow(VEC_ELEM(sinTheta,i),2));
			//if(i==0) std::cout<<"nr "<<nr<<std::endl;
			VEC_ELEM(cosTheta,i) /= nr;
			VEC_ELEM(sinTheta,i) /= nr;
			VEC_ELEM(theta,i) = atan2(VEC_ELEM(sinTheta,i),VEC_ELEM(cosTheta,i));
		}

		nPositive = 0;
		for(int i=1; i<theta.vdim; i++)
			if((VEC_ELEM(theta,i)-VEC_ELEM(theta,i-1))>0)
				nPositive++;

		orientation = (nPositive>L/2)? 1:-1;

		avgCosTheta = 0;
		avgSinTheta = 0;
		FOR_ALL_ELEMENTS_IN_MATRIX1D(theta)
		{
			avgCosTheta += cos(VEC_ELEM(theta,i)-orientation*VEC_ELEM(equiAngles,i));
			avgSinTheta += sin(VEC_ELEM(theta,i)-orientation*VEC_ELEM(equiAngles,i));
		}
		avgCosTheta /= theta.vdim;
		avgSinTheta /= theta.vdim;
		avgTheta = atan2(avgSinTheta, avgCosTheta);
		theta = orientation*equiAngles + avgTheta;

		V.getCol(0,col);
		V.getCol(1,col2);
		col.setRow();
		col2.setRow();

		for(int i=0; i<L; i++)
			PHI2.setRow(((k-1)*L)+i,(cos(VEC_ELEM(theta,i))*col)+(sin(VEC_ELEM(theta,i))*col2));
		//std::cout<<PHI2<<std::endl;
	}
	//std::ofstream fout("/home/xmippuser/Desktop/salida_PHI2.txt");
	//fout << PHI2;
	//fout.close();

}


void ProgAngularGCAR::Q2S2(const Matrix2D<double>& Q, Matrix2D<double>& PR, int NTheta)
{
	Matrix2D<double> R;
	Matrix1D<double> col;
	Matrix1D<double> col2;
	int n;
	double dTheta;

	n = Q.Xdim();
	dTheta = 2*PI/NTheta;

	PR.resizeNoCopy(NTheta*n,3);
	PR.initZeros(NTheta*n,3);

	for(int k=0;k<10;k++)
	{
		Q.getCol(k,col);
		qToRot(col,R);

		R.getRow(0,col);
		R.getRow(1,col2);
		col.setRow();
		col2.setRow();

		std::cout<<"k"<<k<<"e1"<<col<<std::endl;
		std::cout<<"k"<<k<<"e2"<<col2<<std::endl;

		for(int j=0; j<NTheta; j++)
			PR.setRow(k*NTheta+j,(cos(j*dTheta)*col+sin(j*dTheta)*col2));
	}
	//std::cout<<PR<<std::endl;
	std::cout<<"dtheta"<<dTheta<<std::endl;
}

/*void ProgAngularGCAR::cryoSdpunmix(const Matrix2D<double>& PHI, Matrix2D<double>& PHI2)
{
	Matrix2D<double> Q;
	Matrix2D<double> Qb;
	Matrix2D<double> U;
	Matrix2D<double> V;
	Matrix2D<double> D;
	Matrix1D<double> Diag;
	int n;
	int q;
	int q2;
	int idx;
	int m;

	n = PHI.Ydim();
	q = PHI.Xdim();
	q2 = q*q;

	Q.initZeros(n,q2);
	for(j=0;j<n;j++)
	{
		idx = 1;
		for(k=0;k<q;k++)
			for(l=0;l<q;l++)
			{
				MAT_ELEM(Q,j,idx) = MAT_ELEM(PHI,j,k)*MAT_ELEM(PHI,j,l);
				idx++;
			}
	}

	svdcmp(Q.transpose()*Q,U,Diag,V);

	D.initZeros(Diag.vdim,Diag.vdim);
	for(int i=0; i<Diag.vdim; i++)
		MAT_ELEM(D,i,i) = sqrt(VEC_ELEM(Diag,i));

	Qb = U*D*V.transpose();

	m = q2*(q2+1)/2+q2+1;
}*/

void ProgAngularGCAR::registerOrientations(Matrix2D<double>& PHI,const Matrix2D<double>& PHIRef)
{
	Matrix2D<double> T;
	Matrix2D<double> U;
	Matrix2D<double> V;
	Matrix1D<double> D;
	Matrix2D<double> Q;

	T = PHI.inv()*PHIRef;
	svdcmp(T,U,D,V);
	Q = U * V.transpose();

	PHI = PHI*Q;
}

void ProgAngularGCAR::checkOrientations(const Matrix2D<double>& PHI,const Matrix2D<double>& PHIRef)
{
	Matrix2D<double> PHICopia;
	Matrix1D<double> cosAngle;
	Matrix1D<double> err;
	Matrix1D<double> errInAngles;

	PHICopia = PHI;
	registerOrientations(PHICopia,PHIRef);

	cosAngle.initZeros(PHICopia.Ydim());
	err.initZeros(PHICopia.Ydim());
	errInAngles.initZeros(PHICopia.Ydim());

	FOR_ALL_ELEMENTS_IN_MATRIX1D(cosAngle)
	{
		VEC_ELEM(cosAngle,i) = MAT_ELEM(PHICopia,i,1)*MAT_ELEM(PHIRef,i,1) +
							   MAT_ELEM(PHICopia,i,2)*MAT_ELEM(PHIRef,i,2) +
							   MAT_ELEM(PHICopia,i,3)*MAT_ELEM(PHIRef,i,3);
		VEC_ELEM(err,i) = acos(VEC_ELEM(cosAngle,i));
		VEC_ELEM(errInAngles, i) = VEC_ELEM(err,i)*180/PI;
	}
	//Pintar cositas
}

void ProgAngularGCAR::cryoOrientationsSdp(SparseMatrix2D& W, int nEigs, Matrix2D<double> &PHI)
{
	std::vector<EigElement> vec(nEigs);
	Matrix1D<double> vec1;
	Matrix1D<double> vec2;
	Matrix1D<double> vec3;
	Matrix2D<double> VECS;
	int nComplex;
	int complexIndex1;
	int complexIndex2;
	int realIndex;
	double normVec2, normVec3;

	int N=W.nrows();
	ARNonSymStdEig<double, SparseMatrix2D > dpro(N, nEigs, &W, &SparseMatrix2D::multMv);
	dpro.FindEigenvectors();

	PHI.resizeNoCopy(N,3);
	vec1.initZeros(N);
	vec2.initZeros(N);
	vec3.initZeros(N);

	for(int i=0; i<nEigs; i++)
	{
		vec[i].eigenvalue = sqrt(dpro.EigenvalueReal(i)*dpro.EigenvalueReal(i)+dpro.EigenvalueImag(i)*dpro.EigenvalueImag(i));
		vec[i].pos = i;
	}
	sort(vec.begin(),vec.end());

	nComplex = 0;
	complexIndex1 = 0;
	complexIndex2 = 0;
	realIndex = 0;
	for(int j=1; j<4; j++)
	{
		if(abs(dpro.EigenvalueImag(vec[j].pos)) > 2.2204e-06)
		{
			nComplex++;
			if(nComplex==1) complexIndex1 = j;
			else complexIndex2 == j;
		}
		else realIndex = j;
	}

	if(nComplex==2)
	{
		FOR_ALL_ELEMENTS_IN_MATRIX1D(vec1)
		{
			VEC_ELEM(vec1,i) = dpro.RawEigenvector(vec[realIndex].pos)[i];
			VEC_ELEM(vec2,i);//PARTE REAL = dpro.RawEigenvector(vec[complexIndex1].pos)[i];
			VEC_ELEM(vec3,i);//PARTE IMAGINARIA = dpro.RawEigenvector(vec[].pos)[i];
			//CALCULAR LAS NORMAS
		}
	}
	else
		FOR_ALL_ELEMENTS_IN_MATRIX1D(vec1)
		{
			VEC_ELEM(vec1,i) = dpro.RawEigenvector(vec[0].pos)[i];
			VEC_ELEM(vec2,i) = dpro.RawEigenvector(vec[1].pos)[i];
			VEC_ELEM(vec3,i) = dpro.RawEigenvector(vec[2].pos)[i];
		}
	VECS.setCol(0,vec1);
	VECS.setCol(1,vec2);
	VECS.setCol(2,vec3);
	std::cout<<VECS<<std::endl;
}
