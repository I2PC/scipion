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
	/*
	std::vector<SparseElement> vec(4);
	vec.at(0).value = 2;
	vec.at(0).i = 0;
	vec.at(0).j = 3;
	vec.at(1).value = -2;
	vec.at(1).i = 3;
	vec.at(1).j = 0;
	vec.at(2).value = 1;
	vec.at(2).i = 1;
	vec.at(2).j = 2;
	vec.at(3).value = -1;
	vec.at(3).i = 2;
	vec.at(3).j = 1;

	SparseMatrix2D W(vec,4);

	ARNonSymStdEig<double, SparseMatrix2D > dpro(4, 2, &W, &SparseMatrix2D::multMv,"LM");
	dpro.FindEigenvectors();

	std::cout<<"autoval 0:["<<dpro.EigenvalueReal(0)<<"+"<<dpro.EigenvalueImag(0)<<", "<<dpro.EigenvalueReal(1)<<"+"<<dpro.EigenvalueImag(1)<<"]\n";

	std::cout<<"autovec 0: ["<<dpro.EigenvectorReal(0,0)<<"+"<<dpro.EigenvectorImag(0,0)<<+"i,  "<<dpro.EigenvectorReal(0,1)<<"+"<<dpro.EigenvectorImag(0,1)<<+"i,  "
				                 <<dpro.EigenvectorReal(0,2)<<"+"<<dpro.EigenvectorImag(0,2)<<+"i,  "<<dpro.EigenvectorReal(0,3)<<"+"<<dpro.EigenvectorImag(0,3)<<+"i]\n";

	std::cout<<"autovec 0: ["<<dpro.RawEigenvector(0)[0]<<",  "<<dpro.RawEigenvector(0)[1]<<",  "<<dpro.RawEigenvector(0)[2]<<",  "<<dpro.RawEigenvector(0)[3]<<"]\n";
	std::cout<<"autovec 0: ["<<dpro.RawEigenvector(1)[0]<<",  "<<dpro.RawEigenvector(1)[1]<<",  "<<dpro.RawEigenvector(1)[2]<<",  "<<dpro.RawEigenvector(1)[3]<<"]\n";//*/


    Matrix2D<double> q;
    int K = 100;
	int L = 360;

    qrand(q,K);

    Matrix2D<int> clmatrix;
    Matrix2D<int> clcorr;

    clmatrixCheatQ(q,L,clmatrix,clcorr);

    int d = 10;

    SparseMatrix2D W;
    cryo_S2graph(clmatrix,L,d,W);

    int nEigs = 20;
    Matrix2D<double> PHI;
	Matrix2D<double> PHI2;
    Matrix2D<double> PHIRef;

    PHI.resizeNoCopy(K*L,3);

    cryoOrientations(W,nEigs,PHI);

    cryoCosmetify(PHI,PHI2,K,L);

    Q2S2(q,PHIRef,L);

    checkOrientations(PHI2,PHIRef);
    /*std::ofstream fout("/home/xmippuser/Desktop/salida_PHI.txt");
    fout << PHI;
    fout.close();//*/

	//q.initZeros(4,100);
	//q.read("/home/xmippuser/Desktop/Q_Q2S2.txt");
    //Q2S2(q,PHIRef,L);
    //fout.open("/home/xmippuser/Desktop/salida_PHIRef.txt");
    //fout<<PHIRef;
    //fout.close();
}

void ProgAngularGCAR::cryo_S2graph(Matrix2D<int> &clmatrix,int L,int d, SparseMatrix2D &w)
{
    //Description of clmatrix:
	// If k1 and k2 are two projections, then entry (k1,k2) contains the index
	//of their common line in projection k1, and entry (k2,k1) contains the
	//index of their common line in projection k2.

	// Total number of vertices: |V|=N=K*L
	int K=clmatrix.mdimy;
	int N = K*L;
	int daux=2*d+1;

	// The adjacency matrix W is sparse and have only N_entries nonzero entries
	// N_entries is the total number of links
	// first leg of the spider contributes KL(2d+1)
	// there are (K \choose 2) intersections
	// factor of 2 for the antipodal point => 2*(K choose 2)
	// factor of 2 for k1<k2, k2<k1 => 4*(K choose 2)
	// those legs contribute 4*(K choose 2)*(2d+1)
	// overall: 2KLd + KL + 2K(K-1)*(2d+1)
	int N_entries = 2*N*d + N + 2*K*(K-1)*(daux);
	std::vector<IJpair> IJ(N_entries);
	Matrix1D<int> l_qtr(daux,false);
	for (int i=0; i<l_qtr.vdim; i++)
		VEC_ELEM(l_qtr,i)=i-d;
	Matrix1D<int> add_idx(daux,false);
	for (int i=0; i<add_idx.vdim; i++)
	    VEC_ELEM(add_idx,i)=i+1;
	Matrix1D<int> ones_L4(daux,false);
	ones_L4.initConstant(1);
	Matrix1D<int> idx(daux,false);
	idx =  add_idx -daux;

	// the first two legs of the spider - same circle nodes
	int indice;
    for (int k1=1; k1<=K; k1++)
	{
		for (int l1=1;l1 <=L;l1++)
		{
			idx = idx + daux;
			for (int aux=0;aux < daux;aux++)
			{
				indice=VEC_ELEM(idx,aux)-1;
				IJ.at(indice).i=((k1-1)*L+l1)*VEC_ELEM(ones_L4,aux)-1;
				IJ.at(indice).j=(k1-1)*L+ matlab_mod((l1+VEC_ELEM(l_qtr,aux)+L-1), L) ;

			}


		}
	}
    int l1,l2;
    for (int k1 = 1;k1 <= K; k1++)
    {
    	for (int k2 = (k1+1); k2 <= K; k2++ )
    	{
   	        l1 = MAT_ELEM(clmatrix,k1-1,k2-1);
   	        l2 = MAT_ELEM(clmatrix,k2-1,k1-1);
   	        if ((l1 != 0) && (l2 != 0))
   	        {
   	        		// take d points each side
    	        	idx = idx +daux;
    	        	for (int aux=0;aux < daux;aux++)
    	        	{
    	        		indice=VEC_ELEM(idx,aux)-1;
    	        		IJ.at(indice).i=((k1-1)*L+l1)*VEC_ELEM(ones_L4,aux)-1;
    	        		IJ.at(indice).j=(k2-1)*L+matlab_mod((l2+VEC_ELEM(l_qtr,aux)+L-1), L);
    	        	}
    	        	// the antipodal point
    	        	idx = idx + daux;
    	        	for (int aux=0;aux < daux;aux++)
    	        	{
						indice=VEC_ELEM(idx,aux)-1;
						IJ.at(indice).i=((k1-1)*L+ matlab_mod(l1+L/2-1,L) + 1) *VEC_ELEM(ones_L4,aux)-1;
						IJ.at(indice).j=(k2-1)*L+  matlab_mod((l2+VEC_ELEM(l_qtr,aux)+L/2-1),L) ;
    	        	}
      	        }
	    	}
	    }
    // the remaining legs of the spider (k1 > k2)
    for (int k1=2; k1 <= K; k1++)
	    {
	    	for (int k2=1;k2<k1;k2++)
	    	{
	    		l1 = MAT_ELEM(clmatrix,k1-1,k2-1);
	    		l2 = MAT_ELEM(clmatrix,k2-1,k1-1);
	    		if ((l1!=0)&&(l2!=0))
	    		{
	    			idx = idx + daux;
	    			for (int aux=0;aux < daux;aux++)
	    			{
						indice=VEC_ELEM(idx,aux)-1;
						IJ.at(indice).i=((k1-1)*L+l1)*VEC_ELEM(ones_L4,aux)-1;
						IJ.at(indice).j=(k2-1)*L+matlab_mod(l2+VEC_ELEM(l_qtr,aux)+L-1, L) ;
	    			}
	    			idx = idx + daux;
					for (int aux=0;aux < daux;aux++)
					{
						indice=VEC_ELEM(idx,aux)-1;
						IJ.at(indice).i=((k1-1)*L+matlab_mod(l1+L/2-1,L)+1)*VEC_ELEM(ones_L4,aux)-1;
						IJ.at(indice).j=(k2-1)*L+matlab_mod(l2+VEC_ELEM(l_qtr,aux)+L/2-1,L) ;
					}
	    		}
	    	}
	    }

	    int indx=VEC_ELEM(idx,daux-1);
	    std::vector<SparseElement> elems;

	    countElems(elems,IJ,N_entries);

	    elems.resize(elems.size());
	    new (&w) SparseMatrix2D(elems,N);

	    Matrix1D<double> Daux(true);
	    sumRows(w,Daux);
	    MultidimArray<double> D2(N);

	    for (int i=0;i<N;i++)
	    {
	    	DIRECT_MULTIDIM_ELEM(D2,i)=1/VEC_ELEM(Daux,i);
	    }


        w.multMMDiagonal(D2,w);

}

int ProgAngularGCAR::matlab_mod(int a, int b)
{
	int result = a%b;
	if (result<0)
		result = a+b;
	return result;
}

void ProgAngularGCAR::sumRows(SparseMatrix2D &w,Matrix1D<double> &resul)
{

	//auxiliary variables
	double sum=0.0;
	int begin;
	int end;
	// dimension of the matrix
	int N = w.N;
	// resize the result vector
	resul.resizeNoCopy(N);
	// Sum each row
	for (int i=0;i<(N-1);i++)
	{
		begin=DIRECT_MULTIDIM_ELEM((w.iIdx),i)-1;
		end=DIRECT_MULTIDIM_ELEM((w.iIdx),i+1)-1;
	    sum=0.0;
		for (int j=begin;j<end;j++)
			sum += DIRECT_MULTIDIM_ELEM((w.values),j);
		VEC_ELEM(resul,i)=sum;


	}
	// sum the last row
	sum=0.0;
	for (int j=end;j<(w.values).xdim;j++)
		sum += DIRECT_MULTIDIM_ELEM((w.values),j);
	VEC_ELEM(resul,N-1)=sum;
}

void ProgAngularGCAR::countElems(std::vector<SparseElement> &elem,std::vector<IJpair> &IJ,int indx)
{
	sort(IJ.begin(),IJ.end());
    IJpair ij;
    ij.i=IJ.at(0).i;
    ij.j=IJ.at(0).j;
    int contador=1;
    SparseElement e;
    for(int i=1; i<indx; i++)
    {
	    while( ij== IJ.at(i)&&i<indx)
	    {
	    	++contador;
	    	++i;
	    }

	    e.i=ij.i;
	    e.j=ij.j;
	    e.value=contador;
	    elem.push_back(e);
	    ij.i=IJ.at(i).i;
	    ij.j=IJ.at(i).j;
	    contador=1;
    }
}

void ProgAngularGCAR::qrand(Matrix2D<double>& q, int K)
{
	Matrix1D<double> l2Norm;
	l2Norm.initZeros(K);
	q.initZeros(4,K);

	//generar numeros aleatorios segun una normal
	//FOR_ALL_ELLEMENTS_IN_MATRIX2d(q)
	//	q(i,j) = rnd_gaus(double a, double b);
	q.read("matrix_q.txt");

	FOR_ALL_ELEMENTS_IN_MATRIX1D(l2Norm)
	{
		l2Norm(i) = sqrt(pow(q(0,i),2) + pow(q(1,i),2) + pow(q(2,i),2) + pow(q(3,i),2));
		q(0,i)    = q(0,i)/l2Norm(i);
		q(1,i)    = q(1,i)/l2Norm(i);
		q(2,i)    = q(2,i)/l2Norm(i);
		q(3,i)    = q(3,i)/l2Norm(i);
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
	Matrix2D<double> R1, R2, invR1, invR2, XY1, XY2;
	Matrix1D<double> qK1, qK2, X1, Y1, Z1, X2, Y2, Z2, Z3, c1, c2, ev1, ev2;
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
	int N, idx1, idx2;

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
	Matrix2D<double> circ, U, V, PHICopia;
	Matrix1D<double> equiAngles, col, col2, col3, D, diag, cosTheta, sinTheta, theta;
	double avgCosTheta, avgSinTheta, avgTheta, nr;
	int nPositive, orientation, minPos, maxPos;

	circ.initZeros(L,3);

	equiAngles.initZeros(L);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(equiAngles)
	{
		VEC_ELEM(equiAngles,i) = 2*PI*i/L;
	}
	equiAngles.setCol();

	PHI2.resizeNoCopy(K*L,3);
	PHI2.initZeros(K*L,3);


	for(int k=0; k<K; k++)
	{
		PHICopia = PHI;
		//std::cout << "dimx " << PHI.Xdim() << "dimy" << PHI.Ydim() << std::endl;
		//std::cout << "dimx " << PHICopia.Xdim() << "dimy" << PHICopia.Ydim() << std::endl;
		PHICopia.submatrix(k*L,0,((k+1)*L)-1,PHI.Xdim()-1);

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
			PHI2.setRow((k*L)+i,(cos(VEC_ELEM(theta,i))*col)+(sin(VEC_ELEM(theta,i))*col2));
		//std::cout<<PHI2<<std::endl;
	}
	//std::ofstream fout("/home/xmippuser/Desktop/salida_PHI2.txt");
	//fout << PHI2;
	//fout.close();

}


void ProgAngularGCAR::Q2S2(const Matrix2D<double>& Q, Matrix2D<double>& PR, int NTheta)
{
	Matrix2D<double> R;
	Matrix1D<double> col, col2;
	int n;
	double dTheta;

	n = Q.Xdim();
	dTheta = 2*PI/NTheta;

	PR.resizeNoCopy(NTheta*n,3);
	PR.initZeros(NTheta*n,3);

	for(int k=0;k<n;k++)
	{
		Q.getCol(k,col);
		qToRot(col,R);

		R.getRow(0,col);
		R.getRow(1,col2);
		col.setRow();
		col2.setRow();

		//std::cout<<"k"<<k<<"e1"<<col<<std::endl;
		//std::cout<<"k"<<k<<"e2"<<col2<<std::endl;

		for(int j=0; j<NTheta; j++)
			PR.setRow(k*NTheta+j,(cos(j*dTheta)*col+sin(j*dTheta)*col2));
	}
	//std::cout<<PR<<std::endl;
	//std::cout<<"dtheta"<<dTheta<<std::endl;
}


void ProgAngularGCAR::registerOrientations(Matrix2D<double>& PHI,const Matrix2D<double>& PHIRef)
{
	Matrix2D<double> T, U, V, Q;
	Matrix1D<double> D;

	T = PHI.inv()*PHIRef;
	/*std::ofstream fout("/home/xmippuser/Desktop/salida_Tregister.txt");
				fout << T;
				fout.close();//*/
	svdcmp(T,U,D,V);
	Q = U * V.transpose();

	PHI = PHI*Q;
	/*fout.open("/home/xmippuser/Desktop/salida_QOrientations.txt");
			fout << Q;
			fout.close();//*/
}

void ProgAngularGCAR::checkOrientations(const Matrix2D<double>& PHI,const Matrix2D<double>& PHIRef)
{
	Matrix2D<double> PHICopia;
	Matrix1D<double> cosAngle, err, errInAngles;

	PHICopia = PHI;
	registerOrientations(PHICopia,PHIRef);
	/*std::ofstream fout("/home/xmippuser/Desktop/salida_PHI2_registerorientations.txt");
		fout << PHICopia;
		fout.close();//*/
	cosAngle.initZeros(PHICopia.Ydim());
	err.initZeros(PHICopia.Ydim());
	errInAngles.initZeros(PHICopia.Ydim());

	FOR_ALL_ELEMENTS_IN_MATRIX1D(cosAngle)
	{
		VEC_ELEM(cosAngle,i) = MAT_ELEM(PHICopia,i,0)*MAT_ELEM(PHIRef,i,0) +
							   MAT_ELEM(PHICopia,i,1)*MAT_ELEM(PHIRef,i,1) +
							   MAT_ELEM(PHICopia,i,2)*MAT_ELEM(PHIRef,i,2);
/*		if(i<10){
			std::cout << MAT_ELEM(PHICopia,i,0) << "*" << MAT_ELEM(PHIRef,i,0) << "+" <<
						 MAT_ELEM(PHICopia,i,1) << "*" << MAT_ELEM(PHIRef,i,1) << "+" <<
						 MAT_ELEM(PHICopia,i,2) << "*" << MAT_ELEM(PHIRef,i,2) << "\n";
			std::cout << MAT_ELEM(PHICopia,i,0)*MAT_ELEM(PHIRef,i,0) +
					   	 MAT_ELEM(PHICopia,i,1)*MAT_ELEM(PHIRef,i,1) +
					     MAT_ELEM(PHICopia,i,2)*MAT_ELEM(PHIRef,i,2) << "\n";
		}//*/
		if(VEC_ELEM(cosAngle,i)>1) VEC_ELEM(err,i) = 0.0;
		else 					   VEC_ELEM(err,i) = acos(VEC_ELEM(cosAngle,i));

		VEC_ELEM(errInAngles, i) = VEC_ELEM(err,i)*180/PI;
	}
	/*fout.open("/home/xmippuser/Desktop/salida_cosAngle.txt");
	fout << cosAngle;
	fout.close();
	fout.open("/home/xmippuser/Desktop/salida_err.txt");
	fout << err;
	fout.close();
	fout.open("/home/xmippuser/Desktop/salida_errInAngles.txt");
	fout << errInAngles;
	fout.close();
					//*/
	//Paint little things
	double media = errInAngles.sum(true);
	double std = 0;
	double max = 0;
	double min = MAXDOUBLE;
	double elemi;

	FOR_ALL_ELEMENTS_IN_MATRIX1D(errInAngles)
	{
		elemi = VEC_ELEM(errInAngles,i);
		std += (elemi-media)*(elemi-media);
		if(max<elemi) max = elemi;
		if(min>elemi) min = elemi;
	}
	std /= errInAngles.size();
	std = sqrt(std);

	std::cout<< "Mean " << media << std::endl;
	std::cout<< "STD "  << std   << std::endl;
	std::cout<< "Min "  << min   << std::endl;
	std::cout<< "Max "  << max   << std::endl;
}

void ProgAngularGCAR::cryoOrientations(SparseMatrix2D& W, int nEigs, Matrix2D<double> &PHI)
{
	std::vector<EigElement> vec(nEigs);
	Matrix1D<double> vec1, vec2, vec3, col1, col2, ATAVec;
	Matrix2D<double> equations, truncatedEquations, ATA, A, VECS;
	int nComplex, complexIndex1, complexIndex2, realIndex;
	double normVec2, normVec3;

	int N = W.nrows();
	ARNonSymStdEig<double, SparseMatrix2D > dpro(N, nEigs, &W, &SparseMatrix2D::multMv);
	dpro.FindEigenvectors();

	PHI.resizeNoCopy(N,3);
	VECS.resizeNoCopy(N,3);
	vec1.initZeros(N);
	vec2.initZeros(N);
	vec3.initZeros(N);

	for(int i=0; i < nEigs; i++)
	{
		vec[i].eigenvalue = sqrt(dpro.EigenvalueReal(i)*dpro.EigenvalueReal(i)+dpro.EigenvalueImag(i)*dpro.EigenvalueImag(i));
		vec[i].pos = i;
	}
	sort(vec.begin(),vec.end());
	//std::cout<<"primer autovalor "<<vec[1].eigenvalue<<"\nsegundo "<<vec[2].eigenvalue<<"\ntercero  "<<vec[3].eigenvalue<<"\n";
	//std::cout<<"primer autovalor "<<dpro.EigenvalueReal(vec[1].pos)<<"\nsegundo "<<dpro.EigenvalueReal(vec[2].pos)<<"\ntercero  "<<dpro.EigenvalueReal(vec[3].pos)<<"\n";

	nComplex      = 0;
	complexIndex1 = 0;
	complexIndex2 = 0;
	realIndex     = 0;
	for(int j=1; j<4; j++)
	{
		if(abs(dpro.EigenvalueImag(vec[j].pos)) > 2.2204e-06)
		{
			nComplex++;
			if(nComplex==1) complexIndex1 = j;
			else            complexIndex2 = j;
		}
		else realIndex = j;
	}

	if(nComplex==2)
	{
		FOR_ALL_ELEMENTS_IN_MATRIX1D(vec1)
		{
			VEC_ELEM(vec1,i) = dpro.EigenvectorReal(realIndex,i);
			VEC_ELEM(vec2,i) = dpro.EigenvectorReal(complexIndex1,i);
			VEC_ELEM(vec3,i) = dpro.EigenvectorImag(complexIndex1,i);
		}
	}
	else
	{
		FOR_ALL_ELEMENTS_IN_MATRIX1D(vec1)
		{


			VEC_ELEM(vec1,i) = dpro.EigenvectorReal(vec[1].pos,i);
			VEC_ELEM(vec2,i) = dpro.EigenvectorReal(vec[2].pos,i);
			VEC_ELEM(vec3,i) = dpro.EigenvectorReal(vec[3].pos,i);
		}
	}
	VECS.setCol(0,vec1);
	VECS.setCol(1,vec2);
	VECS.setCol(2,vec3);
	/*std::ofstream fout("/home/xmippuser/Desktop/salida_VECS.txt");
	fout << VECS;
	fout.close();//*/
	//std::cout<<VECS<<std::endl;
	equations.initZeros(W.ncols(),9);

	for(int k=0;k<3;k++)
	{
		for(int j=0;j<3;j++)
		{
			VECS.getCol(k,col1);
			VECS.getCol(j,col2);
			equations.setCol(3*k+j,col1*col2);
		}
	}
/*
	fout.open("/home/xmippuser/Desktop/salida_equations.txt");
	    fout << equations;
	    fout.close();
//*/
	truncatedEquations.initZeros(W.ncols(),6);
	equations.getCol(0,col1);
	truncatedEquations.setCol(0,col1);
	equations.getCol(1,col1);
	truncatedEquations.setCol(1,col1);
	equations.getCol(2,col1);
	truncatedEquations.setCol(2,col1);
	equations.getCol(4,col1);
	truncatedEquations.setCol(3,col1);
	equations.getCol(5,col1);
	truncatedEquations.setCol(4,col1);
	equations.getCol(8,col1);
	truncatedEquations.setCol(5,col1);
/*
	fout.open("/home/xmippuser/Desktop/salida_truncated.txt");
	fout << truncatedEquations;
	fout.close();
//*/
	col1.resizeNoCopy(W.ncols());
	col1.initConstant(0.5);

	ATAVec.resizeNoCopy(6);
	ATAVec.setCol();
	ATAVec = truncatedEquations.inv()*col1;
/*
	fout.open("/home/xmippuser/Desktop/salida_ATAVec.txt");
	fout << ATAVec;
	fout.close();
//*/
	//std::cout<<ATAVec<<"\n";

	ATA.resizeNoCopy(3,3);
	ATA.initZeros(3,3);
	int ATAIndex = 0;
	for(int k=0;k<3;k++)
	{
		for(int j=k;j<3;j++)
		{
			MAT_ELEM(ATA,k,j) = VEC_ELEM(ATAVec,ATAIndex);
			++ATAIndex;
		}
	}
	//std::cout<<ATA<<"\n";

	ATA = ATA + ATA.transpose();
/*
	std::cout<<ATA<<"\n";//*/
	//Cholesky
	cholesky(ATA,A);
	/*fout.open("/home/xmippuser/Desktop/salida_ATA.txt");
					fout << ATA;
					fout.close();
					//*/
	//std::cout<<A.transpose()<<"\n";
	PHI = A.transpose()*VECS.transpose();
	PHI = PHI.transpose();
	/*fout.open("/home/xmippuser/Desktop/salida_PHIPHI.txt");
						fout << PHI;
						fout.close();
						//*/
}

