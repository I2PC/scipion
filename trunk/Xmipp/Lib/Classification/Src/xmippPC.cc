/***************************************************************************
 *
 * Authors:     Jorge García de la Nava Ruiz (gdl@ac.uma.es)
 *              Carlos Oscar Sanchez Sorzano
 *
 * Departamento de Arquitectura de Computadores, Universidad de Málaga
 *
 * Copyright (c) 2001 , CSIC/UMA.
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'pascual@cnb.uam.es'
 *
 *****************************************************************************/


#include "../xmippPC.hh"

#include "../xmippDistances.hh"


/**
* Calculate the eigenval/vecs
* @param ts The vectors.
*/
void xmippPC::reset(xmippCTVectors const &ts)
{
	vector<unsigned> dummy;
	for(unsigned i=0;i<ts.size();i++)
		dummy.push_back(i);
	reset(ts,dummy);
}

/**
* Calculate the eigenval/vecs
* @param ts The vectors.
* @param idx The indexes of the vectors to use
*/
void xmippPC::reset(xmippCTVectors const &ts, vector<unsigned> const & idx)
 _THROW {
	vector<xmippVector> a;
	int n=ts.dimension();
	a.resize(n);

	{

    	int verbosity = listener->getVerbosity();
    	if (verbosity)
  		listener->OnReportOperation((string) "Diagonalizing matrix....\n");  
	
		//Get the mean of the given cluster of vectors
		xmippVector mean;
		for(int k=0;k<n;k++){
			a[k].resize(n);
			xmippFeature sum=0.0;
			int l=0;
			for(vector<unsigned>::const_iterator i=idx.begin();i!=idx.end();i++){
				if(finite(ts.itemAt(*i)[k])){
					sum+=ts.itemAt(*i)[k];
					l++;
				}
			}
			mean.push_back(sum/l);
		}

		for(int i=0;i<n;i++){
			for(int j=0;j<=i;j++){
				xmippFeature sum=0.0;
				int l=0;
				for(vector<unsigned>::const_iterator it=idx.begin();it!=idx.end();it++){
					xmippFeature d1=ts.itemAt(*it)[i]-mean[i];
					xmippFeature d2=ts.itemAt(*it)[j]-mean[j];
					if(finite(d1) && finite(d2)){
						sum+=d1*d2;
						l++;
					}
				}
				if(l) a[i][j]=a[j][i]=sum/l;
				else a[i][j]=a[j][i]=0;
			}
		}

//		for(int i=0;i<n;i++)
//			cout << a[i] << endl;

	}


	eigenval.resize(n);
	eigenvec.resize(n);

	xmippVector b;
	b.resize(n);
	xmippVector z;
	z.resize(n);
	xmippVector &d=eigenval;
	vector<xmippVector> &v=eigenvec;

	for(int i=0;i<n;i++){
		v[i].resize(n);
		v[i][i]=1.0;
		b[i]=d[i]=a[i][i];
	}

	int nrot=0;

	//Jacobi method (it=iterationn number)
	for(int it=1;it<=50;it++){
		xmippFeature tresh;
		xmippFeature sm=0.0;
		for (int ip = 0; ip < n - 1; ip++) {
			for (int iq = ip + 1; iq < n; iq++)
				sm += fabs(a[iq][ip]);
		}
		if(sm==0.0){//Done. Sort vectors
			for (int i = 0; i < n - 1; i++) {
				int k=i;
				xmippFeature p = d[i];

				for (int j = i + 1; j < n; j++)
					if (d[j] >= p)
						p = d[k = j];

				if (k != i){//Swap i<->k
					d[k] = d[i];
					d[i] = p;
					xmippVector t=v[i];
					v[i]=v[k];
					v[k]=t;
				}
			}
			return;
		}

		if(it<4) tresh=0.2*sm/(n*n); else tresh=0;

		for(int ip = 0; ip < n - 1; ip++){
			for(int iq = ip + 1; iq < n; iq++){
				xmippFeature g = 100.0 * fabs(a[iq][ip]);

				if (it > 4
				&& fabs(d[ip]) + g == fabs(d[ip])
				&& fabs(d[iq]) + g == fabs(d[iq]))
			    a[iq][ip]=0.0;
				else if (fabs(a[iq][ip])>tresh){
			    xmippFeature tau, t, s, c;
			    xmippFeature h = d[iq] - d[ip];
					if (fabs(h) + g == fabs(h))
						t = a[iq][ip] / h;
					else {
						xmippFeature theta = 0.5 * h / a[iq][ip];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
						if (theta < 0.0)
							t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[iq][ip];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[iq][ip]=0.0;

#define rotate(a,i,j,k,l) \
			g = a[i][j]; \
			h = a[k][l]; \
			a[i][j] = g - s *(h + g*tau); \
			a[k][l] = h + s*(g - h*tau);

					int j;
					for (j = 0; j < ip; j++){
						rotate(a, ip, j, iq, j)
					}
					for (j = ip + 1; j < iq; j++){
						rotate(a, j, ip, iq, j)
					}
					for (j = iq + 1; j < n; j++){
						rotate(a, j, ip, j, iq)
					}
					for (j = 0; j < n; j++){
						rotate(v, ip, j, iq, j)
					}

					nrot += 1;
				}//if
			}//for iq
		}//for ip
		for(int ip = 0; ip < n; ip++) {
	    b[ip] += z[ip];
	    d[ip] = b[ip];
	    z[ip] = 0.0;
		}
	}//for it

	throw Xmipp_error(1,"too many Jacobi iterations");
}

	/**Set identity matrix as eigenvector matrix*/
void xmippPC::setIdentity(int n)
{
	if(n<0) n=0;
	eigenval.resize(n);
	fill(eigenval.begin(),eigenval.end(),1.0);
	eigenvec.resize(n);
	for(int i=0;i<n;i++){
		eigenvec[i].resize(n);
		fill(eigenvec[i].begin(),eigenvec[i].end(),0.0);
		eigenvec[i][i]=1.0;
	}
}



