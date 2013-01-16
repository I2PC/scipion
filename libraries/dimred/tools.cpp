/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2013)
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

#include "generate_data.h"

void GenerateData::generateNewDataset(const String& method, int N, double noise)
{
	label.resizeNoCopy(N);
	X.resizeNoCopy(N,3);
	t.resizeNoCopy(N,2);
	if (noise==0)
		init_random_generator(0);
	if (method=="swiss")
	{
		for (int i=0; i<N; ++i)
		{
			// Generate t
			MAT_ELEM(t,i,0)=(3 * PI / 2) * (1 + 2 * rnd_unif()); // t
			MAT_ELEM(t,i,1)=30*rnd_unif(); // height
			double localT=MAT_ELEM(t,i,0);
			double localHeight=MAT_ELEM(t,i,1);

			// Generate X
			double s,c;
			sincos(localT,&s,&c);
			MAT_ELEM(X,i,0)=localT*c+noise*rnd_gaus();
			MAT_ELEM(X,i,1)=localHeight+noise*rnd_gaus();
			MAT_ELEM(X,i,2)=localT*s+noise*rnd_gaus();

			// Generate label
			VEC_ELEM(label,i)=(unsigned char)(round(localT/2)+round(localHeight/12))%2;
		}
	}
	else if (method=="helix")
	{
		double iN=1.0/N;
		for (int i=0; i<N; ++i)
		{
			// Generate t
			MAT_ELEM(t,i,0)=2 * PI * i*iN;
			MAT_ELEM(t,i,1)=30*rnd_unif(); // height
			double localT=MAT_ELEM(t,i,0);

			// Generate X
			double s,c;
			sincos(localT,&s,&c);
			double s8,c8;
			sincos(8*localT,&s8,&c8);
			MAT_ELEM(X,i,0)=(2 + c8)*c+noise*rnd_gaus();
			MAT_ELEM(X,i,1)=(2 + c8)*s+noise*rnd_gaus();
			MAT_ELEM(X,i,2)=s8+noise*rnd_gaus();

			// Generate label
			VEC_ELEM(label,i)=(unsigned char)(round(localT * 1.5))%2;
		}
	}
	else if (method=="twinpeaks")
	{
        int actualN=round(sqrt(N));
    	X.resizeNoCopy(actualN*actualN,3);
    	t.resizeNoCopy(actualN*actualN,2);
    	label.resizeNoCopy(actualN*actualN);
        for (int ii=0; ii<actualN; ii++)
        {
        	for (int jj=0; jj<actualN; jj++)
        	{
        		int i=ii*actualN+jj;

            	// Generate t
        		double x = 1 - 2 * rnd_unif();
        		double y = 1 - 2 * rnd_unif();
        		MAT_ELEM(t,i,0)=x;
        		MAT_ELEM(t,i,1)=y;

    			// Generate X
    			MAT_ELEM(X,i,0)=x+noise*rnd_gaus();
    			MAT_ELEM(X,i,1)=y+noise*rnd_gaus();
    			double z=10*sin(PI * x) * tanh(3 * y);
    			MAT_ELEM(X,i,2)=z+noise*rnd_gaus();

    			// Generate label
    			VEC_ELEM(label,i)=(unsigned char)(round(0.1*(x+y+z-3)))%2;
        	}
        }
	}
	else if (method=="3d_clusters")
	{
		// Create centers
		std::vector<Matrix1D<double> > centers;
		Matrix1D<double> center(3);
		const int Nclusters=5;
		for (int i=0; i<Nclusters; i++)
		{
			FOR_ALL_ELEMENTS_IN_MATRIX1D(center)
				VEC_ELEM(center,i)=10*rnd_unif();
			centers.push_back(center);
		}

		// Measure the minimum distance between centers
		Matrix1D<double> diff;
		double minDistance=1e38;
		for (int i=0; i<Nclusters-1; ++i)
			for (int j=i+1; j<Nclusters; ++j)
			{
				diff=centers[i]-centers[j];
				double distance=diff.module();
				minDistance=std::min(minDistance,distance);
			}

		// Create clusters
		t.initZeros();
		double sigma=minDistance/sqrt(12);
		for (int n=0; n<N; ++n)
		{
			int i=(Nclusters*n)/N;
			const Matrix1D<double> &center=centers[i];
			MAT_ELEM(X,n,0)=XX(center)+(rnd_unif()-0.5)*sigma+noise*rnd_gaus();
			MAT_ELEM(X,n,1)=YY(center)+(rnd_unif()-0.5)*sigma+noise*rnd_gaus();
			MAT_ELEM(X,n,2)=ZZ(center)+(rnd_unif()-0.5)*sigma+noise*rnd_gaus();
		}
	}
	else if (method=="intersect")
	{
		double iN=1.0/N;
		for (int i=0; i<N; ++i)
		{
			// Generate t
			MAT_ELEM(t,i,0)=2 * PI * i*iN;
			MAT_ELEM(t,i,1)=5*rnd_unif();
			double localT=MAT_ELEM(t,i,0);
			double height=MAT_ELEM(t,i,1);

			// Generate X
			double s,c;
			sincos(localT,&s,&c);
			MAT_ELEM(X,i,0)=c+noise*rnd_gaus();
			MAT_ELEM(X,i,1)=c*s+noise*rnd_gaus();
			MAT_ELEM(X,i,2)=height+noise*rnd_gaus();

			// Generate label
			VEC_ELEM(label,i)=(unsigned char)(round(localT *0.5)+round(height*0.5))%2;
		}
	}
	else
		REPORT_ERROR(ERR_ARG_INCORRECT,"Incorrect method passed to generate data");
}
