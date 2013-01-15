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

	}
	else if (method=="3d_clusters")
	{

	}
	else if (method=="intersect")
	{

	}
	else if (method=="difficult")
	{

	}
	else
		REPORT_ERROR(ERR_ARG_INCORRECT,"Incorrect method passed to generate data");
}

/*
        case 'helix'
        	t = [1:n]' / n;
        	t = t .^ (1.0) * 2 * pi;
			X = [(2 + cos(8 * t)) .* cos(t) (2 + cos(8 * t)) .* sin(t) sin(8 * t)] + noise * randn(n, 3);
        	%labels = uint8(t);
            labels = rem(round(t * 1.5), 2);

        case 'twinpeaks'
            inc = 1.5 / sqrt(n);
            [xx2, yy2] = meshgrid(-1:inc:1);
            xy = 1 - 2 * rand(2, n);
            X = [xy; sin(pi * xy(1,:)) .* tanh(3 * xy(2,:))]' + noise * randn(n, 3);
            X(:,3) = X(:,3) * 10;
            t = xy';
            %labels = uint8(X(:,3));
            labels = rem(sum(round((X + repmat(min(X, [], 1), [size(X, 1) 1])) ./ 10), 2), 2);

        case '3d_clusters'
            numClusters = 5;
            centers = 10 * rand(numClusters, 3);
            D = L2_distance(centers', centers');
            minDistance = min(D(D > 0));
            k = 1;
            n2 = n - (numClusters - 1) * 9;
            X = repmat(0, [n 3]);
            labels = repmat(0, [n 1]);
            for i=1:numClusters
                for j=1:ceil(n2 / numClusters)
                   X(k, 1:3) = centers(i, 1:3) + (rand(1, 3) - 0.5) * minDistance / sqrt(12);
                   labels(k) = i;
                   k = k + 1;
                end
            end
            X = X + noise * randn(size(X, 1), 3);
            t = [];

        case 'intersect'
            t = [1:n]' ./ n .* (2 * pi);
            x = cos(t);
            y = sin(t);
            height = rand(length(x), 1) * 5;
            X = [x x .* y height] + noise * randn(n, 3);
            %labels = uint8(5 * t);
            labels = rem(sum([round(t / 2) round(height / 2)], 2), 2);

        case 'difficult'
            % Generate underlying manifold
            no_dims = 5;
            no_points_per_dim = round(n ^ (1 / no_dims));
            l = linspace(0, 1, no_points_per_dim);
            t = combn(l, no_dims);

            % Generate high-dimensional dataset
            X = [cos(t(:,1)) tanh(3 * t(:,2)) t(:,1) + t(:,3) t(:,4) .* sin(t(:,2)) sin(t(:,1) + t(:,5)) t(:,5) .* cos(t(:,2)) t(:,5) + t(:,4) t(:,2) t(:,3) .* t(:,4) t(:,1)];
            X = X + noise * randn(size(X));

            % Generate labels for dataset (2x2x2x2x2 checkerboard pattern)
            tt = 1 + round(t);
            labels = rem(sum(tt, 2), 2);
 */
