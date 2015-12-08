/***************************************************************************
 * Authors:     AUTHOR_NAME (jlvilas@cnb.csic.es)
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

#include "validation_nontilt.h"
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time
#include <algorithm>
#include "data/sampling.h"


ProgValidationNonTilt::ProgValidationNonTilt()
{
	rank=0;
	Nprocessors=1;
}

void ProgValidationNonTilt::readParams()
{
	fnParticles = getParam("--i");
	fnDir = getParam("--odir");
    fnSym = getParam("--sym");
    fnInit = getParam("--volume");
    useSignificant = checkParam("--useSignificant");
    significance_noise = getDoubleParam("--significance_noise");
}

void ProgValidationNonTilt::defineParams()
{
    addUsageLine("Validate a 3D reconstruction from its projections attending to directionality and spread of the angular assignments from a given significant value");
    addParamsLine("  [--i <md_file=\"\">]         : Metadata file with input projections");
    addParamsLine("  [--volume <md_file=\"\">]    : Volume to validate");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--significance_noise<float=0.95>] : Significance of the alignment with respect the noise");
    addParamsLine("  [--useSignificant]           : Use Significant as alignment method. If not use projection matching");
}

void ProgValidationNonTilt::run()
{
    //Clustering Tendency and Cluster Validity Stephen D. Scott
    randomize_random_generator();
    MetaData md,mdGallery,mdOut,mdOut2;
    MDRow row;

    FileName fnOut,fnOut2, fnGallery;
    fnOut = fnDir+"/clusteringTendency.xmd";
    fnGallery = fnDir+"/gallery.doc";
    fnOut2 = fnDir+"/validation.xmd";
    size_t nSamplesRandom = 500;

    md.read(fnParticles);
    mdGallery.read(fnGallery);

    size_t maxNImg;
    size_t sz = md.size();

    if (useSignificant)
    	md.getValue(MDL_IMAGE_IDX,maxNImg,sz);
    else
    {
    	md.getValue(MDL_ITEM_ID,maxNImg,sz);
    }

    String expression;
    MDRow rowP,row2;
    SymList SL;
    int symmetry, sym_order;
    SL.readSymmetryFile(fnSym.c_str());
    SL.isSymmetryGroup(fnSym.c_str(), symmetry, sym_order);

/*
    double non_reduntant_area_of_sphere = SL.nonRedundantProjectionSphere(symmetry,sym_order);
    double area_of_sphere_no_symmetry = 4.*PI;
    double correction = std::sqrt(non_reduntant_area_of_sphere/area_of_sphere_no_symmetry);
*/
    double correction = 1;
    double validation = 0;
    double num_images = 0;

	MetaData tempMd;
	std::vector<double> sum_u(nSamplesRandom);
	double sum_w=0;
	std::vector<double> H0(nSamplesRandom);
	std::vector<double> H(nSamplesRandom);
	std::vector<double> p(nSamplesRandom);

	if (rank==0)
		init_progress_bar(maxNImg);

	for (size_t idx=0; idx<=maxNImg;idx++)
	{
		if ((idx+1)%Nprocessors==rank)
		{
			if (useSignificant)
				expression = formatString("imageIndex == %lu",idx);
			else
				expression = formatString("itemId == %lu",idx);

			tempMd.importObjects(md, MDExpression(expression));

			if (tempMd.size()==0)
				continue;

			//compute H_0 from noise
			//obtainSumU(tempMd,sum_u,H0);
			obtainSumU_2(mdGallery, tempMd,sum_u,H0);
			//compute H from experimental
			obtainSumW(tempMd,sum_w,sum_u,H,correction);

			std::sort(H0.begin(),H0.end());
			std::sort(H.begin(),H.end());

			double P = 0;
			for(size_t j=0; j<sum_u.size();j++)
			{
				//P += H0.at(j)/H.at(j);
				P += H0.at(size_t((1-significance_noise)*nSamplesRandom))/H.at(j);
				p.at(j) = H0.at(j)/H.at(j);
			}

			P /= (nSamplesRandom);

			if (useSignificant)
				rowP.setValue(MDL_IMAGE_IDX,idx);
			else
				rowP.setValue(MDL_ITEM_ID,idx);

			rowP.setValue(MDL_WEIGHT,P);
			rowP.setValue(MDL_WEIGHT_P,p.at(size_t(significance_noise*nSamplesRandom)));
			mdPartial.addRow(rowP);
			tempMd.clear();

			if (rank==0)
				progress_bar(idx+1);
		}
	}

	if (rank==0)
		progress_bar(maxNImg);

	synchronize();
	gatherClusterability();

	if (rank == 0)
	{
		mdPartial.write(fnOut);
		std::vector<double> P;
		mdPartial.getColumnValues(MDL_WEIGHT,P);

		for (size_t idx=0; idx< P.size();idx++)
		{
			if (P[idx] > 1)
				validation += 1.;
			num_images += 1.;
		}
		validation /= (num_images);
	}

    row2.setValue(MDL_IMAGE,fnInit);
    row2.setValue(MDL_WEIGHT,validation);
    mdOut2.addRow(row2);
    mdOut2.write(fnOut2);
}

void ProgValidationNonTilt::obtainSumU(const MetaData & tempMd,std::vector<double> & sum_u,std::vector<double> & H0)
{

	const size_t tempMdSz= tempMd.size();
    double xRan,yRan,zRan;
    double x,y;
    double sumWRan;
    double * xRanArray = new double[tempMdSz];
    double * yRanArray = new double[tempMdSz];
    double * zRanArray  = new double[tempMdSz];
    std::vector<double> weightV;
    double a;

    for (size_t n=0; n<sum_u.size(); n++)
    {
        sumWRan = 0;
        for (size_t nS=0; nS<tempMd.size(); nS++)
        {
             /*
                x = sin(tilt*PI/180)*cos(rot*PI/180);
        		y = sin(tilt*PI/180)*sin(rot*PI/180);
        		z = std::abs(cos(tilt*PI/180));
             */

        	//http://mathworld.wolfram.com/SpherePointPicking.html
        	x = 2*(double(std::rand())-RAND_MAX/2)/RAND_MAX;
        	y = 2*(double(std::rand())-RAND_MAX/2)/RAND_MAX;

        	while (x*x+y*y >= 1 )
        	{
            	x = 2*(std::rand()-RAND_MAX/2)/RAND_MAX;
            	y = 2*(std::rand()-RAND_MAX/2)/RAND_MAX;
        	}

        	xRan = 2*x*std::sqrt(1-x*x-y*y);
        	yRan = 2*y*std::sqrt(1-x*x-y*y);
        	zRan = std::abs(1-2*(x*x+y*y));

/*        	tilt = (double(std::rand())/RAND_MAX)*(PI);
        	rot  = (std::rand()-RAND_MAX/2)*(2*PI/RAND_MAX);
        	xRan = sin(tilt)*cos(rot);
        	yRan = sin(tilt)*sin(rot);
        	zRan = std::abs(cos(tilt));
*/
        	//std::cout << tilt << " " << rot << std::endl;
        	//std::cout << xRan << " " << yRan << " " << zRan << " " << std::endl;
            //zRan=(std::rand()-RAND_MAX/2);
            //norm = std::sqrt(xRan*xRan+yRan*yRan+zRan*zRan);
            //xRan = (xRan/norm);
            //yRan = (yRan/norm);
            //zRan = std::abs(zRan/norm);

            xRanArray[nS] = xRan;
            yRanArray[nS] = yRan;
            zRanArray[nS] = zRan;
            //tempMd.getColumnValues(MDL_WEIGHT, weightV);
            tempMd.getColumnValues(MDL_MAXCC, weightV);

            std::random_shuffle(weightV.begin(), weightV.end());
        }

        sumWRan = 0;
        double WRan=0;
        double tempWRan, tempW1, tempW2, temp;
        for (size_t nS1=0; nS1<tempMd.size(); nS1++)
        {
            tempWRan = 1e3;
            for (size_t nS2=0; nS2<tempMd.size(); nS2++)
            {
            	temp = xRanArray[nS1]*xRanArray[nS2]+yRanArray[nS1]*yRanArray[nS2]+zRanArray[nS1]*zRanArray[nS2];
                if (temp < 1)
                	a = std::abs(std::acos(temp));
                else
                	a = 0;

                if ( (a<tempWRan) && (a > 0.00001) && (temp<1) )
                {
                    tempWRan = a;
                    tempW2 = weightV[nS2];
                    tempW1 = weightV[nS1];
                    WRan = a*std::exp(std::abs(tempW1-tempW2))*std::exp(-(tempW1+tempW2));
                    if (WRan == 0)
                    	WRan = a;

                }
            }
            sumWRan += WRan;
        }

        if (sumWRan == 0)
        	sumWRan = 0.075*tempMd.size();

        sum_u.at(n)=sumWRan;

    }

    size_t idx = 0;
    while (idx < sum_u.size())
    {
        std::random_shuffle(sum_u.begin(), sum_u.end());

        if(sum_u.at(0) != sum_u.at(1))
        {
            H0[idx] = sum_u.at(0)/(sum_u.at(0)+sum_u.at(1));
            idx += 1;
        }
    }

    delete xRanArray;
    delete yRanArray;
    delete zRanArray;

}

#define _FOR_ALL_OBJECTS_IN_METADATA2(__md) \
        for(MDIterator __iter2(__md); __iter2.hasNext(); __iter2.moveNext())
void ProgValidationNonTilt::obtainSumW(const MetaData & tempMd, double & sum_W, std::vector<double> & sum_u, std::vector<double> & H, const double factor)
{
    double a;
    double rot,tilt,w;
    double x,y,z;
    double xx,yy,zz;
    bool mirror;
    double w2;
    double tempW;
    double W=0;
    double temp;
    double sumW;

    sumW = 0;
    FOR_ALL_OBJECTS_IN_METADATA(tempMd)
    {
        tempMd.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
        tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
        tempMd.getValue(MDL_FLIP,mirror,__iter.objId);
        if (mirror == 1)
        	tilt = tilt + 180;

        //tempMd.getValue(MDL_WEIGHT,w,__iter.objId);
        tempMd.getValue(MDL_MAXCC,w,__iter.objId);
        x = sin(tilt*PI/180.)*cos(rot*PI/180.);
        y = sin(tilt*PI/180.)*sin(rot*PI/180.);
        z = std::abs(cos(tilt*PI/180.));

        tempW = 1e3;
        _FOR_ALL_OBJECTS_IN_METADATA2(tempMd)
        {
            tempMd.getValue(MDL_ANGLE_ROT,rot,__iter2.objId);
            tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter2.objId);
            //tempMd.getValue(MDL_WEIGHT,w2,__iter2.objId);
            tempMd.getValue(MDL_MAXCC,w2,__iter2.objId);
            tempMd.getValue(MDL_FLIP,mirror,__iter2.objId);
            if (mirror == 1)
            	tilt = tilt + 180;

            xx = sin(tilt*PI/180.)*cos(rot*PI/180.);
            yy = sin(tilt*PI/180.)*sin(rot*PI/180.);
            zz = std::abs(cos(tilt*PI/180.));
            temp = x*xx+y*yy+z*zz;

            if (temp < 1)
            	a = std::abs(std::acos(temp));
            else
            	a = 0;

            if ( (a<tempW) && (a > 0.00001) && (temp<1 ))
            {
                W = a*std::exp(std::abs(w-w2))*std::exp(-(w+w2));
                tempW = a;

                if (W == 0)
                	W = a;
            }
        }
        sumW +=  W;
    }

    //Here the problem is that maybe the changes are only in psi angle. We give the best solution assuming an angular sampling of 0.5º
    if (sumW == 0)
    	sumW = 0.075*tempMd.size();

    sum_W = sumW;
    for (size_t n=0; n<sum_u.size(); n++)
    {
        H[n] = sumW/(sumW+factor*sum_u.at(n));
    }
}



void ProgValidationNonTilt::obtainSumU_2(const MetaData & mdGallery, const MetaData & tempMd,std::vector<double> & sum_u,std::vector<double> & H0)
{

	const size_t tempMdSz= tempMd.size();
	const size_t numGallery= mdGallery.size();

    double xRan,yRan,zRan;
    size_t indx;
    double sumWRan;
    double * xRanArray = new double[tempMdSz];
    double * yRanArray = new double[tempMdSz];
    double * zRanArray  = new double[tempMdSz];
    std::vector<double> weightV;
    double a;

    double rot,tilt,w;
    bool mirror;

    for (size_t n=0; n<sum_u.size(); n++)
    {
        sumWRan = 0;

        for (size_t nS=0; nS<tempMd.size(); nS++)
        {
			indx = (size_t) (double( std::rand())*numGallery)/RAND_MAX;

			while ( (indx ==0) || (indx > numGallery) )
				indx = (size_t) (double( std::rand())*numGallery )/RAND_MAX;

        	mdGallery.getValue(MDL_ANGLE_ROT,rot,indx);
        	mdGallery.getValue(MDL_ANGLE_TILT,tilt,indx);
        	mdGallery.getValue(MDL_FLIP,mirror,indx);

            if (mirror == 1)
            	tilt = tilt + 180;

            xRan = sin(tilt*PI/180.)*cos(rot*PI/180.);
            yRan = sin(tilt*PI/180.)*sin(rot*PI/180.);
            zRan = std::abs(cos(tilt*PI/180.));

            xRanArray[nS] = xRan;
            yRanArray[nS] = yRan;
            zRanArray[nS] = zRan;
            //tempMd.getColumnValues(MDL_WEIGHT, weightV);
            tempMd.getColumnValues(MDL_MAXCC, weightV);
            std::random_shuffle(weightV.begin(), weightV.end());
        }


        sumWRan = 0;
        double WRan=0;
        double tempWRan, tempW1, tempW2, temp;
        for (size_t nS1=0; nS1<tempMd.size(); nS1++)
        {
            tempWRan = 1e3;
            for (size_t nS2=0; nS2<tempMd.size(); nS2++)
            {
            	temp = xRanArray[nS1]*xRanArray[nS2]+yRanArray[nS1]*yRanArray[nS2]+zRanArray[nS1]*zRanArray[nS2];
                if (temp < 1)
                	a = std::abs(std::acos(temp));
                else
                	a = 0;

                if ( (a<tempWRan) && (a > 0.00001) && (temp<1) )
                {
                    tempWRan = a;
                    tempW2 = weightV[nS2];
                    tempW1 = weightV[nS1];
                    WRan = a*std::exp(std::abs(tempW1-tempW2))*std::exp(-(tempW1+tempW2));
                    if (WRan == 0)
                    	WRan = a;
                }
            }
            sumWRan += WRan;
        }

        if (sumWRan == 0)
        	sumWRan = 0.075*tempMd.size();

        sum_u.at(n)=sumWRan;

    }

    size_t idx = 0;
    while (idx < sum_u.size())
    {
        std::random_shuffle(sum_u.begin(), sum_u.end());

        if(sum_u.at(0) != sum_u.at(1))
        {
            H0[idx] = sum_u.at(0)/(sum_u.at(0)+sum_u.at(1));
            idx += 1;
        }
    }

    delete xRanArray;
    delete yRanArray;
    delete zRanArray;

}


