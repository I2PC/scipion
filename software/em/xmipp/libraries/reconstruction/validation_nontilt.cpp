/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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
	sampling_rate = 1;
}

void ProgValidationNonTilt::readParams()
{
    fnDir = getParam("--odir");
    fnSym = getParam("--sym");
    fnInit = getParam("--volume");
    sampling_rate = getDoubleParam("--sampling_rate");
}

void ProgValidationNonTilt::defineParams()
{
    //usage
    addUsageLine("Validate a 3D reconstruction from its projections attending to directionality and spread of the angular assignments from a given significant value");
    //params
    addParamsLine("  [--volume <md_file=\"\">]    : Volume to validate");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--sampling_rate <s=1>]      : Sampling rate in A/px");
}

void ProgValidationNonTilt::run()
{
    //Clustering Tendency and Cluster Validity Stephen D. Scott
    randomize_random_generator();
    //char buffer[400];
    //sprintf(buffer, "xmipp_reconstruct_significant -i %s  --initvolumes %s --odir %s --sym  %s --iter 1 --alpha0 %f --angularSampling %f",fnIn.c_str(), fnInit.c_str(),fnDir.c_str(),fnSym.c_str(),alpha0,angularSampling);
    //system(buffer);

    MetaData md,mdOut,mdOut2;
    FileName fnMd,fnOut,fnOut2;
    fnMd = fnDir+"/angles_iter001_00.xmd";
    fnOut = fnDir+"/clusteringTendency.xmd";
    fnOut2 = fnDir+"/validation.xmd";
    size_t nSamplesRandom = 100;

    md.read(fnMd);
    size_t maxNImg;
    size_t sz = md.size();
    md.getValue(MDL_IMAGE_IDX,maxNImg,sz);

    String expression;
    MDRow rowP,row2;
    SymList SL;
    int symmetry, sym_order;
    SL.readSymmetryFile(fnSym.c_str());
    SL.isSymmetryGroup(fnSym.c_str(), symmetry, sym_order);

    double non_reduntant_area_of_sphere = SL.nonRedundantProjectionSphere(symmetry,sym_order);
    double area_of_sphere_no_symmetry = 4.*PI;
    double correction = std::sqrt(non_reduntant_area_of_sphere/area_of_sphere_no_symmetry);
    double validation = 0;

	MetaData tempMd;
	std::vector<double> sum_u(nSamplesRandom);
	std::vector<double> sum_w(nSamplesRandom);
	std::vector<double> H0(nSamplesRandom);
	std::vector<double> H(nSamplesRandom);

	if (rank==0)
		init_progress_bar(maxNImg);

	for (size_t idx=0; idx<=maxNImg;idx++)
	{
		if ((idx+1)%Nprocessors==rank)
		{
			expression = formatString("imageIndex == %lu",idx);
			tempMd.importObjects(md, MDExpression(expression));

			if (tempMd.size()==0)
				continue;

			//compute H_0 from noise
			obtainSumU(tempMd,sum_u,H0);
			//compute H from experimental
			obtainSumW(tempMd,sum_w,sum_u,H,correction);

			std::sort(H0.begin(),H0.end());
			std::sort(H.begin(),H.end());

			double P = 0;
			for(size_t j=0; j<sum_u.size();j++)
				P += H0.at(j)/H.at(j);

			P /= (nSamplesRandom);
			rowP.setValue(MDL_IMAGE_IDX,idx);
			rowP.setValue(MDL_WEIGHT,P);
			mdPartial.addRow(rowP);

			//sum_u.clear();
			//sum_w.clear();
			//H0.clear();
			//H.clear();
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
			if (P[idx] > 1)
				validation++;
	    validation /= (maxNImg+1);

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
            tempMd.getColumnValues(MDL_WEIGHT, weightV);

            std::random_shuffle(weightV.begin(), weightV.end());
        }

        sumWRan = 0;
        double WRan, tempWRan, tempW1, tempW2;
        for (size_t nS1=0; nS1<tempMd.size(); nS1++)
        {
            tempWRan = 1e3;
            for (size_t nS2=0; nS2<tempMd.size(); nS2++)
            {
                a = std::abs(std::acos(xRanArray[nS1]*xRanArray[nS2]+yRanArray[nS1]*yRanArray[nS2]+zRanArray[nS1]*zRanArray[nS2]));
                if ( (a<tempWRan) && (a != 0))
                {
                    tempWRan = a;
                    tempW2 = weightV[nS2];
                    tempW1 = weightV[nS1];
                    WRan = a*std::exp(std::abs(tempW1-tempW2))*std::exp(-(tempW1+tempW2));
                }
            }
            sumWRan += WRan;
        }
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
void ProgValidationNonTilt::obtainSumW(const MetaData & tempMd,std::vector<double> & sum_W,std::vector<double> & sum_u,std::vector<double> & H, const double factor)
{
    double a;
    double rot,tilt,w;
    double x,y,z;
    double xx,yy,zz;
    double w2;
    double tempW;
    double W;
    double sumW;

    sumW = 0;
    FOR_ALL_OBJECTS_IN_METADATA(tempMd)
    {
        tempMd.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
        tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
        tempMd.getValue(MDL_WEIGHT,w,__iter.objId);
        x = sin(tilt*PI/180.)*cos(rot*PI/180.);
        y = sin(tilt*PI/180.)*sin(rot*PI/180.);
        z = std::abs(cos(tilt*PI/180.));

        tempW = 1e3;
        _FOR_ALL_OBJECTS_IN_METADATA2(tempMd)
        {
            tempMd.getValue(MDL_ANGLE_ROT,rot,__iter2.objId);
            tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter2.objId);
            tempMd.getValue(MDL_WEIGHT,w2,__iter2.objId);
            xx = sin(tilt*PI/180.)*cos(rot*PI/180.);
            yy = sin(tilt*PI/180.)*sin(rot*PI/180.);
            zz = std::abs(cos(tilt*PI/180.));
            a = std::abs(std::acos(x*xx+y*yy+z*zz));

            if ( (a<tempW) && (a != 0))
            {
                W = a*std::exp(std::abs(w-w2))*std::exp(-(w+w2));
                tempW = a;
            }
        }
        sumW +=  W;
    }

    for (size_t n=0; n<sum_u.size(); n++)
    {
        H[n] = sumW/(sumW+factor*sum_u.at(n));
    }
}

