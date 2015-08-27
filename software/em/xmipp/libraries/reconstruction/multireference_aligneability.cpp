/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
 * 							(jlvilas@cnb.csic.es)
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

#include "multireference_aligneability.h"
#include "validation_nontilt.h"
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <fstream>
#include "data/sampling.h"
#include "project.h"
#include <string>
#include <data/metadata.h>


void MultireferenceAligneability::readParams()
{
    fnInit = getParam("--volume");
    fin = getParam("--angles_file");
    finRef = getParam("--angles_file_ref");
    fnSym = getParam("--sym");
    fnDir = getParam("--odir");
    donNotUseWeights= checkParam("--dontUseWeights");

}

void MultireferenceAligneability::defineParams()
{
    //usage
    addUsageLine("This function takes a volume and a set of projections with orientations. The volume is projected into the set of projection directions defining the "
    		"the reference projections. Thus, using the projections and references, the calculus of H0k, Hk and Pk is carried out");
    addParamsLine("  [--volume <md_file=\"\">]    : Input volume to be validated");
    addParamsLine("  [--angles_file <file=\".\">]     : Input metadata with projections and orientations");
    addParamsLine("  [--angles_file_ref <file=\".\">] : Input reference metadata with projections and orientations obtained from the projection orientations and the volume ");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections"); //TODO the input will be two doc files one from the exp and the other from refs
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--dontUseWeights]           : Do not use the particle weigths in the clusterability calculation ");

}

void MultireferenceAligneability::run()
{
	//xmipp_multireference_aligneability --volume 1BRD.vol --sym c3 --odir testMultiReference/ --angles_file testMultiReference/angles_iter001_00.xmd --angles_file_ref testMultiReference/gallery_alignment/angles_iter001_00.xmd &
    randomize_random_generator();

    MetaData mdOutCL, mdOutQ;
	MetaData mdExp, mdExpSort, mdProj;
	size_t maxNImg;
	FileName fnOutCL,fnOutQ;
	fnOutCL = fnDir+"/clusteringTendency.xmd";
	fnOutQ = fnDir+"/validation.xmd";

    SymList SL;
    int symmetry, sym_order;
    SL.readSymmetryFile(fnSym.c_str());
    SL.isSymmetryGroup(fnSym.c_str(), symmetry, sym_order);
    double non_reduntant_area_of_sphere = SL.nonRedundantProjectionSphere(symmetry,sym_order);
    double area_of_sphere_no_symmetry = 4.*PI;
    double correction = std::sqrt(non_reduntant_area_of_sphere/area_of_sphere_no_symmetry);

	mdProj.read(finRef);
	mdExp.read(fin);

	mdExpSort.sort(mdExp,MDL_IMAGE_IDX,true,-1,0);
	size_t sz = mdExp.size();
	mdExpSort.getValue(MDL_IMAGE_IDX,maxNImg,sz);

	String expression;
	MDRow row;

	double validation = 0;
	init_progress_bar(maxNImg);

	MetaData tempMdExp, tempMdProj;
	double sum_w_exp;
	double sum_w_proj;
	double sum_noise = 0;

	expression = formatString("imageIndex == %lu",maxNImg);
	tempMdExp.importObjects(mdExp, MDExpression(expression));
	size_t numProjs = tempMdExp.size();
	tempMdExp.clear();

	//Noise
	calc_sumw(numProjs, sum_noise);
	sum_noise *= correction;

	double rank = 0.;
	char hold;
	FileName imagePath;
	for (size_t i=0; i<=maxNImg;i++)
	{
		expression = formatString("imageIndex == %lu",i);
		tempMdExp.importObjects(mdExp, MDExpression(expression));
		tempMdProj.importObjects(mdProj, MDExpression(expression));

		if ( (tempMdExp.size()==0) || (tempMdProj.size()==0))
			continue;

		tempMdExp.getRow(row,1);
		calc_sumu(tempMdExp, sum_w_exp);
		calc_sumu(tempMdProj, sum_w_proj);

		rank = 1/(sum_w_proj-sum_noise)*(sum_w_exp-sum_noise);

		if (rank>1)
			rank=1;
		else if (rank < 0)
			rank = 0;

		validation += (rank>0.5);
		tempMdExp.getValue(MDL_IMAGE,imagePath,1);
		row.setValue(MDL_IMAGE,imagePath);
		row.setValue(MDL_IMAGE_IDX,i);
		row.setValue(MDL_VOLUME_SCORE1,rank);
		row.setValue(MDL_VOLUME_SCORE2,sum_w_proj);
		row.setValue(MDL_VOLUME_SCORE3,sum_w_exp);
		row.setValue(MDL_VOLUME_SCORE4,sum_noise);

		mdOutCL.addRow(row);

		row.clear();
		tempMdExp.clear();
		tempMdProj.clear();
		progress_bar(i+1);
	}

	mdOutCL.write(fnOutCL);
	validation /= (maxNImg+1);
	row.clear();
    row.setValue(MDL_IMAGE,fnInit);
    row.setValue(MDL_WEIGHT,validation);
    mdOutQ.addRow(row);
    mdOutQ.write(fnOutQ);
}


void MultireferenceAligneability::write_projection_file()
{

	String xdim= integerToString(Xdim);
	String ydim= integerToString(Ydim);
	FileName filnam=fnDir+"/params";
	std::ofstream myfile;
	myfile.open(filnam.c_str());
	myfile << "# XMIPP_STAR_1 *\n";
	myfile << "# \n";
	myfile << "data_block1 \n";
	myfile << "_dimensions2D '"+xdim+" "+ydim+"' \n";
	myfile << "_projAngleFile "+fin+" \n";
	myfile << "_noisePixelLevel   '0 0' \n";    //where the variance is the first number and the second one is the mean
	myfile.close();
}

#define _FOR_ALL_OBJECTS_IN_METADATA2(__md) \
        for(MDIterator __iter2(__md); __iter2.hasNext(); __iter2.moveNext())

void MultireferenceAligneability::calc_sumu(const MetaData & tempMd, double & sum_W)
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
    size_t idx = 0;
    FOR_ALL_OBJECTS_IN_METADATA(tempMd)
    {
        tempMd.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
        tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
        //tempMd.getValue(MDL_WEIGHT,w,__iter.objId);
        tempMd.getValue(MDL_MAXCC,w,__iter.objId);

        idx++;
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
            xx = sin(tilt*PI/180.)*cos(rot*PI/180.);
            yy = sin(tilt*PI/180.)*sin(rot*PI/180.);
            zz = std::abs(cos(tilt*PI/180.));
            a = std::abs(std::acos(x*xx+y*yy+z*zz));

            if ( (a<tempW) && (a != 0))
            {
            	if (donNotUseWeights)
            		W = a;
            	else
            		//W = a*std::exp(-2*(w+w2));
            		W = a*std::exp(std::abs(w-w2))*std::exp(-(w+w2));

                tempW = a;
                //std::cout << x << " " << y << " " << z << " " << xx << " " << yy << " " << zz << " " <<  w << " " << w2 << " " << a << " " << W << " " << std::exp(std::abs(w-w2)) << " " << std::exp(-(w+w2)) << std::endl;
            }
        }
        sumW +=  W;
    }

    //calculate the std of weightsDistribution and use it as a weight
/*    double sum=0;
    double mean=0;
    double stdDev=0;
    double tmp=0;
    for (idx = 0; idx < tempMd.size(); idx++)
    	sum += weightsDistribution[idx];

    mean =sum/ tempMd.size();

    for (idx = 0; idx < tempMd.size(); idx++)
    	tmp += (weightsDistribution[idx]-mean)*(weightsDistribution[idx]-mean);

    stdDev = std::sqrt(tmp/(tempMd.size()+1));
    //sum_W = std::exp(-stdDev)*sumW;
 */
    sum_W = sumW;
    //std::cout << " sumW : " << sumW << std::endl;
}

void MultireferenceAligneability::calc_sumw(const size_t num, double & sumw)
{
	size_t trials=200;
    double xRan,yRan,zRan;
    double x,y;
    double sumWRan;
    double * xRanArray = new double[num];
    double * yRanArray = new double[num];
    double * zRanArray  = new double[num];
    double a;
    sumw=0;

    for (size_t n=0; n < trials; n++)
    {
        sumWRan = 0;
        for (size_t nS=0; nS<num; nS++)
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

            xRanArray[nS] = xRan;
            yRanArray[nS] = yRan;
            zRanArray[nS] = zRan;
        }

        sumWRan = 0;
        double WRan, tempWRan, tempW1, tempW2;
        for (size_t nS1=0; nS1<num; nS1++)
        {
            tempWRan = 1e3;
            for (size_t nS2=0; nS2<num; nS2++)
            {
                a = std::abs(std::acos(xRanArray[nS1]*xRanArray[nS2]+yRanArray[nS1]*yRanArray[nS2]+zRanArray[nS1]*zRanArray[nS2]));
                if ( (a<tempWRan) && (a != 0))
                {
                    tempWRan = a;
                    WRan = a;
                }
            }
            sumWRan += WRan;
        }

        sumw += sumWRan;
    }

    sumw /= trials;
}


