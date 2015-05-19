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
}

void MultireferenceAligneability::run()
{
	//xmipp_multireference_aligneability --volume 1BRD.vol --sym c3 --odir testMultiReference/ --angles_file testMultiReference/angles_iter001_00.xmd --angles_file_ref testMultiReference/gallery_alignment/angles_iter001_00.xmd &

	randomize_random_generator();

    MetaData mdOutCL, mdOutQ;
	randomize_random_generator();
	MetaData mdExp, mdExpSort, mdProj;
	size_t nSamplesRandom = 100;
	size_t maxNImg;
	FileName fnOutCL,fnOutQ;
	fnOutCL = fnDir+"/clusteringTendency.xmd";
	fnOutQ = fnDir+"/validation.xmd";

	mdProj.read(finRef);
	mdExp.read(fin);

	mdExpSort.sort(mdExp,MDL_IMAGE_IDX,true,-1,0);
	size_t sz = mdExp.size();
	mdExpSort.getValue(MDL_IMAGE_IDX,maxNImg,sz);

	String expression;
	MDRow row;

	SymList SL;
	int symmetry, sym_order;
	SL.readSymmetryFile(fnSym.c_str());
	SL.isSymmetryGroup(fnSym.c_str(), symmetry, sym_order);

	double non_reduntant_area_of_sphere = SL.nonRedundantProjectionSphere(symmetry,sym_order);
	double area_of_sphere_no_symmetry = 4.*PI;
	double correction = std::sqrt(non_reduntant_area_of_sphere/area_of_sphere_no_symmetry);
	double validation = 0;

	ProgValidationNonTilt PVN;
	init_progress_bar(maxNImg);

	MetaData tempMdExp, tempMdProj;
	std::vector<double> sum_u_exp(nSamplesRandom);
	std::vector<double> sum_u_proj(nSamplesRandom);
	std::vector<double> sum_w_exp(nSamplesRandom);
	std::vector<double> sum_w_proj(nSamplesRandom);
	std::vector<double> H0(nSamplesRandom);
	std::vector<double> H_exp(nSamplesRandom);
	std::vector<double> H_proj(nSamplesRandom);

	for (size_t i=0; i<=maxNImg;i++)
	{
		expression = formatString("imageIndex == %lu",i);
		tempMdExp.importObjects(mdExp, MDExpression(expression));
		tempMdProj.importObjects(mdProj, MDExpression(expression));

		if ( (tempMdExp.size()==0) || (tempMdProj.size()==0))
			continue;

		PVN.obtainSumU(tempMdProj,sum_u_proj,H0);
		PVN.obtainSumU(tempMdExp,sum_u_exp,H0);
		PVN.obtainSumW(tempMdProj,sum_w_proj,sum_u_proj,H_proj,correction);
		PVN.obtainSumW(tempMdExp,sum_w_exp,sum_u_exp,H_exp,correction);

		std::sort(H_proj.begin(),H_proj.end());
		std::sort(H_exp.begin(),H_exp.end());

		double P = 0.;
		for(size_t j=0; j<sum_u_exp.size();j++)  //NOTA, u1 and u2 have 100 elements, then the condition sum_u1.size, does not matter
			P += H_proj.at(j)/H_exp.at(j);

		P /= (nSamplesRandom);
		validation += P;
		row.setValue(MDL_IMAGE_IDX,i);
		row.setValue(MDL_WEIGHT,P);
		mdOutCL.addRow(row);

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


void MultireferenceAligneability::P_calculus(FileName &fnMdProjAngles,FileName &fnMdExpAngles)
{
    MetaData mdOutCL, mdOutQ;
	randomize_random_generator();
	MetaData mdExp, mdExpSort, mdProj;
	size_t nSamplesRandom = 100;
	size_t maxNImg;
	FileName fnOutCL,fnOutQ;
	fnOutCL = fnDir+"/clusteringTendency.xmd";
	fnOutQ = fnDir+"/validation.xmd";

	mdProj.read(fnMdProjAngles);
	mdExp.read(fnMdExpAngles);
	mdExpSort.sort(mdExp,MDL_IMAGE_IDX,true,-1,0);
	size_t sz = mdExp.size();
	mdExpSort.getValue(MDL_IMAGE_IDX,maxNImg,sz);

	String expression;
	MDRow row;

	SymList SL;
	int symmetry, sym_order;
	SL.readSymmetryFile(fnSym.c_str());
	SL.isSymmetryGroup(fnSym.c_str(), symmetry, sym_order);

	double non_reduntant_area_of_sphere = SL.nonRedundantProjectionSphere(symmetry,sym_order);
	double area_of_sphere_no_symmetry = 4.*PI;
	double correction = std::sqrt(non_reduntant_area_of_sphere/area_of_sphere_no_symmetry);
	double validation = 0;

	ProgValidationNonTilt PVN;
	init_progress_bar(maxNImg);

	MetaData tempMdExp, tempMdProj;
	std::vector<double> sum_u_exp(nSamplesRandom);
	std::vector<double> sum_u_proj(nSamplesRandom);
	std::vector<double> sum_w_exp(nSamplesRandom);
	std::vector<double> sum_w_proj(nSamplesRandom);
	std::vector<double> H0(nSamplesRandom);
	std::vector<double> H_exp(nSamplesRandom);
	std::vector<double> H_proj(nSamplesRandom);

	for (size_t i=0; i<=maxNImg;i++)
	{
		expression = formatString("imageIndex == %lu",i);
		tempMdExp.importObjects(mdExp, MDExpression(expression));
		tempMdProj.importObjects(mdProj, MDExpression(expression));

		if ( (tempMdExp.size()==0) || (tempMdProj.size()==0))
			continue;

		PVN.obtainSumU(tempMdProj,sum_u_proj,H0);
		PVN.obtainSumU(tempMdExp,sum_u_exp,H0);
		PVN.obtainSumW(tempMdProj,sum_w_proj,sum_u_proj,H_proj,correction);
		PVN.obtainSumW(tempMdExp,sum_w_exp,sum_u_exp,H_exp,correction);

		std::sort(H_proj.begin(),H_proj.end());
		std::sort(H_exp.begin(),H_exp.end());

		double P = 0.;
		for(size_t j=0; j<sum_u_exp.size();j++)  //NOTA, u1 and u2 have 100 elements, then the condition sum_u1.size, does not matter
			P += H_proj.at(j)/H_exp.at(j);

		P /= (nSamplesRandom);
		validation += P;
		row.setValue(MDL_IMAGE_IDX,i);
		row.setValue(MDL_WEIGHT,P);
		mdOutCL.addRow(row);

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

void MultireferenceAligneability::calc_sumu(MetaData tempMd,std::vector<double> & sum_u)
{
	const size_t tempMdSz= tempMd.size();
    double xRan,yRan,zRan;
    double tilt,rot;
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
	        	rot = 2*PI*rnd_unif(0,1);
	        	tilt = std::acos(2*rnd_unif(0,1)-1);

	        	//tilt =(double(std::rand())/RAND_MAX)*(PI);		[ 0 ,PI]
	        	//rot  =(std::rand()-RAND_MAX/2)*(2*PI/RAND_MAX);   [-PI,PI]
	        	xRan = sin(tilt)*cos(rot);
	        	yRan = sin(tilt)*sin(rot);
	        	zRan = std::abs(cos(tilt));

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
}

