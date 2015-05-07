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

//ejecución del programa
//xmipp_validation_nontilt --volume 4MOA.vol --odir Sig4MOA --sym c1


//xmipp_multireference_aligneability --volume alignment_reference_images.vol  --odir Sig4MOA --sym c1
//xmipp_multireference_aligneability --volume alignment_reference_images.vol  --odir Output_dir --sym c1

//xmipp_multireference_aligneability --stackimages input_volume.stk --volume alignment_reference_images.vol  --odir Output_dir --sym c1  --significance 0.05


#include "multireference_aligneability.h"
#include "validation_nontilt.h"
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time
#include <algorithm>
#include <fstream>
#include "data/sampling.h"
#include "project.h"
#include <string>
#include <data/metadata.h>


void MultireferenceAligneability::readParams()
{
    fnDir = getParam("--odir");
    fnSym = getParam("--sym");
    fnInit = getParam("--volume");
    fin = getParam("--stackimages");
    fsig = getParam("--significance");

}

void MultireferenceAligneability::defineParams()
{
    //usage
    addUsageLine("This function takes a volume and projects it in a set of directions defined in a metadata. Thus, using the projections,"
    		"the calculus of H0 and Hk is carried out");
    //params
    addParamsLine("  [--volume <md_file=\"\">]    : Volume to be projected and from its projections then the calculus of"
    		" H0 and Hk will be performed");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--stackimages <inputDir=\".\">]   : Stack of projections");
    addParamsLine("  [--significance <fsig=0.05>]  : Significance (default signif=0.05");
}

void MultireferenceAligneability::run()
{
    randomize_random_generator();
    char buffer[400];

    //STEP1 From a stack of images adn comparing with a volume, a set of angles (set1) is obtained via significant (alignment)
    	//input_volume.stk
    	//alignment_reference_images.vol Initial volume
    sprintf(buffer, "xmipp_reconstruct_significant  -i %s --initvolumes %s --odir %s --sym %s --iter 1 --alpha0 %s", fin.c_str(), fnInit.c_str(),fnDir.c_str(), fnSym.c_str(),fsig.c_str());
    //system(buffer);

    write_projection_file();
    std::cout << " " << std::endl;
    std::cout << "STEP1 ENDED" << std::endl;
    std::cout << " " << std::endl;

    //STEP2 Using the set of angles, the volume is projected
    MetaData md,mdOut;//,mdOut2,tempMd2; //,mdProjMatch;
    FileName fnVol, fnParams, fnOut;//, fnFSC, fnOut2; //, fnMdProj2;

    fnVol = fnInit;
    fnOut = fnDir+"/projected_images.stk";
    fnParams = fnDir +"/params";

    //Projection in the direction of a set of angles defined in params
    sprintf(buffer, "xmipp_phantom_project -i %s --params %s -o %s", fnVol.c_str() , fnParams.c_str() , fnOut.c_str());
    system(buffer);

    std::cout << " " << std::endl;
    std::cout << "STEP2 ENDED" << std::endl;
    std::cout << " " << std::endl;


    //STEP3 A new set of angles is gotten thought significant
    FileName fnMd_gallery;
    FileName fnOut_gallery;
    ///////////////////////////////////////
    fnOut_gallery = fnDir+"/gallery_alignment";
    ///////////////////////////////////////

    sprintf(buffer,"mkdir %s",fnOut_gallery.c_str());//,fnSym.c_str());
    //system(buffer);// It is necessary to create a new folder for running significant again. Due to, the significant output
    					   // will be the same name that Step1.

    sprintf(buffer, "xmipp_reconstruct_significant -i %s --initvolumes %s --odir %s  --sym %s --iter 1 --alpha0 %s", fnOut.c_str(),fnInit.c_str(),fnOut_gallery.c_str(),fnSym.c_str(),fsig.c_str());
    //system(buffer);

    std::cout << " " << std::endl;
    std::cout << "STEP3 ENDED" << std::endl;
    std::cout << " " << std::endl;


    //STEP4 Calculus of H0 and Hk(using projected data), then P is obtained
    //Files inside gallery_alignment folder
    //FileName fnDir_gallery="gallery_alignment";
    FileName fnDir_gallery="gallery_alignment";
    //FileName fnMd_gallery; //, fnMdProj_signif;
    //fnMd_gallery = fnDir_gallery+"/angles_iter001_00.xmd";
    fnOut_gallery = fnDir+"/"+fnDir_gallery+"/angles_iter001_00.xmd";

    //Calculus of H0 and Hk (using experimental data)
    //standard files
    FileName fnMd;//, fnMdProj;
    fnMd = fnDir+"/angles_iter001_00.xmd";

    MetaData outtput;

    P_calculus(fnOut_gallery,fnMd, outtput);

    std::cout << " " << std::endl;
    std::cout << "FINISHED" << std::endl;
    std::cout << " " << std::endl;

}


void MultireferenceAligneability::P_calculus(FileName &fnMd_gallery,FileName &fnMd,MetaData &mdOut)
{
	randomize_random_generator();

	MetaData md, md_gallery, mdProjMatch, mdProjMatch_gallery;
	size_t nSamplesRandom = 100;
	size_t maxNImg;
	FileName fnOut;
	fnOut = fnDir+"/clusteringTendency.xmd";

	md_gallery.read(fnMd_gallery);
	md.read(fnMd);

	size_t sz_gallery = md_gallery.size();
	size_t sz = md.size();

	md_gallery.getValue(MDL_IMAGE_IDX,maxNImg,sz_gallery);
	md.getValue(MDL_IMAGE_IDX,maxNImg,sz);

	String expression, expression2;
	MDRow row, row_gallery, row2, row2_gallery;

	init_progress_bar(maxNImg);

	MDIterator __iter(mdProjMatch);
	//MDIterator __iter(mdProjMatch_signif);

	SymList SL;
	int symmetry, sym_order;
	SL.readSymmetryFile(fnSym.c_str());
	SL.isSymmetryGroup(fnSym.c_str(), symmetry, sym_order);

	double non_reduntant_area_of_sphere = SL.nonRedundantProjectionSphere(symmetry,sym_order);
	double area_of_sphere_no_symmetry = 4.*PI;
	double correction = std::sqrt(non_reduntant_area_of_sphere/area_of_sphere_no_symmetry);
	double validation = 0;

	ProgValidationNonTilt PVN;

	for (size_t i=0; i<=maxNImg;i++)
	{
		std::cout << "Iteration" << i << std::endl;
		MetaData tempMd, tempMd_gallery;
		std::vector<double> sum_u1(nSamplesRandom);
		std::vector<double> sum_u2(nSamplesRandom);
		std::vector<double> sum_w1(nSamplesRandom);
		std::vector<double> sum_w2(nSamplesRandom);
		std::vector<double> H0(nSamplesRandom);
		std::vector<double> H(nSamplesRandom);

		expression = formatString("imageIndex == %lu",i);
		tempMd.importObjects(md, MDExpression(expression));
		tempMd_gallery.importObjects(md_gallery, MDExpression(expression));

		if ( (tempMd.size()==0) || (tempMd_gallery.size()==0))
			continue;

		calc_sumu(tempMd_gallery,sum_u1);
		calc_sumu(tempMd,sum_u2);
		PVN.obtainSumW(tempMd_gallery,sum_w1,sum_u1,H0,correction);
		PVN.obtainSumW(tempMd,sum_w2,sum_u2,H,correction);

		std::sort(H0.begin(),H0.end());
		std::sort(H.begin(),H.end());

		double P = 0.;
		for(size_t j=0; j<sum_u1.size();j++)  //NOTA, u1 and u2 have 100 elements, then the condition sum_u1.size, does not matter
			P += H0.at(j)/H.at(j);

		P /= (nSamplesRandom);
		row.setValue(MDL_IMAGE_IDX,i);
		row.setValue(MDL_WEIGHT,P);
		mdOut.addRow(row);

		std::cout << mdOut << std::endl;

		sum_u1.clear();
		sum_u2.clear();
		sum_w1.clear();
		sum_w2.clear();
		H0.clear();
		H.clear();
		tempMd.clear();
		progress_bar(i+1);
		__iter.moveNext();
	}

	mdOut.write(fnOut);
}


void MultireferenceAligneability::write_projection_file()
{
	FileName filnam=fnDir+"/params";
	std::ofstream myfile;
	myfile.open(filnam.c_str());
	myfile << "# XMIPP_STAR_1 *\n";
	myfile << "# \n";
	myfile << "data_block1 \n";
	myfile << "_dimensions2D   '101 101' \n";
	myfile << "_projAngleFile "+fnDir+"/images_iter001_00.xmd \n";
	myfile << "_noisePixelLevel   '0 4' \n";    //where the variance is the first number and the second one is the mean
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

