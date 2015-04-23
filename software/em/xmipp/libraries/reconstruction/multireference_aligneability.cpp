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

//ejecución del programa
//xmipp_validation_nontilt --volume 4MOA.vol --odir Sig4MOA --sym c1


//xmipp_multireference_aligneability --volume 4MOA.vol  --odir Sig4MOA


#include "multireference_aligneability.h"
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
}

void MultireferenceAligneability::run()
{
    //Clustering Tendency and Cluster Validity Stephen D. Scott
    randomize_random_generator();
    char buffer[400];

    //Step1 From a volume, a set of angles (set1) is obtained via sigfinicant
    	//plates.stk Volumen
    	//4MOA.vol Volumen inicial
    sprintf(buffer, "xmipp_reconstruct_significant  -i plates.stk --initvolumes 4MOA.vol --odir %s  --sym c1 --iter 1 --alpha0 0.05", fnDir.c_str());
    system(buffer);

    writeparams();
    std::cout << " " << std::endl;
    std::cout << "STEP1 ENDED" << std::endl;
    std::cout << " " << std::endl;

    //Step2 Using the set of angles, set1, the volume is projected
    MetaData md,mdOut,mdOut2,tempMd2,mdProjMatch;
    FileName fnVol, fnParams, fnMdProj2, fnOut, fnFSC, fnOut2;

    fnVol = fnInit;//+"/4MOA.vol";
    fnOut = fnDir+"/plates2.stk";
    fnParams = fnDir +"/params.rst";


    //Projection in the direction of a set of angles defined in params
    sprintf(buffer, "xmipp_phantom_project -i %s --params %s -o %s", fnVol.c_str() , fnParams.c_str() , fnOut.c_str());
    system(buffer);

    std::cout << " " << std::endl;
    std::cout << "STEP2 ENDED" << std::endl;
    std::cout << " " << std::endl;

    //Step3 A new set of angles is gotten thought significant
    system("mkdir signif");// It is necessary to create a new folder for running significant again. Due to, the significant output
    					   // will be the same name that Step1.

    sprintf(buffer, "xmipp_reconstruct_significant  -i %s --initvolumes 4MOA.vol --odir signif  --sym c1 --iter 1 --alpha0 0.05", fnOut.c_str());
    system(buffer);

    std::cout << " " << std::endl;
    std::cout << "STEP3 ENDED" << std::endl;
    std::cout << " " << std::endl;

    //Step4 Calculus of H0 and Hk(using projected data)

    //Files inside signif folder
    FileName fnDir_signif="signif";
    FileName fnMd_signif, fnMdProj_signif;
    fnMd_signif = fnDir_signif+"/angles_iter001_00.xmd";
    fnMdProj_signif = fnDir_signif+"/images_iter001_00.xmd";

    std::cout << " " << std::endl;
    std::cout << "STEP4 ENDED" << std::endl;
    std::cout << " " << std::endl;

    //Step5 Calculus of H0 and Hk (using experimental data)
    //standard files
    FileName fnMd, fnMdProj;
    fnMd = fnDir+"/angles_iter001_00.xmd";
    fnMdProj = fnDir+"/images_iter001_00.xmd";

    MetaData outtput;
    H0andHk_calculus(fnMd_signif,fnMdProj_signif,fnMd,fnMdProj, outtput);

    std::cout << " " << std::endl;
    std::cout << "FINISHED" << std::endl;
    std::cout << " " << std::endl;

}



void MultireferenceAligneability::H0andHk_calculus(FileName &fnMd_signif,FileName &fnMdProj_signif,
		FileName &fnMd,FileName &fnMdProj, MetaData &mdOut)
{
			randomize_random_generator();

			MetaData md, md_signif, mdProjMatch, mdProjMatch_signif;
			size_t nSamplesRandom = 100;
			size_t maxNImg;


	        md_signif.read(fnMd_signif);
	        md.read(fnMd);

	        size_t sz_signif = md_signif.size();
	        size_t sz = md.size();

	        md_signif.getValue(MDL_IMAGE_IDX,maxNImg,sz_signif);
	        md.getValue(MDL_IMAGE_IDX,maxNImg,sz);

	        mdProjMatch_signif.read(fnMdProj_signif);
	        mdProjMatch.read(fnMdProj);

	        String expression, expression_signif;
	        MDRow row, row_signif, row2, row2_signif;

	        mdProjMatch_signif.addLabel(MDL_ENABLED,0);
	        mdProjMatch.addLabel(MDL_ENABLED,0);

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

	        std::ofstream myfile2;
	        myfile2.open ("Pvalues.rst");

	        for (size_t i=0; i<=maxNImg;i++)
	        {
	        	MetaData tempMd, tempMd_signif;
	            std::vector<double> sum_u(nSamplesRandom);
	            std::vector<double> sum_w(nSamplesRandom);
	            std::vector<double> H0(nSamplesRandom);
	            std::vector<double> H(nSamplesRandom);

	            expression = formatString("imageIndex == %lu",i);
	            tempMd.importObjects(md_signif, MDExpression(expression));

	            if (tempMd.size()==0)
	                continue;

	            if (tempMd_signif.size()==0)
	            	continue;

	            obtainSumU(tempMd_signif,sum_u,H0);
	            obtainSumW(tempMd,sum_w,sum_u,H);

	            std::sort(H0.begin(),H0.end());
	            std::sort(H.begin(),H.end());

	            double P = 0.;
	            for(size_t j=0; j<sum_u.size();j++)
	                  	P += H0.at(j)/H.at(j);

	            P /= (nSamplesRandom/correction);
	            row.setValue(MDL_IMAGE_IDX,i);
	            row.setValue(MDL_WEIGHT,P);
	            mdOut.addRow(row);

	            std::cout << mdOut << std::endl;

	            myfile2 << mdOut << "\n";

	            if ( P > 1)
	               	validation++;
	            else
	               	mdProjMatch.setValue(MDL_ENABLED,0,__iter.objId);


	            sum_u.clear();
	            sum_w.clear();
	            H0.clear();
	            H.clear();
	            tempMd.clear();
	            progress_bar(i+1);
	            __iter.moveNext();
	        }
	    	myfile2.close();
}


void MultireferenceAligneability::writeparams()
{
	std::ofstream myfile;
	myfile.open ("Sig4MOA/params.rst");
	myfile << "# XMIPP_STAR_1 *\n";
	myfile << "# \n";
	myfile << "data_block1 \n";
	myfile << "_dimensions2D   '101 101' \n";
	myfile << "_projAngleFile Sig4MOA/images_iter001_00.xmd \n";
	myfile << "_noisePixelLevel   '0 0' \n";

	myfile.close();

}




void MultireferenceAligneability::obtainSumU(const MetaData & tempMd,std::vector<double> & sum_u,std::vector<double> & H0)
{//Determine H0

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

        	tilt =2*PI*rnd_unif(0,1);			// Defines the tilt between 0 and 2pi
        	rot  =std::acos(2*rnd_unif(0,1)-1);     // Defines the rot, between 0 and pi. See reference:
        											//http://mathworld.wolfram.com/SpherePointPicking.html
        	xRan = sin(tilt)*cos(rot);
        	yRan = sin(tilt)*sin(rot);
        	zRan = std::abs(cos(tilt));							// Projection coordinates

            xRanArray[nS] = xRan;
            yRanArray[nS] = yRan;
            zRanArray[nS] = zRan;
            tempMd.getColumnValues(MDL_WEIGHT, weightV);  //selects the weights column

            std::random_shuffle(weightV.begin(), weightV.end()); // Arranges in a random way
        }

        sumWRan = 0;
        double WRan, tempWRan, tempW1, tempW2;
        for (size_t nS1=0; nS1<tempMd.size(); nS1++)
        {
            tempWRan = 1e3;
            for (size_t nS2=0; nS2<tempMd.size(); nS2++)
            {
                a = std::abs(std::acos(xRanArray[nS1]*xRanArray[nS2]+yRanArray[nS1]*yRanArray[nS2]+zRanArray[nS1]*zRanArray[nS2]));
                	//Calculates all possible distances (without weights)
                if ( (a<tempWRan) && (a != 0)) //searches the closest neighbor
                {
                    tempWRan = a;
                    tempW2 = weightV[nS2];
                    tempW1 = weightV[nS1];
                    WRan = a*std::exp(std::abs(tempW1-tempW2))*std::exp(-(tempW1+tempW2)); //distances weighted
                }
            }
            sumWRan += WRan; //sum of alfa_i
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
void MultireferenceAligneability::obtainSumW(const MetaData & tempMd,std::vector<double> & sum_W,std::vector<double> & sum_u,std::vector<double> & H)
{//Determines Hk
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
                W = a*std::exp(std::abs(w-w2))*std::exp(-(w+w2));  //adds weights to the distances
                tempW = a;
            }
        }
        sumW +=  W;
    }

    for (size_t n=0; n<sum_u.size(); n++)
    {
        H[n] = sumW/(sumW+sum_u.at(n));
    }
}







