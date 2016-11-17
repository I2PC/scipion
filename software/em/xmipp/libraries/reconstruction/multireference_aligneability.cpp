/***************************************************************************
 * Authors:     Javier Vargas (jvargas@cnb.csic.es) (2016)
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


MultireferenceAligneability::MultireferenceAligneability()
{
	rank=0;
	Nprocessors=1;
}

void MultireferenceAligneability::readParams()
{
	fnParticles = getParam("-i");
	fnParticlesRef = getParam("-i2");
	fnInit = getParam("--volume");
    fin = getParam("--angles_file");
    finRef = getParam("--angles_file_ref");
    fnGallery = getParam("--gallery");
    fnSym = getParam("--sym");
    fnDir = getParam("--odir");
    donNotUseWeights= checkParam("--dontUseWeights");
    check_mirror = false;
    check_mirror = checkParam("--check_mirrors");
}

void MultireferenceAligneability::defineParams()
{
    addUsageLine("This function takes a volume and a set of projections with orientations. The volume is projected into the set of projection directions defining the "
    		"the reference projections. Thus, using the projections and references, calculation of the particle alignment precision and accuracy is carried out");
    addParamsLine("  [ -i <md_file> ]  : Input metadata with images and their angles to be validated against the volume");
    addParamsLine("  [ -i2 <md_file> ] : Input metadata of reference particles obtained projecting the volume at the same orientations that the input particles");
    addParamsLine("  [--volume <md_file=\"\">]    : Input volume to be validated");
    addParamsLine("  [--angles_file <file=\".\">]     : Input metadata with particle projections at several orientations from (usually) Significant");
    addParamsLine("  [--angles_file_ref <file=\".\">] : Input reference metadata with projections and orientations obtained projecting the volume at the several orientations from (usually) Significant");
    addParamsLine("  [--gallery <file=\".\">]		  : Reference Gallery of projections ");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections"); //TODO the input will be two doc files one from the exp and the other from refs
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--check_mirrors] 			  : Correct for mirrors in the alignment precision and accuracy estimation. In this case the precision of the projection axis without direction is computed");
    addParamsLine("  [--dontUseWeights]           : Do not use the particle weigths in the clusterability calculation ");

}

void MultireferenceAligneability::run()
{
	//xmipp_multireference_aligneability --volume 1BRD.vol --sym c3 --odir testMultiReference/ --angles_file testMultiReference/angles_iter001_00.xmd --angles_file_ref testMultiReference/gallery_alignment/angles_iter001_00.xmd &
    randomize_random_generator();

    MetaData mdOutQ;
	MetaData mdExp, mdExpSort, mdProj, mdGallery, mdInputParticles, mdInputParticlesRef;
	size_t maxNImg;
	FileName fnOutCL, fnOutQ;
	fnOutCL = fnDir+"/pruned_particles_alignability.xmd";
	fnOutQ = fnDir+"/validationAlignability.xmd";

    int symmetry, sym_order;
    SL.readSymmetryFile(fnSym.c_str());

	mdInputParticles.read(fnParticles);
	mdInputParticlesRef.read(fnParticles);
	mdProj.read(finRef);
	mdExp.read(fin);
	mdGallery.read(fnGallery);

	mdExpSort.sort(mdExp,MDL_IMAGE_IDX,true,-1,0);
	size_t sz = mdExp.size();
	mdExpSort.getValue(MDL_IMAGE_IDX,maxNImg,sz);

	String expression;
	MDRow row,rowInput,rowInputRef;

	double validationAlignabilityPrecision, validationAlignabilityAccuracy, validationAlignability, validationMirror;
	validationAlignabilityPrecision = 0;
	validationAlignabilityAccuracy = 0;
	validationAlignability = 0;
	validationMirror = 0;

	if (rank==0)
		init_progress_bar(maxNImg);

	MetaData tempMdExp, tempMdProj;
	double sum_w_exp;
	double sum_w_proj;
	double sum_noise;
	double error_mirror_exp;
	double error_mirror_proj;

	sum_w_exp = 0;
	sum_w_proj = 0;
	sum_noise = 0;
	error_mirror_exp = 0;
	error_mirror_proj = 0;

	expression = formatString("imageIndex == %lu",maxNImg);
	tempMdExp.importObjects(mdExp, MDExpression(expression));
	size_t numProjs = tempMdExp.size();
	tempMdExp.clear();

	//Noise
	calc_sumw2(numProjs, sum_noise, mdGallery);

	double rankPrec = 0.;
	double rankAccNoMirror = 0.;
	double rankAccMirror = 0.;

	double accuracy = 0.;
	double accuracyRef = 0.;
	double accuracyMirror = 0.;
	double accuracyMirrorRef = 0.;

	char hold;
	FileName imagePath;

	for (size_t i=0; i<=maxNImg;i++)
	{
		if ((i+1)%Nprocessors==rank)
		{
			expression = formatString("imageIndex == %lu",i);
			tempMdExp.importObjects(mdExp, MDExpression(expression));
			tempMdProj.importObjects(mdProj, MDExpression(expression));
			mdInputParticles.getRow(rowInput,i+1);
			mdInputParticles.getRow(rowInputRef,i+1);


			if ( (tempMdExp.size()==0) || (tempMdProj.size()==0))
				continue;

			calc_sumu(tempMdExp, sum_w_exp, error_mirror_exp);
			calc_sumu(tempMdProj, sum_w_proj, error_mirror_proj);

			obtainAngularAccuracy(tempMdExp, rowInput, accuracy, accuracyMirror);
			obtainAngularAccuracy(tempMdProj, rowInputRef, accuracyRef, accuracyMirrorRef);

#ifdef DEBUG

			std::cout << " " << std::endl;
			std::cout << "accuracy "    << accuracy << std::endl;
			std::cout << "accuracyRef " << accuracyRef << std::endl;
			std::cout << "accuracy2f " << accuracyRef/accuracy << std::endl;



#endif
#undef DEBUG

			rankPrec = 1/(sum_w_proj-sum_noise)*(sum_w_exp-sum_noise);
			rankAccMirror = 1/(accuracyRef-sum_noise)*(accuracy-sum_noise);
			rankAccNoMirror = 1/(accuracyMirrorRef-sum_noise)*(accuracyMirror-sum_noise);

			tempMdExp.getValue(MDL_IMAGE,imagePath,1);
			rowInput.setValue(MDL_IMAGE,imagePath);
			rowInput.setValue(MDL_IMAGE_IDX,i);
			rowInput.setValue(MDL_SCORE_BY_ALIGNABILITY_PRECISION, rankPrec);
			rowInput.setValue(MDL_SCORE_BY_ALIGNABILITY_ACCURACY, rankAccMirror);
			rowInput.setValue(MDL_SCORE_BY_MIRROR, rankAccNoMirror);
			rowInput.setValue(MDL_SCORE_BY_ALIGNABILITY_PRECISION_EXP,sum_w_exp);
			rowInput.setValue(MDL_SCORE_BY_ALIGNABILITY_PRECISION_REF,sum_w_proj);
			rowInput.setValue(MDL_SCORE_BY_ALIGNABILITY_ACCURACY_EXP,accuracy);
			rowInput.setValue(MDL_SCORE_BY_ALIGNABILITY_ACCURACY_REF,accuracyRef);
			rowInput.setValue(MDL_SCORE_BY_ALIGNABILITY_NOISE,sum_noise);

			mdPartialParticles.addRow(rowInput);

			rowInput.clear();
			rowInputRef.clear();
			tempMdExp.clear();
			tempMdProj.clear();

			if (rank==0)
				progress_bar(i+1);
		}
	}

	synchronize();
	gatherResults();


	if (rank == 0)
	{
		mdPartialParticles.write(fnOutCL);
		row.clear();

		FOR_ALL_OBJECTS_IN_METADATA(mdPartialParticles)
		{
			mdPartialParticles.getValue(MDL_SCORE_BY_ALIGNABILITY_PRECISION,rankPrec,__iter.objId);
			mdPartialParticles.getValue(MDL_SCORE_BY_ALIGNABILITY_ACCURACY,rankAccMirror,__iter.objId);
			mdPartialParticles.getValue(MDL_SCORE_BY_MIRROR,rankAccNoMirror,__iter.objId);

			validationAlignabilityPrecision += (rankPrec>0.5);
			validationAlignabilityAccuracy += (rankAccMirror > 0.5);
			validationAlignability += ( (rankAccMirror > 0.5) && (rankPrec>0.5));
			validationMirror += (rankAccNoMirror> 0.5);

		}

		validationAlignabilityPrecision /= (maxNImg+1);
		validationAlignabilityAccuracy /= (maxNImg+1);
		validationAlignability /= (maxNImg+1);
		validationMirror /= (maxNImg+1);

		row.setValue(MDL_IMAGE,fnInit);
		row.setValue(MDL_WEIGHT_PRECISION_ALIGNABILITY,validationAlignabilityPrecision);
		row.setValue(MDL_WEIGHT_ACCURACY_ALIGNABILITY,validationAlignabilityAccuracy);
		row.setValue(MDL_WEIGHT_ALIGNABILITY,validationAlignability);
		row.setValue(MDL_WEIGHT_PRECISION_MIRROR,validationMirror);

		mdOutQ.addRow(row);
	    mdOutQ.write(fnOutQ);

	}

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

void MultireferenceAligneability::calc_sumu(const MetaData & tempMd, double & sum_W, double & mirrorProb)
{
    double a;
    double rotRef,tiltRef,psiRef, wRef;
    double rot,tilt,psi;
    double w2;
    double W;
    double sumW;

    W = 0;
    sumW = 0;
    sum_W = 0;
    mirrorProb = 0;

    FOR_ALL_OBJECTS_IN_METADATA(tempMd)
    {
        tempMd.getValue(MDL_ANGLE_ROT,rotRef,__iter.objId);
        tempMd.getValue(MDL_ANGLE_TILT,tiltRef,__iter.objId);
        tempMd.getValue(MDL_ANGLE_PSI,psiRef,__iter.objId);
        tempMd.getValue(MDL_MAXCC,wRef,__iter.objId);

        _FOR_ALL_OBJECTS_IN_METADATA2(tempMd)
        {
        	tempMd.getValue(MDL_ANGLE_ROT,rot,__iter2.objId);
        	tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter2.objId);
            tempMd.getValue(MDL_ANGLE_PSI,psi,__iter2.objId);
        	tempMd.getValue(MDL_MAXCC,w2,__iter2.objId);


#ifdef DEBUG
{

        std::cout << xx << " " << yy << " " << zz << " " << mirror <<  std::endl;

}
#endif
#undef DEBUG

			a = SL.computeDistance(rotRef, tiltRef, psiRef,
								   rot, tilt, psi, true,check_mirror, false);

        	if (donNotUseWeights)
        	{
        		W += a;
        		sumW += 1;
        	}

        	else
        	{
        		W += a*(wRef*w2);
        		sumW += (wRef*w2);
        	}
        }
    }

    sum_W = (W / sumW);

#ifdef DEBUG
{
        std::cout << sum_W << std::endl;
		char c;
		std::getchar();
}
#endif
#undef DEBUG

}

void MultireferenceAligneability::calc_sumw(const size_t num, double & sumw)
{
	size_t trials=500;
    double xRan,yRan,zRan;
    double x,y;
    double sumWRan;
    double * xRanArray = new double[num];
    double * yRanArray = new double[num];
    double * zRanArray  = new double[num];
    double a;
    sumw=0;
    a = 0;

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

        double WRan, tempWRan, tempW1, tempW2;
        WRan = 0;
        for (size_t nS1=0; nS1<num; nS1++)
        {
            for (size_t nS2=0; nS2<num; nS2++)
            {
                a = std::abs(std::acos(xRanArray[nS1]*xRanArray[nS2]+yRanArray[nS1]*yRanArray[nS2]+zRanArray[nS1]*zRanArray[nS2]));
                if ( (a != 0))
                    WRan += a;
            }
        }

        sumw += WRan;
    }

    sumw /= trials;
}



void MultireferenceAligneability::calc_sumw2(const size_t num, double & sumw, const MetaData & mdGallery)
{
	//In the distance calculation the distance of a point with itself gives no distance.
	const size_t numGallery= mdGallery.size()+1;
	const double trials = 500;
    double xRan,yRan,zRan;
    size_t indx;
    size_t * indxArray = new size_t[numGallery];
    double sumWRan;
    double * rotArray = new double[num];
    double * tiltArray = new double[num];
    double * psiArray  = new double[num];
    double a;
    double rot,tilt,psi,w;
    bool mirror;
    sumWRan = 0;

    if (numGallery < num)
        REPORT_ERROR(ERR_ARG_INCORRECT, "The gallery size is smaller than the number of orientations per particle. Increase the angular sampling of the gallery");

    for (size_t n=0; n<trials; n++)
    {
    	size_t currentIndx = 0;

        for (size_t i=1; i<numGallery; i++)
        	indxArray[i]=i;

    	while ( currentIndx < num )
    	{
			indx = (size_t) (double( std::rand())*(numGallery-1))/RAND_MAX+1;

			while ( (indx == 0) || (indxArray[indx] == -1) )
				indx = (size_t) (double( std::rand())*(numGallery-1))/RAND_MAX;

			indxArray[indx] = -1;

			mdGallery.getValue(MDL_ANGLE_ROT,rot,indx);
        	mdGallery.getValue(MDL_ANGLE_TILT,tilt,indx);
        	mdGallery.getValue(MDL_ANGLE_PSI,psi,indx);

            rotArray[currentIndx]  = rot;
            tiltArray[currentIndx] = tilt;
            psiArray[currentIndx]  = psi;
            currentIndx++;
    	}


#ifdef DEBUG
{
		for (int var = 0; var < num; var++)
		{
	        std::cout << rotArray[var] << " " << tiltArray[var] << " " << psiArray[var] << std::endl;
		}

		char c;
		std::getchar();
}
#endif
#undef DEBUG

		double tempW1, tempW2, temp;
        a = 0;
        for (size_t nS1=0; nS1<num; nS1++)
        {
            for (size_t nS2=0; nS2<num; nS2++)
            {
    			a += SL.computeDistance(rotArray[nS1],tiltArray[nS1],psiArray[nS1],rotArray[nS2],tiltArray[nS2],psiArray[nS2], true, check_mirror, false);
            }

        }

        sumWRan += a;
    }

    sumw = ( sumWRan / (trials*(num-1)*(num-1)) );


#ifdef DEBUG
        std::cout << "   " << std::endl;
        std::cout << " sumw  : " << sumw << std::endl;
        std::cout << "   " << std::endl;
		std::getchar();
#endif
#undef DEBUG

    delete tiltArray;
    delete rotArray;
    delete psiArray;
    delete indxArray;


}

void MultireferenceAligneability::obtainAngularAccuracy(const MetaData & tempMd, const MDRow & row, double & accuracy, double & accuracyMirror)
{
	double rot, tilt, psi, w;
	double rotAux, tiltAux, psiAux;
	double rotRef, tiltRef, psiRef, wRef;
    bool mirror;;
	double sumOfW;
	double tempAccuracy, tempAccuracyMirror;


	//First for the reference:
    row.getValue(MDL_ANGLE_ROT,rotRef);
    row.getValue(MDL_ANGLE_TILT,tiltRef);
    row.getValue(MDL_ANGLE_PSI,psiRef);
    row.getValue(MDL_MAXCC,wRef);
    row.getValue(MDL_FLIP,mirror);

    accuracyMirror = 0;
    accuracy = 0;
	tempAccuracy = 0;
	tempAccuracyMirror = 0;
	sumOfW = 0;

    FOR_ALL_OBJECTS_IN_METADATA(tempMd)
    {
        tempMd.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
        tempMd.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
        tempMd.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
        tempMd.getValue(MDL_MAXCC,w,__iter.objId);
        tempMd.getValue(MDL_FLIP,mirror,__iter.objId);

        rotAux = rot;
        tiltAux = tilt;
        psiAux = psi;

        sumOfW += w;



//        tempAccuracy = SL.computeDistance(rotRef, tiltRef, psiRef,
//        		                   rotAux, tiltAux, psiAux, false,true, false);


        tempAccuracy = SL.computeDistance(rotRef, tiltRef, psiRef,
        		rotAux, tiltAux, psiAux, true, check_mirror, false);

        tempAccuracyMirror = SL.computeDistance(rotRef, tiltRef, psiRef,
        		                   rot, tilt, psi, true,true, false);

    	accuracy += tempAccuracy*w;
        accuracyMirror += tempAccuracyMirror*w;


#ifdef DEBUG
        std::cout << " tempAccuracy : " << " " << tempAccuracy << " " << std::endl;
        std::cout << " tempAccuracyMirror: " << tempAccuracyMirror  <<  std::endl;

#endif
#undef DEBUG


#ifdef DEBUG
        std::cout << tilt << " " << rot << " " << std::endl;
        std::cout << "weight: " << w << " " <<   "arco: " << tempAccuracy/(w) << " " << std::endl;

#endif

    }

    accuracy /= (sumOfW);
    accuracyMirror /= (sumOfW);

#ifdef DEBUG
        std::cout << " accuracy : " << " " << accuracy << " " << std::endl;
        std::cout << " accuracyMirror: " << accuracyMirror  <<  std::endl;
        std::cout << "-------------"  <<  std::endl;


#endif
#undef DEBUG


#ifdef DEBUG

    std::cout << " " << std::endl;
    std::cout << accuracy << std::endl;
    std::cout << "-----------" << std::endl;

#endif

}


