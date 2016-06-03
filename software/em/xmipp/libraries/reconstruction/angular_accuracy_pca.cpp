/***************************************************************************
 * Authors:     Javier Vargas (jvargas@cnb.csic.es)
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

#include "angular_accuracy_pca.h"
#include <math.h>
#include"data/transformations.h"

ProgAngularAccuracyPCA::ProgAngularAccuracyPCA()
{
	rank=0;
	Nprocessors=1;
}

void ProgAngularAccuracyPCA::readParams()
{
	fnPhantom = getParam("-i");
	fnNeighbours = getParam("--i2");
    fnOut = getParam("-o");
	fnOutQ = fnOut.getDir()+"/validationAlignabilityAccuracy.xmd";
    newXdim = getIntParam("--dim");
    newYdim = newXdim;

	std::cout <<  newXdim << std::endl;

}

void ProgAngularAccuracyPCA::defineParams()
{
    addUsageLine("Determine the angular determination accuracy of a set of particles and a 3D reconstruction ");
    addParamsLine("  [ -i <volume_file> ]      	: Voxel volume");
    addParamsLine("  [--i2 <md_file=\"\">]    	: Metadata file with neighbour projections");
    addParamsLine("  [ -o <md_file=\"\">]    	: Metadata file with obtained weights");
    addParamsLine("  [--dim <d=-1>]             : Scale images to this size if they are larger.");
    addParamsLine("                             : Set to -1 for no rescaling");
}

void ProgAngularAccuracyPCA::run()
{
	MetaData md;
    StringVector blocks;
    getBlocksInMetaDataFile(fnNeighbours, blocks);

    phantomVol.read(fnPhantom);
    phantomVol().setXmippOrigin();

    size_t numPCAs;

	if (rank==0)
		init_progress_bar(blocks.size());

    for (size_t i = 0; i < blocks.size(); ++i)
    {
    	if ((i+1)%Nprocessors==rank)
    	{
    		md.read((String) blocks[i].c_str()+'@'+fnNeighbours);

    		if (md.size() <= 1)
    			continue;

    		else if ( (md.size() > 20) )
    			numPCAs = 3;
    		else if ( (md.size() >= 5) & (md.size() < 20) )
    			numPCAs = 2;
    		else
    			numPCAs = 1;

    		obtainPCAs(md,numPCAs);

    		MDRow row;
    		FOR_ALL_OBJECTS_IN_METADATA(md)
    		{
    			md.getRow(row,__iter.objId);
    			mdPartial.addRow(row);
    		}

			if (rank==0)
				progress_bar(i+1);
    	}
    }

	synchronize();
	gatherResults();

	if (rank == 0)
	{
		double pcaResidualProj,pcaResidualExp,pcaResidual,Zscore,temp,qResidualProj,qResidualExp,qZscore;

		qResidualProj = 0;
		qResidualExp = 0;
		qZscore = 0;

		String expression;
		size_t maxIdx;
		MetaData MDSort, tempMd,MDOut, MDOutQ;
		MDSort.sort(mdPartial,MDL_ITEM_ID,true,-1,0);
		MDSort.getValue(MDL_ITEM_ID,maxIdx,MDSort.lastObject());

		MDRow row;

		for (size_t i=0; i<=maxIdx;i++)
		{
			expression = formatString("itemId == %lu",i);
			tempMd.importObjects(MDSort, MDExpression(expression));

			if ( (tempMd.size() <= 0))
				continue;

			pcaResidualProj = -1e3;
			pcaResidualExp = -1e3;
			pcaResidual = -1e3;
			Zscore = -1e3;

			tempMd.getRow(row,tempMd.firstObject());

			FOR_ALL_OBJECTS_IN_METADATA(tempMd)
			{
				tempMd.getValue(MDL_SCORE_BY_PCA_RESIDUAL_PROJ,temp,__iter.objId);
				if (temp > pcaResidualProj)
					pcaResidualProj=temp;

				tempMd.getValue(MDL_SCORE_BY_PCA_RESIDUAL_EXP,temp,__iter.objId);
				if (temp > pcaResidualExp)
					pcaResidualExp=temp;

				tempMd.getValue(MDL_SCORE_BY_PCA_RESIDUAL,temp,__iter.objId);
				if (temp > pcaResidual)
					pcaResidual=temp;

				tempMd.getValue(MDL_SCORE_BY_ZSCORE,temp,__iter.objId);
				if (temp > Zscore)
					Zscore=temp;
			}

			qResidualProj += pcaResidualProj;
			qResidualExp  += pcaResidualExp;
			qZscore       += Zscore;

			row.setValue(MDL_SCORE_BY_PCA_RESIDUAL_PROJ,pcaResidualProj);
			row.setValue(MDL_SCORE_BY_PCA_RESIDUAL_EXP,pcaResidualExp);
			row.setValue(MDL_SCORE_BY_PCA_RESIDUAL,pcaResidual);
			row.setValue(MDL_SCORE_BY_ZSCORE,Zscore);
			MDOut.addRow(row);
			row.clear();
		}

		MDOut.write(fnOut);

		qResidualProj /= MDOut.size();
		qResidualExp /= MDOut.size();
		qZscore /= MDOut.size();

	    row.setValue(MDL_IMAGE,fnPhantom);
	    row.setValue(MDL_SCORE_BY_PCA_RESIDUAL_PROJ,qResidualProj);
	    row.setValue(MDL_SCORE_BY_PCA_RESIDUAL_EXP,qResidualExp);
	    row.setValue(MDL_SCORE_BY_ZSCORE,qZscore);

	    MDOutQ.addRow(row);
	    MDOutQ.write(fnOutQ);
		progress_bar(blocks.size());
	}
}

void ProgAngularAccuracyPCA::obtainPCAs(MetaData &SF, size_t numPCAs)
{
	size_t numIter = 200;
	pca.clear();
	size_t imgno;
	Image<double> img;
	double rot, tilt, psi;
	double shiftX, shiftY;
	bool mirror;
	size_t  Xdim, Ydim, Zdim, Ndim;
	phantomVol().getDimensions(Xdim,Ydim,Zdim,Ndim);

	if (  (newXdim == -1) )
	{
		newXdim = Xdim;
		newYdim = Ydim;
	}

	Matrix2D<double> proj;;
	imgno = 0;
	Projection P;
	FileName image;
	MultidimArray<float> temp;
	MultidimArray<double> avg;
	Matrix2D<double> E, Trans(3,3);
	Matrix1D<double> opt_offsets(2);

#ifdef DEBUG
	FOR_ALL_OBJECTS_IN_METADATA(SF)
	{
		int enabled;
		SF.getValue(MDL_ENABLED,enabled,__iter.objId);
		if ( (enabled==-1)  )
		{
			imgno++;
			continue;
		}

		SF.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
		SF.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
		SF.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
		SF.getValue(MDL_FLIP,mirror,__iter.objId);

		if (mirror)
		{
			double newrot;
			double newtilt;
			double newpsi;
			Euler_mirrorY(rot,tilt,psi,newrot,newtilt,newpsi);
			rot = newrot;
			tilt = newtilt;
			psi = newpsi;
		}

		projectVolume(phantomVol(), P, Ydim, Xdim, rot, tilt, psi);

		Euler_angles2matrix(rot, tilt, psi, E, false);
		double angle = atan2(MAT_ELEM(E,0,1),MAT_ELEM(E,0,0));
		selfRotate(LINEAR, P(),-(angle*180)/3.14159 , WRAP);
		typeCast(P(), temp);
		selfScaleToSize(LINEAR,temp,newXdim,newYdim,1);
		temp.resize(newXdim*newYdim);
		temp.statisticsAdjust(0,1);
		pca.addVector(temp);
		imgno++;

		#ifdef DEBUG
		{
			{
				size_t val;
				SF.getValue(MDL_ITEM_ID,val,__iter.objId);
				{
					std::cout << E << std::endl;
					std::cout << (angle*180)/3.14159 << std::endl;
					P.write("kk_proj.tif");
					SF.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
					std::cout << rot << " " << tilt << " " << psi << std::endl;
					char c;
					std::getchar();
				}
			}
		}
		#endif

	}
	pca.subtractAvg();
	avg = pca.avg;
	pca.learnPCABasis(numPCAs,numIter);
	//pca.projectOnPCABasis(projRef);
	pca.v.clear();

#endif

	imgno = 0;
	FileName f;

	FOR_ALL_OBJECTS_IN_METADATA(SF)
	{
		int enabled;
		SF.getValue(MDL_ENABLED,enabled,__iter.objId);
		SF.getValue(MDL_SHIFT_X,shiftX,__iter.objId);
		SF.getValue(MDL_SHIFT_Y,shiftY,__iter.objId);

		if ( (enabled==-1)  )
		{
			imgno++;
			continue;
		}

		//geo2TransformationMatrix(input, Trans);
		//ApplyGeoParams params;
		//params.only_apply_shifts = true;
		//img.readApplyGeo(SF,__iter.objId,params);

		SF.getValue(MDL_IMAGE,f,__iter.objId);
		img.read(f);
		SF.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
		SF.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
		SF.getValue(MDL_ANGLE_PSI,psi,__iter.objId);

		if (mirror)
		{
			double newrot;
			double newtilt;
			double newpsi;
			Euler_mirrorY(rot,tilt,psi,newrot,newtilt,newpsi);
			rot = newrot;
			tilt = newtilt;
			psi = newpsi;
		}

		Euler_angles2matrix(rot, tilt, psi, E, false);
		double angle = atan2(MAT_ELEM(E,0,1),MAT_ELEM(E,0,0));
		angle=-(angle*180)/3.14159;
		rotation2DMatrix(angle, Trans, true);
	    dMij(Trans, 0, 2) = shiftX;
	    dMij(Trans, 1, 2) = shiftY;

	    selfApplyGeometry(BSPLINE3, img(), Trans, IS_NOT_INV, WRAP);

#ifdef DEBUG
		{
			std::cout <<  MAT_ELEM(Trans,0,0) << " " << MAT_ELEM(Trans,0,1) << " " << MAT_ELEM(Trans,0,2) << " " << MAT_ELEM(Trans,1,0) << " " <<  MAT_ELEM(Trans,1,1) << " " << MAT_ELEM(Trans,1,2) << " " << MAT_ELEM(Trans,2,0) << " " << MAT_ELEM(Trans,2,1) << " " << MAT_ELEM(Trans,2,2) << " " << std::endl;
			std::cout <<  MAT_ELEM(E,0,0) << " " << MAT_ELEM(E,0,1) << " " << MAT_ELEM(E,0,2) << " " << MAT_ELEM(E,1,0) << " " <<  MAT_ELEM(E,1,1) << " " << MAT_ELEM(E,1,2) << " " << MAT_ELEM(E,2,0) << " " << MAT_ELEM(E,2,1) << " " << MAT_ELEM(E,2,2) << " " << std::endl;
			size_t val;
			SF.getValue(MDL_ITEM_ID,val,__iter.objId);
			if (true)
			{
				SF.getValue(MDL_IMAGE,f,__iter.objId);
				std::cout << f << std::endl;
				img.write("kk_exp.tif");

				char c;
				std::getchar();

			}
		}
#endif

		typeCast(img(), temp);
		selfScaleToSize(LINEAR,temp,newXdim,newYdim,1);
		temp.resize(newXdim*newYdim);
		temp.statisticsAdjust(0,1);
		pca.addVector(temp);
		imgno++;
	}

	pca.subtractAvg();
	pca.learnPCABasis(numPCAs,numIter);
	pca.projectOnPCABasis(proj);

	std::vector< MultidimArray<float> > recons(SF.size());
	for(int n=0; n<SF.size(); n++)
		recons[n] = MultidimArray<float>(newXdim*newYdim);

	pca.reconsFromPCA(proj,recons);
	pca.evaluateZScore(numPCAs,numIter, false);

	imgno = 0;
	Image<float> imgRes;
	double R2_Proj,R2_Exp;
	R2_Proj=0;
	R2_Exp=0;

	MultidimArray <int> ROI;
    ROI.resizeNoCopy(newYdim,newXdim);
    ROI.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( temp < (newXdim/2))
            A2D_ELEM(ROI,i,j)= 1;
        else
            A2D_ELEM(ROI,i,j)= 0;
    }

    ROI.resize(newYdim*newXdim);
    ROI.setXmippOrigin();

	FOR_ALL_OBJECTS_IN_METADATA(SF)
	{
		int enabled;
		SF.getValue(MDL_ENABLED,enabled,__iter.objId);
		if ( (enabled==-1)  )
		{
			imgno++;
			continue;
		}

		//Projected Image
		SF.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
		SF.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
		SF.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
		SF.getValue(MDL_FLIP,mirror,__iter.objId);

		if (mirror)
		{
			double newrot;
			double newtilt;
			double newpsi;
			Euler_mirrorY(rot,tilt,psi,newrot,newtilt,newpsi);
			rot = newrot;
			tilt = newtilt;
			psi = newpsi;
		}

		projectVolume(phantomVol(), P, Ydim, Xdim, rot, tilt, psi);
		Euler_angles2matrix(rot, tilt, psi, E, false);
		double angle = atan2(MAT_ELEM(E,0,1),MAT_ELEM(E,0,0));
		angle = -(angle*180)/3.14159;
		selfRotate(LINEAR, P(),angle , WRAP);
		typeCast(P(), temp);
		selfScaleToSize(LINEAR,temp,newXdim,newYdim,1);
		temp.resize(newXdim*newYdim);
		temp.statisticsAdjust(0,1);
		temp.setXmippOrigin();

		//Reconstructed Image
		recons[imgno].statisticsAdjust(0,1);
		recons[imgno].resize(newYdim*newXdim);
		recons[imgno].setXmippOrigin();

		R2_Proj = correlationIndex(temp,recons[imgno],&ROI);

		SF.getValue(MDL_SHIFT_X,shiftX,__iter.objId);
		SF.getValue(MDL_SHIFT_Y,shiftY,__iter.objId);

		SF.getValue(MDL_IMAGE,f,__iter.objId);
		img.read(f);

		rotation2DMatrix(angle, Trans, true);
		dMij(Trans, 0, 2) = shiftX;
		dMij(Trans, 1, 2) = shiftY;
		selfApplyGeometry(BSPLINE3, img(), Trans, IS_NOT_INV, WRAP);

		typeCast(img(), temp);
		selfScaleToSize(LINEAR,temp,newXdim,newYdim,1);
		temp.resize(newXdim*newYdim);
		temp.statisticsAdjust(0,1);
		temp.setXmippOrigin();

		R2_Exp = correlationIndex(temp,recons[imgno],&ROI);

		SF.setValue(MDL_SCORE_BY_PCA_RESIDUAL_PROJ,R2_Proj,__iter.objId);
		SF.setValue(MDL_SCORE_BY_PCA_RESIDUAL_EXP,R2_Exp,__iter.objId);
		SF.setValue(MDL_SCORE_BY_PCA_RESIDUAL,R2_Proj*R2_Exp,__iter.objId);
		SF.setValue(MDL_SCORE_BY_ZSCORE, exp(-A1D_ELEM(pca.Zscore,imgno)),__iter.objId);

		#ifdef DEBUG
		{
			size_t val;
			SF.getValue(MDL_ITEM_ID,val,__iter.objId);
			if (true)
			{
			Image<float>  imgRecons;
			Image<double> imgAvg;
			SF.getValue(MDL_IMAGE,f,__iter.objId);

			img.write("kk_exp0.tif");
			apply_binary_mask(ROI,temp,imgRecons());
			imgRecons().resize(newYdim,newXdim);
			imgRecons.write("kk_exp.tif");

			recons[imgno].statisticsAdjust(0,1);
			apply_binary_mask(ROI,recons[imgno],imgRecons());
			imgRecons().resize(newYdim,newXdim);
			imgRecons.write("kk_reconstructed.tif");

			imgAvg()=pca.avg;
			imgAvg().resize(newYdim,newXdim);
			imgAvg.write("kk_average.tif");

			//Projected Image
			SF.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
			SF.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
			SF.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
			SF.getValue(MDL_FLIP,mirror,__iter.objId);

			if (mirror)
			{
				double newrot;
				double newtilt;
				double newpsi;
				Euler_mirrorY(rot,tilt,psi,newrot,newtilt,newpsi);
				rot = newrot;
				tilt = newtilt;
				psi = newpsi;
			}

			projectVolume(phantomVol(), P, Ydim, Xdim, rot, tilt, psi);
			Euler_angles2matrix(rot, tilt, psi, E, false);
			double angle = atan2(MAT_ELEM(E,0,1),MAT_ELEM(E,0,0));
			selfRotate(LINEAR, P(),-(angle*180)/3.14159 , WRAP);
			P.write("kk_proj.tif");

			//std::cout <<  exp(-R2_Proj) << " " << exp(-R2_Exp) << std::endl;
			std::cout <<  R2_Proj << " " << R2_Exp <<  " " << R2_Proj*R2_Exp << std::endl;

			for(int i=0; i<numPCAs;i++)
			{
				std::cout << "proj " << MAT_ELEM(proj,i,imgno) <<  "   " <<  "projRef " <<  MAT_ELEM(proj,i,imgno) << std::endl;

			}

			char c;
			std::getchar();

			}

		}

#endif

		imgno++;
	}

	recons.clear();
	img.clear();
	temp.clear();
	imgRes.clear();
	pca.clear();

}
