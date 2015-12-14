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

#include "angular_accuracy_pca.h"
#include <math.h>

ProgAngularAccuracyPCA::ProgAngularAccuracyPCA()
{
	rank=0;
	Nprocessors=1;
}

void ProgAngularAccuracyPCA::readParams()
{
	fnPhantom = getParam("-i");
	fnNeighbours = getParam("--i2");
    fnResiduals = getParam("--oresiduals");
    fnReconstructed = getParam("--oreconstructed");
    fnOut = getParam("-o");

    newXdim = 60;
    newYdim = 60;

}

void ProgAngularAccuracyPCA::defineParams()
{
    addUsageLine("Determine the angular determination accuracy of a set of particles and a 3D reconstruction ");
    addParamsLine("  [ -i <volume_file> ]      	: Voxel volume");
    addParamsLine("  [--i2 <md_file=\"\">]    	: Metadata file with neighbour projections");
    addParamsLine("  [ -o <md_file=\"\">]    	: Metadata file with obtained weights");
    addParamsLine("  [--oresiduals <stack=\"\">]  : Suffix of the output stack for the residuals");
    addParamsLine("  [--oreconstructed <stack=\"\">]  : Suffix of output stack for the reconstructed projections");
}

void ProgAngularAccuracyPCA::run()
{
	MetaData md;
    StringVector blocks;
    getBlocksInMetaDataFile(fnNeighbours, blocks);

    phantomVol.read(fnPhantom);
    phantomVol().setXmippOrigin();

    String fnTempResiduals, fnTempReconstructed;
    fnTempResiduals = "";
    fnTempReconstructed = "";
    size_t numPCAs;

	if (rank==0)
		init_progress_bar(blocks.size());

    for (size_t i = 0; i < blocks.size(); ++i)
    {
    	if ((i+1)%Nprocessors==rank)
    	{
    		md.read((String) blocks[i].c_str()+'@'+fnNeighbours);

    		if (fnResiduals!="")
    		{
    			fnTempResiduals =  blocks[i].c_str()+fnResiduals;
    			createEmptyFile(fnTempResiduals, newXdim, newYdim, 1, md.size(), true, WRITE_OVERWRITE);
    		}

    		if (fnReconstructed!="")
    		{
    			fnTempReconstructed =  blocks[i].c_str()+fnReconstructed;
    			createEmptyFile(fnTempReconstructed, newXdim, newYdim, 1, md.size(), true, WRITE_OVERWRITE);
    		}

    		if (md.size() < 10)
    			numPCAs = md.size()/2;
    		else
    			numPCAs = 2;

    		obtainPCAs(md,fnTempResiduals,fnTempReconstructed,numPCAs);

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
		mdPartial.write(fnOut);
		progress_bar(blocks.size());
	}
}

void ProgAngularAccuracyPCA::obtainPCAs(MetaData &SF, String fnTempResiduals, String fnTempReconstructed, size_t numPCAs)
{
	size_t numIter = 200;

	pca.clear();
	size_t imgno;
	Image<double> img;
	double rot, tilt, psi;
	bool mirror;
	size_t  Xdim, Ydim, Zdim, Ndim;
	phantomVol().getDimensions(Xdim,Ydim,Zdim,Ndim);
	Matrix2D<double> proj, projRef;
	imgno = 1;
	Projection P;
	FileName image;
	MultidimArray<float> temp;
	Matrix2D<double> E;

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

		pca.addVector(temp);
		imgno++;

		#ifdef DEBUG
		{
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
		#endif
	}

	pca.subtractAvg();
	pca.learnPCABasis(numPCAs,numIter);
	pca.projectOnPCABasis(projRef);
	pca.v.clear();

	imgno = 0;
	FileName f;
	FOR_ALL_OBJECTS_IN_METADATA(SF)
	{

		int enabled;
		SF.getValue(MDL_ENABLED,enabled,__iter.objId);

		if ( (enabled==-1)  )
		{
			imgno++;
			continue;
		}

		ApplyGeoParams params;
		params.only_apply_shifts = true;
		img.readApplyGeo(SF,__iter.objId,params);
		Matrix2D<double> E;
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
		selfRotate(LINEAR, img(),-(angle*180)/3.14159 , WRAP);

#ifdef DEBUG
		{
			std::cout << E << std::endl;
			std::cout << (angle*180)/3.14159 << std::endl;
			std::cout << rot << " " << tilt << " " << psi << std::endl;
			char c;
			std::getchar();
			SF.getValue(MDL_IMAGE,f,__iter.objId);
			std::cout << f << std::endl;
			img.write("kk_exp.tif");
		}
#endif

		typeCast(img(), temp);
		selfScaleToSize(LINEAR,temp,newXdim,newYdim,1);
		temp.resize(newXdim*newYdim);
		pca.addVector(temp);
		imgno++;
	}


	pca.subtractAvg();
	pca.projectOnPCABasis(proj);
	std::vector< MultidimArray<float> > v;
	v.reserve(pca.v.size());
	for (size_t i=0; i<pca.v.size(); i++)
		v.push_back(pca.v[0]);
	pca.reconsFromPCA(proj,v);

	pca.evaluateZScore(numPCAs,numIter, false);

	size_t idx = 0;
	imgno=1;
	Image<float> res;
	FOR_ALL_OBJECTS_IN_METADATA(SF)
	{
		int enabled;
		SF.getValue(MDL_ENABLED,enabled,__iter.objId);

		if ( (enabled==-1)  )
		{
			imgno++;
			continue;
		}

		ApplyGeoParams params;
		params.only_apply_shifts = true;
		img.readApplyGeo(SF,__iter.objId,params);
		SF.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
		SF.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
		SF.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
		SF.getValue(MDL_IMAGE,image,__iter.objId);

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
		selfRotate(LINEAR, img(),-(angle*180)/3.14159 , WRAP);
		typeCast(img(), temp);
		temp.resize(newXdim*newYdim);

		res() = v[idx]-temp;
		res().resize(newXdim,newYdim);
		double stdRes=(res().computeStddev());
		double stdTemp=temp.computeStddev();
		double R2 = 1-(stdRes/stdTemp)*((stdRes/stdTemp));
		SF.setValue(MDL_SCORE_BY_PCA_RESIDUAL,R2,__iter.objId);
		SF.setValue(MDL_SCORE_BY_ZSCORE, exp(-A1D_ELEM(pca.Zscore,idx)/3.),__iter.objId);

#ifdef DEBUG
{
	*temp = v[idx];
	temp->resize(newXdim,newYdim);
	res()= (*temp);
	std::cout << E << std::endl;
	std::cout << (angle*180)/3.14159 << std::endl;
	std::cout << rot << " " << tilt << " " << psi << std::endl;
	SF.getValue(MDL_IMAGE,f,__iter.objId);
	std::cout << f << std::endl;
	res.write("kk_exp.tif");
	img.write("kk_proj.tif");
	std::cout <<R2 <<std::endl;
	char c;
	std::getchar();
}
#endif

		if (fnTempResiduals!="")
		{
			FileName fnResidual;
			fnResidual.compose(imgno,fnTempResiduals);
			SF.setValue(MDL_IMAGE_RESIDUAL,fnResidual,__iter.objId);
			res.write(fnResidual);
		}

		if (fnTempReconstructed !="")
		{
			res() = v[idx];
			res().resize(newXdim,newYdim);
			selfRotate(LINEAR, res(),(angle*180)/3.14159 , WRAP);
			FileName fnReconstruted;
			fnReconstruted.compose(imgno,fnTempReconstructed);
			SF.setValue(MDL_IMAGE_REF,fnReconstruted,__iter.objId);
			res.write(fnReconstruted);
		}

		res.clear();
		temp.clear();
		idx++;
		imgno++;
	}

#ifdef DEBUG

	std::cout << proj << std::endl;

	Image<float> tt;
	Image<double> tt2;
	pca.v[0].resize(newXdim,newYdim);

	tt() = v[0];
	tt.write("kk_phan.tif");

	tt2() = pca.PCAbasis[0];
	tt2().resize(newXdim,newYdim);
	tt2.write("kk_pca0.tif");

	tt2() = pca.PCAbasis[1];
	tt2().resize(newXdim,newYdim);
	tt2.write("kk_pca1.tif");


#endif

	v.clear();
	img.clear();
	pca.clear();

}
