#include "unitCell.h"

// 60 is the angle between the vectors
// _5f_to_2f and _5f_2fp
//that is hte vector that joing a vertex with the two closes 2-fold
// symmetry axis
#define tg60   1.73205081//tg60 = -tg120
#define sin60  0.8660254//sin60=sin120
//#define DEBUG
#ifdef DEBUG
	#define scale 300//780-use to scale chimera bild files so they fit your actual 3Dmap
#endif
UnitCell::UnitCell(String sym, double rmin, double rmax, double expanded,
		double offset) {
	this->rmin = rmin; //delete voxels closer to the center than this radius
	this->rmax = rmax; //delete voxels more far away to the center than this radius
	this->offset = offset; //rotate unit cell this degrees
	SL.isSymmetryGroup(sym, symmetry, sym_order); //parse symmetry string
	if (symmetry == pg_I2) {
		// golden = (1 + sqrt(5)) / 2
		// 5f  = [+1.000, 0, golden] //one vertex
		// 5f' = [-1.000, 0, golden] //another vertex
		// 5f''= [ 0.,    -golden, +1.000]//the third one needed to create a face
		// 2f = (5f + 5f') / 2 // 2 fold symmetry axis
		// 2f' = (5f'' + 5f) /2 // another2 fold symmetry axis
		// 3f = (5f + 5f' +5f'')/3 // 3 fold symmetry axis
		////unit cell is enclosed by planes
		////Points defining plane 1, center 5f 2f
		////Points defining plane 2, center 2f 3f
		////Points defining plane 3, center 3f 2f'
		////Points defining plane 4, center 2f', 5f
		//All this may be precomputed but
		// It is easier to understand this way
		double golden = (1. + sqrt(5.)) / 2.;
		Matrix1D<double> _5f = vectorR3(1., 0., golden);
		_5f.selfNormalize();
		Matrix1D<double> _5fp = vectorR3(0., golden, +1.);
		_5fp.selfNormalize();
		Matrix1D<double> _5fpp = vectorR3(-1., 0., golden);
		_5fpp.selfNormalize();
		Matrix1D<double> _2f = (_5f + _5fp) / 2.;
		Matrix1D<double> _2fp = (_5fpp + _5f) / 2.;
		Matrix1D<double> _3f = (_5f + _5fp + _5fpp) / 3.;
		//vectors that join a symmetry axis with the next only.
		//all the points defined by the vectors are in the same
		// plane
		//The plane is defined by _5f_to_2f x _2f_to_5f
		// where x is the vector product
		Matrix1D<double> _5f_to_2f = _2f - _5f;
		_5f_to_2f.selfNormalize();
		Matrix1D<double> _2f_to_3f = _3f - _2f;
		_2f_to_3f.selfNormalize();
		Matrix1D<double> _3f_to_2fp = _2fp - _3f;
		_3f_to_2fp.selfNormalize();
		Matrix1D<double> _2fp_to_5f = _5f - _2fp;
		_2fp_to_5f.selfNormalize();

		// vector perpendicular to the triangle face.
		Matrix1D<double> planeVector = vectorProduct(_5f_to_2f, _2fp_to_5f);

		//vectors perpendicular to the unit cell edges
		//the "positive version" are the same vectors
		//but pointing in a direction that makes
		//easier to know if a vector is enclosed by the poly-
		//hedron that defines the unit cell
		vectExpansion.push_back(
				expanded * vectorProduct(planeVector, _5f_to_2f));
		vectExpansion.push_back(
				expanded * vectorProduct(planeVector, _2f_to_3f));
		vectExpansion.push_back(
				expanded * vectorProduct(planeVector, _3f_to_2fp));
		vectExpansion.push_back(
				expanded * vectorProduct(planeVector, _2fp_to_5f));

		//vertex for an expanded unitCell
		//60 is the angle between _5f_to_2f and _2fp_to_5f
		// sin60 * expanded ==  vectExpansion[0].module()
		//vectExpansion[0].module()==vectExpansion[3].module()
		//1) it takes a while but this is the posiction of the
		//   expanded _5f
		//1/t60*vectExpansion[0].module()*_5f_to_2f
		//1/s60*vectExpansion[3].module()*_5f_to_2f
		double mod = (1. / tg60 + 1. / sin60) * sin60 * expanded;
		expandedUnitCell.push_back(_5f + vectExpansion[0] + (-_5f_to_2f) * mod);
		// the next vertex is easy because the angle between the
		// expanded vector and the edge is 90º
		expandedUnitCell.push_back(_2f + vectExpansion[0] + vectExpansion[1]);
		mod = ((-1.) / tg60 + 1. / sin60) * sin60 * expanded;
		// 3fold symmetry axis
		expandedUnitCell.push_back(_3f + vectExpansion[1] + (_2f_to_3f) * mod);
		// _2fp symmetry axis
		expandedUnitCell.push_back(_2fp + vectExpansion[2] + vectExpansion[3]);
#include "chimeraTesterI2.txt" //draws all vectors using chimera

		//normalize expansion vector to mod_1
		for (std::vector<Matrix1D<double> >::iterator it =
				vectExpansion.begin(); it != vectExpansion.end(); ++it) {
			(*it).selfNormalize();
		}
		//vectors normal to faces of the expanded polyhedra
		for (int i = 0; i < 4; i++) {
			planeVectors.push_back(
					vectorProduct(expandedUnitCell[i],
							expandedUnitCell[(i + 1) % 4]));
		}
	} else
		REPORT_ERROR(ERR_ARG_INCORRECT, "Symmetry not implemented");
#ifdef DEBUG
	{
		std::ofstream testFile;
		testFile.open("ico1.bild");

		testFile << ".color blue\n";
		Matrix1D<double> t;
		Matrix1D<double> tt;
		for (int i = 0; i < 4; i++) {
			t = expandedUnitCell[i];
			tt = expandedUnitCell[i] + planeVectors[i];
			t *= scale;
			tt *= scale;
			testFile << ".arrow " << t(0) << " " << t(1) << " " << t(2) << " "
					<< tt(0) << " " << tt(1) << " " << tt(2) << " " << .011*scale <<"\n";
		}
		testFile.close();
	}
#endif
}

void UnitCell::maskUnitCell(ImageGeneric & in3Dmap,
		ImageGeneric & out3DDmap) {
	//1) get dimensions
	size_t xDim, yDim, zDim;
	in3Dmap.getDimensions(xDim, yDim, zDim);
	if (rmax == 0)
		rmax = (double) xDim / 2.;
#ifdef DEBUG1
	{
		std::ofstream testFile;
		testFile.open("ico2.bild");

		Matrix1D<double> t;
		double th = 1.;
		testFile << ".scale 1\n.color blue\n";
		Matrix1D<double> tt;
		for (int i = 0; i < 4; i++) {
			t = expandedUnitCell[i] * scale ;
			tt = expandedUnitCell[i] * scale + planeVectors[i] * scale;
			testFile << ".v " << t(0) << " " << t(1) << " " << t(2) << " "
					<< tt(0) << " " << tt(1) << " " << tt(2) << " 1\n";
		}
		testFile.close();
	}
#endif

	//2) get expanded unitcell enclosing box
	if (symmetry == pg_I2) {
		double minX = rmax;
		double minY = rmax;
		double minZ = rmax;
		double maxX = rmin;
		double maxY = rmin;
		double maxZ = rmin;
		Matrix1D<double> minVector, maxVector;
		for (std::vector<Matrix1D<double> >::iterator it =
				expandedUnitCell.begin(); it != expandedUnitCell.end(); ++it) {
			(*it).selfNormalize();
			minVector = (*it) * rmin;
			maxVector = (*it) * std::min(rmax, (double) xDim / 2.);
			expandedUnitCellMin.push_back(minVector);
			expandedUnitCellMax.push_back(maxVector);

			minX = std::min(minX, minVector(0));
			minX = std::min(minX, maxVector(0));

			minY = std::min(minY, minVector(1));
			minY = std::min(minY, maxVector(1));

			minZ = std::min(minZ, minVector(2));
			minZ = std::min(minZ, maxVector(2));

			maxX = std::max(maxX, minVector(0));
			maxX = std::max(maxX, maxVector(0));

			maxY = std::max(maxY, minVector(1));
			maxY = std::max(maxY, maxVector(1));

			maxZ = std::max(maxZ, minVector(2));
			maxZ = std::max(maxZ, maxVector(2));
		}
#ifdef DEBUG
		//draw real unitcell
		{
#include <iostream>
#include <fstream>
			//Debug chimera file
			std::ofstream testFile;

			//Reference points: 5f, 3f and 2f
			testFile.open("ico2.bild");

			Matrix1D<double> t, tt;
			double th = 1.;
			testFile << ".scale 1\n.color cyan\n";
			for (int i = 0; i < 4; i++) {
				t = expandedUnitCellMin[i];
				tt = expandedUnitCellMax[i];
				testFile << ".cylinder " << t(0) << " " << t(1) << " " << t(2)
						<< " " << tt(0) << " " << tt(1) << " " << tt(2) << " "
						<< th << "\n";
			}
			for (int i = 0; i < 3; i++) {
				t = expandedUnitCellMax[i];
				tt = expandedUnitCellMax[i + 1];
				testFile << ".cylinder " << t(0) << " " << t(1) << " " << t(2)
						<< " " << tt(0) << " " << tt(1) << " " << tt(2) << " "
						<< th << "\n";
			}
			t = expandedUnitCellMax[3];
			tt = expandedUnitCellMax[0];
			testFile << ".cylinder " << t(0) << " " << t(1) << " " << t(2)
					<< " " << tt(0) << " " << tt(1) << " " << tt(2) << " " << th
					<< "\n";
			for (int i = 0; i < 3; i++) {
				t = expandedUnitCellMin[i];
				tt = expandedUnitCellMin[i + 1];
				testFile << ".cylinder " << t(0) << " " << t(1) << " " << t(2)
						<< " " << tt(0) << " " << tt(1) << " " << tt(2) << " "
						<< th << "\n";
			}
			t = expandedUnitCellMin[3];
			tt = expandedUnitCellMin[0];
			testFile << ".cylinder " << t(0) << " " << t(1) << " " << t(2)
					<< " " << tt(0) << " " << tt(1) << " " << tt(2) << " " << th
					<< "\n";
			testFile.close();
			//draw box
			testFile.open("ico3.bild");
			t = vectorR3(minX, minY, minZ);
			tt = vectorR3(maxX, maxY, maxZ);
			testFile << ".color red\n";
			testFile << ".box " << t(0) << " " << t(1) << " " << t(2) << " "
					<< tt(0) << " " << tt(1) << " " << tt(2) << "\n";
			testFile.close();
		}
#endif

		in3Dmap.data->im;
		MultidimArray<float> * map;
		MultidimArray<float> * imageMap2;
		in3Dmap().getMultidimArrayPointer(map);
		out3DDmap.setDatatype(DT_Float);
		out3DDmap.resize(xDim, yDim, zDim, 1);
		out3DDmap().getMultidimArrayPointer(imageMap2);
		imageMap2->setXmippOrigin();
		int r2;
		int rmin2 = rmin * rmin;
		int rmax2 = rmax * rmax;
		int iMinZ, iMinY, iMinX, iMaxZ, iMaxY, iMaxX;
		iMinZ = (int) ceil(minZ);
		iMaxZ = (int) floor(maxZ);
		iMinY = (int) ceil(minY);
		iMaxY = (int) floor(maxY);
		iMinX = (int) ceil(minX);
		iMaxX = (int) floor(maxX);
		double dotproduct;
		bool doIt = false;

		//for all points in the enclosing box
		//check if they are inside the expanded unit cell
		for (int k = iMinZ; k <= iMaxZ; ++k)
			for (int i = iMinY; i <= iMaxY; ++i)
				for (int j = iMinX; j <= iMaxX; ++j) {
					r2 = (i * i + j * j + k * k);
					if (r2 > rmin2 && r2 < rmax2) {
						doIt = true;
						for (std::vector<Matrix1D<double> >::iterator it1 =
								planeVectors.begin(); it1 != planeVectors.end();
								++it1) {
							dotproduct = dotProduct(*it1,
									vectorR3((double) j, (double) i,
											(double) k));
							if (dotproduct < 0) {
								doIt = false;
								break;
							}
						}
						if (doIt)
							A3D_ELEM(*imageMap2,k,i,j) = A3D_ELEM(*map, k, i, j);
					}
				}
		imageMap2->selfWindow(iMinZ, iMinY, iMinX,
				iMaxZ, iMaxY, iMaxX, 0.);
		MDRow MD;
		out3DDmap.setDataMode(_DATA_ALL);
		MD.setValue(MDL_SHIFT_X, -(double)iMinX);
		MD.setValue(MDL_SHIFT_Y, -(double)iMinY);
		MD.setValue(MDL_SHIFT_Z, -(double)iMinZ);
		MD.setValue(MDL_SHIFT_X, -(double)iMinX);
		MD.setValue(MDL_SHIFT_Y, -(double)iMinY);
		MD.setValue(MDL_SHIFT_Z, -(double)iMinZ);
		double sampling;
		in3Dmap.image->MDMainHeader.getValue(MDL_SAMPLINGRATE_X, sampling);
		out3DDmap.image->MDMainHeader.setValue(MDL_SAMPLINGRATE_X, sampling);

		in3Dmap.image->MDMainHeader.getValue(MDL_SAMPLINGRATE_Y, sampling);
		out3DDmap.image->MDMainHeader.setValue(MDL_SAMPLINGRATE_Y, sampling);

		in3Dmap.image->MDMainHeader.getValue(MDL_SAMPLINGRATE_Z, sampling);
		out3DDmap.image->MDMainHeader.setValue(MDL_SAMPLINGRATE_Z, sampling);

		out3DDmap.image->setGeo(MD, 0);
	} else
		REPORT_ERROR(ERR_ARG_INCORRECT, "Symmetry not implemented");

}
