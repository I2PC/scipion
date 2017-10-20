#include "unitCell.h"
#include <math.h>
// 60 is the angle between the vectors
// icosahedron _5f_to_2f and _5f_2fp: vector that joing a vertex with the two closes 2-fold symmetry axis
#define tg60   1.73205081
#define sin60  0.8660254
#define DEBUG
#ifdef DEBUG
#define scale 1.//780-use to scale chimera bild files so they fit your actual 3Dmap
#endif

void UnitCell::doCyclic(int order) { //check offset is in radians
	Matrix1D<double> _centre, _2f, _2fp;
	_centre = vectorR3(0., 0., 0.);
	_2f = vectorR3(cos(offset), sin(offset), 0.) * rmax;
	_2fp = vectorR3(cos(TWOPI / order + offset), sin(TWOPI / order + offset),
			0.) * rmax;
	_minZ = -rmax;
	_maxZ = +rmax;

	cyclicSymmetry(_centre, _2f, _2fp, order, expanded, offset);
}

void UnitCell::doDihedral(int order) {
	Matrix1D<double> _centre, _2f, _2fp;
	_centre = vectorR3(0., 0., 0.);
	_2f = vectorR3(cos(offset), sin(offset), 0.) * rmax;
	_2fp = vectorR3(cos(TWOPI / order + offset), sin(TWOPI / order + offset),
			0.) * rmax;
	_minZ = -rmax;
	_maxZ = +rmax;

	dihedralSymmetry(_centre, _2f, _2fp, order, expanded, offset);
}

void UnitCell::doTetrahedral(int symmetry) {
	Matrix1D<double> _centroid, _3f, _3fp, _3fpp;
	_centroid = vectorR3(0., 0., 0.);
	_3f = vectorR3(0., 0., 1.) * rmax;
	_3fp = vectorR3(0., 0.94281, -0.33333) * rmax;
	_3fpp = vectorR3(0.81650, -0.47140, -0.33333) * rmax;
	_minZ = -rmax;
	_maxZ = +rmax;

	tetrahedralSymmetry(_centroid, _3f, _3fp, _3fpp, expanded);
}

void UnitCell::doOctahedral(int symmetry) {
	Matrix1D<double> _centre, _4f, _4fp, _4fpp;
	_centre = vectorR3(0., 0., 0.);
	_4f = vectorR3(1., 0., 0.) * rmax;
	_4fp = vectorR3(0., 1., 0.) * rmax;
	_4fpp = vectorR3(0., 0., 1.) * rmax;
	_minZ = -rmax;
	_maxZ = +rmax;

	octahedralSymmetry(_centre, _4f, _4fp, _4fpp, expanded);
}

void UnitCell::doIcosahedral(int symmetry) {
	double phi = (1. + sqrt(5.)) / 2.;
	Matrix1D<double> _centre, _5f, _5fp, _5fpp;
	_centre = vectorR3(0., 0., 0.);
	if (symmetry == pg_I1) {
		_5f = vectorR3(0., -0.52573111, 0.85065081);
		_5fp = vectorR3(0., 0.52573111, 0.85065081);
		_5fpp = vectorR3(-8.50650808e-01, -5.55111512e-17, 5.25731112e-01);
	} else if (symmetry == pg_I2) {
		_5f = vectorR3(-0.52573111, 0., 0.85065081);
		_5fp = vectorR3(0.52573111, 0., 0.85065081);
		_5fpp = vectorR3(-5.55111512e-17, 8.50650808e-01, 5.25731112e-01);
	} else if (symmetry == pg_I3) {
		_5f = vectorR3(0., 0., -1.);
		_5fp = vectorR3(0.89442719, 0., -0.4472136);
		_5fpp = vectorR3(0.2763932, -0.85065081, -0.4472136);
	} else if (symmetry == pg_I4) {
		_5f = vectorR3(0., 0., 1.);
		_5fp = vectorR3(0.89442719, 0., 0.4472136);
		_5fpp = vectorR3(0.2763932, 0.85065081, 0.4472136);
	} else
		REPORT_ERROR(ERR_ARG_INCORRECT, "Symmetry not implemented");
	_minZ = rmax;
	_maxZ = rmin;

	icoSymmetry(_centre, _5f, _5fp, _5fpp, expanded);
}

UnitCell::UnitCell(String sym, double rmin, double rmax, double expanded,
		double offset, double sampling, double x_origin, double y_origin, double z_origin) {
	// list of numbers smaller than 10000 which have prime decomposition
	// that does not contain prime number greater than 19
	// This is the greates prime number that can handle ccp4 fft routine

	this->rmin = rmin; //delete voxels closer to the center than this radius
	this->rmax = rmax; //delete voxels more far away to the center than this radius
	this->offset = offset; //rotate unit cell this degrees
	this->expanded = expanded; //coefficient of expansion
	this->sampling = sampling; // Angstrom/pixels relation 
	this->x_origin = x_origin; // origin x coordinate introduced with the input volume
	this->y_origin = y_origin; // origin y coordinate introduced with the input volume
	this->z_origin = z_origin; // origin z coordinate introduced with the input volume
	this->newOriginAfterExpansion = vectorR3(0., 0., 0.);
	SL.isSymmetryGroup(sym, symmetry, sym_order); //parse symmetry string
	std::cerr << "symmetry: " << symmetry << std::endl;
	if (symmetry == pg_CN && sym_order != 1)
		doCyclic(sym_order);
	else if (symmetry == pg_DN)
		doDihedral(sym_order);
	else if (symmetry == pg_T)
		doTetrahedral();
	else if (symmetry == pg_O)
		doOctahedral();
	else if (symmetry >= pg_I1 && symmetry <= pg_I4)
		doIcosahedral(symmetry);

#ifdef DEBUG
	{
		std::ofstream testFile;
		testFile.open("planeVectors_in_expandedUnitCell.bild");

		testFile << ".color blue\n";
		Matrix1D<double> t;
		Matrix1D<double> tt;
		for (int i = 0; i < expandedUnitCell.size(); i++) {
			if (symmetry == pg_CN || symmetry == pg_T || symmetry == pg_O) {
				if (i == 0)
					continue;
				planeVectors[i - 1].selfNormalize();
				t = expandedUnitCell[i];
				tt = expandedUnitCell[i] + (planeVectors[i - 1] * 10.);
			} else if (symmetry == pg_DN) {
				if (i == 0) {
					t = expandedUnitCell[i];
					planeVectors[expandedUnitCell.size() - 1].selfNormalize();
					tt = expandedUnitCell[i]
							+ (planeVectors[expandedUnitCell.size() - 1] * 10.);
				} else if (i != 0) {
					t = expandedUnitCell[i];
					planeVectors[i - 1].selfNormalize();
					tt = expandedUnitCell[i] + (planeVectors[i - 1] * 10.);
				}
			} else if (symmetry >= pg_I1 || symmetry <= pg_I4) {
				t = expandedUnitCell[i] * 74;
				planeVectors[i].selfNormalize();
				tt = (expandedUnitCell[i] + (planeVectors[i] * 0.1)) * 74;
			}
			t *= scale;
			tt *= scale;
			testFile << ".arrow " << t(0) << " " << t(1) << " " << t(2) << " "
					<< tt(0) << " " << tt(1) << " " << tt(2) << " "
					<< .011 * scale * 10 << "\n";
		}
		testFile.close();
	}
#endif

}

//FUNCTION cyclicSymmetry fills out the variable planeVectors
//(normal vectors to the two flat faces that define a unit cell)
void UnitCell::cyclicSymmetry(const Matrix1D<double> & _centre,
		const Matrix1D<double> & _2f, const Matrix1D<double> & _2fp, int order,
		double expanded, double offset) {
	Matrix1D<double> _centre_to_2f = _2f - _centre;
	Matrix1D<double> _centre_to_2fp = _2fp - _centre;
	Matrix1D<double> _2f_to_2fp = _2fp - _2f;
	_2f_to_2fp.selfNormalize();
	// vector perpendicular to the triangle face, oriented in the positive sense of z axis
	Matrix1D<double> planeVector;
	if (order == 2) {
		planeVector = vectorR3(0., 0., 1.);
	} else if (order > 2) {
		planeVector = vectorProduct(_centre_to_2f, _2f_to_2fp);
	}
	planeVector.selfNormalize();

	//vectExpansion: vectors perpendicular to the unit cell faces
	Matrix1D<double> vp_0 = vectorProduct(_centre_to_2f, planeVector);
	vp_0.selfNormalize();
	vectExpansion.push_back(expanded * vp_0 * rmax / 2);
	Matrix1D<double> vp_1 = vectorProduct(planeVector, _centre_to_2fp);
	vp_1.selfNormalize();
	vectExpansion.push_back(expanded * vp_1 * rmax / 2);

	//computation of coordinates for expandedUnitCell vertices
	Matrix1D<double> d_2f = (expanded / sin(TWOPI / order + offset)) * (-1)
			* _2fp;
	Matrix1D<double> d_2fp = (expanded / sin(TWOPI / order + offset)) * (-1)
			* _2f;

	if (expanded == 0) {
		expandedUnitCell.push_back(_centre);
		expandedUnitCell.push_back(_2f);
		expandedUnitCell.push_back(_2fp);
	} else if (expanded >= 0) {
		//_centre expands to expandedUnitCell[0], which is equivalent to newOriginAfterExpansion
		if (order == 2) {
			expandedUnitCell.push_back(_centre + vectExpansion[0]);
		} else if (order > 2) {
			expandedUnitCell.push_back(_centre + d_2f + d_2fp);
		}
		//_2f expands to expandedUnitCell[1]
		expandedUnitCell.push_back(_2f + vectExpansion[0]);
		//_2fp expands to expandedUnitCell[2]
		expandedUnitCell.push_back(_2fp + vectExpansion[1]);
	}
	newOriginAfterExpansion = expandedUnitCell[0];

	//computation of coordinates from expandedUnitCells regarding to the newOriginAfterExpansion
	Matrix1D<double> new_centre_to_new_2f = expandedUnitCell[1]
			- newOriginAfterExpansion;
	Matrix1D<double> new_centre_to_new_2fp = expandedUnitCell[2]
			- newOriginAfterExpansion;
	//vectors normal to faces of the expanded unit cell
	planeVectors.push_back(
			(-1) * vectorProduct(new_centre_to_new_2f, planeVector));
	planeVectors.push_back(
			(-1) * vectorProduct(planeVector, new_centre_to_new_2fp));
#include "chimeraTesterC.txt" //draws the (expanded) unit cell and directions of vectExpansion vectors using chimera
}
//FUNCTION dihedralSymmetry fills out the variable planeVectors
//(normal vectors to the three flat faces that define a unit cell)
void UnitCell::dihedralSymmetry(const Matrix1D<double> & _centre,
		const Matrix1D<double> & _2f, const Matrix1D<double> & _2fp, int order,
		double expanded, double offset) {
	Matrix1D<double> _centre_to_2f = _2f - _centre;
	Matrix1D<double> _centre_to_2fp = _2fp - _centre;
	Matrix1D<double> _2f_to_2fp = _2fp - _2f;
	_2f_to_2fp.selfNormalize();
	// vector perpendicular to the triangle face, oriented in the positive sense of z axis
	Matrix1D<double> planeVector;
	if (order == 2) {
		planeVector = vectorR3(0., 0., 1.);
	} else if (order > 2) {
		planeVector = vectorProduct(_centre_to_2f, _2f_to_2fp);
	}
	planeVector.selfNormalize();

	// vector perpendicular to the triangle face, oriented in the negative sense of z axis
	Matrix1D<double> planeVector_down;
	if (order == 2) {
		planeVector_down = vectorR3(0., 0., -1.);
	} else if (order > 2) {
		planeVector_down = (-1) * vectorProduct(_centre_to_2f, _2f_to_2fp);
	}
	planeVector_down.selfNormalize();

	//vectExpansion: vectors perpendicular to the unit cell faces
	Matrix1D<double> vp_0 = vectorProduct(_centre_to_2f, planeVector);
	vp_0.selfNormalize();
	vectExpansion.push_back(expanded * vp_0 * rmax / 2);
	Matrix1D<double> vp_1 = vectorProduct(planeVector, _centre_to_2fp);
	vp_1.selfNormalize();
	vectExpansion.push_back(expanded * vp_1 * rmax / 2);
	Matrix1D<double> vp_2 = planeVector_down;
	vectExpansion.push_back(expanded * vp_2 * rmax / 2);

	//computation of coordinates for expandedUnitCell vertices
	Matrix1D<double> d_2f = (expanded / sin(TWOPI / order + offset)) * (-1)
			* _2fp;
	Matrix1D<double> d_2fp = (expanded / sin(TWOPI / order + offset)) * (-1)
			* _2f;

	if (expanded == 0) {
		expandedUnitCell.push_back(_centre);
		expandedUnitCell.push_back(_2f);
		expandedUnitCell.push_back(_2fp);
	} else if (expanded >= 0) {
		//_centre expands to expandedUnitCell[0], which is equivalent to newOriginAfterExpansion
		if (order == 2) {
			expandedUnitCell.push_back(
					_centre + vectExpansion[0] + vectExpansion[2]);
		} else if (order > 2) {
			expandedUnitCell.push_back(
					_centre + d_2f + d_2fp + vectExpansion[2]);
		}
		//_2f expands to expandedUnitCell[1]
		expandedUnitCell.push_back(_2f + vectExpansion[0] + vectExpansion[2]);
		//_2fp expands to expandedUnitCell[2]
		expandedUnitCell.push_back(_2fp + vectExpansion[1] + vectExpansion[2]);
	}
	newOriginAfterExpansion = expandedUnitCell[0];

	//computation of coordinates from expandedUnitCells regarding to the newOriginAfterExpansion
	Matrix1D<double> new_centre_to_new_2f = expandedUnitCell[1]
			- newOriginAfterExpansion;
	Matrix1D<double> new_centre_to_new_2fp = expandedUnitCell[2]
			- newOriginAfterExpansion;
	//vectors normal to faces of the expanded unit cell
	planeVectors.push_back(
			(-1) * vectorProduct(new_centre_to_new_2f, planeVector));
	planeVectors.push_back(
			(-1) * vectorProduct(planeVector, new_centre_to_new_2fp));
	if (order == 2) {
		planeVectors.push_back(planeVector);
	} else if (order > 2) {
		planeVectors.push_back(
				vectorProduct(new_centre_to_new_2f, new_centre_to_new_2fp));
	}
#include "chimeraTesterD.txt" //draws the (expanded) unit cell and directions of vectExpansion vectors using chimera
}
//FUNCTION tetrahedralSymmetry fills out the variable planeVectors
// (normal vectors to the polyhedron that define a unit cell)
void UnitCell::tetrahedralSymmetry(const Matrix1D<double> & _centroid,
		const Matrix1D<double> & _3f, const Matrix1D<double> & _3fp,
		const Matrix1D<double> & _3fpp, const double expanded) {
	//amplified_expanded is the resulting product of multiplying the expanded selected value by the factor 25, which could be changed
	double amplified_expanded = expanded * 25;

	//vectExpansion: vectors perpendicular to the unit cell edges
	Matrix1D<double> vectorPlane_0 = vectorProduct(_3f, _3fp);
	vectorPlane_0.selfNormalize();
	vectExpansion.push_back(amplified_expanded * vectorPlane_0);
	Matrix1D<double> vectorPlane_1 = vectorProduct(_3fp, _3fpp);
	vectorPlane_1.selfNormalize();
	vectExpansion.push_back(amplified_expanded * vectorPlane_1);
	Matrix1D<double> vectorPlane_2 = vectorProduct(_3fpp, _3f);
	vectorPlane_2.selfNormalize();
	vectExpansion.push_back(amplified_expanded * vectorPlane_2);

	//computation of coordinates for expandedUnitCell vertices
	//computation of first_expandedUnitCells considering vectExpansion vectors only
	Matrix1D<double> first_expandedUnitCell_0 = _3f + vectExpansion[0]
			+ vectExpansion[2];
	Matrix1D<double> first_expandedUnitCell_1 = _3fp + vectExpansion[0]
			+ vectExpansion[1];
	Matrix1D<double> first_expandedUnitCell_2 = _3fpp + vectExpansion[1]
			+ vectExpansion[2];

	//_centroid expands to expandedUnitCell[0], which is equivalent to newOriginAfterExpansion
	//intersection of the three planes of coordinates
	if (expanded == 0) {
		expandedUnitCell.push_back(_centroid);
		expandedUnitCell.push_back(_3f);
		expandedUnitCell.push_back(_3fp);
		expandedUnitCell.push_back(_3fpp);
	} else if (expanded >= 0) {
		double x;
		x = first_expandedUnitCell_0(0);
		double y;
		y = first_expandedUnitCell_2(1)
				- (vectExpansion[2](0) / vectExpansion[2](1))
						* (x - first_expandedUnitCell_2(0));
		double z;
		z = first_expandedUnitCell_1(2)
				- (vectExpansion[1](1) / vectExpansion[1](2))
						* (y - first_expandedUnitCell_1(1));
		expandedUnitCell.push_back(_centroid + vectorR3(x, y, z));
		//_3f expands to expandedUnitCell[1]
		expandedUnitCell.push_back(first_expandedUnitCell_0);
		//computation of the module m of elongation of first_expandedUnitCells 1 and 2
		double m1 = _3fp(1) - first_expandedUnitCell_1(1);
		double m2 = (1 / sin60) * (_3fpp(0) - first_expandedUnitCell_2(0));
		double m = std::max(m1, m2);
		//_3fp expands to expandedUnitCell[2]
		Matrix1D<double> v1 = first_expandedUnitCell_1
				- newOriginAfterExpansion;
		v1.selfNormalize();
		expandedUnitCell.push_back(first_expandedUnitCell_1 + (m * v1));
		//_3fpp expands to expandUnitCell[3]
		Matrix1D<double> v2 = first_expandedUnitCell_2
				- newOriginAfterExpansion;
		v2.selfNormalize();
		expandedUnitCell.push_back(first_expandedUnitCell_2 + (m * v2));
	}
	newOriginAfterExpansion = expandedUnitCell[0];

	//computation of coordinates from expandedUnitCells regarding to the newOriginAfterExpansion
	Matrix1D<double> new_centroid_to_new_3f = expandedUnitCell[1]
			- newOriginAfterExpansion;
	Matrix1D<double> new_centroid_to_new_3fp = expandedUnitCell[2]
			- newOriginAfterExpansion;
	Matrix1D<double> new_centroid_to_new_3fpp = expandedUnitCell[3]
			- newOriginAfterExpansion;
	//vectors normal to faces of the expanded unit cell
	planeVectors.push_back(
			(-1)
					* vectorProduct(new_centroid_to_new_3f,
							new_centroid_to_new_3fp));
	planeVectors.push_back(
			(-1)
					* vectorProduct(new_centroid_to_new_3fp,
							new_centroid_to_new_3fpp));
	planeVectors.push_back(
			(-1)
					* vectorProduct(new_centroid_to_new_3fpp,
							new_centroid_to_new_3f));
#include "chimeraTesterT.txt" //draws the (expanded) unit cell and directions of vectExpansion vectors using chimera
}
//FUNCTION octahedralSymmetry fills out the variable planeVectors
// (normal vectors to the polyhedron that define a unit cell)
void UnitCell::octahedralSymmetry(const Matrix1D<double> & _centre,
		const Matrix1D<double> & _4f, const Matrix1D<double> & _4fp,
		const Matrix1D<double> & _4fpp, double expanded) {
	//amplified_expanded is the resulting product of multiplying the expanded selected value by the factor 25, which could be changed
	double amplified_expanded = expanded * 25;
	Matrix1D<double> _2f = (_4f + _4fp) / 2.;
	Matrix1D<double> _2fp = (_4fpp + _4f) / 2.;
	Matrix1D<double> _3f = (_4f + _4fp + _4fpp) / 3.;
	//vectors that join a symmetry axis with the next only.
	Matrix1D<double> _4f_to_2f = _2f - _4f;
	Matrix1D<double> _2f_to_3f = _3f - _2f;
	Matrix1D<double> _3f_to_2fp = _2fp - _3f;
	Matrix1D<double> _2fp_to_4f = _4f - _2fp;
	// vector perpendicular to the triangle face.
	Matrix1D<double> planeVector = vectorProduct(_4f_to_2f, _2fp_to_4f);
	planeVector.selfNormalize();
	//vectExpansion: vectors perpendicular to the unit cell edges
	Matrix1D<double> v0 = vectorProduct(planeVector, _4f_to_2f);
	v0.selfNormalize();
	Matrix1D<double> v1 = vectorProduct(planeVector, _2f_to_3f);
	v1.selfNormalize();
	Matrix1D<double> v2 = vectorProduct(planeVector, _3f_to_2fp);
	v2.selfNormalize();
	Matrix1D<double> v3 = vectorProduct(planeVector, _2fp_to_4f);
	v3.selfNormalize();

	vectExpansion.push_back(amplified_expanded * planeVector);
	vectExpansion.push_back(amplified_expanded * v0);
	vectExpansion.push_back(amplified_expanded * v1);
	vectExpansion.push_back(amplified_expanded * v2);
	vectExpansion.push_back(amplified_expanded * v3);

	//computation of coordinates for expandedUnitCell vertices
	_4f_to_2f.selfNormalize();
	_2fp_to_4f.selfNormalize();
	Matrix1D<double> first_expandedUnitCell_0 = _4f
			+ (amplified_expanded / sin60) * (-1) * (_4f_to_2f)
			+ (amplified_expanded / sin60) * (_2fp_to_4f);
	Matrix1D<double> first_expandedUnitCell_1 = _2f + vectExpansion[1]
			+ vectExpansion[2];
	Matrix1D<double> v4 = (_2f_to_3f - _3f_to_2fp);
	v4.selfNormalize();
	Matrix1D<double> first_expandedUnitCell_2 = _3f
			+ v4 * amplified_expanded / sin60;
	//_centre expands to expandedUnitCell[0], which is equivalent to newOriginAfterExpansion
	//intersection of the three planes of coordinates
	if (expanded == 0) {
		expandedUnitCell.push_back(_centre);
		expandedUnitCell.push_back(_4f);
		expandedUnitCell.push_back(_2f);
		expandedUnitCell.push_back(_3f);
		expandedUnitCell.push_back(_2fp);
	} else if (expanded >= 0) {
		double x;
		x = first_expandedUnitCell_1(0)
				+ (first_expandedUnitCell_0(1) - first_expandedUnitCell_1(1));
		double y;
		y = first_expandedUnitCell_0(1);
		double z;
		z = first_expandedUnitCell_0(2);
		expandedUnitCell.push_back(_centre + vectorR3(x, y, z));
		//_4f expands to expandedUnitCell[1]
		expandedUnitCell.push_back(first_expandedUnitCell_0);
		//_2f expands to expandedUnitCell[2]
		expandedUnitCell.push_back(first_expandedUnitCell_1);
		//_3f expands to expandeUnitCell[3]
		expandedUnitCell.push_back(first_expandedUnitCell_2);
		//_2fp expands to expandedUnitCell[4]
		expandedUnitCell.push_back(_2fp + vectExpansion[3] + vectExpansion[4]);
	}
	newOriginAfterExpansion = expandedUnitCell[0];
	//computation of coordinates from expandedUnitCells regarding to the newOriginAfterExpansion
	Matrix1D<double> new_expandedUnitCell_1 = expandedUnitCell[1]
			- newOriginAfterExpansion;
	Matrix1D<double> new_expandedUnitCell_2 = expandedUnitCell[2]
			- newOriginAfterExpansion;
	Matrix1D<double> new_expandedUnitCell_3 = expandedUnitCell[3]
			- newOriginAfterExpansion;
	Matrix1D<double> new_expandedUnitCell_4 = expandedUnitCell[4]
			- newOriginAfterExpansion;
	//vectors normal to faces of the expanded unit cell
	planeVectors.push_back(vectorProduct(new_expandedUnitCell_1,
							new_expandedUnitCell_2));
	planeVectors.push_back(vectorProduct(new_expandedUnitCell_2,
							new_expandedUnitCell_3));
	planeVectors.push_back(vectorProduct(new_expandedUnitCell_3,
							new_expandedUnitCell_4));
	planeVectors.push_back(vectorProduct(new_expandedUnitCell_4,
							new_expandedUnitCell_1));
#include "chimeraTesterO.txt" //draws the (expanded) unit cell 
	//and directions of vectExpansion vectors using chimera
}//FUNCTION icoSymmetry fills out the variable planeVectors
// (normal vectors to the polyhedron that define a unit cell)
void UnitCell::icoSymmetry(const Matrix1D<double> & _centre,
		const Matrix1D<double> & _5f, const Matrix1D<double> & _5fp,
		const Matrix1D<double> & _5fpp, double expanded) {
	Matrix1D<double> _2f = (_5f + _5fp) / 2.;
	Matrix1D<double> _2fp = (_5fpp + _5f) / 2.;
	Matrix1D<double> _3f = (_5f + _5fp + _5fpp) / 3.;

	//vectors that join a symmetry axis with the next only.
	Matrix1D<double> _5f_to_2f = _2f - _5f;
	Matrix1D<double> _2f_to_3f = _3f - _2f;
	Matrix1D<double> _3f_to_2fp = _2fp - _3f;
	Matrix1D<double> _2fp_to_5f = _5f - _2fp;

	// vector perpendicular to the triangle face.
	Matrix1D<double> planeVector = vectorProduct(_5f_to_2f, _2fp_to_5f);

	//vectExpansion: vectors perpendicular to the unit cell edges
	Matrix1D<double> v0 = vectorProduct(planeVector, _5f_to_2f);
	v0.selfNormalize();
	vectExpansion.push_back(expanded * v0);
	Matrix1D<double> v1 = vectorProduct(planeVector, _2f_to_3f);
	v1.selfNormalize();
	vectExpansion.push_back(expanded * v1);
	Matrix1D<double> v2 = vectorProduct(planeVector, _3f_to_2fp);
	v2.selfNormalize();
	vectExpansion.push_back(expanded * v2);
	Matrix1D<double> v3 = vectorProduct(planeVector, _2fp_to_5f);
	v3.selfNormalize();
	vectExpansion.push_back(expanded * v3);

	//computation of coordinates for expandedUnitCell vertices
	//_centre remains constant and it is equivalent to newOriginAfterExpansion
	newOriginAfterExpansion = _centre;
	if (expanded == 0) {
		expandedUnitCell.push_back(_5f);
		expandedUnitCell.push_back(_2f);
		expandedUnitCell.push_back(_3f);
		expandedUnitCell.push_back(_2fp);
	} else if (expanded >= 0) {
		//60 is the angle between _5f_to_2f and _2fp_to_5f
		//mod: proyection (module) of expanded in the negative direction of the _5f_to_2f vector, and
		//proyection (module) of expanded in the direction of the _2fp_to_5f vector
		//double mod = tg60 * (sin60 * expanded);
		double mod = expanded / sin60;
		_2fp_to_5f.selfNormalize();
		_5f_to_2f.selfNormalize();
		//_5f expands to expandedUnitCell[0]
		expandedUnitCell.push_back(
				_5f + mod * ((_2fp_to_5f) + ((-1) * _5f_to_2f)));
		//_2f expands to expandedUnitCell[1]
		expandedUnitCell.push_back(_2f + vectExpansion[0] + vectExpansion[1]);
		//mod: proyection (module) of expanded in the negative direction of the 3f_to_2fp vector
		mod = (1 / tg60) * expanded;
		_3f_to_2fp.selfNormalize();
		//_3f expands to expandedUnitCell[2]
		expandedUnitCell.push_back(
				_3f + vectExpansion[2] + ((-1) * _3f_to_2fp * mod));
		//_2fp expands to expandedUnitCell[3]
		expandedUnitCell.push_back(_2fp + vectExpansion[2] + vectExpansion[3]);
	}
	//vectors normal to faces of the expanded polyhedra
	for (int i = 0; i < 4; i++) {
		planeVectors.push_back(
				vectorProduct(expandedUnitCell[i],
						expandedUnitCell[(i + 1) % 4]));
	}
#include "chimeraTesterI2.txt" //draws the (expanded) unit cell and directions of vectExpansion vectors using chimera
}void UnitCell::maskUnitCell(ImageGeneric & in3Dmap, ImageGeneric & out3DDmap) {
	//1) get dimensions
	size_t xDim, yDim, zDim;
	in3Dmap.getDimensions(xDim, yDim, zDim);
	double x_offset = x_origin - xDim/2.; //difference between default and user introduced x_origin coordinate
	double y_offset = y_origin - yDim/2.; //difference between default and user introduced y_origin coordinate
	double z_offset = z_origin - zDim/2; //difference between default and user introduced z_origin coordinate
	
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
			t = expandedUnitCell[i] * scale;
			tt = expandedUnitCell[i] * scale + planeVectors[i] * scale;
			testFile << ".v " << t(0) << " " << t(1) << " " << t(2) << " "
			<< tt(0) << " " << tt(1) << " " << tt(2) << " 1\n";
		}
		testFile.close();
	}
#endif

	//2) get expanded unitcell enclosing box
	{
		double minX = rmin;
		double minY = rmin;
		double minZ = rmin;
		double maxX = rmax;
		double maxY = rmax;
		double maxZ = rmax;

		if (sym_order != 1) {
			 minX = rmax;
			 minY = rmax;
			 minZ = _minZ;
			 maxX = rmin;
			 maxY;
			if (symmetry == pg_CN || symmetry == pg_DN) {
				maxY = rmax;
			} else if (symmetry != pg_CN && symmetry != pg_DN) {
				maxY = rmin;
			}
			 maxZ = _maxZ;

			Matrix1D<double> minVector, maxVector;
			for (std::vector<Matrix1D<double> >::iterator it =
					expandedUnitCell.begin(); it != expandedUnitCell.end();
					++it) {
				if (symmetry == pg_CN || symmetry == pg_DN || symmetry == pg_T
						|| symmetry == pg_O) {
					minVector = (*it);
					maxVector = (*it);
				} else if (symmetry >= pg_I1 || symmetry <= pg_I4) {
					(*it).selfNormalize();
					minVector = (*it) * rmin;
					maxVector = (*it) * std::min(rmax, (double) xDim / 2.);
				}
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

		}

		else if(sym_order == 1 && symmetry == pg_CN) {

			minX = - rmax;
			minY = - rmax;
			minZ = - rmax;
			maxX = rmax;
			maxY = rmax;
			maxZ = rmax;
		}
		std::cout << "minX " << minX << std::endl;
		std::cout << "minY " << minY << std::endl;
		std::cout << "minZ " << minZ << std::endl;
		std::cout << "maxX " << maxX << std::endl;
		std::cout << "maxY " << maxY << std::endl;
		std::cout << "maxZ " << maxZ << std::endl;
#ifdef DEBUG
		//draw real unitcell
		{
#include <iostream>
#include <fstream>
			//Debug chimera file
			std::ofstream testFile;

			//draw box
			testFile.open("box.bild");
			Matrix1D<double> t, tt;
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

		//check that the output 3Dmap has a valid size

		const size_t N = 1168;
		int a[N] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
					17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 30, 32, 33, 34,
					35, 36, 38, 39, 40, 42, 44, 45, 48, 49, 50, 51, 52, 54, 55,
					56, 57, 60, 63, 64, 65, 66, 68, 70, 72, 75, 76, 77, 78, 80,
					81, 84, 85, 88, 90, 91, 95, 96, 98, 99, 100, 102, 104, 105,
					108, 110, 112, 114, 117, 119, 120, 121, 125, 126, 128, 130,
					132, 133, 135, 136, 140, 143, 144, 147, 150, 152, 153, 154,
					156, 160, 162, 165, 168, 169, 170, 171, 175, 176, 180, 182,
					187, 189, 190, 192, 195, 196, 198, 200, 204, 208, 209, 210,
					216, 220, 221, 224, 225, 228, 231, 234, 238, 240, 242, 243,
					245, 247, 250, 252, 255, 256, 260, 264, 266, 270, 272, 273,
					275, 280, 285, 286, 288, 289, 294, 297, 300, 304, 306, 308,
					312, 315, 320, 323, 324, 325, 330, 336, 338, 340, 342, 343,
					350, 351, 352, 357, 360, 361, 363, 364, 374, 375, 378, 380,
					384, 385, 390, 392, 396, 399, 400, 405, 408, 416, 418, 420,
					425, 429, 432, 440, 441, 442, 448, 450, 455, 456, 459, 462,
					468, 475, 476, 480, 484, 486, 490, 494, 495, 500, 504, 507,
					510, 512, 513, 520, 525, 528, 532, 539, 540, 544, 546, 550,
					560, 561, 567, 570, 572, 576, 578, 585, 588, 594, 595, 600,
					605, 608, 612, 616, 624, 625, 627, 630, 637, 640, 646, 648,
					650, 660, 663, 665, 672, 675, 676, 680, 684, 686, 693, 700,
					702, 704, 714, 715, 720, 722, 726, 728, 729, 735, 741, 748,
					750, 756, 760, 765, 768, 770, 780, 784, 792, 798, 800, 810,
					816, 819, 825, 832, 833, 836, 840, 845, 847, 850, 855, 858,
					864, 867, 875, 880, 882, 884, 891, 896, 900, 910, 912, 918,
					924, 931, 935, 936, 945, 950, 952, 960, 968, 969, 972, 975,
					980, 988, 990, 1000, 1001, 1008, 1014, 1020, 1024, 1026,
					1029, 1040, 1045, 1050, 1053, 1056, 1064, 1071, 1078, 1080,
					1083, 1088, 1089, 1092, 1100, 1105, 1120, 1122, 1125, 1134,
					1140, 1144, 1152, 1155, 1156, 1170, 1176, 1183, 1188, 1190,
					1197, 1200, 1210, 1215, 1216, 1224, 1225, 1232, 1235, 1248,
					1250, 1254, 1260, 1274, 1275, 1280, 1287, 1292, 1296, 1300,
					1309, 1320, 1323, 1326, 1330, 1331, 1344, 1350, 1352, 1360,
					1365, 1368, 1372, 1375, 1377, 1386, 1400, 1404, 1408, 1425,
					1428, 1430, 1440, 1444, 1445, 1452, 1456, 1458, 1463, 1470,
					1482, 1485, 1496, 1500, 1512, 1520, 1521, 1530, 1536, 1539,
					1540, 1547, 1560, 1568, 1573, 1575, 1584, 1596, 1600, 1615,
					1617, 1620, 1625, 1632, 1638, 1650, 1664, 1666, 1672, 1680,
					1683, 1690, 1694, 1700, 1701, 1710, 1715, 1716, 1728, 1729,
					1734, 1750, 1755, 1760, 1764, 1768, 1782, 1785, 1792, 1800,
					1805, 1815, 1820, 1824, 1836, 1848, 1859, 1862, 1870, 1872,
					1875, 1881, 1890, 1900, 1904, 1911, 1920, 1925, 1936, 1938,
					1944, 1950, 1960, 1976, 1980, 1989, 1995, 2000, 2002, 2016,
					2023, 2025, 2028, 2040, 2048, 2052, 2057, 2058, 2079, 2080,
					2090, 2100, 2106, 2112, 2125, 2128, 2142, 2145, 2156, 2160,
					2166, 2176, 2178, 2184, 2187, 2197, 2200, 2205, 2210, 2223,
					2240, 2244, 2250, 2261, 2268, 2275, 2280, 2288, 2295, 2299,
					2304, 2310, 2312, 2340, 2352, 2366, 2375, 2376, 2380, 2394,
					2400, 2401, 2420, 2430, 2431, 2432, 2448, 2450, 2457, 2464,
					2470, 2475, 2496, 2499, 2500, 2508, 2520, 2527, 2535, 2541,
					2548, 2550, 2560, 2565, 2574, 2584, 2592, 2600, 2601, 2618,
					2625, 2640, 2646, 2652, 2660, 2662, 2673, 2688, 2695, 2700,
					2704, 2717, 2720, 2730, 2736, 2744, 2750, 2754, 2772, 2793,
					2800, 2805, 2808, 2816, 2835, 2850, 2856, 2860, 2873, 2880,
					2888, 2890, 2904, 2907, 2912, 2916, 2925, 2926, 2940, 2964,
					2970, 2975, 2992, 3000, 3003, 3024, 3025, 3040, 3042, 3060,
					3072, 3078, 3080, 3087, 3094, 3120, 3125, 3135, 3136, 3146,
					3150, 3159, 3168, 3179, 3185, 3192, 3200, 3211, 3213, 3230,
					3234, 3240, 3249, 3250, 3264, 3267, 3276, 3300, 3315, 3325,
					3328, 3332, 3344, 3360, 3366, 3375, 3380, 3388, 3400, 3402,
					3420, 3430, 3432, 3456, 3458, 3465, 3468, 3500, 3510, 3520,
					3528, 3536, 3549, 3553, 3564, 3570, 3575, 3584, 3591, 3600,
					3610, 3630, 3640, 3645, 3648, 3672, 3675, 3696, 3705, 3718,
					3724, 3740, 3744, 3750, 3757, 3762, 3773, 3780, 3800, 3808,
					3822, 3825, 3840, 3850, 3861, 3872, 3876, 3888, 3900, 3920,
					3927, 3952, 3960, 3969, 3971, 3978, 3990, 3993, 4000, 4004,
					4032, 4046, 4050, 4056, 4080, 4095, 4096, 4104, 4114, 4116,
					4125, 4131, 4158, 4160, 4165, 4180, 4199, 4200, 4212, 4224,
					4225, 4235, 4250, 4256, 4275, 4284, 4290, 4312, 4320, 4332,
					4335, 4352, 4356, 4368, 4374, 4375, 4389, 4394, 4400, 4410,
					4420, 4446, 4455, 4459, 4480, 4488, 4500, 4522, 4536, 4550,
					4560, 4563, 4576, 4590, 4598, 4608, 4617, 4620, 4624, 4641,
					4655, 4675, 4680, 4693, 4704, 4719, 4725, 4732, 4750, 4752,
					4760, 4788, 4800, 4802, 4840, 4845, 4851, 4860, 4862, 4864,
					4875, 4896, 4900, 4913, 4914, 4928, 4940, 4950, 4992, 4998,
					5000, 5005, 5016, 5040, 5049, 5054, 5070, 5082, 5096, 5100,
					5103, 5120, 5130, 5145, 5148, 5168, 5184, 5187, 5200, 5202,
					5225, 5236, 5250, 5265, 5280, 5292, 5304, 5320, 5324, 5346,
					5355, 5376, 5390, 5400, 5408, 5415, 5434, 5440, 5445, 5460,
					5472, 5488, 5491, 5500, 5508, 5525, 5544, 5577, 5586, 5600,
					5610, 5616, 5625, 5632, 5643, 5670, 5700, 5712, 5720, 5733,
					5746, 5760, 5775, 5776, 5780, 5808, 5814, 5824, 5831, 5832,
					5850, 5852, 5880, 5915, 5928, 5929, 5940, 5950, 5967, 5984,
					5985, 6000, 6006, 6048, 6050, 6069, 6075, 6080, 6084, 6120,
					6125, 6137, 6144, 6156, 6160, 6171, 6174, 6175, 6188, 6237,
					6240, 6250, 6270, 6272, 6292, 6300, 6318, 6336, 6358, 6370,
					6375, 6384, 6400, 6422, 6426, 6435, 6460, 6468, 6480, 6498,
					6500, 6517, 6528, 6534, 6545, 6552, 6561, 6591, 6600, 6615,
					6630, 6650, 6655, 6656, 6664, 6669, 6688, 6720, 6732, 6750,
					6760, 6776, 6783, 6800, 6804, 6825, 6840, 6859, 6860, 6864,
					6875, 6885, 6897, 6912, 6916, 6930, 6936, 7000, 7007, 7020,
					7040, 7056, 7072, 7098, 7106, 7125, 7128, 7140, 7150, 7168,
					7182, 7200, 7203, 7220, 7225, 7260, 7280, 7290, 7293, 7296,
					7315, 7344, 7350, 7371, 7392, 7410, 7425, 7436, 7448, 7480,
					7488, 7497, 7500, 7514, 7524, 7546, 7560, 7581, 7600, 7605,
					7616, 7623, 7644, 7650, 7680, 7695, 7700, 7722, 7735, 7744,
					7752, 7776, 7800, 7803, 7840, 7854, 7865, 7875, 7904, 7920,
					7938, 7942, 7956, 7980, 7986, 8000, 8008, 8019, 8064, 8075,
					8085, 8092, 8100, 8112, 8125, 8151, 8160, 8190, 8192, 8208,
					8228, 8232, 8250, 8262, 8281, 8316, 8320, 8330, 8360, 8379,
					8398, 8400, 8415, 8424, 8448, 8450, 8470, 8500, 8505, 8512,
					8550, 8568, 8575, 8580, 8619, 8624, 8640, 8645, 8664, 8670,
					8704, 8712, 8721, 8736, 8748, 8750, 8775, 8778, 8788, 8800,
					8820, 8840, 8892, 8910, 8918, 8925, 8960, 8976, 9000, 9009,
					9025, 9044, 9072, 9075, 9100, 9120, 9126, 9152, 9163, 9180,
					9196, 9216, 9234, 9240, 9248, 9261, 9282, 9295, 9310, 9317,
					9350, 9360, 9375, 9386, 9405, 9408, 9438, 9450, 9464, 9477,
					9500, 9504, 9520, 9537, 9555, 9576, 9600, 9604, 9625, 9633,
					9639, 9680, 9690, 9702, 9720, 9724, 9728, 9747, 9750, 9792,
					9800, 9801, 9826, 9828, 9856, 9880, 9900, 9945, 9975, 9984,
					9996 };

		std::vector<int> validSizes(a, a + sizeof(a) / sizeof(int));
		std::vector<int>::iterator low; // iterator to search in validSizes
		low = std::lower_bound(validSizes.begin(), validSizes.end(), iMaxX);
		if (validSizes[low - validSizes.begin()] != iMaxX)
			iMaxX = validSizes[low - validSizes.begin()];
		low = std::lower_bound(validSizes.begin(), validSizes.end(), iMaxY);
		if (validSizes[low - validSizes.begin()] != iMaxY)
			iMaxY = validSizes[low - validSizes.begin()];
		low = std::lower_bound(validSizes.begin(), validSizes.end(), iMaxZ);
		if (validSizes[low - validSizes.begin()] != iMaxZ)
			iMaxZ = validSizes[low - validSizes.begin()];

		double dotproduct;
		bool doIt = false;

		//for all points in the enclosing box
		//check if they are inside the expanded unit cell
		int ii, jj, kk; //expandedUnitCell
		for (int k = iMinZ; k <= iMaxZ; ++k)
			for (int i = iMinY; i <= iMaxY; ++i)
				for (int j = iMinX; j <= iMaxX; ++j) {
					r2 = (i * i + j * j + k * k);
					if (r2 >= rmin2 && r2 < rmax2) {
						doIt = true;
						if (sym_order != 1) {
							for (std::vector<Matrix1D<double> >::iterator it1 =
									planeVectors.begin();
									it1 != planeVectors.end(); ++it1) {

								dotproduct =
										dotProduct(*it1,
												vectorR3(
														((double) j)
																- newOriginAfterExpansion(
																		0),
														((double) i)
																- newOriginAfterExpansion(
																		1),
														((double) k)
																- newOriginAfterExpansion(
																		2)));

								if (dotproduct < 0) {
									doIt = false;
									break;
								}
							}
						}

						if (doIt) {
							if (x_offset == 0.0 && y_offset == 0.0 && z_offset == 0.0){
								A3D_ELEM(*imageMap2,k,i,j) = A3D_ELEM(*map, k, i, j);
							}else {
								int k1 = k + z_offset;
								int i1 = i + y_offset;
								int j1 = j + x_offset;
								A3D_ELEM(*imageMap2,k1,i1,j1) = A3D_ELEM(*map, k1, i1, j1);
							}
						}
					}
				}
		if (x_offset != 0.0 || y_offset != 0.0 || z_offset != 0.0){
			iMinZ = iMinZ + z_offset;
			iMinY = iMinY + y_offset;
			iMinX = iMinX + x_offset;
			iMaxZ = iMaxZ + z_offset;
			iMaxY = iMaxY + y_offset;
			iMaxX = iMaxX + x_offset;
		}
		imageMap2->selfWindow(iMinZ, iMinY, iMinX, iMaxZ, iMaxY, iMaxX, 0.);
		MDRow MD;
		out3DDmap.setDataMode(_DATA_ALL);
		MD.setValue(MDL_SHIFT_X, -(double) iMinX);
		MD.setValue(MDL_SHIFT_Y, -(double) iMinY);
		MD.setValue(MDL_SHIFT_Z, -(double) iMinZ);

		out3DDmap.image->MDMainHeader.setValue(MDL_SAMPLINGRATE_X, sampling);

		out3DDmap.image->MDMainHeader.setValue(MDL_SAMPLINGRATE_Y, sampling);

		out3DDmap.image->MDMainHeader.setValue(MDL_SAMPLINGRATE_Z, sampling);

		out3DDmap.image->setGeo(MD, 0);
	}
}
