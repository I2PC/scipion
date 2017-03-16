#ifndef _UNITCELL_HH
#define _UNITCELL_HH
/*
#include "matrix1d.h"
#include "matrix2d.h"
#include "xmipp_funcs.h"
#include "args.h"
#include "grids.h"
*/
#include "symmetries.h"
#include <vector>
#include <iterator>
#include "xmipp_image_generic.h"
#include "data/xmipp_image.h"

class UnitCell{

/** Typename to contain a list of plane_vector */

private:
	//point group symmetry
	int symmetry;
	//pointgroup symmetry order
	int sym_order;
	//symmetry object
	SymList SL;
	//min max radii for masking crown
	double rmin, rmax;
	//offset for unitcell
	double offset;//, expand;
	//vectors defining expansion direction
	std::vector<Matrix1D<double> > vectExpansion;
	// normal vector to the polyhedron that define a unit cell.
	std::vector<Matrix1D<double> > planeVectors;
	//vectors with corner of expanded unit cell
	std::vector<Matrix1D<double> > expandedUnitCell;
	//above vector multiplied by rmin
	std::vector<Matrix1D<double> > expandedUnitCellMin;
	//above vector multiplied by rmax
	std::vector<Matrix1D<double> > expandedUnitCellMax;
public:
	/** set to zero everything that is not the unit cell */
	void maskUnitCell(ImageGeneric & in3DDmap,
			          ImageGeneric & out3DDmap);
	/** close debug file */
	//~UnitCell();
	/** unit auxiliary vectors */
	UnitCell(String sym, double rmin, double rmax, double expanded,
			double offset);
};//end unitcell class

#endif
