#ifndef _NPE
#define _NPE

#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include "dimred_tools.h"

class NPE: public DimRedAlgorithm
{
public:
	int k;
	double tol;
public:
	/// Set specific parameters
	void setSpecificParameters(int k=12);

	/// Reduce dimensionality
	void reduceDimensionality();
};
//@}
#endif
