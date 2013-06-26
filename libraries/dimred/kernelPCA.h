/***************************************************************************
 *
 * Authors:    Albrecht Wolf
 *
 ***************************************************************************/
#ifndef _KERNEL_PCA
#define _KERNEL_PCA

#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include "dimred_tools.h"

/**@defgroup KernelPCA Kernel PCA
   @ingroup DimRedLibrary */
//@{
/** Class for making a Kernel PCA dimensionality reduction */
class KernelPCA: public DimRedAlgorithm
{
public:
	double sigma;
public:
	/// Set specific parameters
	void setSpecificParameters(double sigma=1.0);

	/// Reduce dimensionality
	void reduceDimensionality();
};
//@}
#endif
