/*
 * matrix2d.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: scheres
 */
/* Is diagonal ------------------------------------------------------------- */
#include "matrix2d.h"

template <>
bool Matrix2D< std::complex<double> >::isDiagonal() const
{
    if (mdimx != mdimy)
        return false;
    for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; i++)
        {
        	if (i != j && abs((*this)(i, j)) > XMIPP_EQUAL_ACCURACY)
        	            return false;
        }
    return true;
}

/* Is Scalar --------------------------------------------------------------- */
template <>
bool Matrix2D< std::complex<double> >::isScalar() const
{
    if (!isDiagonal())
        return false;
    for (int i = 1; i < mdimy; i++)
        if (abs((*this)(i, i) - (*this)(0, 0)) >
            XMIPP_EQUAL_ACCURACY)
            return false;
    return true;
}


/* Interface to numerical recipes: svbksb ---------------------------------- */
void svbksb(Matrix2D<double> &u, Matrix1D<double> &w, Matrix2D<double> &v,
            Matrix1D<double> &b, Matrix1D<double> &x)
{
    // Call to the numerical recipes routine. Results will be stored in X
    svbksb(u.adaptForNumericalRecipes2(),
           w.adaptForNumericalRecipes(),
           v.adaptForNumericalRecipes2(),
           u.mdimy, u.mdimx,
           b.adaptForNumericalRecipes(),
           x.adaptForNumericalRecipes());
}

