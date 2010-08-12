/*
 * matrix2d.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: scheres
 */
/* Is diagonal ------------------------------------------------------------- */
#include "matrix2d.h"

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

