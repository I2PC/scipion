/*
 * matrix2d.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: scheres
 */
/* Is diagonal ------------------------------------------------------------- */
template <>
bool Matrix2D< std::complex<double> >::isDiagonal() const
{
    if (XSIZE(*this) != YSIZE(*this))
        return false;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        if (i != j && abs(DIRECT_MAT_ELEM(*this, i, j)) > XMIPP_EQUAL_ACCURACY)
            return false;
    return true;
}

/* Is Scalar --------------------------------------------------------------- */
template <>
bool Matrix2D< std::complex<double> >::isScalar() const
{
    if (!isDiagonal())
        return false;
    for (int i = 1; i < YSIZE(*this); i++)
        if (abs(DIRECT_MAT_ELEM(*this, i, i) - DIRECT_MAT_ELEM(*this, 0, 0)) >
            XMIPP_EQUAL_ACCURACY)
            return false;
    return true;
}


/* Interface to numerical recipes: svbksb ---------------------------------- */
void svbksb(Matrix2D<double> &u, Matrix1D<double> &w, Matrix2D<double> &v,
            Matrix1D<double> &b, Matrix1D<double> &x)
{
    // Call to the numerical recipes routine. Results will be stored in X
    svbksb(u.adaptForNumericalRecipes22D(),
           w.adaptForNumericalRecipes1D(),
           v.adaptForNumericalRecipes22D(),
           u.ydim, u.xdim,
           b.adaptForNumericalRecipes1D(),
           x.adaptForNumericalRecipes1D());
}

