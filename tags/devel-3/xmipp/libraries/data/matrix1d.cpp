/*
 * matrix1d.cpp
 *
 *  Created on: Apr 13, 2010
 *      Author: scheres
 */

/* Vector R2 and R3 -------------------------------------------------------- */
Matrix1D<double> vectorR2(double x, double y)
{
    Matrix1D<double> result(2);
    result( 0) = x;
    result( 1) = y;
    return result;
}

Matrix1D<double> vectorR3(double x, double y, double z)
{
    Matrix1D<double> result(3);
    result( 0) = x;
    result( 1) = y;
    result( 2) = z;
    return result;
}

Matrix1D<int> vectorR3(int x, int y, int z)
{
    Matrix1D<int> result(3);
    result( 0) = x;
    result( 1) = y;
    result( 2) = z;
    return result;
}
