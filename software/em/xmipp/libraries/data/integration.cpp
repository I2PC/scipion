/***************************************************************************
 *
 * Authors:     Javier Rodriguez Falces (jrodriguez@cnb.csic.es)
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

#include "integration.h"

/* Integrate --------------------------------------------------------------- */
double integrateNewtonCotes(double(*f)(double),
                            double a, double b, int N)
{
    if (N < 2 || N > 9)
        REPORT_ERROR(ERR_VALUE_INCORRECT, "integrateNewtonCotes: N must be greater than 1");
    double h = (b - a) / (N - 1);
    Matrix1D<double> fx(N);
    for (int i = 0; i < N; i++)
        fx(i) = (*f)(a + i * h);
    switch (N)
    {
    case 2:
        return h / 2*(fx(0) + fx(1));
    case 3:
        return h / 3*(fx(0) + 4*fx(1) + fx(2));
    case 4:
        return h*3.0 / 8.0*((fx(0) + fx(3)) + 3*(fx(1) + fx(2)));
    case 5:
        return h*2.0 / 45.0*(7*(fx(0) + fx(4)) + 32*(fx(1) + fx(3)) + 12*fx(2));
    case 6:
        return h*5.0 / 288.0*(19*(fx(0) + fx(5)) + 75*(fx(1) + fx(4)) + 50*(fx(2) + fx(3)));
    case 7:
        return h / 140.0*(41*(fx(0) + fx(6)) + 216*(fx(1) + fx(5)) + 27*(fx(2) + fx(4)) +
                          272*fx(3));
    case 8:
        return h*7.0 / 17280.0*(751*(fx(0) + fx(7)) + 3577*(fx(1) + fx(6)) + 1323*(fx(2) + fx(5)) +
                                2989*(fx(3) + fx(4)));
    case 9:
        return 4.0 / 14175.0*h*(989*(fx(0) + fx(8)) + 5888*(fx(1) + fx(7)) +
                                -928*(fx(2) + fx(6)) + 10496*(fx(3) + fx(5)) - 4540*fx(4));
    }
    REPORT_ERROR(ERR_ARG_INCORRECT,"Number of points is too high");
}

//**********************************************************
// Implementation of the integral using the Trapeze method
//**********************************************************

double Trapeze::operator()()
{   //adapted from qtrap
    double s,olds;
    int j;
    olds = -1.0e30;
    for (j = 1;j <= JMAX;j++)
    {
        s = Trap(j);    //changed; Trap is integrating fcn
        if (fabs(s - olds) <= EPS*fabs(olds))
            return s;
        if (s == 0.0 && olds == 0.0 && j > 6)
            return s;
        olds = s;
    }
    printf("Too many steps in routine qtrap_y\n");
    exit(1);
    return 0.0;
}


double Trapeze::Trap(int n)
{ //adapted from trapzd
    double tnm, sum, del;
    int j, it;
    if (n == 1)
    {
        it = 1;
        x = a;
        s = func();  //changed
        x = b;
        s += func(); //changed
        return (s *= 0.5 * (b - a));
    }
    else
    {
        for (it = 1, j = 1;j < n - 1;j++)
            it <<= 1;
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 1;j <= it;j++, x += del)
            sum += func(); //changed
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
}

//**********************************************************
// Implementation of the integral using the Romberg method
//**********************************************************mask.cpp to adapt
#define JMAXP 30
#define K 5

double Romberg::operator()()
{  //adapted from qromb
    int j;
    double ss,dss, h[JMAXP+2], s[JMAXP+2];
    h[1] = 1.0;
    for (j = 1;j <= JMAXP;j++)
    {
        s[j] = midpnt(j); //changed; midpnt is integrating
        if (j >= K)
        {     //function
            polint(&h[j-K], &s[j-K], K, 0.0, ss, dss);
            if (fabs(dss) <= EPS*fabs(ss))
                return ss;
        }
        s[j+1] = s[j];
        h[j+1] = h[j] / 9.0;
    }
    REPORT_ERROR(ERR_NUMERICAL,"Too many steps in routine Romberg");
    return 0.0;
}

//*
// The midpnt function is used in the Romberg::operator only
//*
double Romberg::midpnt(int n)
{   //adapted from midpoint
    double tnm, sum, del, ddel;
    int it, j;
    if (n == 1)
    {
        x = 0.5 * (a + b);     //changed; set x
        return (s = (b - a) * func());  //changed; evaluate func
    }
    else
    {
        for (it = 1, j = 1;j < n - 1;j++)
            it *= 3;
        tnm = it;
        del = (b - a) / (3.0 * tnm);
        ddel = del + del;
        x = a + 0.5 * del;
        sum = 0.0;
        for (j = 1;j <= it;j++)
        {
            sum += func();   //changed; evaluate func
            x += ddel;    //changed; set x
            sum += func();   //changed; evaluate func
            x += del;    //changed; set x
        }
        s = (s + (b - a) * sum / tnm) / 3.0;
        return s;
    }
}

