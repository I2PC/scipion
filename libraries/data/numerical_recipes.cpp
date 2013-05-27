/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "numerical_recipes.h"


/* NUMERICAL UTILITIES ----------------------------------------------------- */
void nrerror(const char error_text[])
{
    fprintf(stderr, "Numerical Recipes run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* RANDOM NUMBERS ---------------------------------------------------------- */
#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

/* Chapter 7 Section 1: UNIFORM RANDOM NUMBERS */
double ran1(int *idum)
{
    static long ix1, ix2, ix3;
    static double r[98];
    double temp;
    static int iff = 0;
    int j;

    if (*idum < 0 || iff == 0)
    {
        iff = 1;
        ix1 = (IC1 - (*idum)) % M1;
        ix1 = (IA1 * ix1 + IC1) % M1;
        ix2 = ix1 % M2;
        ix1 = (IA1 * ix1 + IC1) % M1;
        ix3 = ix1 % M3;
        for (j = 1;j <= 97;j++)
        {
            ix1 = (IA1 * ix1 + IC1) % M1;
            ix2 = (IA2 * ix2 + IC2) % M2;
            r[j] = (ix1 + ix2 * RM2) * RM1;
        }
        *idum = 1;
    }
    ix1 = (IA1 * ix1 + IC1) % M1;
    ix2 = (IA2 * ix2 + IC2) % M2;
    ix3 = (IA3 * ix3 + IC3) % M3;
    j = 1 + ((97 * ix3) / M3);
    if (j > 97 || j < 1)
        nrerror("RAN1: This cannot happen.");
    temp = r[j];
    r[j] = (ix1 + ix2 * RM2) * RM1;
    return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

/* Chapter 7 Section 3: GAUSSIAN RANDOM NUMBERS */
double gasdev(int *idum)
{
    static int iset = 0;
    static double gset;
    double fac, r, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            r = v1 * v1 + v2 * v2;
        }
        while (r >= 1.0);
        fac = sqrt(-2.0 * log(r) / r);
        gset = v1 * fac;
        iset = 1;
        return v2*fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}

// t-distribution (nor Numerical Recipes, but Mathematics of Computation, vol. 62, 779-781.
// I downloaded sem-code from http://ideas.repec.org/c/ega/comcod/200703.html
// Sjors May 2008
double tdev(double nu, int *idum)
{
    static int iset = 0;
    static double gset;
    double fac, r, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            r = v1 * v1 + v2 * v2;
        }
        while (r >= 1.0);
        fac = sqrt(nu*(pow(r,-2.0/nu) -1.0)/r);
        gset = v1 * fac;
        iset = 1;
        return v2*fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}


// Kolmogorov-Smirnov test
void ksone(double data[], int n, double(*func)(double), double * d, double * prob)
{
    qcksrt(n, data);
    double fn, ff, en, dt, fo=0.;
    en = (double)n;
    *d = 0.;
    for (int j=1; j<=n; j++)
    {
        fn = j / en;
        ff = (*func)(data[j]);
        dt = XMIPP_MAX(fabs(fo - ff), fabs(fn - ff));
        if (dt> *d)
            *d = dt;
        fo = fn;
    }
    *prob = probks(sqrt(en)*(*d));
}

// Calculate KS-confidence level
double probks(double alam)
{
    int j;
    double a2, fac=2.0, sum=0.0, term, termbf=0.0;
    double EPS1=0.001, EPS2=1.0e-8;

    a2 = -2.0 * alam * alam;
    for (j = 1; j<= 100; j++)
    {
        term = fac * exp(a2*j*j);
        sum += term;
        if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum)
            return sum;
        fac = -fac;
        termbf = fabs(term);
    }
    return 1.0;

}

/* SORTING ----------------------------------------------------------------- */
/* Chapter 8, Section 3: Indexing */
void indexx(int n, double arrin[], int indx[])
{
    int l, j, ir, indxt, i;
    double q;

    for (j = 1;j <= n;j++)
        indx[j] = j;
    l = (n >> 1) + 1;
    ir = n;
    for (;;)
    {
        if (l > 1)
            q = arrin[(indxt=indx[--l])];
        else
        {
            q = arrin[(indxt=indx[ir])];
            indx[ir] = indx[1];
            if (--ir == 1)
            {
                indx[1] = indxt;
                return;
            }
        }
        i = l;
        j = l << 1;
        while (j <= ir)
        {
            if (j < ir && arrin[indx[j]] < arrin[indx[j+1]])
                j++;
            if (q < arrin[indx[j]])
            {
                indx[i] = indx[j];
                j += (i = j);
            }
            else
                j = ir + 1;
        }
        indx[i] = indxt;
    }
}

/* Chapter 8, Section 4: Quicksort */
#define M 7
#define NSTACK 50
#define FM 7875
#define FA 211
#define FC 1663

void qcksrt(int n, double arr[])
{
    int l = 1, jstack = 0, j, ir, iq, i;
    int istack[NSTACK+1];
    long int fx = 0L;
    double a;

    ir = n;
    for (;;)
    {
        if (ir - l < M)
        {
            for (j = l + 1;j <= ir;j++)
            {
                a = arr[j];
                for (i = j - 1;arr[i] > a && i > 0;i--)
                    arr[i+1] = arr[i];
                arr[i+1] = a;
            }
            if (jstack == 0)
                return;
            ir = istack[jstack--];
            l = istack[jstack--];
        }
        else
        {
            i = l;
            j = ir;
            fx = (fx * FA + FC) % FM;
            iq = l + ((ir - l + 1) * fx) / FM;
            a = arr[iq];
            arr[iq] = arr[l];
            for (;;)
            {
                while (j > 0 && a < arr[j])
                    j--;
                if (j <= i)
                {
                    arr[i] = a;
                    break;
                }
                arr[i++] = arr[j];
                while (a > arr[i] && i <= n)
                    i++;
                if (j <= i)
                {
                    arr[(i=j)] = a;
                    break;
                }
                arr[j--] = arr[i];
            }
            if (ir - i >= i - l)
            {
                istack[++jstack] = i + 1;
                istack[++jstack] = ir;
                ir = i - 1;
            }
            else
            {
                istack[++jstack] = l;
                istack[++jstack] = i - 1;
                l = i + 1;
            }
            if (jstack > NSTACK)
                nrerror("NSTACK too small in QCKSRT");
        }
    }
}

#undef M
#undef NSTACK
#undef FM
#undef FA
#undef FC

/* BESSEL FUNCTIONS -------------------------------------------------------- */
/* CO: They may not come in the numerical recipes but it is not a bad
   idea to put them here, in fact they come from Gabor's group in Feb'84     */
double bessj0(double x)
{
    double ax, z;
    double xx, y, ans, ans1, ans2;

    if ((ax = fabs(x)) < 8.0)
    {
        y = x * x;
        ans1 = 57568490574.0 + y * (-13362590354.0 +
                                    y * (651619640.7
                                         + y * (-11214424.18 +
                                                y * (77392.33017 +
                                                     y * (-184.9052456)))));
        ans2 = 57568490411.0 + y * (1029532985.0 +
                                    y * (9494680.718
                                         + y * (59272.64853 +
                                                y * (267.8532712 +
                                                     y * 1.0))));
        ans = ans1 / ans2;
    }
    else
    {
        z = 8.0 / ax;
        y = z * z;
        xx = ax - 0.785398164;
        ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                          + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                                       + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                                                                      - y * 0.934935152e-7)));
        ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
    }
    return ans;
}

/*............................................................................*/
double bessi0(double x)
{
    double y, ax, ans;
    if ((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                                          + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
    }
    else
    {
        y = 3.75 / ax;
        ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
                                      + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
                                                                + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
                                                                                            + y * 0.392377e-2))))))));
    }
    return ans;
}

/*............................................................................*/
double bessi1(double x)
{
    double ax, ans;
    double y;
    if ((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                               + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
    }
    else
    {
        y = 3.75 / ax;
        ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
                                  - y * 0.420059e-2));
        ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
                                + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
        ans *= (exp(ax) / sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
}

/* General Bessel functions ------------------------------------------------ */
double chebev(double a, double b, double c[], int m, double x)
{
    double d = 0.0, dd = 0.0, sv, y, y2;
    int j;

    if ((x - a)*(x - b) > 0.0)
        nrerror("x not in range in routine chebev");
    y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a));
    for (j = m - 1;j >= 1;j--)
    {
        sv = d;
        d = y2 * d - dd + c[j];
        dd = sv;
    }
    return y*d - dd + 0.5*c[0];
}
#define NUSE1 5
#define NUSE2 5

void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
    double xx;
    static double c1[] =
        {
            -1.142022680371172e0, 6.516511267076e-3,
            3.08709017308e-4, -3.470626964e-6, 6.943764e-9,
            3.6780e-11, -1.36e-13
        };
    static double c2[] =
        {
            1.843740587300906e0, -0.076852840844786e0,
            1.271927136655e-3, -4.971736704e-6, -3.3126120e-8,
            2.42310e-10, -1.70e-13, -1.0e-15
        };

    xx = 8.0 * x * x - 1.0;
    *gam1 = chebev(-1.0, 1.0, c1, NUSE1, xx);
    *gam2 = chebev(-1.0, 1.0, c2, NUSE2, xx);
    *gampl = *gam2 - x * (*gam1);
    *gammi = *gam2 + x * (*gam1);
}

#undef NUSE1
#undef NUSE2

#define EPS 1.0e-16
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
void bessjy(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp)
{
    int i, isign, l, nl;
    double a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2,
    fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl,
    rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
    temp, w, x2, xi, xi2, xmu, xmu2;

    if (x <= 0.0 || xnu < 0.0)
        nrerror("bad arguments in bessjy");
    nl = (x < XMIN ? (int)(xnu + 0.5) : XMIPP_MAX(0, (int)(xnu - x + 1.5)));
    xmu = xnu - nl;
    xmu2 = xmu * xmu;
    xi = 1.0 / x;
    xi2 = 2.0 * xi;
    w = xi2 / PI;
    isign = 1;
    h = xnu * xi;
    if (h < FPMIN)
        h = FPMIN;
    b = xi2 * xnu;
    d = 0.0;
    c = h;
    for (i = 1;i <= MAXIT;i++)
    {
        b += xi2;
        d = b - d;
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = b - 1.0 / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        del = c * d;
        h = del * h;
        if (d < 0.0)
            isign = -isign;
        if (fabs(del - 1.0) < EPS)
            break;
    }
    if (i > MAXIT)
        nrerror("x too large in bessjy; try asymptotic expansion");
    rjl = isign * FPMIN;
    rjpl = h * rjl;
    rjl1 = rjl;
    rjp1 = rjpl;
    fact = xnu * xi;
    for (l = nl;l >= 1;l--)
    {
        rjtemp = fact * rjl + rjpl;
        fact -= xi;
        rjpl = fact * rjtemp - rjl;
        rjl = rjtemp;
    }
    if (rjl == 0.0)
        rjl = EPS;
    f = rjpl / rjl;
    if (x < XMIN)
    {
        x2 = 0.5 * x;
        pimu = PI * xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu / sin(pimu));
        d = -log(x2);
        e = xmu * d;
        fact2 = (fabs(e) < EPS ? 1.0 : sinh(e) / e);
        beschb(xmu, &gam1, &gam2, &gampl, &gammi);
        ff = 2.0 / PI * fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        e = exp(e);
        p = e / (gampl * PI);
        q = 1.0 / (e * PI * gammi);
        pimu2 = 0.5 * pimu;
        fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2) / pimu2);
        r = PI * pimu2 * fact3 * fact3;
        c = 1.0;
        d = -x2 * x2;
        sum = ff + r * q;
        sum1 = p;
        for (i = 1;i <= MAXIT;i++)
        {
            ff = (i * ff + p + q) / (i * i - xmu2);
            c *= (d / i);
            p /= (i - xmu);
            q /= (i + xmu);
            del = c * (ff + r * q);
            sum += del;
            del1 = c * p - i * del;
            sum1 += del1;
            if (fabs(del) < (1.0 + fabs(sum))*EPS)
                break;
        }
        if (i > MAXIT)
            nrerror("bessy series failed to converge");
        rymu = -sum;
        ry1 = -sum1 * xi2;
        rymup = xmu * xi * rymu - ry1;
        rjmu = w / (rymup - f * rymu);
    }
    else
    {
        a = 0.25 - xmu2;
        p = -0.5 * xi;
        q = 1.0;
        br = 2.0 * x;
        bi = 2.0;
        fact = a * xi / (p * p + q * q);
        cr = br + q * fact;
        ci = bi + p * fact;
        den = br * br + bi * bi;
        dr = br / den;
        di = -bi / den;
        dlr = cr * dr - ci * di;
        dli = cr * di + ci * dr;
        temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        for (i = 2;i <= MAXIT;i++)
        {
            a += 2 * (i - 1);
            bi += 2.0;
            dr = a * dr + br;
            di = a * di + bi;
            if (fabs(dr) + fabs(di) < FPMIN)
                dr = FPMIN;
            fact = a / (cr * cr + ci * ci);
            cr = br + cr * fact;
            ci = bi - ci * fact;
            if (fabs(cr) + fabs(ci) < FPMIN)
                cr = FPMIN;
            den = dr * dr + di * di;
            dr /= den;
            di /= -den;
            dlr = cr * dr - ci * di;
            dli = cr * di + ci * dr;
            temp = p * dlr - q * dli;
            q = p * dli + q * dlr;
            p = temp;
            if (fabs(dlr - 1.0) + fabs(dli) < EPS)
                break;
        }
        if (i > MAXIT)
            nrerror("cf2 failed in bessjy");
        gam = (p - f) / q;
        rjmu = sqrt(w / ((p - f) * gam + q));
        rjmu = NRSIGN(rjmu, rjl);
        rymu = rjmu * gam;
        rymup = rymu * (p + q / gam);
        ry1 = xmu * xi * rymu - rymup;
    }
    fact = rjmu / rjl;
    *rj = rjl1 * fact;
    *rjp = rjp1 * fact;
    for (i = 1;i <= nl;i++)
    {
        rytemp = (xmu + i) * xi2 * ry1 - rymu;
        rymu = ry1;
        ry1 = rytemp;
    }
    *ry = rymu;
    *ryp = xnu * xi * rymu - ry1;
}
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN

/*............................................................................*/
double bessi0_5(double x)
{
    return (x == 0) ? 0 : sqrt(2 / (PI*x))*sinh(x);
}
double bessi1_5(double x)
{
    return (x == 0) ? 0 : sqrt(2 / (PI*x))*(cosh(x) - sinh(x) / x);
}
double bessi2(double x)
{
    return (x == 0) ? 0 : bessi0(x) - ((2*1) / x) * bessi1(x);
}
double bessi2_5(double x)
{
    return (x == 0) ? 0 : bessi0_5(x) - ((2*1.5) / x) * bessi1_5(x);
}
double bessi3(double x)
{
    return (x == 0) ? 0 : bessi1(x) - ((2*2) / x) * bessi2(x);
}
double bessi3_5(double x)
{
    return (x == 0) ? 0 : bessi1_5(x) - ((2*2.5) / x) * bessi2_5(x);
}
double bessi4(double x)
{
    return (x == 0) ? 0 : bessi2(x) - ((2*3) / x) * bessi3(x);
}
double bessj1_5(double x)
{
    double rj, ry, rjp, ryp;
    bessjy(x, 1.5, &rj, &ry, &rjp, &ryp);
    return rj;
}
double bessj3_5(double x)
{
    double rj, ry, rjp, ryp;
    bessjy(x, 3.5, &rj, &ry, &rjp, &ryp);
    return rj;
}

/* Special functions ------------------------------------------------------- */
double gammln(double xx)
{
    double x, tmp, ser;
    static double cof[6] =
        {
            76.18009173, -86.50532033, 24.01409822,
            -1.231739516, 0.120858003e-2, -0.536382e-5
        };
    int j;

    x = xx - 1.0;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.0;
    for (j = 0;j <= 5;j++)
    {
        x += 1.0;
        ser += cof[j] / x;
    }
    return -tmp + log(2.50662827465*ser);
}


double betai(double a, double b, double x)
{
    double bt;
    if (x < 0.0 || x > 1.0)
        nrerror("Bad x in routine BETAI");
    if (x == 0.0 || x == 1.0)
        bt = 0.0;
    else
        bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt*betacf(a, b, x) / a;
    else
        return 1.0 -bt*betacf(b, a, 1.0 - x) / b;

}

#define ITMAX 100
#define EPS 3.0e-7
double betacf(double a, double b, double x)
{
    double qap, qam, qab, em, tem, d;
    double bz, bm = 1.0, bp, bpp;
    double az = 1.0, am = 1.0, ap, app, aold;
    int m;

    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    bz = 1.0 - qab * x / qap;
    for (m = 1;m <= ITMAX;m++)
    {
        em = (double) m;
        tem = em + em;
        d = em * (b - em) * x / ((qam + tem) * (a + tem));
        ap = az + d * am;
        bp = bz + d * bm;
        d = -(a + em) * (qab + em) * x / ((qap + tem) * (a + tem));
        app = ap + d * az;
        bpp = bp + d * bz;
        aold = az;
        am = ap / bpp;
        bm = bp / bpp;
        az = app / bpp;
        bz = 1.0;
        if (fabs(az - aold) < (EPS*fabs(az)))
            return az;
    }
    nrerror("a or b too big, or ITMAX too small in BETACF");
    return 0;
}
#undef ITMAX
#undef EPS

void instantiate_recipes()
{
    double **DD1;
    double *D1;

    double **FF1;
    double *F1;

    int **II1;
    int *I1;
    int i1 = 0;

    char *C1;

    ask_Tvector(D1, i1, i1);
    free_Tvector(D1, i1, i1);
    ask_Tvector(F1, i1, i1);
    free_Tvector(F1, i1, i1);
    ask_Tvector(I1, i1, i1);
    free_Tvector(I1, i1, i1);
    ask_Tvector(C1, i1, i1);
    free_Tvector(C1, i1, i1);

    ask_Tmatrix(DD1, i1, i1, i1, i1);
    free_Tmatrix(DD1, i1, i1, i1, i1);
    ask_Tmatrix(FF1, i1, i1, i1, i1);
    free_Tmatrix(FF1, i1, i1, i1, i1);
    ask_Tmatrix(II1, i1, i1, i1, i1);
    free_Tmatrix(II1, i1, i1, i1, i1);
}

/* Optimization ------------------------------------------------------------ */
#undef MAX
#undef SIGN
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define F1DIM(x,f) {\
    for (int j = 1; j<=ncom; j++) \
        xt[j] = pcom[j] + x * xicom[j]; \
    f = (*func)(xt,prm);}

void mnbrak(double *ax, double *bx, double *cx,
            double *fa, double *fb, double *fc, double(*func)(double *, void*),
            void *prm, int ncom, double *pcom, double *xicom)
{
    double ulim, u, r, q, fu, dum;
    double *xt=NULL;
    ask_Tvector(xt, 1, ncom);

    F1DIM(*ax,*fa);
    F1DIM(*bx,*fb);
    if (*fb > *fa)
    {
        SHFT(dum, *ax, *bx, dum)
        SHFT(dum, *fb, *fa, dum)
    }
    *cx = (*bx) + GOLD * (*bx - *ax);
    F1DIM(*cx,*fc);
    while (*fb > *fc)
    {
        r = (*bx - *ax) * (*fb - *fc);
        q = (*bx - *cx) * (*fb - *fa);
        u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
            (2.0 * SIGN(MAX(fabs(q - r), TINY), q - r));
        ulim = (*bx) + GLIMIT * (*cx - *bx);
        if ((*bx - u)*(u - *cx) > 0.0)
        {
            F1DIM(u,fu);
            if (fu < *fc)
            {
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                return;
            }
            else if (fu > *fb)
            {
                *cx = u;
                *fc = fu;
                return;
            }
            u = (*cx) + GOLD * (*cx - *bx);
            F1DIM(u,fu);
        }
        else if ((*cx - u)*(u - ulim) > 0.0)
        {
            F1DIM(u,fu);
            if (fu < *fc)
            {
                SHFT(*bx, *cx, u, *cx + GOLD*(*cx - *bx))
                double aux;
                F1DIM(u,aux);
                SHFT(*fb, *fc, fu, aux)
            }
        }
        else if ((u - ulim)*(ulim - *cx) >= 0.0)
        {
            u = ulim;
            F1DIM(u,fu);
        }
        else
        {
            u = (*cx) + GOLD * (*cx - *bx);
            F1DIM(u,fu);
        }
        SHFT(*ax, *bx, *cx, u)
        SHFT(*fa, *fb, *fc, fu)
    }
    free_Tvector(xt, 1, ncom);
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
double brent(double ax, double bx, double cx, double(*func)(double *,void*),
             void *prm, double tol, double *xmin,
             int ncom, double *pcom, double *xicom)
{
    int iter;
    double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double e = 0.0;
    double *xt=NULL;
    ask_Tvector(xt, 1, ncom);

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;
    F1DIM(x,fx);
    fw = fv = fx;
    for (iter = 1;iter <= ITMAX;iter++)
    {
        xm = 0.5 * (a + b);
        tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
        if (fabs(x - xm) <= (tol2 - 0.5*(b - a)))
        {
            *xmin = x;
            free_Tvector(xt, 1, ncom);
            return fx;
        }
        if (fabs(e) > tol1)
        {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
                p = -p;
            q = fabs(q);
            etemp = e;
            e = d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a - x) || p >= q*(b - x))
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            else
            {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = SIGN(tol1, xm - x);
            }
        }
        else
        {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }
        u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
        F1DIM(u,fu);
        if (fu <= fx)
        {
            if (u >= x)
                a = x;
            else
                b = x;
            SHFT(v, w, x, u)
            SHFT(fv, fw, fx, fu)
        }
        else
        {
            if (u < x)
                a = u;
            else
                b = u;
            if (fu <= fw || w == x)
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }
    nrerror("Too many iterations in brent");
    *xmin = x;
    free_Tvector(xt, 1, ncom);
    return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef F1DIM

#define TOL 2.0e-4
void linmin(double *p, double *xi, int n, double &fret,
            double(*func)(double *, void*), void *prm)
{
    int j;
    double xx, xmin, fx, fb, fa, bx, ax;

    int ncom = n;
    double *pcom=NULL;
    double *xicom=NULL;
    ask_Tvector(pcom, 1, n);
    ask_Tvector(xicom, 1, n);
    for (j = 1;j <= n;j++)
    {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }
    ax = 0.0;
    xx = 1.0;
    bx = 2.0;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, func, prm, ncom, pcom, xicom);
    fret = brent(ax, xx, bx, func, prm, TOL, &xmin, ncom, pcom, xicom);
    for (j = 1;j <= n;j++)
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_Tvector(xicom, 1, n);
    free_Tvector(pcom, 1, n);
}
#undef TOL

#define ITMAX 200
void powell(double *p, double *xi, int n, double ftol, int &iter,
            double &fret, double(*func)(double *, void *), void *prm,
            bool show)
{
    int i, ibig, j;
    double t, fptt, fp, del;
    double *pt, *ptt, *xit;
    bool   different_from_0;

    ask_Tvector(pt, 1, n);
    ask_Tvector(ptt, 1, n);
    ask_Tvector(xit, 1, n);
    fret = (*func)(p,prm);
    for (j = 1;j <= n;j++)
        pt[j] = p[j];

    for (iter = 1;;(iter)++)
    {
        /* By coss ----- */
        if (show)
        {
            std::cout << iter << " (" << p[1];
            for (int co = 2; co <= n; co++)
                std::cout << "," << p[co];
            std::cout << ")--->" << fret << std::endl;
        }
        /* ------------- */

        fp = fret;
        ibig = 0;
        del = 0.0;
        for (i = 1;i <= n;i++)
        {
            different_from_0 = false; // CO
            for (j = 1;j <= n;j++)
            {
                xit[j] = xi[j*n+i];
                if (xit[j] != 0)
                    different_from_0 = true;
            }
            if (different_from_0)
            {
                fptt = fret;
                linmin(p, xit, n, fret, func, prm);
                if (fabs(fptt - fret) > del)
                {
                    del = fabs(fptt - fret);
                    ibig = i;
                }
                /* By coss ----- */
                if (show)
                {
                    std::cout << "   (";
                    if (i == 1)
                        std::cout << "***";
                    std::cout << p[1];
                    for (int co = 2; co <= n; co++)
                    {
                        std::cout << ",";
                        if (co == i)
                            std::cout << "***";
                        std::cout << p[co];
                    }
                    std::cout << ")--->" << fret << std::endl;
                }
                /* ------------- */
            }
        }
        if (2.0*fabs(fp - fret) <= ftol*(fabs(fp) + fabs(fret)) || n==1)
        {
            free_Tvector(xit, 1, n);
            free_Tvector(ptt, 1, n);
            free_Tvector(pt, 1, n);
            return;
        }
        if (iter == ITMAX)
            nrerror("Too many iterations in routine POWELL");
        for (j = 1;j <= n;j++)
        {
            ptt[j] = 2.0 * p[j] - pt[j];
            xit[j] = p[j] - pt[j];
            pt[j] = p[j];
        }
        fptt = (*func)(ptt,prm);
        if (fptt < fp)
        {
#define SQR(a) ((a)*(a))
            t = 2.0 * (fp - 2.0 * fret + fptt) * SQR(fp - fret - del) - del * SQR(fp - fptt);
            if (t < 0.0)
            {
                linmin(p, xit, n, fret, func, prm);
                for (j = 1;j <= n;j++)
                    xi[j*n+ibig] = xit[j];
            }
        }
    }
}
#undef ITMAX
#undef SQR

/* Non linear least squares ------------------------------------------------ */
// These routines have been taken from
// http://users.utu.fi/vesoik/userdocs/programs/libpet
// and they implement an algorithm of Lawson-Hanson of
// nonnegative least squares

/* Example of use:
   double a[]={ 5, 0, -2,
               0, 3,  0,
               1, 1, -1,
              -1, 1, -1,
               9, 9, -9};
   double b[]={1, 9, -1};
   double x[5];
   double rnorm;
   int i;

   int success=nnls(a,3,5,b,x,&rnorm,NULL,NULL,NULL);
   printf("success=%d\n",success);
   printf("rnorm=%d\n",rnorm);
   for (i=0; i<5; i++)
      printf("%f\n",x[i]);

   In this case: x=0 2.666 0 0 0.111

   This program resolves A^t*x=b subject to x>=0.
   In terms of basis vectors, the rows of A are the basis axes, b is
   the vector we want to represent in the subspace spanned by the rows of A
   and x are the nonnegative coordinates of the representation of b in A.
*/

/*****************************************************************************
 *
 *  Compute orthogonal rotation matrix:
 *    (C, S) so that (C, S)(A) = (sqrt(A**2+B**2))
 *    (-S,C)         (-S,C)(B)   (   0          )
 *  Compute sig = sqrt(A**2+B**2):
 *    sig is computed last to allow for the possibility that sig may be in
 *    the same location as A or B.
 */
void _nnls_g1(double a, double b, double *cterm, double *sterm, double *sig)
{
    double d1, xr, yr;

    if (fabs(a) > fabs(b))
    {
        xr = b / a;
        d1 = xr;
        yr = sqrt(d1 * d1 + 1.);
        d1 = 1. / yr;
        *cterm = (a >= 0.0 ? fabs(d1) : -fabs(d1));
        *sterm = (*cterm) * xr;
        *sig = fabs(a) * yr;
    }
    else if (b != 0.)
    {
        xr = a / b;
        d1 = xr;
        yr = sqrt(d1 * d1 + 1.);
        d1 = 1. / yr;
        *sterm = (b >= 0.0 ? fabs(d1) : -fabs(d1));
        *cterm = (*sterm) * xr;
        *sig = fabs(b) * yr;
    }
    else
    {
        *sig = 0.;
        *cterm = 0.;
        *sterm = 1.;
    }
} /* _nnls_g1 */
/****************************************************************************/

/*****************************************************************************
 *
 *  Construction and/or application of a single Householder transformation:
 *           Q = I + U*(U**T)/B
 *
 *  Function returns 0 if succesful, or >0 in case of erroneous parameters.
 *
 */
int _nnls_h12(
    int mode,
    /* mode=1 to construct and apply a Householder transformation, or
       mode=2 to apply a previously constructed transformation */
    int lpivot,     /* Index of the pivot element */
    int l1, int m,
    /* Transformation is constructed to zero elements indexed from l1 to M */
    double *u, int u_dim1, double *up,
    /* With mode=1: On entry, u[] must contain the pivot vector.
       On exit, u[] and up contain quantities defining the vector u[] of
       the Householder transformation. */
    /* With mode=2: On entry, u[] and up should contain quantities previously
       computed with mode=1. These will not be modified. */
    /* u_dim1 is the storage increment between elements. */
    double *cm,
    /* On entry, cm[] must contain the matrix (set of vectors) to which the
       Householder transformation is to be applied. On exit, cm[] will contain
       the set of transformed vectors */
    int ice,        /* Storage increment between elements of vectors in cm[] */
    int icv,        /* Storage increment between vectors in cm[] */
    int ncv         /* Nr of vectors in cm[] to be transformed;
                                     if ncv<=0, then no operations will be done on cm[] */
)
{
    double d1, d2, b, clinv, cl, sm;
    int incr, k, j, i2, i3, i4;

    /* Check parameters */
    if (mode != 1 && mode != 2)
        return(1);
    if (m < 1 || u == NULL || u_dim1 < 1 || cm == NULL)
        return(2);
    if (lpivot < 0 || lpivot >= l1 || l1 >= m)
        return(0);
    /* Function Body */
    cl = (d1 = u[lpivot*u_dim1], fabs(d1));
    if (mode == 2)
    { /* Apply transformation I+U*(U**T)/B to cm[] */
        if (cl <= 0.)
            return(0);
    }
    else
    { /* Construct the transformation */
        for (j = l1; j < m; j++)
        { /* Computing MAX */
            d2 = (d1 = u[j*u_dim1], fabs(d1));
            if (d2 > cl)
                cl = d2;
        }
        if (cl <= 0.)
            return(0);
        clinv = 1.0 / cl;
        /* Computing 2nd power */
        d1 = u[lpivot*u_dim1] * clinv;
        sm = d1 * d1;
        for (j = l1; j < m; j++)
        {
            d1 = u[j*u_dim1] * clinv;
            sm += d1 * d1;
        }
        cl *= sqrt(sm);
        if (u[lpivot*u_dim1] > 0.)
            cl = -cl;
        *up = u[lpivot*u_dim1] - cl;
        u[lpivot*u_dim1] = cl;
    }
    if (ncv <= 0)
        return(0);
    b = (*up) * u[lpivot*u_dim1];
    /* b must be nonpositive here; if b>=0., then return */
    if (b >= 0.)
        return(0);
    b = 1.0 / b;
    i2 = 1 - icv + ice * lpivot;
    incr = ice * (l1 - lpivot);
    for (j = 0; j < ncv; j++)
    {
        i2 += icv;
        i3 = i2 + incr;
        i4 = i3;
        sm = cm[i2-1] * (*up);
        for (k = l1; k < m; k++)
        {
            sm += cm[i3-1] * u[k*u_dim1];
            i3 += ice;
        }
        if (sm != 0.0)
        {
            sm *= b;
            cm[i2-1] += sm * (*up);
            for (k = l1; k < m; k++)
            {
                cm[i4-1] += sm * u[k*u_dim1];
                i4 += ice;
            }
        }
    }
    return(0);
} /* _nnls_h12 */

/*****************************************************************************
 *  Algorithm NNLS (Non-negative least-squares)
 *
 *  Given an m by n matrix A, and an m-vector B, computes an n-vector X,
 *  that solves the least squares problem
 *      A * X = B   , subject to X>=0
 *
 *  Function returns 0 if succesful, 1, if iteration count exceeded 3*N,
 *  or 2 in case of invalid problem dimensions or memory allocation error.
 *
 *  Instead of pointers for working space, NULL can be given to let this
 *  function to allocate and free the required memory.
 */
int nnls(
    double *a, int m, int n,
    /* On entry, a[n][m] contains the m by n matrix A. On exit, a[][] contains
       the product matrix Q*A, where Q is an m by n orthogonal matrix generated
       implicitly by this function.*/
    double *b,
    /* On entry, b[] must contain the m-vector B.
       On exit, b[] contains Q*B */
    double *x,
    /* On exit, x[] will contain the solution vector */
    double *rnorm,
    /* On exit, rnorm contains the Euclidean norm of the residual vector */
    double *wp,  /* An n-array of working space, w[]. */
    /* On exit, w[] will contain the dual solution vector.
       w[i]=0.0 for all i in set p and w[i]<=0.0 for all i in set z. */
    double *zzp, /* An m-array of working space, zz[]. */
    int *indexp  /* An n-array of working space, index[]. */
)
{
    int pfeas, ret = 0, iz, jz, iz1, iz2, npp1, *index;
    double d1, d2, sm, up, ss, *w, *zz;
    int iter, k, j = 0, l, itmax, izmax = 0, nsetp, ii, jj = 0, ip;
    double temp, wmax, t, alpha, asave, dummy, unorm, ztest, cc;


    /* Check the parameters and data */
    if (m <= 0 || n <= 0 || a == NULL || b == NULL || x == NULL)
        return(2);
    /* Allocate memory for working space, if required */
    if (wp != NULL)
        w = wp;
    else
        w = (double*)calloc(n, sizeof(double));
    if (zzp != NULL)
        zz = zzp;
    else
        zz = (double*)calloc(m, sizeof(double));
    if (indexp != NULL)
        index = indexp;
    else
        index = (int*)calloc(n, sizeof(int));
    if (w == NULL || zz == NULL || index == NULL)
        return(2);

    /* Initialize the arrays INDEX[] and X[] */
    for (k = 0; k < n; k++)
    {
        x[k] = 0.;
        index[k] = k;
    }
    iz2 = n - 1;
    iz1 = 0;
    nsetp = 0;
    npp1 = 0;

    /* Main loop; quit if all coeffs are already in the solution or */
    /* if M cols of A have been triangularized */
    iter = 0;
    itmax = n * 3;
    while (iz1 <= iz2 && nsetp < m)
    {
        /* Compute components of the dual (negative gradient) vector W[] */
        for (iz = iz1; iz <= iz2; iz++)
        {
            j = index[iz];
            sm = 0.;
            for (l = npp1; l < m; l++)
                sm += a[j*m+l] * b[l];
            w[j] = sm;
        }

        while (1)
        {
            /* Find largest positive W[j] */
            for (wmax = 0., iz = iz1; iz <= iz2; iz++)
            {
                j = index[iz];
                if (w[j] > wmax)
                {
                    wmax = w[j];
                    izmax = iz;
                }
            }

            /* Terminate if wmax<=0.; */
            /* it indicates satisfaction of the Kuhn-Tucker conditions */
            if (wmax <= 0.0)
                break;
            iz = izmax;
            j = index[iz];

            /* The sign of W[j] is ok for j to be moved to set P. */
            /* Begin the transformation and check new diagonal element to avoid */
            /* near linear dependence. */
            asave = a[j*m+npp1];
            _nnls_h12(1, npp1, npp1 + 1, m, &a[j*m+0], 1, &up, &dummy, 1, 1, 0);
            unorm = 0.;
            if (nsetp != 0)
                for (l = 0; l < nsetp; l++)
                {
                    d1 = a[j*m+l];
                    unorm += d1 * d1;
                }
            unorm = sqrt(unorm);
            d2 = unorm + (d1 = a[j*m+npp1], fabs(d1)) * 0.01;
            if ((d2 - unorm) > 0.)
            {
                /* Col j is sufficiently independent. Copy B into ZZ, update ZZ */
                /* and solve for ztest ( = proposed new value for X[j] ) */
                for (l = 0; l < m; l++)
                    zz[l] = b[l];
                _nnls_h12(2, npp1, npp1 + 1, m, &a[j*m+0], 1, &up, zz, 1, 1, 1);
                ztest = zz[npp1] / a[j*m+npp1];
                /* See if ztest is positive */
                if (ztest > 0.)
                    break;
            }

            /* Reject j as a candidate to be moved from set Z to set P. Restore */
            /* A[npp1,j], set W[j]=0., and loop back to test dual coeffs again */
            a[j*m+npp1] = asave;
            w[j] = 0.;
        } /* while(1) */
        if (wmax <= 0.0)
            break;

        /* Index j=INDEX[iz] has been selected to be moved from set Z to set P. */
        /* Update B and indices, apply householder transformations to cols in */
        /* new set Z, zero subdiagonal elts in col j, set W[j]=0. */
        for (l = 0; l < m; ++l)
            b[l] = zz[l];
        index[iz] = index[iz1];
        index[iz1] = j;
        iz1++;
        nsetp = npp1 + 1;
        npp1++;
        if (iz1 <= iz2)
            for (jz = iz1; jz <= iz2; jz++)
            {
                jj = index[jz];
                _nnls_h12(2, nsetp - 1, npp1, m, &a[j*m+0], 1, &up,
                          &a[jj*m+0], 1, m, 1);
            }
        if (nsetp != m)
            for (l = npp1; l < m; l++)
                a[j*m+l] = 0.;
        w[j] = 0.;
        /* Solve the triangular system; store the solution temporarily in Z[] */
        for (l = 0; l < nsetp; l++)
        {
            ip = nsetp - (l + 1);
            if (l != 0)
                for (ii = 0; ii <= ip; ii++)
                    zz[ii] -= a[jj*m+ii] * zz[ip+1];
            jj = index[ip];
            zz[ip] /= a[jj*m+ip];
        }

        /* Secondary loop begins here */
        while (++iter < itmax)
        {
            /* See if all new constrained coeffs are feasible; if not, compute alpha */
            for (alpha = 2.0, ip = 0; ip < nsetp; ip++)
            {
                l = index[ip];
                if (zz[ip] <= 0.)
                {
                    t = -x[l] / (zz[ip] - x[l]);
                    if (alpha > t)
                    {
                        alpha = t;
                        jj = ip - 1;
                    }
                }
            }

            /* If all new constrained coeffs are feasible then still alpha==2. */
            /* If so, then exit from the secondary loop to main loop */
            if (alpha == 2.0)
                break;
            /* Use alpha (0.<alpha<1.) to interpolate between old X and new ZZ */
            for (ip = 0; ip < nsetp; ip++)
            {
                l = index[ip];
                x[l] += alpha * (zz[ip] - x[l]);
            }

            /* Modify A and B and the INDEX arrays to move coefficient i */
            /* from set P to set Z. */
            k = index[jj+1];
            pfeas = 1;
            do
            {
                x[k] = 0.;
                if (jj != (nsetp - 1))
                {
                    jj++;
                    for (j = jj + 1; j < nsetp; j++)
                    {
                        ii = index[j];
                        index[j-1] = ii;
                        _nnls_g1(a[ii*m+j-1], a[ii*m+j], &cc, &ss, &a[ii*m+j-1]);
                        for (a[ii*m+j] = 0., l = 0; l < n; l++)
                            if (l != ii)
                            {
                                /* Apply procedure G2 (CC,SS,A(J-1,L),A(J,L)) */
                                temp = a[l*m+j-1];
                                a[l*m+j-1] = cc * temp + ss * a[l*m+j];
                                a[l*m+j] = -ss * temp + cc * a[l*m+j];
                            }
                        /* Apply procedure G2 (CC,SS,B(J-1),B(J)) */
                        temp = b[j-1];
                        b[j-1] = cc * temp + ss * b[j];
                        b[j] = -ss * temp + cc * b[j];
                    }
                }
                npp1 = nsetp - 1;
                nsetp--;
                iz1--;
                index[iz1] = k;

                /* See if the remaining coeffs in set P are feasible; they should */
                /* be because of the way alpha was determined. If any are */
                /* infeasible it is due to round-off error. Any that are */
                /* nonpositive will be set to zero and moved from set P to set Z */
                for (jj = 0; jj < nsetp; jj++)
                {
                    k = index[jj];
                    if (x[k] <= 0.)
                    {
                        pfeas = 0;
                        break;
                    }
                }
            }
            while (pfeas == 0);

            /* Copy B[] into zz[], then solve again and loop back */
            for (k = 0; k < m; k++)
                zz[k] = b[k];
            for (l = 0; l < nsetp; l++)
            {
                ip = nsetp - (l + 1);
                if (l != 0)
                    for (ii = 0; ii <= ip; ii++)
                        zz[ii] -= a[jj*m+ii] * zz[ip+1];
                jj = index[ip];
                zz[ip] /= a[jj*m+ip];
            }
        } /* end of secondary loop */
        if (iter > itmax)
        {
            ret = 1;
            break;
        }
        for (ip = 0; ip < nsetp; ip++)
        {
            k = index[ip];
            x[k] = zz[ip];
        }
    } /* end of main loop */
    /* Compute the norm of the final residual vector */
    sm = 0.;
    if (npp1 < m)
        for (k = npp1; k < m; k++)
            sm += (b[k] * b[k]);
    else
        for (j = 0; j < n; j++)
            w[j] = 0.;
    *rnorm = sqrt(sm);
    /* Free working space, if it was allocated here */
    if (wp == NULL)
        free(w);
    if (zzp == NULL)
        free(zz);
    if (indexp == NULL)
        free(index);
    return(ret);
} /* nnls_ */
/****************************************************************************/
/****************************************************************************/
/*
  nnlsWght()

  Algorithm for weighting the problem that is given to nnls-algorithm.
  Square roots of weights are used because in nnls the difference
  w*A-w*b is squared.
  Algorithm returns zero if successful, 1 if arguments are inappropriate.

*/
int nnlsWght(int N, int M, double *A, double *b, double *weight)
{
    int n, m;
    double *w;

    /* Check the arguments */
    if (N < 1 || M < 1 || A == NULL || b == NULL || weight == NULL)
        return(1);

    /* Allocate memory */
    w = (double*)malloc(M * sizeof(double));
    if (w == NULL)
        return(2);

    /* Check that weights are not zero and get the square roots of them to w[] */
    for (m = 0; m < M; m++)
    {
        if (weight[m] <= 1.0e-20)
            w[m] = 0.0;
        else
            w[m] = sqrt(weight[m]);
    }

    /* Multiply rows of matrix A and elements of vector b with weights*/
    for (m = 0; m < M; m++)
    {
        for (n = 0; n < N; n++)
        {
            A[n*M+m] *= w[m];
        }
        b[m] *= w[m];
    }

    free(w);
    return(0);
}
/****************************************************************************/

/* Singular value descomposition ------------------------------------------- */
/* Copied from Bilib library (linearalgebra.h) */
double Pythag(double a, double b)
{
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absb < absa)
        return(absa * sqrt(1.0 + absb * absb / (absa * absa)));
    else
        return((absb == 0.0) ? (0.0) : (absb * sqrt(1.0 + absa * absa / (absb * absb))));
}

#define SVDMAXITER 1000
void svdcmp(double *U, int Lines, int Columns, double *W, double *V)
{
    double *rv1 = (double *)NULL;
    double Norm, Scale;
    double c, f, g, h, s;
    double x, y, z;
    long i, its, j, jj, k, l = 0L, nm = 0L;
    bool    Flag;
    int     MaxIterations = SVDMAXITER;

    ask_Tvector(rv1, 0, Columns*Columns - 1);
    g = Scale = Norm = 0.0;
    for (i = 0L; (i < Columns); i++)
    {
        l = i + 1L;
        rv1[i] = Scale * g;
        g = s = Scale = 0.0;
        if (i < Lines)
        {
            for (k = i; (k < Lines); k++)
            {
                Scale += fabs(U[k * Columns + i]);
            }
            if (Scale != 0.0)
            {
                for (k = i; (k < Lines); k++)
                {
                    U[k * Columns + i] /= Scale;
                    s += U[k * Columns + i] * U[k * Columns + i];
                }
                f = U[i * Columns + i];
                g = (0.0 <= f) ? (-sqrt(s)) : (sqrt(s));
                h = f * g - s;
                U[i * Columns + i] = f - g;
                for (j = l; (j < Columns); j++)
                {
                    for (s = 0.0, k = i; (k < Lines); k++)
                    {
                        s += U[k * Columns + i] * U[k * Columns + j];
                    }
                    f = s / h;
                    for (k = i; (k < Lines); k++)
                    {
                        U[k * Columns + j] += f * U[k * Columns + i];
                    }
                }
                for (k = i; (k < Lines); k++)
                {
                    U[k * Columns + i] *= Scale;
                }
            }
        }
        W[i] = Scale * g;
        g = s = Scale = 0.0;
        if ((i < Lines) && (i != (Columns - 1L)))
        {
            for (k = l; (k < Columns); k++)
            {
                Scale += fabs(U[i * Columns + k]);
            }
            if (Scale != 0.0)
            {
                for (k = l; (k < Columns); k++)
                {
                    U[i * Columns + k] /= Scale;
                    s += U[i * Columns + k] * U[i * Columns + k];
                }
                f = U[i * Columns + l];
                g = (0.0 <= f) ? (-sqrt(s)) : (sqrt(s));
                h = f * g - s;
                U[i * Columns + l] = f - g;
                for (k = l; (k < Columns); k++)
                {
                    rv1[k] = U[i * Columns + k] / h;
                }
                for (j = l; (j < Lines); j++)
                {
                    for (s = 0.0, k = l; (k < Columns); k++)
                    {
                        s += U[j * Columns + k] * U[i * Columns + k];
                    }
                    for (k = l; (k < Columns); k++)
                    {
                        U[j * Columns + k] += s * rv1[k];
                    }
                }
                for (k = l; (k < Columns); k++)
                {
                    U[i * Columns + k] *= Scale;
                }
            }
        }
        Norm = ((fabs(W[i]) + fabs(rv1[i])) < Norm) ? (Norm) : (fabs(W[i]) + fabs(rv1[i]));
    }
    for (i = Columns - 1L; (0L <= i); i--)
    {
        if (i < (Columns - 1L))
        {
            if (g != 0.0)
            {
                for (j = l; (j < Columns); j++)
                {
                    V[j * Columns + i] = U[i * Columns + j] / (U[i * Columns + l] * g);
                }
                for (j = l; (j < Columns); j++)
                {
                    for (s = 0.0, k = l; (k < Columns); k++)
                    {
                        s += U[i * Columns + k] * V[k * Columns + j];
                    }
                    for (k = l; (k < Columns); k++)
                    {
                        if (s != 0.0)
                        {
                            V[k * Columns + j] += s * V[k * Columns + i];
                        }
                    }
                }
            }
            for (j = l; (j < Columns); j++)
            {
                V[i * Columns + j] = V[j * Columns + i] = 0.0;
            }
        }
        V[i * Columns + i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for (i = (Lines < Columns) ? (Lines - 1L) : (Columns - 1L); (0L <= i); i--)
    {
        l = i + 1L;
        g = W[i];
        for (j = l; (j < Columns); j++)
        {
            U[i * Columns + j] = 0.0;
        }
        if (g != 0.0)
        {
            g = 1.0 / g;
            for (j = l; (j < Columns); j++)
            {
                for (s = 0.0, k = l; (k < Lines); k++)
                {
                    s += U[k * Columns + i] * U[k * Columns + j];
                }
                f = s * g / U[i * Columns + i];
                for (k = i; (k < Lines); k++)
                {
                    if (f != 0.0)
                    {
                        U[k * Columns + j] += f * U[k * Columns + i];
                    }
                }
            }
            for (j = i; (j < Lines); j++)
            {
                U[j * Columns + i] *= g;
            }
        }
        else
        {
            for (j = i; (j < Lines); j++)
            {
                U[j * Columns + i] = 0.0;
            }
        }
        U[i * Columns + i] += 1.0;
    }
    for (k = Columns - 1L; (0L <= k); k--)
    {
        for (its = 1L; (its <= MaxIterations); its++)
        {
            Flag = true;
            for (l = k; (0L <= l); l--)
            {
                nm = l - 1L;
                if ((fabs(rv1[l]) + Norm) == Norm)
                {
                    Flag = false;
                    break;
                }
                if ((fabs(W[nm]) + Norm) == Norm)
                {
                    break;
                }
            }
            if (Flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; (i <= k); i++)
                {
                    f = s * rv1[i];
                    rv1[i] *= c;
                    if ((fabs(f) + Norm) == Norm)
                    {
                        break;
                    }
                    g = W[i];
                    h = Pythag(f, g);
                    W[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 0L; (j < Lines); j++)
                    {
                        y = U[j * Columns + nm];
                        z = U[j * Columns + i];
                        U[j * Columns + nm] = y * c + z * s;
                        U[j * Columns + i] = z * c - y * s;
                    }
                }
            }
            z = W[k];
            if (l == k)
            {
                if (z < 0.0)
                {
                    W[k] = -z;
                    for (j = 0L; (j < Columns); j++)
                    {
                        V[j * Columns + k] = -V[j * Columns + k];
                    }
                }
                break;
            }
            if (its == MaxIterations)
            {
                free_Tvector(rv1, 0, Columns*Columns - 1);
                return;
            }
            x = W[l];
            nm = k - 1L;
            y = W[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = Pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + ((0.0 <= f) ? (fabs(g))
                                                : (-fabs(g))))) - h)) / x;
            c = s = 1.0;
            for (j = l; (j <= nm); j++)
            {
                i = j + 1L;
                g = rv1[i];
                y = W[i];
                h = s * g;
                g = c * g;
                z = Pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 0L; (jj < Columns); jj++)
                {
                    x = V[jj * Columns + j];
                    z = V[jj * Columns + i];
                    V[jj * Columns + j] = x * c + z * s;
                    V[jj * Columns + i] = z * c - x * s;
                }
                z = Pythag(f, h);
                W[j] = z;
                if (z != 0.0)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 0L; (jj < Lines); jj++)
                {
                    y = U[jj * Columns + j];
                    z = U[jj * Columns + i];
                    U[jj * Columns + j] = y * c + z * s;
                    U[jj * Columns + i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            W[k] = x;
        }
    }
    free_Tvector(rv1, 0, Columns*Columns - 1);
}

void svbksb(double *u, double *w, double *v, int m, int n, double *b, double *x)
{
    int jj, j, i;
    double s, *tmp;

    ask_Tvector(tmp, 1, n);
    for (j = 1;j <= n;j++)
    {
        s = 0.0;
        if (w[j])
        {
            for (i = 1;i <= m;i++)
                s += u[i*n+j] * b[i];
            s /= w[j];
        }
        tmp[j] = s;
    }
    for (j = 1;j <= n;j++)
    {
        s = 0.0;
        for (jj = 1;jj <= n;jj++)
            s += v[j*n+jj] * tmp[jj];
        x[j] = s;
    }
    free_Tvector(tmp, 1, n);
}

/* Savitzky-Golay filter coefficients. ------------------------------------- */
// This routine is used in Xmipp to perform Numerical Derivatives of equally spaced
// data
void savgol(double *c, int np, int nl, int nr, int ld, int m)
/* Returns in c[1..np], in wrap-around order a set of Savitzky-Golay
 filter coeficients. nl is the number of leftward (past) data points used,
 while nr is the number of rightward (future) data points,
 making the total number of data points used nl +nr +1.
 ld is the order of the derivative desired (e.g., ld = 0 for smoothed function).
 m is the order of the smoothing polynomial, also equal to the highest conserved
 moment; usual values are m = 2or m = 4. */
{
    int imj, ipj, j, k, kk, mm, *indx;
    double d, fac, sum, *a, *b;

    if (np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m)
        nrerror("SAVGOL: bad arguments");
    ask_Tvector(indx, 1, m + 1);
    ask_Tvector(a, 1, (m + 1)*(m + 1));
    ask_Tvector(b, 1, m + 1);
    // Set up the normal equations of the desired least-squares fit.
    for (ipj = 0;ipj <= (m << 1);ipj++)
    {
        sum = (ipj ? 0.0 : 1.0);
        for (k = 1;k <= nr;k++)
            sum += pow((double)k, (double)ipj);
        for (k = 1;k <= nl;k++)
            sum += pow((double) - k, (double)ipj);
        mm = XMIPP_MIN(ipj, 2 * m - ipj);
        for (imj = -mm;imj <= mm;imj += 2)
            a[(1+(ipj+imj)/2)*(m+1)+1+(ipj-imj)/2] = sum;
    }
    // Solve them: LU decomposition.
    ludcmp(a, m + 1, indx, &d);

    for (j = 1;j <= m + 1;j++)
        b[j] = 0.0;
    b[ld+1] = 1.0;
    // Right-hand side vector is unit vector, depending on which derivative we want.
    lubksb(a, m + 1, indx, b); // Get one row of the inverse matrix.
    for (kk = 1;kk <= np;kk++)
        c[kk] = 0.0; // Zero the output array (it may be bigger than number of coeficients).
    for (k = -nl;k <= nr;k++)
    {
        sum = b[1]; // Each Savitzky-Golay coeficient is the dot product of powers of an integer with the inverse matrix row.
        fac = 1.0;
        for (mm = 1;mm <= m;mm++)
            sum += b[mm+1] * (fac *= k);
        kk = ((np - k) % np) + 1; // Store in wrap-around order.
        c[kk] = sum;
    }

    free_Tvector(b, 1, m + 1);
    free_Tvector(a, 1, (m + 1)*(m + 1));
    free_Tvector(indx, 1, m + 1);
}

// CFSQP -------------------------------------------------------------------

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
int x_is_new = TRUE;

/* Declare and initialize user-accessible stopping criterion */
double objeps = -1.e0;
double objrep = -1.e0;
double gLgeps = -1.e0;
extern int nstop;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                     !!!! NOTICE !!!!

1. The routines contained in this file are due to Prof. K.Schittkowski
    of the University of Bayreuth, Germany (modification of routines
    due to Prof. MJD Powell at the University of Cambridge).  They can
    be freely distributed.

2. A few minor modifications were performed at the University of
    Maryland. They are marked in the code by "umd".

                                      A.L. Tits, J.L. Zhou, and
          Craig Lawrence
                                      University of Maryland

 ***********************************************************************



             SOLUTION OF QUADRATIC PROGRAMMING PROBLEMS



   QL0001 SOLVES THE QUADRATIC PROGRAMMING PROBLEM

   MINIMIZE        .5*X'*C*X + D'*X
   SUBJECT TO      A(J)*X  +  B(J)   =  0  ,  J=1,...,ME
                   A(J)*X  +  B(J)  >=  0  ,  J=ME+1,...,M
                   XL  <=  X  <=  XU

HERE C MUST BE AN N BY N SYMMETRIC AND POSITIVE MATRIX, D AN N-DIMENSIONAL
VECTOR, A AN M BY N MATRIX AND B AN M-DIMENSIONAL VECTOR. THE ABOVE
SITUATION IS INDICATED BY IWAR(1)=1. ALTERNATIVELY, I.E. IF IWAR(1)=0,
THE OBJECTIVE FUNCTION MATRIX CAN ALSO BE PROVIDED IN FACTORIZED FORM.
IN THIS CASE, C IS AN UPPER TRIANGULAR MATRIX.

THE SUBROUTINE REORGANIZES SOME DATA SO THAT THE PROBLEM CAN BE SOLVED
BY A MODIFICATION OF AN ALGORITHM PROPOSED BY POWELL (1983).


USAGE:

      QL0001(M,ME,MMAX,N,NMAX,MNN,C,D,A,B,XL,XU,X,U,IOUT,IFAIL,IPRINT,
             WAR,LWAR,IWAR,LIWAR)


   DEFINITION OF THE PARAMETERS:

   M :        TOTAL NUMBER OF CONSTRAINTS.
   ME :       NUMBER OF EQUALITY CONSTRAINTS.
   MMAX :     ROW DIMENSION OF A. MMAX MUST BE AT LEAST ONE AND GREATER
              THAN M.
   N :        NUMBER OF VARIABLES.
   NMAX :     ROW DIMENSION OF C. NMAX MUST BE GREATER OR EQUAL TO N.
   MNN :      MUST BE EQUAL TO M + N + N.
   C(NMAX,NMAX): OBJECTIVE FUNCTION MATRIX WHICH SHOULD BE SYMMETRIC AND
              POSITIVE DEFINITE. IF IWAR(1) = 0, C IS SUPPOSED TO BE THE
              CHOLESKEY-FACTOR OF ANOTHER MATRIX, I.E. C IS UPPER
              TRIANGULAR.
   D(NMAX) :  CONTAINS THE CONSTANT VECTOR OF THE OBJECTIVE FUNCTION.
   A(MMAX,NMAX): CONTAINS THE DATA MATRIX OF THE LINEAR CONSTRAINTS.
   B(MMAX) :  CONTAINS THE CONSTANT DATA OF THE LINEAR CONSTRAINTS.
   XL(N),XU(N): CONTAIN THE LOWER AND UPPER BOUNDS FOR THE VARIABLES.
   X(N) :     ON RETURN, X CONTAINS THE OPTIMAL SOLUTION VECTOR.
   U(MNN) :   ON RETURN, U CONTAINS THE LAGRANGE MULTIPLIERS. THE FIRST
              M POSITIONS ARE RESERVED FOR THE MULTIPLIERS OF THE M
              LINEAR CONSTRAINTS AND THE SUBSEQUENT ONES FOR THE
              MULTIPLIERS OF THE LOWER AND UPPER BOUNDS. ON SUCCESSFUL
              TERMINATION, ALL VALUES OF U WITH RESPECT TO INEQUALITIES
              AND BOUNDS SHOULD BE GREATER OR EQUAL TO ZERO.
   IOUT :     INTEGER INDICATING THE DESIRED OUTPUT UNIT NUMBER, I.E.
              ALL WRITE-STATEMENTS START WITH 'WRITE(IOUT,... '.
   IFAIL :    SHOWS THE TERMINATION REASON.
      IFAIL = 0 :   SUCCESSFUL RETURN.
      IFAIL = 1 :   TOO MANY ITERATIONS (MORE THAN 40*(N+M)).
      IFAIL = 2 :   ACCURACY INSUFFICIENT TO SATISFY CONVERGENCE
                    CRITERION.
      IFAIL = 5 :   LENGTH OF A WORKING ARRAY IS TOO SHORT.
      IFAIL > 10 :  THE CONSTRAINTS ARE INCONSISTENT.
   IPRINT :   OUTPUT CONTROL.
      IPRINT = 0 :  NO OUTPUT OF QL0001.
      IPRINT > 0 :  BRIEF OUTPUT IN ERROR CASES.
   WAR(LWAR) : REAL WORKING ARRAY. THE LENGTH LWAR SHOULD BE GRATER THAN
               3*NMAX*NMAX/2 + 10*NMAX + 2*MMAX.
   IWAR(LIWAR): INTEGER WORKING ARRAY. THE LENGTH LIWAR SHOULD BE AT
              LEAST N.
              IF IWAR(1)=0 INITIALLY, THEN THE CHOLESKY DECOMPOSITION
              WHICH IS REQUIRED BY THE DUAL ALGORITHM TO GET THE FIRST
              UNCONSTRAINED MINIMUM OF THE OBJECTIVE FUNCTION, IS
              PERFORMED INTERNALLY. OTHERWISE, I.E. IF IWAR(1)=1, THEN
              IT IS ASSUMED THAT THE USER PROVIDES THE INITIAL FAC-
              TORIZATION BY HIMSELF AND STORES IT IN THE UPPER TRIAN-
              GULAR PART OF THE ARRAY C.

   A NAMED COMMON-BLOCK  /CMACHE/EPS   MUST BE PROVIDED BY THE USER,
   WHERE EPS DEFINES A GUESS FOR THE UNDERLYING MACHINE PRECISION.


   AUTHOR:    K. SCHITTKOWSKI,
              MATHEMATISCHES INSTITUT,
              UNIVERSITAET BAYREUTH,
              8580 BAYREUTH,
              GERMANY, F.R.


   VERSION:   1.4  (MARCH, 1987)
*/
/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

 - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef int integer;
typedef char *address;
typedef short int shortint;
typedef float cfsqpreal;
typedef double doublereal;
typedef struct
{
    cfsqpreal r, i;
}
cfsqpcomplex;
typedef struct
{
    doublereal r, i;
}
doublecomplex;
typedef long int logical;
typedef short int shortlogical;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long flag;
typedef long ftnlen;
typedef long ftnint;
#endif

/*external read, write*/
typedef struct
{
    flag cierr;
    ftnint ciunit;
    flag ciend;
    char *cifmt;
    ftnint cirec;
}
cilist;

/*internal read, write*/
typedef struct
{
    flag icierr;
    char *iciunit;
    flag iciend;
    char *icifmt;
    ftnint icirlen;
    ftnint icirnum;
}
icilist;

/*open*/
typedef struct
{
    flag oerr;
    ftnint ounit;
    char *ofnm;
    ftnlen ofnmlen;
    char *osta;
    char *oacc;
    char *ofm;
    ftnint orl;
    char *oblnk;
}
olist;

/*close*/
typedef struct
{
    flag cerr;
    ftnint cunit;
    char *csta;
}
cllist;

/*rewind, backspace, endfile*/
typedef struct
{
    flag aerr;
    ftnint aunit;
}
alist;

/* inquire */
typedef struct
{
    flag inerr;
    ftnint inunit;
    char *infile;
    ftnlen infilen;
    ftnint *inex; /*parameters in standard's order*/
    ftnint *inopen;
    ftnint *innum;
    ftnint *innamed;
    char *inname;
    ftnlen innamlen;
    char *inacc;
    ftnlen inacclen;
    char *inseq;
    ftnlen inseqlen;
    char  *indir;
    ftnlen indirlen;
    char *infmt;
    ftnlen infmtlen;
    char *inform;
    ftnint informlen;
    char *inunf;
    ftnlen inunflen;
    ftnint *inrecl;
    ftnint *innrec;
    char *inblank;
    ftnlen inblanklen;
}
inlist;

#define VOID void

union Multitype { /* for multiple entry points */
    shortint h;
    integer i;
    cfsqpreal r;
    doublereal d;
    cfsqpcomplex cc;
    doublecomplex z;
};

typedef union Multitype Multitype;

typedef long Long;

struct Vardesc
{ /* for Namelist */
    char *name;
    char *addr;
    Long *dims;
    int  type;
};
typedef struct Vardesc Vardesc;

struct Namelist
{
    char *name;
    Vardesc **vars;
    int nvars;
};
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint(*J_fp)(...);
typedef integer(*I_fp)(...);
typedef cfsqpreal(*R_fp)(...);
typedef doublereal(*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID(*C_fp)(...);
typedef /* Double Complex */ VOID(*Z_fp)(...);
typedef logical(*L_fp)(...);
typedef shortlogical(*K_fp)(...);
typedef /* Character */ VOID(*H_fp)(...);
typedef /* Subroutine */ int(*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint(*J_fp)();
typedef integer(*I_fp)();
typedef cfsqpreal(*R_fp)();
typedef doublereal(*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID(*C_fp)();
typedef /* Double Complex */ VOID(*Z_fp)();
typedef logical(*L_fp)();
typedef shortlogical(*K_fp)();
typedef /* Character */ VOID(*H_fp)();
typedef /* Subroutine */ int(*S_fp)();
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f; /* complex function */
typedef VOID H_f; /* character function */
typedef VOID Z_f; /* double complex function */
typedef doublereal E_f; /* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

/* asolano: confuses Portland Group compiler, doesn't seem to affect anything */
#define Skip_f2c_Undefs

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif
#endif



/* Common Block Declarations */

struct Tcmache
{
    doublereal eps;
}
cmache_;

#define cmache_1 cmache_

/* umd */
/*
ql0002_ is declared here to provide ANSI C compliance.
(Thanks got to Martin Wauchope for providing this correction)
*/
#ifdef __STDC__

int ql0002_(integer *n, integer *m, integer *meq, integer *mmax,
            integer *mn, integer *mnn, integer *nmax,
            logical *lql,
            doublereal *a, doublereal *b, doublereal *grad,
            doublereal *g, doublereal *xl, doublereal *xu, doublereal *x,
            integer *nact, integer *iact, integer *maxit,
            doublereal *vsmall,
            integer *info,
            doublereal *diag, doublereal *w,
            integer *lw);
#else
int ql0002_();
#endif
/* umd */
/*
When the fortran code was f2c converted, the use of fortran COMMON
blocks was no longer available. Thus an additional variable, eps1,
was added to the parameter list to account for this.
*/
/* umd */
/*
Two alternative definitions are provided in order to give ANSI
compliance.
*/
#ifdef __STDC__
int ql0001_(int *m, int *me, int *mmax, int *n, int *nmax, int *mnn,
            double *c, double *d, double *a, double *b, double *xl,
            double *xu, double *x, double *u, int *iout, int *ifail,
            int *iprint, double *war, int *lwar, int *iwar, int *liwar,
            double *eps1)
#else
/* Subroutine */
int ql0001_(m, me, mmax, n, nmax, mnn, c, d, a, b, xl, xu, x,
            u, iout, ifail, iprint, war, lwar, iwar, liwar, eps1)
integer *m, *me, *mmax, *n, *nmax, *mnn;
doublereal *c, *d, *a, *b, *xl, *xu, *x, *u;
integer *iout, *ifail, *iprint;
doublereal *war;
integer *lwar, *iwar, *liwar;
doublereal *eps1;
#endif
{
    /* System generated locals */
    integer c_dim1, c_offset, a_dim1, a_offset, i__1;

    /* Builtin functions */
    /*    integer s_wsfe(), do_fio(), e_wsfe(); */

    /* Local variables */
    static doublereal diag;
    /* extern int ql0002_(); */
    static integer nact, info;
    static doublereal zero;
    static integer i, j, maxit;
    static doublereal qpeps;
    static integer in, mn, lw;
    static logical lql;
    static integer inw1, inw2;

    /*     INTRINSIC FUNCTIONS:  DSQRT */

    /* Parameter adjustments */
    --iwar;
    --war;
    --u;
    --x;
    --xu;
    --xl;
    --b;
    a_dim1 = *mmax;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --d;
    c_dim1 = *nmax;
    c_offset = c_dim1 + 1;
    c -= c_offset;

    /* Function Body */
    cmache_1.eps = *eps1;

    /*     CONSTANT DATA */

    /* ################################################################# */

    if (fabs(c[*nmax + *nmax * c_dim1]) == 0.e0)
    {
        c[*nmax + *nmax * c_dim1] = cmache_1.eps;
    }

    /* umd */
    /*  This prevents a subsequent more major modification of the Hessian */
    /*  matrix in the important case when a minmax problem (yielding a */
    /*  singular Hessian matrix) is being solved. */
    /*                                 ----UMCP, April 1991, Jian L. Zhou */
    /* ################################################################# */

    lql = FALSE_;
    if (iwar[1] == 1)
    {
        lql = TRUE_;
    }
    zero = 0.;
    maxit = (*m + *n) * 40;
    qpeps = cmache_1.eps;
    inw1 = 1;
    inw2 = inw1 + *mmax;

    /*     PREPARE PROBLEM DATA FOR EXECUTION */

    if (*m <= 0)
    {
        goto L20;
    }
    in = inw1;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j)
    {
        war[in] = -b[j];
        /* L10: */
        ++in;
    }
L20:
    lw = *nmax * 3 * *nmax / 2 + *nmax * 10 + *m;
    if (inw2 + lw > *lwar)
    {
        goto L80;
    }
    if (*liwar < *n)
    {
        goto L81;
    }
    if (*mnn < *m + *n + *n)
    {
        goto L82;
    }
    mn = *m + *n;

    /*     CALL OF QL0002 */

    ql0002_(n, m, me, mmax, &mn, mnn, nmax, &lql, &a[a_offset], &war[inw1], &
            d[1], &c[c_offset], &xl[1], &xu[1], &x[1], &nact, &iwar[1], &
            maxit, &qpeps, &info, &diag, &war[inw2], &lw);

    /*     TEST OF MATRIX CORRECTIONS */

    *ifail = 0;
    if (info == 1)
    {
        goto L40;
    }
    if (info == 2)
    {
        goto L90;
    }
    if (info < 0)
    {
        goto L70;
    }

    /*     REORDER MULTIPLIER */

    i__1 = *mnn;
    for (j = 1; j <= i__1; ++j)
    {
        /* L50: */
        u[j] = zero;
    }
    in = inw2 - 1;
    if (nact == 0)
    {
        goto L30;
    }
    i__1 = nact;
    for (i = 1; i <= i__1; ++i)
    {
        j = iwar[i];
        u[j] = war[in + i];
        /* L60: */
    }
L30:
    return 0;

    /*     ERROR MESSAGES */

L70:
    *ifail = -info + 10;
    /*
        if (*iprint > 0 && nact > 0) {
     io___18.ciunit = *iout;
     s_wsfe(&io___18);
     i__1 = -info;
     do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
     i__2 = nact;
     for (i = 1; i <= i__2; ++i) {
         do_fio(&c__1, (char *)&iwar[i], (ftnlen)sizeof(integer));
     }
     e_wsfe();
        }
    */
    return 0;
L80:
    *ifail = 5;
    /*
        if (*iprint > 0) {
     io___19.ciunit = *iout;
     s_wsfe(&io___19);
     e_wsfe();
        }
    */
    return 0;
L81:
    *ifail = 5;
    /*
        if (*iprint > 0) {
     io___20.ciunit = *iout;
     s_wsfe(&io___20);
     e_wsfe();
        }
    */
    return 0;
L82:
    *ifail = 5;
    /*
        if (*iprint > 0) {
     io___21.ciunit = *iout;
     s_wsfe(&io___21);
     e_wsfe();
        }
    */
    return 0;
L40:
    *ifail = 1;
    /*
        if (*iprint > 0) {
     io___22.ciunit = *iout;
     s_wsfe(&io___22);
     do_fio(&c__1, (char *)&maxit, (ftnlen)sizeof(integer));
     e_wsfe();
        }
    */
    return 0;
L90:
    *ifail = 2;
    /*
        if (*iprint > 0) {
     io___23.ciunit = *iout;
     s_wsfe(&io___23);
     e_wsfe();
        }
    */
    return 0;

    /*     FORMAT-INSTRUCTIONS */

} /* ql0001_ */


/* umd
Two alternative definitions are provided in order to give ANSI
compliance.
(Thanks got to Martin Wauchope for providing this correction)
*/
#ifdef __STDC__
int ql0002_(integer *n, integer *m, integer *meq, integer *mmax,
            integer *mn, integer *mnn, integer *nmax,
            logical *lql,
            doublereal *a, doublereal *b, doublereal *grad,
            doublereal *g, doublereal *xl, doublereal *xu, doublereal *x,
            integer *nact, integer *iact, integer *maxit,
            doublereal *vsmall,
            integer *info,
            doublereal *diag, doublereal *w,
            integer *lw)
#else
/* Subroutine */ int ql0002_(n, m, meq, mmax, mn, mnn, nmax, lql, a, b, grad,
                             g, xl, xu, x, nact, iact, maxit, vsmall, info, diag, w, lw)
integer *n, *m, *meq, *mmax, *mn, *mnn, *nmax;
logical *lql;
doublereal *a, *b, *grad, *g, *xl, *xu, *x;
integer *nact, *iact, *maxit;
doublereal *vsmall;
integer *info;
doublereal *diag, *w;
integer *lw;
#endif
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    /* umd */
    /* double sqrt();    */

    /* Local variables */
    static doublereal onha, xmag, suma, sumb, sumc, temp, step, zero;
    static integer iwwn;
    static doublereal sumx, sumy;
    static integer i, j, k;
    static doublereal fdiff;
    static integer iflag, jflag, kflag, lflag;
    static doublereal diagr;
    static integer ifinc, kfinc, jfinc, mflag, nflag;
    static doublereal vfact, tempa;
    static integer iterc, itref;
    static doublereal cvmax, ratio, xmagr;
    static integer kdrop;
    static logical lower;
    static integer knext, k1;
    static doublereal ga, gb;
    static integer ia, id;
    static doublereal fdiffa;
    static integer ii, il, kk, jl, ir, nm, is, iu, iw, ju, ix, iz, nu, iy;

    static doublereal parinc, parnew;
    static integer ira, irb, iwa;
    static doublereal one;
    static integer iwd, iza;
    static doublereal res;
    static integer iwr, iws;
    static doublereal sum;
    static integer iww, iwx, iwy;
    static doublereal two;
    static integer iwz;


    /*       WHETHER THE CONSTRAINT IS ACTIVE. */


    /*   AUTHOR:    K. SCHITTKOWSKI, */
    /*              MATHEMATISCHES INSTITUT, */
    /*              UNIVERSITAET BAYREUTH, */
    /*              8580 BAYREUTH, */
    /*              GERMANY, F.R. */

    /*   AUTHOR OF ORIGINAL VERSION: */
    /*              M.J.D. POWELL, DAMTP, */
    /*              UNIVERSITY OF CAMBRIDGE, SILVER STREET */
    /*              CAMBRIDGE, */
    /*              ENGLAND */


    /*   REFERENCE: M.J.D. POWELL: ZQPCVX, A FORTRAN SUBROUTINE FOR CONVEX */
    /*              PROGRAMMING, REPORT DAMTP/1983/NA17, UNIVERSITY OF */
    /*              CAMBRIDGE, ENGLAND, 1983. */


    /*   VERSION :  2.0 (MARCH, 1987) */


    /************************************************************************
    ***/


    /*   INTRINSIC FUNCTIONS:   DMAX1,DSQRT,DABS,DMIN1 */


    /*   INITIAL ADDRESSES */

    /* Parameter adjustments */
    --w;
    --iact;
    --x;
    --xu;
    --xl;
    g_dim1 = *nmax;
    g_offset = g_dim1 + 1;
    g -= g_offset;
    --grad;
    --b;
    a_dim1 = *mmax;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    iwz = *nmax;
    iwr = iwz + *nmax * *nmax;
    iww = iwr + *nmax * (*nmax + 3) / 2;
    iwd = iww + *nmax;
    iwx = iwd + *nmax;
    iwa = iwx + *nmax;

    /*     SET SOME CONSTANTS. */

    zero = 0.;
    one = 1.;
    two = 2.;
    onha = 1.5;
    vfact = 1.;

    /*     SET SOME PARAMETERS. */
    /*     NUMBER LESS THAN VSMALL ARE ASSUMED TO BE NEGLIGIBLE. */
    /*     THE MULTIPLE OF I THAT IS ADDED TO G IS AT MOST DIAGR TIMES */
    /*       THE LEAST MULTIPLE OF I THAT GIVES POSITIVE DEFINITENESS. */
    /*     X IS RE-INITIALISED IF ITS MAGNITUDE IS REDUCED BY THE */
    /*       FACTOR XMAGR. */
    /*     A CHECK IS MADE FOR AN INCREASE IN F EVERY IFINC ITERATIONS, */
    /*       AFTER KFINC ITERATIONS ARE COMPLETED. */

    diagr = two;
    xmagr = .01;
    ifinc = 3;
    kfinc = max(10, *n);

    /*     FIND THE RECIPROCALS OF THE LENGTHS OF THE CONSTRAINT NORMALS. */
    /*     RETURN IF A CONSTRAINT IS INFEASIBLE DUE TO A ZERO NORMAL. */

    *nact = 0;
    if (*m <= 0)
    {
        goto L45;
    }
    i__1 = *m;
    for (k = 1; k <= i__1; ++k)
    {
        sum = zero;
        i__2 = *n;
        for (i = 1; i <= i__2; ++i)
        {
            /* L10: */
            /* Computing 2nd power */
            d__1 = a[k + i * a_dim1];
            sum += d__1 * d__1;
        }
        if (sum > zero)
        {
            goto L20;
        }
        if (b[k] == zero)
        {
            goto L30;
        }
        *info = -k;
        if (k <= *meq)
        {
            goto L730;
        }
        if (b[k] <= 0.)
        {
            goto L30;
        }
        else
        {
            goto L730;
        }
L20:
        sum = one / sqrt(sum);
L30:
        ia = iwa + k;
        /* L40: */
        w[ia] = sum;
    }
L45:
    i__1 = *n;
    for (k = 1; k <= i__1; ++k)
    {
        ia = iwa + *m + k;
        /* L50: */
        w[ia] = one;
    }

    /*     IF NECESSARY INCREASE THE DIAGONAL ELEMENTS OF G. */

    if (!(*lql))
    {
        goto L165;
    }
    *diag = zero;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        id = iwd + i;
        w[id] = g[i + i * g_dim1];
        /* Computing MAX */
        d__1 = *diag, d__2 = *vsmall - w[id];
        *diag = max(d__1, d__2);
        if (i == *n)
        {
            goto L60;
        }
        ii = i + 1;
        i__2 = *n;
        for (j = ii; j <= i__2; ++j)
        {
            /* Computing MIN */
            d__1 = w[id], d__2 = g[j + j * g_dim1];
            ga = -min(d__1, d__2);
            gb = (d__1 = w[id] - g[j + j * g_dim1], abs(d__1)) + (d__2 = g[i
                    + j * g_dim1], abs(d__2));
            if (gb > zero)
            {
                /* Computing 2nd power */
                d__1 = g[i + j * g_dim1];
                ga += d__1 * d__1 / gb;
            }
            /* L55: */
            *diag = max(*diag, ga);
        }
L60:
        ;
    }
    if (*diag <= zero)
    {
        goto L90;
    }
L70:
    *diag = diagr * *diag;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        id = iwd + i;
        /* L80: */
        g[i + i * g_dim1] = *diag + w[id];
    }

    /*     FORM THE CHOLESKY FACTORISATION OF G. THE TRANSPOSE */
    /*     OF THE FACTOR WILL BE PLACED IN THE R-PARTITION OF W. */

L90:
    ir = iwr;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j)
    {
        ira = iwr;
        irb = ir + 1;
        i__2 = j;
        for (i = 1; i <= i__2; ++i)
        {
            temp = g[i + j * g_dim1];
            if (i == 1)
            {
                goto L110;
            }
            i__3 = ir;
            for (k = irb; k <= i__3; ++k)
            {
                ++ira;
                /* L100: */
                temp -= w[k] * w[ira];
            }
L110:
            ++ir;
            ++ira;
            if (i < j)
            {
                w[ir] = temp / w[ira];
            }
            /* L120: */
        }
        if (temp < *vsmall)
        {
            goto L140;
        }
        /* L130: */
        w[ir] = sqrt(temp);
    }
    goto L170;

    /*     INCREASE FURTHER THE DIAGONAL ELEMENT OF G. */

L140:
    w[j] = one;
    sumx = one;
    k = j;
L150:
    sum = zero;
    ira = ir - 1;
    i__1 = j;
    for (i = k; i <= i__1; ++i)
    {
        sum -= w[ira] * w[i];
        /* L160: */
        ira += i;
    }
    ir -= k;
    --k;
    w[k] = sum / w[ir];
    /* Computing 2nd power */
    d__1 = w[k];
    sumx += d__1 * d__1;
    if (k >= 2)
    {
        goto L150;
    }
    *diag = *diag + *vsmall - temp / sumx;
    goto L70;

    /*     STORE THE CHOLESKY FACTORISATION IN THE R-PARTITION */
    /*     OF W. */

L165:
    ir = iwr;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i)
    {
        i__2 = i;
        for (j = 1; j <= i__2; ++j)
        {
            ++ir;
            /* L166: */
            w[ir] = g[j + i * g_dim1];
        }
    }

    /*     SET Z THE INVERSE OF THE MATRIX IN R. */

L170:
    nm = *n - 1;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        iz = iwz + i;
        if (i == 1)
        {
            goto L190;
        }
        i__1 = i;
        for (j = 2; j <= i__1; ++j)
        {
            w[iz] = zero;
            /* L180: */
            iz += *n;
        }
L190:
        ir = iwr + (i + i * i) / 2;
        w[iz] = one / w[ir];
        if (i == *n)
        {
            goto L220;
        }
        iza = iz;
        i__1 = nm;
        for (j = i; j <= i__1; ++j)
        {
            ir += i;
            sum = zero;
            i__3 = iz;
            i__4 = *n;
            for (k = iza; i__4 < 0 ? k >= i__3 : k <= i__3; k += i__4)
            {
                sum += w[k] * w[ir];
                /* L200: */
                ++ir;
            }
            iz += *n;
            /* L210: */
            w[iz] = -sum / w[ir];
        }
L220:
        ;
    }

    /*     SET THE INITIAL VALUES OF SOME VARIABLES. */
    /*     ITERC COUNTS THE NUMBER OF ITERATIONS. */
    /*     ITREF IS SET TO ONE WHEN ITERATIVE REFINEMENT IS REQUIRED. */
    /*     JFINC INDICATES WHEN TO TEST FOR AN INCREASE IN F. */

    iterc = 1;
    itref = 0;
    jfinc = -kfinc;

    /*     SET X TO ZERO AND SET THE CORRESPONDING RESIDUALS OF THE */
    /*     KUHN-TUCKER CONDITIONS. */

L230:
    iflag = 1;
    iws = iww - *n;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        x[i] = zero;
        iw = iww + i;
        w[iw] = grad[i];
        if (i > *nact)
        {
            goto L240;
        }
        w[i] = zero;
        is = iws + i;
        k = iact[i];
        if (k <= *m)
        {
            goto L235;
        }
        if (k > *mn)
        {
            goto L234;
        }
        k1 = k - *m;
        w[is] = xl[k1];
        goto L240;
L234:
        k1 = k - *mn;
        w[is] = -xu[k1];
        goto L240;
L235:
        w[is] = b[k];
L240:
        ;
    }
    xmag = zero;
    vfact = 1.;
    if (*nact <= 0)
    {
        goto L340;
    }
    else
    {
        goto L280;
    }

    /*     SET THE RESIDUALS OF THE KUHN-TUCKER CONDITIONS FOR GENERAL X. */

L250:
    iflag = 2;
    iws = iww - *n;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        iw = iww + i;
        w[iw] = grad[i];
        if (*lql)
        {
            goto L259;
        }
        id = iwd + i;
        w[id] = zero;
        i__1 = *n;
        for (j = i; j <= i__1; ++j)
        {
            /* L251: */
            w[id] += g[i + j * g_dim1] * x[j];
        }
        i__1 = i;
        for (j = 1; j <= i__1; ++j)
        {
            id = iwd + j;
            /* L252: */
            w[iw] += g[j + i * g_dim1] * w[id];
        }
        goto L260;
L259:
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
            /* L261: */
            w[iw] += g[i + j * g_dim1] * x[j];
        }
L260:
        ;
    }
    if (*nact == 0)
    {
        goto L340;
    }
    i__2 = *nact;
    for (k = 1; k <= i__2; ++k)
    {
        kk = iact[k];
        is = iws + k;
        if (kk > *m)
        {
            goto L265;
        }
        w[is] = b[kk];
        i__1 = *n;
        for (i = 1; i <= i__1; ++i)
        {
            iw = iww + i;
            w[iw] -= w[k] * a[kk + i * a_dim1];
            /* L264: */
            w[is] -= x[i] * a[kk + i * a_dim1];
        }
        goto L270;
L265:
        if (kk > *mn)
        {
            goto L266;
        }
        k1 = kk - *m;
        iw = iww + k1;
        w[iw] -= w[k];
        w[is] = xl[k1] - x[k1];
        goto L270;
L266:
        k1 = kk - *mn;
        iw = iww + k1;
        w[iw] += w[k];
        w[is] = -xu[k1] + x[k1];
L270:
        ;
    }

    /*     PRE-MULTIPLY THE VECTOR IN THE S-PARTITION OF W BY THE */
    /*     INVERS OF R TRANSPOSE. */

L280:
    ir = iwr;
    il = iws + 1;
    iu = iws + *nact;
    i__2 = iu;
    for (i = il; i <= i__2; ++i)
    {
        sum = zero;
        if (i == il)
        {
            goto L300;
        }
        ju = i - 1;
        i__1 = ju;
        for (j = il; j <= i__1; ++j)
        {
            ++ir;
            /* L290: */
            sum += w[ir] * w[j];
        }
L300:
        ++ir;
        /* L310: */
        w[i] = (w[i] - sum) / w[ir];
    }

    /*     SHIFT X TO SATISFY THE ACTIVE CONSTRAINTS AND MAKE THE */
    /*     CORRESPONDING CHANGE TO THE GRADIENT RESIDUALS. */

    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        iz = iwz + i;
        sum = zero;
        i__1 = iu;
        for (j = il; j <= i__1; ++j)
        {
            sum += w[j] * w[iz];
            /* L320: */
            iz += *n;
        }
        x[i] += sum;
        if (*lql)
        {
            goto L329;
        }
        id = iwd + i;
        w[id] = zero;
        i__1 = *n;
        for (j = i; j <= i__1; ++j)
        {
            /* L321: */
            w[id] += g[i + j * g_dim1] * sum;
        }
        iw = iww + i;
        i__1 = i;
        for (j = 1; j <= i__1; ++j)
        {
            id = iwd + j;
            /* L322: */
            w[iw] += g[j + i * g_dim1] * w[id];
        }
        goto L330;
L329:
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
            iw = iww + j;
            /* L331: */
            w[iw] += sum * g[i + j * g_dim1];
        }
L330:
        ;
    }

    /*     FORM THE SCALAR PRODUCT OF THE CURRENT GRADIENT RESIDUALS */
    /*     WITH EACH COLUMN OF Z. */

L340:
    kflag = 1;
    goto L930;
L350:
    if (*nact == *n)
    {
        goto L380;
    }

    /*     SHIFT X SO THAT IT SATISFIES THE REMAINING KUHN-TUCKER */
    /*     CONDITIONS. */

    il = iws + *nact + 1;
    iza = iwz + *nact * *n;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        sum = zero;
        iz = iza + i;
        i__1 = iww;
        for (j = il; j <= i__1; ++j)
        {
            sum += w[iz] * w[j];
            /* L360: */
            iz += *n;
        }
        /* L370: */
        x[i] -= sum;
    }
    *info = 0;
    if (*nact == 0)
    {
        goto L410;
    }

    /*     UPDATE THE LAGRANGE MULTIPLIERS. */

L380:
    lflag = 3;
    goto L740;
L390:
    i__2 = *nact;
    for (k = 1; k <= i__2; ++k)
    {
        iw = iww + k;
        /* L400: */
        w[k] += w[iw];
    }

    /*     REVISE THE VALUES OF XMAG. */
    /*     BRANCH IF ITERATIVE REFINEMENT IS REQUIRED. */

L410:
    jflag = 1;
    goto L910;
L420:
    if (iflag == itref)
    {
        goto L250;
    }

    /*     DELETE A CONSTRAINT IF A LAGRANGE MULTIPLIER OF AN */
    /*     INEQUALITY CONSTRAINT IS NEGATIVE. */

    kdrop = 0;
    goto L440;
L430:
    ++kdrop;
    if (w[kdrop] >= zero)
    {
        goto L440;
    }
    if (iact[kdrop] <= *meq)
    {
        goto L440;
    }
    nu = *nact;
    mflag = 1;
    goto L800;
L440:
    if (kdrop < *nact)
    {
        goto L430;
    }

    /*     SEEK THE GREATEAST NORMALISED CONSTRAINT VIOLATION, DISREGARDING */

    /*     ANY THAT MAY BE DUE TO COMPUTER ROUNDING ERRORS. */

L450:
    cvmax = zero;
    if (*m <= 0)
    {
        goto L481;
    }
    i__2 = *m;
    for (k = 1; k <= i__2; ++k)
    {
        ia = iwa + k;
        if (w[ia] <= zero)
        {
            goto L480;
        }
        sum = -b[k];
        i__1 = *n;
        for (i = 1; i <= i__1; ++i)
        {
            /* L460: */
            sum += x[i] * a[k + i * a_dim1];
        }
        sumx = -sum * w[ia];
        if (k <= *meq)
        {
            sumx = abs(sumx);
        }
        if (sumx <= cvmax)
        {
            goto L480;
        }
        temp = (d__1 = b[k], abs(d__1));
        i__1 = *n;
        for (i = 1; i <= i__1; ++i)
        {
            /* L470: */
            temp += (d__1 = x[i] * a[k + i * a_dim1], abs(d__1));
        }
        tempa = temp + abs(sum);
        if (tempa <= temp)
        {
            goto L480;
        }
        temp += onha * abs(sum);
        if (temp <= tempa)
        {
            goto L480;
        }
        cvmax = sumx;
        res = sum;
        knext = k;
L480:
        ;
    }
L481:
    i__2 = *n;
    for (k = 1; k <= i__2; ++k)
    {
        lower = TRUE_;
        ia = iwa + *m + k;
        if (w[ia] <= zero)
        {
            goto L485;
        }
        sum = xl[k] - x[k];
        if (sum < 0.)
        {
            goto L482;
        }
        else if (sum == 0)
        {
            goto L485;
        }
        else
        {
            goto L483;
        }
L482:
        sum = x[k] - xu[k];
        lower = FALSE_;
L483:
        if (sum <= cvmax)
        {
            goto L485;
        }
        cvmax = sum;
        res = -sum;
        knext = k + *m;
        if (lower)
        {
            goto L485;
        }
        knext = k + *mn;
L485:
        ;
    }

    /*     TEST FOR CONVERGENCE */

    *info = 0;
    if (cvmax <= *vsmall)
    {
        goto L700;
    }

    /*     RETURN IF, DUE TO ROUNDING ERRORS, THE ACTUAL CHANGE IN */
    /*     X MAY NOT INCREASE THE OBJECTIVE FUNCTION */

    ++jfinc;
    if (jfinc == 0)
    {
        goto L510;
    }
    if (jfinc != ifinc)
    {
        goto L530;
    }
    fdiff = zero;
    fdiffa = zero;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        sum = two * grad[i];
        sumx = abs(sum);
        if (*lql)
        {
            goto L489;
        }
        id = iwd + i;
        w[id] = zero;
        i__1 = *n;
        for (j = i; j <= i__1; ++j)
        {
            ix = iwx + j;
            /* L486: */
            w[id] += g[i + j * g_dim1] * (w[ix] + x[j]);
        }
        i__1 = i;
        for (j = 1; j <= i__1; ++j)
        {
            id = iwd + j;
            temp = g[j + i * g_dim1] * w[id];
            sum += temp;
            /* L487: */
            sumx += abs(temp);
        }
        goto L495;
L489:
        i__1 = *n;
        for (j = 1; j <= i__1; ++j)
        {
            ix = iwx + j;
            temp = g[i + j * g_dim1] * (w[ix] + x[j]);
            sum += temp;
            /* L490: */
            sumx += abs(temp);
        }
L495:
        ix = iwx + i;
        fdiff += sum * (x[i] - w[ix]);
        /* L500: */
        fdiffa += sumx * (d__1 = x[i] - w[ix], abs(d__1));
    }
    *info = 2;
    sum = fdiffa + fdiff;
    if (sum <= fdiffa)
    {
        goto L700;
    }
    temp = fdiffa + onha * fdiff;
    if (temp <= sum)
    {
        goto L700;
    }
    jfinc = 0;
    *info = 0;
L510:
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        ix = iwx + i;
        /* L520: */
        w[ix] = x[i];
    }

    /*     FORM THE SCALAR PRODUCT OF THE NEW CONSTRAINT NORMAL WITH EACH */
    /*     COLUMN OF Z. PARNEW WILL BECOME THE LAGRANGE MULTIPLIER OF */
    /*     THE NEW CONSTRAINT. */

L530:
    ++iterc;
    if (iterc <= *maxit)
    {
        goto L531;
    }
    *info = 1;
    goto L710;
L531:
    iws = iwr + (*nact + *nact * *nact) / 2;
    if (knext > *m)
    {
        goto L541;
    }
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        iw = iww + i;
        /* L540: */
        w[iw] = a[knext + i * a_dim1];
    }
    goto L549;
L541:
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        iw = iww + i;
        /* L542: */
        w[iw] = zero;
    }
    k1 = knext - *m;
    if (k1 > *n)
    {
        goto L545;
    }
    iw = iww + k1;
    w[iw] = one;
    iz = iwz + k1;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        is = iws + i;
        w[is] = w[iz];
        /* L543: */
        iz += *n;
    }
    goto L550;
L545:
    k1 = knext - *mn;
    iw = iww + k1;
    w[iw] = -one;
    iz = iwz + k1;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        is = iws + i;
        w[is] = -w[iz];
        /* L546: */
        iz += *n;
    }
    goto L550;
L549:
    kflag = 2;
    goto L930;
L550:
    parnew = zero;

    /*     APPLY GIVENS ROTATIONS TO MAKE THE LAST (N-NACT-2) SCALAR */
    /*     PRODUCTS EQUAL TO ZERO. */

    if (*nact == *n)
    {
        goto L570;
    }
    nu = *n;
    nflag = 1;
    goto L860;

    /*     BRANCH IF THERE IS NO NEED TO DELETE A CONSTRAINT. */

L560:
    is = iws + *nact;
    if (*nact == 0)
    {
        goto L640;
    }
    suma = zero;
    sumb = zero;
    sumc = zero;
    iz = iwz + *nact * *n;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        ++iz;
        iw = iww + i;
        suma += w[iw] * w[iz];
        sumb += (d__1 = w[iw] * w[iz], abs(d__1));
        /* L563: */
        /* Computing 2nd power */
        d__1 = w[iz];
        sumc += d__1 * d__1;
    }
    temp = sumb + abs(suma) * .1;
    tempa = sumb + abs(suma) * .2;
    if (temp <= sumb)
    {
        goto L570;
    }
    if (tempa <= temp)
    {
        goto L570;
    }
    if (sumb > *vsmall)
    {
        goto L5;
    }
    goto L570;
L5:
    sumc = sqrt(sumc);
    ia = iwa + knext;
    if (knext <= *m)
    {
        sumc /= w[ia];
    }
    temp = sumc + abs(suma) * .1;
    tempa = sumc + abs(suma) * .2;
    if (temp <= sumc)
    {
        goto L567;
    }
    if (tempa <= temp)
    {
        goto L567;
    }
    goto L640;

    /*     CALCULATE THE MULTIPLIERS FOR THE NEW CONSTRAINT NORMAL */
    /*     EXPRESSED IN TERMS OF THE ACTIVE CONSTRAINT NORMALS. */
    /*     THEN WORK OUT WHICH CONTRAINT TO DROP. */

L567:
    lflag = 4;
    goto L740;
L570:
    lflag = 1;
    goto L740;

    /*     COMPLETE THE TEST FOR LINEARLY DEPENDENT CONSTRAINTS. */

L571:
    if (knext > *m)
    {
        goto L574;
    }
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        suma = a[knext + i * a_dim1];
        sumb = abs(suma);
        if (*nact == 0)
        {
            goto L581;
        }
        i__1 = *nact;
        for (k = 1; k <= i__1; ++k)
        {
            kk = iact[k];
            if (kk <= *m)
            {
                goto L568;
            }
            kk -= *m;
            temp = zero;
            if (kk == i)
            {
                temp = w[iww + kk];
            }
            kk -= *n;
            if (kk == i)
            {
                temp = -w[iww + kk];
            }
            goto L569;
L568:
            iw = iww + k;
            temp = w[iw] * a[kk + i * a_dim1];
L569:
            suma -= temp;
            /* L572: */
            sumb += abs(temp);
        }
L581:
        if (suma <= *vsmall)
        {
            goto L573;
        }
        temp = sumb + abs(suma) * .1;
        tempa = sumb + abs(suma) * .2;
        if (temp <= sumb)
        {
            goto L573;
        }
        if (tempa <= temp)
        {
            goto L573;
        }
        goto L630;
L573:
        ;
    }
    lflag = 1;
    goto L775;
L574:
    k1 = knext - *m;
    if (k1 > *n)
    {
        k1 -= *n;
    }
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        suma = zero;
        if (i != k1)
        {
            goto L575;
        }
        suma = one;
        if (knext > *mn)
        {
            suma = -one;
        }
L575:
        sumb = abs(suma);
        if (*nact == 0)
        {
            goto L582;
        }
        i__1 = *nact;
        for (k = 1; k <= i__1; ++k)
        {
            kk = iact[k];
            if (kk <= *m)
            {
                goto L579;
            }
            kk -= *m;
            temp = zero;
            if (kk == i)
            {
                temp = w[iww + kk];
            }
            kk -= *n;
            if (kk == i)
            {
                temp = -w[iww + kk];
            }
            goto L576;
L579:
            iw = iww + k;
            temp = w[iw] * a[kk + i * a_dim1];
L576:
            suma -= temp;
            /* L577: */
            sumb += abs(temp);
        }
L582:
        temp = sumb + abs(suma) * .1;
        tempa = sumb + abs(suma) * .2;
        if (temp <= sumb)
        {
            goto L578;
        }
        if (tempa <= temp)
        {
            goto L578;
        }
        goto L630;
L578:
        ;
    }
    lflag = 1;
    goto L775;

    /*     BRANCH IF THE CONTRAINTS ARE INCONSISTENT. */

L580:
    *info = -knext;
    if (kdrop == 0)
    {
        goto L700;
    }
    parinc = ratio;
    parnew = parinc;

    /*     REVISE THE LAGRANGE MULTIPLIERS OF THE ACTIVE CONSTRAINTS. */

L590:
    if (*nact == 0)
    {
        goto L601;
    }
    i__2 = *nact;
    for (k = 1; k <= i__2; ++k)
    {
        iw = iww + k;
        w[k] -= parinc * w[iw];
        if (iact[k] > *meq)
        {
            /* Computing MAX */
            d__1 = zero, d__2 = w[k];
            w[k] = max(d__1, d__2);
        }
        /* L600: */
    }
L601:
    if (kdrop == 0)
    {
        goto L680;
    }

    /*     DELETE THE CONSTRAINT TO BE DROPPED. */
    /*     SHIFT THE VECTOR OF SCALAR PRODUCTS. */
    /*     THEN, IF APPROPRIATE, MAKE ONE MORE SCALAR PRODUCT ZERO. */

    nu = *nact + 1;
    mflag = 2;
    goto L800;
L610:
    iws = iws - *nact - 1;
    nu = min(*n, nu);
    i__2 = nu;
    for (i = 1; i <= i__2; ++i)
    {
        is = iws + i;
        j = is + *nact;
        /* L620: */
        w[is] = w[j + 1];
    }
    nflag = 2;
    goto L860;

    /*     CALCULATE THE STEP TO THE VIOLATED CONSTRAINT. */

L630:
    is = iws + *nact;
L640:
    sumy = w[is + 1];
    step = -res / sumy;
    parinc = step / sumy;
    if (*nact == 0)
    {
        goto L660;
    }

    /*     CALCULATE THE CHANGES TO THE LAGRANGE MULTIPLIERS, AND REDUCE */
    /*     THE STEP ALONG THE NEW SEARCH DIRECTION IF NECESSARY. */

    lflag = 2;
    goto L740;
L650:
    if (kdrop == 0)
    {
        goto L660;
    }
    temp = one - ratio / parinc;
    if (temp <= zero)
    {
        kdrop = 0;
    }
    if (kdrop == 0)
    {
        goto L660;
    }
    step = ratio * sumy;
    parinc = ratio;
    res = temp * res;

    /*     UPDATE X AND THE LAGRANGE MULTIPIERS. */
    /*     DROP A CONSTRAINT IF THE FULL STEP IS NOT TAKEN. */

L660:
    iwy = iwz + *nact * *n;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        iy = iwy + i;
        /* L670: */
        x[i] += step * w[iy];
    }
    parnew += parinc;
    if (*nact >= 1)
    {
        goto L590;
    }

    /*     ADD THE NEW CONSTRAINT TO THE ACTIVE SET. */

L680:
    ++(*nact);
    w[*nact] = parnew;
    iact[*nact] = knext;
    ia = iwa + knext;
    if (knext > *mn)
    {
        ia -= *n;
    }
    w[ia] = -w[ia];

    /*     ESTIMATE THE MAGNITUDE OF X. THEN BEGIN A NEW ITERATION, */
    /*     RE-INITILISING X IF THIS MAGNITUDE IS SMALL. */

    jflag = 2;
    goto L910;
L690:
    if (sum < xmagr * xmag)
    {
        goto L230;
    }
    if (itref <= 0)
    {
        goto L450;
    }
    else
    {
        goto L250;
    }

    /*     INITIATE ITERATIVE REFINEMENT IF IT HAS NOT YET BEEN USED, */
    /*     OR RETURN AFTER RESTORING THE DIAGONAL ELEMENTS OF G. */

L700:
    if (iterc == 0)
    {
        goto L710;
    }
    ++itref;
    jfinc = -1;
    if (itref == 1)
    {
        goto L250;
    }
L710:
    if (!(*lql))
    {
        return 0;
    }
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        id = iwd + i;
        /* L720: */
        g[i + i * g_dim1] = w[id];
    }
L730:
    return 0;


    /*     THE REMAINIG INSTRUCTIONS ARE USED AS SUBROUTINES. */


    /* ******************************************************************** */



    /*     CALCULATE THE LAGRANGE MULTIPLIERS BY PRE-MULTIPLYING THE */
    /*     VECTOR IN THE S-PARTITION OF W BY THE INVERSE OF R. */

L740:
    ir = iwr + (*nact + *nact * *nact) / 2;
    i = *nact;
    sum = zero;
    goto L770;
L750:
    ira = ir - 1;
    sum = zero;
    if (*nact == 0)
    {
        goto L761;
    }
    i__2 = *nact;
    for (j = i; j <= i__2; ++j)
    {
        iw = iww + j;
        sum += w[ira] * w[iw];
        /* L760: */
        ira += j;
    }
L761:
    ir -= i;
    --i;
L770:
    iw = iww + i;
    is = iws + i;
    w[iw] = (w[is] - sum) / w[ir];
    if (i > 1)
    {
        goto L750;
    }
    if (lflag == 3)
    {
        goto L390;
    }
    if (lflag == 4)
    {
        goto L571;
    }

    /*     CALCULATE THE NEXT CONSTRAINT TO DROP. */

L775:
    kdrop = 0;
    if (*nact == 0)
    {
        goto L791;
    }
    i__2 = *nact;
    for (k = 1; k <= i__2; ++k)
    {
        if (iact[k] <= *meq)
        {
            goto L790;
        }
        iw = iww + k;
        if (res * w[iw] >= zero)
        {
            goto L790;
        }
        temp = w[k] / w[iw];
        if (kdrop == 0)
        {
            goto L780;
        }
        if (abs(temp) >= abs(ratio))
        {
            goto L790;
        }
L780:
        kdrop = k;
        ratio = temp;
L790:
        ;
    }
L791:
    switch ((int)lflag)
    {
    case 1:
        goto L580;
    case 2:
        goto L650;
    }


    /* ******************************************************************** */



    /*     DROP THE CONSTRAINT IN POSITION KDROP IN THE ACTIVE SET. */

L800:
    ia = iwa + iact[kdrop];
    if (iact[kdrop] > *mn)
    {
        ia -= *n;
    }
    w[ia] = -w[ia];
    if (kdrop == *nact)
    {
        goto L850;
    }

    /*     SET SOME INDICES AND CALCULATE THE ELEMENTS OF THE NEXT */
    /*     GIVENS ROTATION. */

    iz = iwz + kdrop * *n;
    ir = iwr + (kdrop + kdrop * kdrop) / 2;
L810:
    ira = ir;
    ir = ir + kdrop + 1;
    /* Computing MAX */
    d__3 = (d__1 = w[ir - 1], abs(d__1)), d__4 = (d__2 = w[ir], abs(d__2));
    temp = max(d__3, d__4);
    /* Computing 2nd power */
    d__1 = w[ir - 1] / temp;
    /* Computing 2nd power */
    d__2 = w[ir] / temp;
    sum = temp * sqrt(d__1 * d__1 + d__2 * d__2);
    ga = w[ir - 1] / sum;
    gb = w[ir] / sum;

    /*     EXCHANGE THE COLUMNS OF R. */

    i__2 = kdrop;
    for (i = 1; i <= i__2; ++i)
    {
        ++ira;
        j = ira - kdrop;
        temp = w[ira];
        w[ira] = w[j];
        /* L820: */
        w[j] = temp;
    }
    w[ir] = zero;

    /*     APPLY THE ROTATION TO THE ROWS OF R. */

    w[j] = sum;
    ++kdrop;
    i__2 = nu;
    for (i = kdrop; i <= i__2; ++i)
    {
        temp = ga * w[ira] + gb * w[ira + 1];
        w[ira + 1] = ga * w[ira + 1] - gb * w[ira];
        w[ira] = temp;
        /* L830: */
        ira += i;
    }

    /*     APPLY THE ROTATION TO THE COLUMNS OF Z. */

    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        ++iz;
        j = iz - *n;
        temp = ga * w[j] + gb * w[iz];
        w[iz] = ga * w[iz] - gb * w[j];
        /* L840: */
        w[j] = temp;
    }

    /*     REVISE IACT AND THE LAGRANGE MULTIPLIERS. */

    iact[kdrop - 1] = iact[kdrop];
    w[kdrop - 1] = w[kdrop];
    if (kdrop < *nact)
    {
        goto L810;
    }
L850:
    --(*nact);
    switch ((int)mflag)
    {
    case 1:
        goto L250;
    case 2:
        goto L610;
    }


    /* ******************************************************************** */



    /*     APPLY GIVENS ROTATION TO REDUCE SOME OF THE SCALAR */
    /*     PRODUCTS IN THE S-PARTITION OF W TO ZERO. */

L860:
    iz = iwz + nu * *n;
L870:
    iz -= *n;
L880:
    is = iws + nu;
    --nu;
    if (nu == *nact)
    {
        goto L900;
    }
    if (w[is] == zero)
    {
        goto L870;
    }
    /* Computing MAX */
    d__3 = (d__1 = w[is - 1], abs(d__1)), d__4 = (d__2 = w[is], abs(d__2));
    temp = max(d__3, d__4);
    /* Computing 2nd power */
    d__1 = w[is - 1] / temp;
    /* Computing 2nd power */
    d__2 = w[is] / temp;
    sum = temp * sqrt(d__1 * d__1 + d__2 * d__2);
    ga = w[is - 1] / sum;
    gb = w[is] / sum;
    w[is - 1] = sum;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        k = iz + *n;
        temp = ga * w[iz] + gb * w[k];
        w[k] = ga * w[k] - gb * w[iz];
        w[iz] = temp;
        /* L890: */
        --iz;
    }
    goto L880;
L900:
    switch ((int)nflag)
    {
    case 1:
        goto L560;
    case 2:
        goto L630;
    }


    /* ******************************************************************** */



    /*     CALCULATE THE MAGNITUDE OF X AN REVISE XMAG. */

L910:
    sum = zero;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        sum += (d__1 = x[i], abs(d__1)) * vfact * ((d__2 = grad[i], abs(d__2))
                + (d__3 = g[i + i * g_dim1] * x[i], abs(d__3)));
        if (*lql)
        {
            goto L920;
        }
        if (sum < 1e-30)
        {
            goto L920;
        }
        vfact *= 1e-10;
        sum *= 1e-10;
        xmag *= 1e-10;
L920:
        ;
    }
    /* L925: */
    xmag = max(xmag, sum);
    switch ((int)jflag)
    {
    case 1:
        goto L420;
    case 2:
        goto L690;
    }


    /* ******************************************************************** */



    /*     PRE-MULTIPLY THE VECTOR IN THE W-PARTITION OF W BY Z TRANSPOSE. */

L930:
    jl = iww + 1;
    iz = iwz;
    i__2 = *n;
    for (i = 1; i <= i__2; ++i)
    {
        is = iws + i;
        w[is] = zero;
        iwwn = iww + *n;
        i__1 = iwwn;
        for (j = jl; j <= i__1; ++j)
        {
            ++iz;
            /* L940: */
            w[is] += w[iz] * w[j];
        }
    }
    switch ((int)kflag)
    {
    case 1:
        goto L350;
    case 2:
        goto L550;
    }
    return 0;
} /* ql0002_ */

#ifdef uNdEfInEd
comments from the converter:
(stderr from f2c)
ql0001:
ql0002:
#endif


#define DMAX1(a, b) ((a) > (b) ? (a) : (b))
#define DMIN1(a, b) ((a) < (b) ? (a) : (b))
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#define NONE 0
#define OBJECT 1
#define CONSTR 2

/***************************************************************/
/*     Global Variables and Data Structures                 */
/***************************************************************/

struct _objective
{
    double val;
    double *grad;
    double mult;
    double mult_L; /* mode A=1 */
    int act_sip;   /* SIP      */
};

struct _constraint
{
    double val;
    double *grad;
    double mult;
    int act_sip;   /* SIP      */
    int d1bind;    /* SR constraints  */
};

struct _parameter
{
    double *x;
    double *bl;
    double *bu;
    double *mult;
    void *cd;      /* Client data pointer */
};

struct _violation
{    /* SIP      */
    int type;
    int index;
};

double  bgbnd, tolfea;
int  nstop, maxit;

struct Tglob_info
{
    int nnineq, M, ncallg, ncallf, mode, modec;
    int tot_actf_sip, tot_actg_sip, nfsip, ncsipl, ncsipn; /* SIP */
}
glob_info;

struct Tglob_prnt
{
    int iprint, info, ipd, iter, initvl, iter_mod;
    FILE *io;
}
glob_prnt;

struct Tglob_grd
{
    double epsmac, rteps, udelta, valnom;
}
glob_grd;

struct Tglob_log
{
    int dlfeas, local, update, first, rhol_is1, d0_is0, get_ne_mult;
}
glob_log;

/* User-accessible stopping criterion (see cfsqpusr.h)          */
extern double objeps;
extern double objrep;
extern double gLgeps;
extern int x_is_new;

/* Workspace        */
int     *iw;
double *w;
int  lenw, leniw;

/***************************************************************/
/*     Memory Utilities                                        */
/***************************************************************/

#ifdef __STDC__
static int  *make_iv(int);
static double  *make_dv(int);
static double  **make_dm(int, int);
static void  free_iv(int *);
static void  free_dv(double *);
static void  free_dm(double **, int);
static double  *convert(double **, int, int);
#else
static int  *make_iv();
static double  *make_dv();
static double  **make_dm();
static void  free_iv();
static void  free_dv();
static void  free_dm();
static double  *convert();
#endif

/***************************************************************/
/*     Utility Subroutines                                     */
/***************************************************************/

#ifdef __STDC__
int
ql0001_(int *, int *, int *, int *, int *, int *, double *, double *,
        double *, double *, double *, double *, double *, double *,
        int *, int *, int *, double *, int *, int *, int *, double *);
static void  diagnl(int, double, double **);
static void  error(const char string[], int *);
static void
estlam(int, int, int *, double, double **, double *, double *, double *,
       struct _constraint *, double *, double *, double *, double *);
static double  *colvec(double **, int, int);
static double  scaprd(int, double *, double *);
static double  smallNumber();
static int  fuscmp(double, double);
static int  indexs(int, int);
static void  matrcp(int, double **, int, double **);
static void  matrvc(int, int, double **, double *, double *);
static void  nullvc(int, double *);
static void
resign(int, int, double *, double *, double *, struct _constraint *,
       double *, int, int);
static void  sbout1(FILE *, int, const char *, double, double *, int, int);
static void  sbout2(FILE *, int, int, const char *, const char *, double *);
static void  shift(int, int, int *);
static double
slope(int, int, int, int, int, struct _objective *, double *, double *,
      double *, double, double, int, double *, int);
static int element(int *, int, int);
#else
int   ql0001_();      /* QLD Subroutine */
static void  diagnl();
static void  error();
static void  estlam();
static double  *colvec();
static double  scaprd();
static double  smallNumber();
static int  fuscmp();
static int  indexs();
static void  matrcp();
static void  matrvc();
static void  nullvc();
static void  resign();
static void  sbout1();
static void  sbout2();
static void  shift();
static double  slope();
static int element();
#endif

/**************************************************************/
/*     Gradients - Finite Difference                          */
/**************************************************************/

#ifdef __STDC__
void  grobfd(int, int, double *, double *, void(*)(int, int,
             double *, double *, void *), void *);
void  grcnfd(int, int, double *, double *, void(*)(int, int,
             double *, double *, void *), void *);
#else
void  grobfd();
void  grcnfd();
#endif

/**************************************************************/
/*     Main routines for optimization -                       */
/**************************************************************/

#ifdef __STDC__
static void
cfsqp1(int, int, int, int, int, int, int, int, int, int, int *, int,
       int, int, int, double, double, int *, int *, struct _parameter *,
       struct _constraint *, struct _objective *, double *,
       void(*)(int, int, double *, double *, void *),
       void(*)(int, int, double *, double *, void *),
       void(*)(int, int, double *, double *,
               void(*)(int, int, double *, double *, void *), void *),
       void(*)(int, int, double *, double *,
               void(*)(int, int, double *, double *, void *), void *));
static void
check(int, int, int, int *, int, int, int, int, int, int, int, int *, double,
      double, struct _parameter *);
static void
initpt(int, int, int, int, int, int, int, struct _parameter *,
       struct _constraint *, void(*)(int, int, double *, double *, void *),
       void(*)(int, int, double *, double *,
               void(*)(int, int, double *, double *, void *), void *));
static void
dir(int, int, int, int, int, int, int, int, int, int, int, int, double *,
    double, double, double *, double *, double, double *, double *, int *,
    int *, int *, int *, int *, int *, struct _parameter *, double *,
    double *, struct _constraint *, struct _objective *, double *,
    double *, double *, double *, double *, double *, double **, double *,
    double *, double *, double *, double **, double **, double *,
    double *, struct _violation *, void(*)(int, int, double *, double *,
                                           void *), void(*)(int, int, double *, double *, void *));
static void
step1(int, int, int, int, int, int, int, int, int, int, int, int *, int *, int *,
      int *, int *, int *, int *, int *, int, double, struct _objective *,
      double *, double *, double *, double *, double *, double *, double *,
      double *, double *, double *, double *, struct _constraint *,
      double *, double *, struct _violation *viol,
      void(*)(int, int, double *, double *, void *),
      void(*)(int, int, double *, double *, void *), void *);
static void
hessian(int, int, int, int, int, int, int, int, int, int, int, int *, int,
        double *, struct _parameter *, struct _objective *,
        double, double *, double *, double *, double *, double *,
        struct _constraint *, double *, int *, int *, double *,
        double *, double *, double **, double *, double, int *,
        double *, double *, void(*)(int, int, double *, double *, void *),
        void(*)(int, int, double *, double *, void *),
        void(*)(int, int, double *, double *,
                void(*)(int, int, double *, double *, void *), void *),
        void(*)(int, int, double *, double *,
                void(*)(int, int, double *, double *, void *), void *),
        double **, double *, double *, struct _violation *);
static void
out(int, int, int, int, int, int, int, int, int, int, int, int *, double *,
    struct _constraint *, struct _objective *, double,
    double, double, double, double, int);
static void
update_omega(int, int, int, int *, int, int, int, int, double, double,
             struct _constraint *, struct _objective *, double *,
             struct _violation *, void(*)(int, int, double *, double *,
                                          void *), void(*)(int, int, double *, double *, void *),
             void(*)(int, int, double *, double *,
                     void(*)(int, int, double *, double *, void *), void *),
             void(*)(int, int, double *, double *,
                     void(*)(int, int, double *, double *, void *), void *),
             void *, int);
#else
static void  cfsqp1();
static void  check();
static void initpt();
static void  dir();
static void  step1();
static void  hessian();
static void  out();
static void update_omega();
#endif

#ifdef __STDC__
static void
dealloc(int, int, double *, int *, int *, struct _constraint *cs,
        struct _parameter *);
#else
static void dealloc();
#endif

#ifdef __STDC__
void
cfsqp(int nparam, int nf, int nfsr, int nineqn, int nineq, int neqn,
      int neq, int ncsrl, int ncsrn, int *mesh_pts,
      int mode, int iprint, int miter, int *inform, double bigbnd,
      double eps, double epseqn, double udelta, double *bl, double *bu,
      double *x, double *f, double *g, double *lambda,
      void(*obj)(int, int, double *, double *, void *),
      void(*constr)(int, int, double *, double *, void *),
      void(*gradob)(int, int, double *, double *,
                    void(*)(int, int, double *, double *, void *), void *),
      void(*gradcn)(int, int, double *, double *,
                    void(*)(int, int, double *, double *, void *), void *),
      void *cd)
#else
void
cfsqp(nparam, nf, nfsr, nineqn, nineq, neqn, neq, ncsrl, ncsrn, mesh_pts,
      mode, iprint, miter, inform, bigbnd, eps, epseqn, udelta, bl, bu, x,
      f, g, lambda, obj, constr, gradob, gradcn, cd)
int  nparam, nf, nfsr, neqn, nineqn, nineq, neq, ncsrl, ncsrn, mode,
iprint, miter, *mesh_pts, *inform;
double bigbnd, eps, epseqn, udelta;
double *bl, *bu, *x, *f, *g, *lambda;
void(* obj)(), (* constr)(), (* gradob)(), (* gradcn)();
void    *cd;
#endif

/*---------------------------------------------------------------------
* Brief specification of various arrays and parameters in the calling
* sequence. See manual for a more detailed description.
*
* nparam : number of variables
* nf     : number of objective functions (count each set of sequentially
*          related objective functions once)
* nfsr   : number of sets of sequentially related objectives (possibly
*          zero)
* nineqn : number of nonlinear inequality constraints
* nineq  : total number of inequality constraints
* neqn   : number of nonlinear equality constraints
* neq    : total number of equality constraints
* ncsrl  : number of sets of linear sequentially related inequality
*          constraints
* ncsrn  : number of sets of nonlinear sequentially related inequality
*          constraints
* mesh_pts : array of integers giving the number of actual objectives/
*            constraints in each sequentially related objective or
*      constraint set. The order is as follows:
*      (i) objective sets, (ii) nonlinear constraint sets,
*            (iii) linear constraint sets. If one or no sequentially
*            related constraint or objectives sets are present, the
*            user may simply pass the address of an integer variable
*            containing the appropriate number (possibly zero).
* mode   : mode=CBA specifies job options as described below:
*          A = 0 : ordinary minimax problems
*            = 1 : ordinary minimax problems with each individual
*                  function replaced by its absolute value, ie,
*                  an L_infty problem
*          B = 0 : monotone decrease of objective function
*                  after each iteration
*            = 1 : monotone decrease of objective function after
*                  at most four iterations
*          C = 1 : default operation.
*            = 2 : requires that constraints always be evaluated
*                  before objectives during the line search.
* iprint : print level indicator with the following options-
*          iprint=0: no normal output, only error information
*                    (this option is imposed during phase 1)
*          iprint=1: a final printout at a local solution
*          iprint=2: a brief printout at the end of each iteration
*          iprint=3: detailed infomation is printed out at the end
*                    of each iteration (for debugging purposes)
*          For iprint=2 or 3, the information may be printed at
*          iterations that are multiples of 10, instead of every
*          iteration. This may be done by adding the desired number
*          of iterations to skip printing to the desired iprint value
*          as specified above. e.g., sending iprint=23 would give
*          the iprint=3 information once every 20 iterations.
* miter  : maximum number of iterations allowed by the user to solve
*          the problem
* inform : status report at the end of execution
*          inform= 0:normal termination
*          inform= 1:no feasible point found for linear constraints
*          inform= 2:no feasible point found for nonlinear constraints
*          inform= 3:no solution has been found in miter iterations
*          inform= 4:stepsize smaller than machine precision before
*                    a successful new iterate is found
*          inform= 5:failure in attempting to construct d0
*          inform= 6:failure in attempting to construct d1
*          inform= 7:inconsistent input data
*    inform= 8:new iterate essentially identical to previous
*       iterate, though stopping criterion not satisfied.
*          inform= 9:penalty parameter too large, unable to satisfy
*                    nonlinear equality constraint
* bigbnd : plus infinity
* eps    : stopping criterion. Execution stopped when the norm of the
*          Newton direction vector is smaller than eps
* epseqn : tolerance of the violation of nonlinear equality constraints
*          allowed by the user at an optimal solution
* udelta : perturbation size in computing gradients by finite
*          difference. The actual perturbation is determined by
*          sign(x_i) X max{udelta, rteps X max{1, |x_i|}} for each
*          component of x, where rteps is the square root of machine
*          precision.
* bl     : array of dimension nparam,containing lower bound of x
* bu     : array of dimension nparam,containing upper bound of x
* x      : array of dimension nparam,containing initial guess in input
*          and final iterate at the end of execution
* f      : array of dimension sufficient enough to hold the value of
*          all regular objective functions and the value of all
*          members of the sequentially related objective sets.
*          (dimension must be at least 1)
* g      : array of dimension sufficient enough to hold the value of
*          all regular constraint functions and the value of all
*          members of the sequentially related constraint sets.
*          (dimension must be at least 1)
* lambda : array of dimension nparam+dim(f)+dim(g), containing
*          Lagrange multiplier values at x in output. (A concerns the
*          mode, see above). The first nparam positions contain the
*          multipliers associated with the simple bounds, the next
*          dim(g) positions contain the multipliers associated with
*          the constraints. The final dim(f) positions contain the
*          multipliers associated with the objective functions. The
*          multipliers are in the order they were specified in the
*          user-defined objective and constraint functions.
* obj    : Pointer to function that returns the value of objective
*          functions, one upon each call
* constr : Pointer to function that returns the value of constraints
*          one upon each call
* gradob : Pointer to function that computes gradients of f,
*          alternatively it can be replaced by grobfd to compute
*          finite difference approximations
* gradcn : Pointer to function that computes gradients of g,
*          alternatively it can be replaced by grcnfd to compute
*          finite difference approximations
* cd     : Void pointer that may be used by the user for the passing of
*          "client data" (untouched by CFSQP)
*
*----------------------------------------------------------------------
*
*
*                       CFSQP  Version 2.5b
*
*                  Craig Lawrence, Jian L. Zhou
*                         and Andre Tits
*                  Institute for Systems Research
*                               and
*                Electrical Engineering Department
*                     University of Maryland
*                     College Park, Md 20742
*
*                            June, 1997
*
*
*  The purpose of CFSQP is to solve general nonlinear constrained
*  minimax optimization problems of the form
*
*   (A=0 in mode)     minimize    max_i f_i(x)   for i=1,...,n_f
*                        or
*   (A=1 in mode)     minimize    max_j |f_i(x)|   for i=1,...,n_f
*                       s.t.      bl   <= x <=  bu
*                                 g_j(x) <= 0,   for j=1,...,nineqn
*                                 A_1 x - B_1 <= 0
*
*                                 h_i(x)  = 0,   for i=1,...,neqn
*                                 A_2 x - B_2  = 0
*
* CFSQP is also able to efficiently handle problems with large sets of
* sequentially related objectives or constraints, see the manual for
* details.
*
*
*                  Conditions for External Use
*                  ===========================
*
*   1. The CFSQP routines may not be distributed to third parties.
*      Interested parties should contact the authors directly.
*   2. If modifications are performed on the routines, these
*      modifications will be communicated to the authors.  The
*      modified routines will remain the sole property of the authors.
*   3. Due acknowledgment  must be  made of the  use of the CFSQP
*      routines in research  reports  or  publications.  Whenever
*      such reports are released for public access, a copy should
*      be forwarded to the authors.
*   4. The CFSQP  routines  may  only  be  used  for research and
*      development, unless it has been agreed  otherwise with the
*      authors in writing.
*
* Copyright (c) 1993-1997 by Craig T. Lawrence, Jian L. Zhou, and
*                         Andre L. Tits
* All Rights Reserved.
*
*
* Enquiries should be directed to:
*
*      Prof. Andre L. Tits
*      Electrical Engineering Dept.
*      and Systems Research Center
*      University of Maryland
*      College Park, Md 20742
*      U. S. A.
*
*      Phone : 301-405-3669
*      Fax   : 301-405-6707
*      E-mail: andre@eng.umd.edu
*
*  References:
*  [1] E. Panier and A. Tits, `On Combining Feasibility, Descent and
*      Superlinear Convergence In Inequality Constrained Optimization',
*      Mathematical Programming, Vol. 59(1993), 261-276.
*  [2] J. F. Bonnans, E. Panier, A. Tits and J. Zhou, `Avoiding the
*      Maratos Effect by Means of a Nonmonotone Line search: II.
*      Inequality Problems - Feasible Iterates', SIAM Journal on
*      Numerical Analysis, Vol. 29, No. 4, 1992, pp. 1187-1202.
*  [3] J.L. Zhou and A. Tits, `Nonmonotone Line Search for Minimax
*      Problems', Journal of Optimization Theory and Applications,
*      Vol. 76, No. 3, 1993, pp. 455-476.
*  [4] C.T. Lawrence, J.L. Zhou and A. Tits, `User's Guide for CFSQP
*      Version 2.5: A C Code for Solving (Large Scale) Constrained
*      Nonlinear (Minimax) Optimization Problems, Generating Iterates
*      Satisfying All Inequality Constraints,' Institute for
*      Systems Research, University of Maryland,Technical Report
*      TR-94-16r1, College Park, MD 20742, 1997.
*  [5] C.T. Lawrence and A.L. Tits, `Nonlinear Equality Constraints
*      in Feasible Sequential Quadratic Programming,' Optimization
*      Methods and Software, Vol. 6, March, 1996, pp. 265-282.
*  [6] J.L. Zhou and A.L. Tits, `An SQP Algorithm for Finely
*      Discretized Continuous Minimax Problems and Other Minimax
*      Problems With Many Objective Functions,' SIAM Journal on
*      Optimization, Vol. 6, No. 2, May, 1996, pp. 461--487.
*  [7] C. T. Lawrence and A. L. Tits, `Feasible Sequential Quadratic
*      Programming for Finely Discretized Problems from SIP,'
*      To appear in R. Reemtsen, J.-J. Ruckmann (eds.): Semi-Infinite
*      Programming, in the series Nonconcex Optimization and its
*      Applications. Kluwer Academic Publishers, 1997.
*
***********************************************************************
*/
{
    int i, ipp, j, ncnstr, nclin, nctotl, nob, nobL, modem, nn,
    nppram, nrowa, ncsipl1, ncsipn1, nfsip1;
    int  feasbl, feasb, prnt, Linfty;
    int *indxob, *indxcn, *mesh_pts1;
    double *signeq;
    double xi, gi, gmax, dummy, epskt;
    struct _constraint *cs;      /* pointer to array of constraints */
    struct _objective  *ob;      /* pointer to array of objectives  */
    struct _parameter  *param;   /* pointer to parameter structure  */
    struct _parameter  _param;

    /*     Make adjustments to parameters for SIP constraints       */
    glob_info.tot_actf_sip = glob_info.tot_actg_sip = 0;
    mesh_pts = mesh_pts - 1;
    glob_info.nfsip = nfsr;
    glob_info.ncsipl = ncsrl;
    glob_info.ncsipn = ncsrn;
    nf = nf - nfsr;
    nfsip1 = nfsr;
    nfsr = 0;
    for (i = 1; i <= nfsip1; i++)
        nfsr = nfsr + mesh_pts[i];
    nf = nf + nfsr;
    nineqn = nineqn - ncsrn;
    nineq = nineq - ncsrl - ncsrn;
    ncsipl1 = ncsrl;
    ncsipn1 = ncsrn;
    ncsrl = 0;
    ncsrn = 0;
    if (ncsipn1)
        for (i = 1; i <= ncsipn1; i++)
            ncsrn = ncsrn + mesh_pts[nfsip1+i];
    if (ncsipl1)
        for (i = 1; i <= ncsipl1; i++)
            ncsrl = ncsrl + mesh_pts[nfsip1+ncsipn1+i];
    nineqn = nineqn + ncsrn;
    nineq = nineq + ncsrn + ncsrl;
    /* Create array of constraint structures    */
    cs = (struct _constraint *)calloc(nineq + neq + 1,
                                      sizeof(struct _constraint));
    for (i = 1; i <= nineq + neq; i++)
    {
        cs[i].grad = make_dv(nparam);
        cs[i].act_sip = FALSE;
        cs[i].d1bind = FALSE;
    }
    /* Create parameter structure     */
    _param.x = make_dv(nparam + 1);
    _param.bl = make_dv(nparam);
    _param.bu = make_dv(nparam);
    _param.mult = make_dv(nparam + 1);
    param = &_param;

    /*   Initialize, compute the machine precision, etc.   */
    bl = bl - 1;
    bu = bu - 1;
    x = x - 1;
    for (i = 1; i <= nparam; i++)
    {
        param->x[i] = x[i];
        param->bl[i] = bl[i];
        param->bu[i] = bu[i];
    }
    param->cd = cd;    /* Initialize client data */
    dummy = 0.e0;
    f = f - 1;
    g = g - 1;
    lambda = lambda - 1;
    glob_prnt.iter = 0;
    nstop = 1;
    nn = nineqn + neqn;
    glob_grd.epsmac = smallNumber();
    tolfea = glob_grd.epsmac * 1.e2;
    bgbnd = bigbnd;
    glob_grd.rteps = sqrt(glob_grd.epsmac);
    glob_grd.udelta = udelta;
    glob_log.rhol_is1 = FALSE;
    glob_log.get_ne_mult = FALSE;
    signeq = make_dv(neqn);

    nob = 0;
    gmax = -bgbnd;
    glob_prnt.info = 0;
    glob_prnt.iprint = iprint % 10;
    ipp = iprint;
    glob_prnt.iter_mod = DMAX1(iprint - iprint % 10, 1);
    glob_prnt.io = stdout;
    ncnstr = nineq + neq;
    glob_info.nnineq = nineq;
    if (glob_prnt.iprint > 0)
    {
        fprintf(glob_prnt.io,
                "\n\n    CFSQP Version 2.5b (Released June 1997) \n");
        fprintf(glob_prnt.io,
                "          Copyright (c) 1993 --- 1997       \n");
        fprintf(glob_prnt.io,
                "           C.T. Lawrence, J.L. Zhou         \n");
        fprintf(glob_prnt.io,
                "                and A.L. Tits               \n");
        fprintf(glob_prnt.io,
                "             All Rights Reserved            \n\n");
    }
    /*-----------------------------------------------------*/
    /*   Check the input data      */
    /*-----------------------------------------------------*/
    check(nparam, nf, nfsr, &Linfty, nineq, nineqn, neq, neqn,
          ncsrl, ncsrn, mode, &modem, eps, bgbnd, param);
    if (glob_prnt.info == 7)
    {
        *inform = glob_prnt.info;
        return;
    }

    maxit = DMAX1(DMAX1(miter, 10 * DMAX1(nparam, ncnstr)), 1000);
    feasb = TRUE;
    feasbl = TRUE;
    prnt = FALSE;
    nppram = nparam + 1;

    /*-----------------------------------------------------*/
    /*   Check whether x is within bounds    */
    /*-----------------------------------------------------*/
    for (i = 1; i <= nparam; i++)
    {
        xi = param->x[i];
        if (param->bl[i] <= xi && param->bu[i] >= xi)
            continue;
        feasbl = FALSE;
        break;
    }
    nclin = ncnstr - nn;
    /*-----------------------------------------------------*/
    /*   Check whether linear constraints are feasbile     */
    /*-----------------------------------------------------*/
    if (nclin != 0)
    {
        for (i = 1; i <= nclin; i++)
        {
            j = i + nineqn;
            if (j <= nineq)
            {
                constr(nparam, j, (param->x) + 1, &gi, param->cd);
                if (gi > glob_grd.epsmac)
                    feasbl = FALSE;
            }
            else
            {
                constr(nparam, j + neqn, (param->x) + 1, &gi, param->cd);
                if (fabs(gi) > glob_grd.epsmac)
                    feasbl = FALSE;
            }
            cs[j].val = gi;
        }
    }
    /*-------------------------------------------------------*/
    /*   Generate a new point if infeasible      */
    /*-------------------------------------------------------*/
    if (!feasbl)
    {
        if (glob_prnt.iprint > 0)
        {
            fprintf(glob_prnt.io,
                    " The given initial point is infeasible for inequality\n");
            fprintf(glob_prnt.io,
                    " constraints and linear equality constraints:\n");
            sbout1(glob_prnt.io, nparam, "                    ", dummy,
                   param->x, 2, 1);
            prnt = TRUE;
        }
        nctotl = nparam + nclin;
        lenw = 2 * nparam * nparam + 10 * nparam + 2 * nctotl + 1;
        leniw = DMAX1(2 * nparam + 2 * nctotl + 3, 2 * nclin + 2 * nparam + 6);
        /*-----------------------------------------------------*/
        /*   Attempt to generate a point satisfying all linear */
        /*   constraints.          */
        /*-----------------------------------------------------*/
        nrowa = DMAX1(nclin, 1);
        iw = make_iv(leniw);
        w = make_dv(lenw);
        initpt(nparam, nineqn, neq, neqn, nclin, nctotl, nrowa, param,
               &cs[nineqn], constr, gradcn);
        free_iv(iw);
        free_dv(w);
        if (glob_prnt.info != 0)
        {
            *inform = glob_prnt.info;
            return;
        }
    }
    indxob = make_iv(DMAX1(nineq + neq, nf));
    indxcn = make_iv(nineq + neq);
L510:
    if (glob_prnt.info != -1)
    {
        for (i = 1; i <= nineqn; i++)
        {
            constr(nparam, i, (param->x) + 1, &(cs[i].val), param->cd);
            if (cs[i].val > 0.e0)
                feasb = FALSE;
        }
        glob_info.ncallg = nineqn;
        if (!feasb)
        {
            /* Create array of objective structures for Phase 1  */
            ob = (struct _objective *)calloc(nineqn + 1,
                                             sizeof(struct _objective));
            for (i = 1; i <= nineqn; i++)
            {
                ob[i].grad = make_dv(nparam);
                ob[i].act_sip = FALSE;
            }
            for (i = 1; i <= nineqn; i++)
            {
                nob++;
                indxob[nob] = i;
                ob[nob].val = cs[i].val;
                gmax = DMAX1(gmax, ob[nob].val);
            }
            for (i = 1; i <= nineq - nineqn; i++)
                indxcn[i] = nineqn + i;
            for (i = 1; i <= neq - neqn; i++)
                indxcn[i+nineq-nineqn] = nineq + neqn + i;
            goto L605;
        }
    }

    /* Create array of objective structures for Phase 2 and      */
    /* initialize.      */
    ob = (struct _objective *)calloc(nf + 1, sizeof(struct _objective));
    for (i = 1; i <= nf; i++)
    {
        ob[i].grad = make_dv(nparam);
        ob[i].act_sip = FALSE;
    }
    for (i = 1; i <= nineqn; i++)
    {
        indxcn[i] = i;
    }
    for (i = 1; i <= neq - neqn; i++)
        cs[i+nineq+neqn].val = cs[i+nineq].val;
    for (i = 1; i <= neqn; i++)
    {
        j = i + nineq;
        constr(nparam, j, (param->x) + 1, &(cs[j].val), param->cd);
        indxcn[nineqn+i] = j;
    }
    for (i = 1; i <= nineq - nineqn; i++)
        indxcn[i+nn] = nineqn + i;
    for (i = 1; i <= neq - neqn; i++)
        indxcn[i+nineq+neqn] = nineq + neqn + i;
    glob_info.ncallg += neqn;

    L605:
    if (glob_prnt.iprint > 0 && feasb && !prnt)
    {
        fprintf(glob_prnt.io,
                "The given initial point is feasible for inequality\n");
        fprintf(glob_prnt.io,
                "         constraints and linear equality constraints:\n");
        sbout1(glob_prnt.io, nparam, "                    ", dummy,
               param->x, 2, 1);
        prnt = TRUE;
    }
    if (nob == 0)
    {
        if (glob_prnt.iprint > 0)
        {
            if (glob_prnt.info != 0)
            {
                fprintf(glob_prnt.io,
                        "To generate a feasible point for nonlinear inequality\n");
                fprintf(glob_prnt.io,
                        "constraints and linear equality constraints, ");
                fprintf(glob_prnt.io, "ncallg = %10d\n", glob_info.ncallg);
                if (ipp == 0)
                    fprintf(glob_prnt.io, " iteration           %26d\n",
                            glob_prnt.iter);
                if (ipp > 0)
                    fprintf(glob_prnt.io, " iteration           %26d\n",
                            glob_prnt.iter - 1);
                if (ipp == 0)
                    glob_prnt.iter++;
            }
            if (feasb && !feasbl)
            {
                fprintf(glob_prnt.io,
                        "Starting from the generated point feasible for\n");
                fprintf(glob_prnt.io,
                        "inequality constraints and linear equality constraints:\n");
                sbout1(glob_prnt.io, nparam, "                    ",
                       dummy, param->x, 2, 1);

            }
            if (glob_prnt.info != 0 || !prnt || !feasb)
            {
                fprintf(glob_prnt.io,
                        "Starting from the generated point feasible for\n");
                fprintf(glob_prnt.io,
                        "inequality constraints and linear equality constraints:\n");
                sbout1(glob_prnt.io, nparam, "                    ",
                       dummy, param->x, 2, 1);
            }
        }
        feasb = TRUE;
        feasbl = TRUE;
    }
    if (ipp > 0 && !feasb && !prnt)
    {
        fprintf(glob_prnt.io,
                " The given initial point is infeasible for inequality\n");
        fprintf(glob_prnt.io,
                " constraints and linear equality constraints:\n");
        sbout1(glob_prnt.io, nparam, "                    ", dummy,
               param->x, 2, 1);
        prnt = TRUE;
    }
    if (nob == 0)
        nob = 1;
    if (feasb)
    {
        nob = nf;
        glob_prnt.info = 0;
        glob_prnt.iprint = iprint % 10;
        ipp = iprint;
        glob_prnt.iter_mod = DMAX1(iprint - iprint % 10, 1);
        glob_info.mode = modem;
        epskt = eps;
        if (Linfty)
            nobL = 2 * nob;
        else
            nobL = nob;
        if (nob != 0 || neqn != 0)
            goto L910;
        fprintf(glob_prnt.io,
                "current feasible iterate with no objective specified\n");
        *inform = glob_prnt.info;
        for (i = 1; i <= nineq + neq; i++)
            g[i] = cs[i].val;
        dealloc(nineq, neq, signeq, indxcn, indxob, cs, param);
        free((char *) ob);
        return;
    }
    ipp = 0;
    glob_info.mode = 0;
    nobL = nob;
    glob_prnt.info = -1;
    epskt = 1.e-10;
    L910:
    nctotl = nppram + ncnstr + DMAX1(nobL, 1);
    leniw = 2 * (ncnstr + DMAX1(nobL, 1)) + 2 * nppram + 6;
    lenw = 2 * nppram * nppram + 10 * nppram + 6 * (ncnstr + DMAX1(nobL, 1) + 1);
    glob_info.M = 4;
    if (modem == 1 && nn == 0)
        glob_info.M = 3;

    param->x[nparam+1] = gmax;
    if (feasb)
    {
        for (i = 1; i <= neqn; i++)
        {
            if (cs[i+nineq].val > 0.e0)
                signeq[i] = -1.e0;
            else
                signeq[i] = 1.e0;
        }
    }
    if (!feasb)
    {
        ncsipl1 = ncsrl;
        ncsipn1 = 0;
        nfsip1 = ncsrn;
        mesh_pts1 = &mesh_pts[nfsr];
    }
    else
    {
        ncsipl1 = ncsrl;
        ncsipn1 = ncsrn;
        nfsip1 = nfsr;
        mesh_pts1 = mesh_pts;
    }
    /*---------------------------------------------------------------*/
    /*    either attempt to generate a point satisfying all      */
    /*    constraints or try to solve the original problem           */
    /*---------------------------------------------------------------*/
    nrowa = DMAX1(ncnstr + DMAX1(nobL, 1), 1);
    w = make_dv(lenw);
    iw = make_iv(leniw);

    cfsqp1(miter, nparam, nob, nobL, nfsip1, nineqn, neq, neqn, ncsipl1, ncsipn1,
           mesh_pts1, ncnstr, nctotl, nrowa, feasb, epskt, epseqn, indxob,
           indxcn, param, cs, ob, signeq, obj, constr, gradob, gradcn);

    free_iv(iw);
    free_dv(w);
    if (glob_prnt.info == -1)
    { /* Successful phase 1 termination  */
        for (i = 1; i <= nob; i++)
            cs[i].val = ob[i].val;
        nob = 0;
        for (i = 1; i <= nineqn; i++)
            free_dv(ob[i].grad);
        free((char *) ob);
        goto L510;
    }
    if (glob_prnt.info != 0)
    {
        if (feasb)
        {
            for (i = 1; i <= nparam; i++)
                x[i] = param->x[i];
            for (i = 1; i <= nineq + neq; i++)
                g[i] = cs[i].val;
            *inform = glob_prnt.info;
            dealloc(nineq, neq, signeq, indxcn, indxob, cs, param);
            for (i = 1; i <= nf; i++)
            {
                f[i] = ob[i].val;
                free_dv(ob[i].grad);
            }
            free((char *) ob);
            return;
        }
        glob_prnt.info = 2;
        fprintf(glob_prnt.io,
                "Error: No feasible point is found for nonlinear inequality\n");
        fprintf(glob_prnt.io,
                "constraints and linear equality constraints\n");
        *inform = glob_prnt.info;
        dealloc(nineq, neq, signeq, indxcn, indxob, cs, param);
        for (i = 1; i <= nineqn; i++)
            free_dv(ob[i].grad);
        free((char *) ob);
        return;
    }
    /* Successful phase 2 termination    */
    *inform = glob_prnt.info;
    for (i = 1; i <= nparam; i++)
    {
        x[i] = param->x[i];
        lambda[i] = param->mult[i];
    }
    for (i = 1; i <= nineq + neq; i++)
    {
        g[i] = cs[i].val;
        lambda[i+nparam] = cs[i].mult;
    }
    for (i = 1; i <= nf; i++)
    {
        f[i] = ob[i].val;
        lambda[i+nparam+nineq+neq] = ob[i].mult;
        free_dv(ob[i].grad);
    }
    /* If just one objective, set multiplier=1 */
    if (nf == 1)
        lambda[1+nparam+nineq+neq] = 1.e0;
    free((char *) ob);
    dealloc(nineq, neq, signeq, indxcn, indxob, cs, param);
    return;
}

/***************************************************************/
/*     Free allocated memory           */
/***************************************************************/

#ifdef __STDC__
static void
dealloc(int nineq, int neq, double *signeq, int *indxob,
        int *indxcn, struct _constraint *cs, struct _parameter *param)
#else
static void
dealloc(nineq, neq, signeq, indxob, indxcn, cs, param)
int nineq, neq;
double *signeq;
int    *indxob, *indxcn;
struct _constraint *cs;
struct _parameter  *param;
#endif
{
    int i;

    free_dv(param->x);
    free_dv(param->bl);
    free_dv(param->bu);
    free_dv(param->mult);
    free_dv(signeq);
    free_iv(indxob);
    free_iv(indxcn);
    for (i = 1; i <= nineq + neq; i++)
        free_dv(cs[i].grad);
    free((char *) cs);
}
/************************************************************/
/*   CFSQP : Main routine                                   */
/************************************************************/


#ifdef __STDC__
static void
dealloc1(int, int, double **, double **, double **, double *, double *,
         double *, double *, double *, double *, double *, double *,
         double *, double *, double *, double *, int *, int *, int *);
#else
static void dealloc1();
#endif

#ifdef __STDC__
static void
cfsqp1(int miter, int nparam, int nob, int nobL, int nfsip, int nineqn,
       int neq, int neqn, int ncsipl, int ncsipn, int *mesh_pts, int ncnstr,
       int nctotl, int nrowa, int feasb, double epskt, double epseqn,
       int *indxob, int *indxcn, struct _parameter *param,
       struct _constraint *cs, struct _objective *ob,
       double *signeq, void(*obj)(int, int, double *, double *, void *),
       void(*constr)(int, int, double *, double *, void *),
       void(*gradob)(int, int, double *, double *,
                     void(*)(int, int, double *, double *, void *), void *),
       void(*gradcn)(int, int, double *, double *,
                     void(*)(int, int, double *, double *, void *), void *))
#else
static void
cfsqp1(miter, nparam, nob, nobL, nfsip, nineqn, neq, neqn, ncsipl, ncsipn,
       mesh_pts, ncnstr, nctotl, nrowa, feasb, epskt, epseqn, indxob,
       indxcn, param, cs, ob, signeq, obj, constr, gradob, gradcn)
int miter, nparam, nob, nobL, nfsip, nineqn, neq, neqn, ncnstr,
nctotl, nrowa, feasb, ncsipl, ncsipn, *mesh_pts;
int *indxob, *indxcn;
double  epskt, epseqn;
double  *signeq;
struct _constraint *cs;
struct _objective  *ob;
struct _parameter  *param;
void(* obj)(), (* constr)(), (* gradob)(), (* gradcn)();
#endif
{
    int   i, iskp, nfs, ncf, ncg, nn, nstart, nrst, ncnst1;
    int   *iact, *iskip, *istore;
    double Cbar, Ck, dbar, fmax, fM, fMp, steps, d0nm, dummy,
    sktnom, scvneq, grdftd, psf;
    double *di, *d, *gm, *grdpsf, *penp, *bl, *bu, *clamda,
    *cvec, *psmu, *span, *backup;
    double **hess, **hess1, **a;
    double *tempv;
    struct _violation *viol;
    struct _violation _viol;

    /*   Allocate memory                              */

    hess = make_dm(nparam, nparam);
    hess1 = make_dm(nparam + 1, nparam + 1);
    a = make_dm(nrowa, nparam + 2);
    di = make_dv(nparam + 1);
    d = make_dv(nparam + 1);
    gm = make_dv(4 * neqn);
    grdpsf = make_dv(nparam);
    penp = make_dv(neqn);
    bl = make_dv(nctotl);
    bu = make_dv(nctotl);
    clamda = make_dv(nctotl + nparam + 1);
    cvec = make_dv(nparam + 1);
    psmu = make_dv(neqn);
    span = make_dv(4);
    backup = make_dv(nob + ncnstr);
    iact = make_iv(nob + nineqn + neqn);
    iskip = make_iv(glob_info.nnineq + 1);
    istore = make_iv(nineqn + nob + neqn);

    viol = &_viol;
    viol->index = 0;
    viol->type = NONE;

    glob_prnt.initvl = 1;
    glob_log.first = TRUE;
    nrst = glob_prnt.ipd = 0;
    dummy = 0.e0;
    scvneq = 0.e0;
    steps = 0.e0;
    sktnom = 0.e0;
    d0nm = 0.e0;
    if (glob_prnt.iter == 0)
        diagnl(nparam, 1.e0, hess);
    if (feasb)
    {
        glob_log.first = TRUE;
        if (glob_prnt.iter > 0)
            glob_prnt.iter--;
        if (glob_prnt.iter != 0)
            diagnl(nparam, 1.e0, hess);
    }
    Ck = Cbar = 1.e-2;
    dbar = 5.e0;
    nstart = 1;
    glob_info.ncallf = 0;
    nstop = 1;
    nfs = 0;
    if (glob_info.mode != 0)
        nfs = glob_info.M;
    if (feasb)
    {
        nn = nineqn + neqn;
        ncnst1 = ncnstr;
    }
    else
    {
        nn = 0;
        ncnst1 = ncnstr - nineqn - neqn;
    }
    scvneq = 0.e0;
    for (i = 1; i <= ncnst1; i++)
    {
        glob_grd.valnom = cs[indxcn[i]].val;
        backup[i] = glob_grd.valnom;
        if (feasb && i > nineqn && i <= nn)
        {
            gm[i-nineqn] = glob_grd.valnom * signeq[i-nineqn];
            scvneq = scvneq + fabs(glob_grd.valnom);
        }
        if (feasb && i <= nn)
        {
            iact[i] = indxcn[i];
            istore[i] = 0;
            if (i > nineqn)
                penp[i-nineqn] = 2.e0;
        }
        gradcn(nparam, indxcn[i], (param->x) + 1, (cs[indxcn[i]].grad) + 1,
               constr, param->cd);
    }
    nullvc(nparam, grdpsf);
    psf = 0.e0;
    if (feasb && neqn != 0)
        resign(nparam, neqn, &psf, grdpsf, penp, cs, signeq, 12, 12);
    fmax = -bgbnd;
    for (i = 1; i <= nob; i++)
    {
        if (!feasb)
        {
            glob_grd.valnom = ob[i].val;
            iact[i] = i;
            istore[i] = 0;
            gradcn(nparam, indxob[i], (param->x) + 1, (ob[i].grad) + 1, constr,
                   param->cd);
        }
        else
        {
            iact[nn+i] = i;
            istore[nn+i] = 0;
            obj(nparam, i, (param->x) + 1, &(ob[i].val), param->cd);
            glob_grd.valnom = ob[i].val;
            backup[i+ncnst1] = glob_grd.valnom;
            gradob(nparam, i, (param->x) + 1, (ob[i].grad) + 1, obj, param->cd);
            glob_info.ncallf++;
            if (nobL != nob)
                fmax = DMAX1(fmax, -ob[i].val);
        }
        fmax = DMAX1(fmax, ob[i].val);
    }
    if (feasb && nob == 0)
        fmax = 0.e0;
    fM = fmax;
    fMp = fmax - psf;
    span[1] = fM;

    if (glob_prnt.iprint >= 3 && glob_log.first)
    {
        for (i = 1; i <= nob; i++)
        {
            if (feasb)
            {
                if (nob > 1)
                {
                    tempv = ob[i].grad;
                    sbout2(glob_prnt.io, nparam, i, "gradf(j,", ")", tempv);
                }
                if (nob == 1)
                {
                    tempv = ob[1].grad;
                    sbout1(glob_prnt.io, nparam, "gradf(j)            ",
                           dummy, tempv, 2, 2);
                }
                continue;
            }
            tempv = ob[i].grad;
            sbout2(glob_prnt.io, nparam, indxob[i], "gradg(j,", ")", tempv);
        }
        if (ncnstr != 0)
        {
            for (i = 1; i <= ncnst1; i++)
            {
                tempv = cs[indxcn[i]].grad;
                sbout2(glob_prnt.io, nparam, indxcn[i], "gradg(j,", ")", tempv);
            }
            if (neqn != 0)
            {
                sbout1(glob_prnt.io, nparam, "grdpsf(j)           ", dummy,
                       grdpsf, 2, 2);
                sbout1(glob_prnt.io, neqn, "P                   ", dummy,
                       penp, 2, 2);
            }
        }
        for (i = 1; i <= nparam; i++)
        {
            tempv = colvec(hess, i, nparam);
            sbout2(glob_prnt.io, nparam, i, "hess (j,", ")", tempv);
            free_dv(tempv);
        }
    }

    /*----------------------------------------------------------*
     *              Main loop of the algorithm                  *
     *----------------------------------------------------------*/

    nstop = 1;
    for (;;)
    {
        out(miter, nparam, nob, nobL, nfsip, nineqn, nn, nineqn, ncnst1,
            ncsipl, ncsipn, mesh_pts, param->x, cs, ob, fM, fmax, steps,
            sktnom, d0nm, feasb);
        if (nstop == 0)
        {
            if (!feasb)
            {
                dealloc1(nparam, nrowa, a, hess, hess1, di, d, gm,
                         grdpsf, penp, bl, bu, clamda, cvec, psmu, span, backup,
                         iact, iskip, istore);
                return;
            }
            for (i = 1; i <= ncnst1; i++)
                cs[i].val = backup[i];
            for (i = 1; i <= nob; i++)
                ob[i].val = backup[i+ncnst1];
            for (i = 1; i <= neqn; i++)
                cs[glob_info.nnineq+i].mult = signeq[i] * psmu[i];
            dealloc1(nparam, nrowa, a, hess, hess1, di, d, gm,
                     grdpsf, penp, bl, bu, clamda, cvec, psmu, span, backup,
                     iact, iskip, istore);
            return;
        }
        if (!feasb && glob_prnt.iprint == 0)
            glob_prnt.iter++;
        /*   Update the SIP constraint set Omega_k  */
        if ((ncsipl + ncsipn) != 0 || nfsip)
            update_omega(nparam, ncsipl, ncsipn, mesh_pts, nineqn, nob, nobL,
                         nfsip, steps, fmax, cs, ob, param->x, viol,
                         constr, obj, gradob, gradcn, param->cd, feasb);
        /*   Compute search direction               */
        dir(nparam, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn,
            ncnst1, feasb, &steps, epskt, epseqn, &sktnom, &scvneq, Ck, &d0nm,
            &grdftd, indxob, indxcn, iact, &iskp, iskip, istore, param, di, d,
            cs, ob, &fM, &fMp, &fmax, &psf, grdpsf, penp, a, bl, bu, clamda, cvec,
            hess, hess1, backup, signeq, viol, obj, constr);
        if (nstop == 0 && !glob_log.get_ne_mult)
            continue;
        glob_log.first = FALSE;
        if (!glob_log.update && !glob_log.d0_is0)
        {
            /*   Determine step length                                */
            step1(nparam, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn,
                  ncnst1, &ncg, &ncf, indxob, indxcn, iact, &iskp, iskip, istore,
                  feasb, grdftd, ob, &fM, &fMp, &fmax, &psf, penp, &steps, &scvneq,
                  bu, param->x, di, d, cs, backup, signeq, viol, obj, constr,
                  param->cd);
            if (nstop == 0)
                continue;
        }
        /*   Update the Hessian                                      */
        hessian(nparam, nob, nfsip, nobL, nineqn, neq, neqn, nn, ncsipn, ncnst1,
                nfs, &nstart, feasb, bu, param, ob, fmax, &fM, &fMp, &psf, grdpsf,
                penp, cs, gm, indxob, indxcn, bl, clamda, di, hess, d, steps, &nrst,
                signeq, span, obj, constr, gradob, gradcn, hess1, cvec, psmu, viol);
        if (nstop == 0 || glob_info.mode == 0)
            continue;
        if (d0nm > dbar)
            Ck = DMAX1(0.5e0 * Ck, Cbar);
        if (d0nm <= dbar && glob_log.dlfeas)
            Ck = Ck;
        if (d0nm <= dbar && !glob_log.dlfeas &&
            !glob_log.rhol_is1)
            Ck = 10.e0 * Ck;
    }
}

/*******************************************************************/
/*    Free up memory used by CFSQP1       */
/*******************************************************************/

#ifdef __STDC__
static void
dealloc1(int nparam, int nrowa, double **a, double **hess, double **hess1,
         double *di, double *d, double *gm, double *grdpsf, double *penp,
         double *bl, double *bu, double *clamda, double *cvec, double *psmu,
         double *span, double *backup, int *iact, int *iskip, int *istore)
#else
static void
dealloc1(nparam, nrowa, a, hess, hess1, di, d, gm, grdpsf, penp, bl, bu, clamda,
         cvec, psmu, span, backup, iact, iskip, istore)
int nparam, nrowa;
double **a, **hess, **hess1;
double  *di, *d, *gm, *grdpsf, *penp, *bl, *bu, *clamda, *cvec, *psmu, *span,
*backup;
int *iact, *iskip, *istore;
#endif
{
    free_dm(a, nrowa);
    free_dm(hess, nparam);
    free_dm(hess1, nparam + 1);
    free_dv(di);
    free_dv(d);
    free_dv(gm);
    free_dv(grdpsf);
    free_dv(penp);
    free_dv(bl);
    free_dv(bu);
    free_dv(clamda);
    free_dv(cvec);
    free_dv(psmu);
    free_dv(span);
    free_dv(backup);
    free_iv(iact);
    free_iv(iskip);
    free_iv(istore);
}
/************************************************************/
/*   CFSQP - Check the input data                           */
/************************************************************/


#ifdef __STDC__
static void
check(int nparam, int nf, int nfsip, int *Linfty, int nineq,
      int nnl, int neq, int neqn, int ncsipl, int ncsipn, int mode,
      int *modem, double eps, double bigbnd, struct _parameter *param)
#else
static void
check(nparam, nf, nfsip, Linfty, nineq, nnl, neq, neqn, ncsipl, ncsipn,
      mode, modem, eps, bigbnd, param)
int nparam, nf, nfsip, nineq, nnl, neq, neqn, ncsipl, ncsipn, mode, *modem,
*Linfty;
double  bigbnd, eps;
struct  _parameter *param;
#endif
{
    int i;
    double bli, bui;

    if (nparam <= 0)
        error("nparam should be positive!                ",
              &glob_prnt.info);
    if (nf < 0)
        error("nf should not be negative!                ",
              &glob_prnt.info);
    if (nineq < 0)
        error("nineq should not be negative!             ",
              &glob_prnt.info);
    if (nineq >= 0 && nnl >= 0 && nineq < nnl)
        error("nineq should be no smaller then nnl!     ",
              &glob_prnt.info);
    if (neqn < 0)
        error("neqn should not be negative!              ",
              &glob_prnt.info);
    if (neq < neqn)
        error("neq should not be smaller then neqn      ",
              &glob_prnt.info);
    if (nf < nfsip)
        error("nf should not be smaller then nfsip      ",
              &glob_prnt.info);
    if (nineq < ncsipn + ncsipl)
        error("ncsrl+ncsrn should not be larger then nineq",
              &glob_prnt.info);
    if (nparam <= neq - neqn)
        error("Must have nparam > number of linear equalities",
              &glob_prnt.info);
    if (glob_prnt.iprint < 0 || glob_prnt.iprint > 3)
        error("iprint mod 10 should be 0,1,2 or 3!       ",
              &glob_prnt.info);
    if (eps <= glob_grd.epsmac)
    {
        error("eps should be bigger than epsmac!         ",
              &glob_prnt.info);
        fprintf(glob_prnt.io,
                "epsmac = %22.14e which is machine dependent.\n",
                glob_grd.epsmac);
    }
    if (!(mode == 100 || mode == 101 || mode == 110 || mode == 111
          || mode == 200 || mode == 201 || mode == 210 || mode == 211))
        error("mode is not properly specified!           ",
              &glob_prnt.info);
    if (glob_prnt.info != 0)
    {
        fprintf(glob_prnt.io,
                "Error: Input parameters are not consistent.\n");
        return;
    }
    for (i = 1; i <= nparam; i++)
    {
        bli = param->bl[i];
        bui = param->bu[i];
        if (bli > bui)
        {
            fprintf(glob_prnt.io,
                    "lower bounds should be smaller than upper bounds\n");
            glob_prnt.info = 7;
        }
        if (glob_prnt.info != 0)
            return;
        if (bli < (-bigbnd))
            param->bl[i] = -bigbnd;
        if (bui > bigbnd)
            param->bu[i] = bigbnd;
    }
    if (mode >= 200)
    {
        i = mode - 200;
        glob_info.modec = 2;
    }
    else
    {
        i = mode - 100;
        glob_info.modec = 1;
    }
    if (i < 10)
        *modem = 0;
    else
    {
        *modem = 1;
        i -= 10;
    }
    if (!i)
        *Linfty = FALSE;
    else
        *Linfty = TRUE;
}
/****************************************************************/
/*    CFSQP : Generate a feasible point satisfying simple       */
/*            bounds and linear constraints.                    */
/****************************************************************/


#ifdef __STDC__
static void
initpt(int nparam, int nnl, int neq, int neqn, int nclin, int nctotl,
       int nrowa, struct _parameter *param, struct _constraint *cs,
       void(*constr)(int, int, double *, double *, void *),
       void(*gradcn)(int, int, double *, double *,
                     void(*)(int, int, double *, double *, void *), void *))
#else
static void
initpt(nparam, nnl, neq, neqn, nclin, nctotl, nrowa, param, cs,
       constr, gradcn)
int nparam, nnl, neq, neqn, nclin, nctotl, nrowa;
struct _constraint *cs;
struct _parameter  *param;
void(* constr)(), (* gradcn)();
#endif
{
    int i, j, infoql, mnn, temp1, iout, zero;
    double x0i, *atemp, *htemp;
    double *x, *bl, *bu, *cvec, *clamda, *bj;
    double **a, **hess;

    hess = make_dm(nparam, nparam);
    a = make_dm(nrowa, nparam);
    x = make_dv(nparam);
    bl = make_dv(nctotl);
    bu = make_dv(nctotl);
    cvec = make_dv(nparam);
    clamda = make_dv(nctotl + nparam + 1);
    bj = make_dv(nclin);

    glob_prnt.info = 1;
    for (i = 1; i <= nclin; i++)
    {
        glob_grd.valnom = cs[i].val;
        j = i + nnl;
        if (j <= glob_info.nnineq)
            gradcn(nparam, j, (param->x) + 1, cs[i].grad + 1, constr, param->cd);
        else
            gradcn(nparam, j + neqn, (param->x) + 1, cs[i].grad + 1, constr,
                   param->cd);
    }
    for (i = 1; i <= nparam; i++)
    {
        x0i = param->x[i];
        bl[i] = param->bl[i] - x0i;
        bu[i] = param->bu[i] - x0i;
        cvec[i] = 0.e0;
    }
    for (i = nclin; i >= 1; i--)
        bj[nclin-i+1] = -cs[i].val;
    for (i = nclin; i >= 1; i--)
        for (j = 1; j <= nparam; j++)
            a[nclin-i+1][j] = -cs[i].grad[j];
    diagnl(nparam, 1.e0, hess);
    nullvc(nparam, x);

    iout = 6;
    zero = 0;
    mnn = nrowa + 2 * nparam;
    iw[1] = 1;
    temp1 = neq - neqn;
    htemp = convert(hess, nparam, nparam);
    atemp = convert(a, nrowa, nparam);

    ql0001_(&nclin, &temp1, &nrowa, &nparam, &nparam, &mnn, (htemp + 1),
            (cvec + 1), (atemp + 1), (bj + 1), (bl + 1), (bu + 1), (x + 1), (clamda + 1),
            &iout, &infoql, &zero, (w + 1), &lenw, (iw + 1), &leniw,
            &glob_grd.epsmac);

    free_dv(htemp);
    free_dv(atemp);
    if (infoql == 0)
    {
        for (i = 1; i <= nparam; i++)
            param->x[i] = param->x[i] + x[i];
        x_is_new = TRUE;
        for (i = 1; i <= nclin; i++)
        {
            j = i + nnl;
            if (j <= glob_info.nnineq)
                constr(nparam, j, (param->x) + 1,
                       &(cs[i].val), param->cd);
            else
                constr(nparam, j + neqn, (param->x) + 1, &(cs[i].val), param->cd);
        }
        glob_prnt.info = 0;
    }
    if (glob_prnt.info == 1 && glob_prnt.iprint != 0)
    {
        fprintf(glob_prnt.io,
                "\n Error: No feasible point is found for the");
        fprintf(glob_prnt.io, " linear constraints.\n");
    }
    free_dm(a, nrowa);
    free_dm(hess, nparam);
    free_dv(x);
    free_dv(bl);
    free_dv(bu);
    free_dv(cvec);
    free_dv(clamda);
    free_dv(bj);
    return;
}
/****************************************************************/
/*   CFSQP : Update the SIP "active" objective and constraint   */
/*           sets Omega_k and Xi_k.                         */
/****************************************************************/


#ifdef __STDC__
static void
update_omega(int nparam, int ncsipl, int ncsipn, int *mesh_pts,
             int nineqn, int nob, int nobL, int nfsip, double steps,
             double fmax, struct _constraint *cs, struct _objective *ob,
             double *x, struct _violation *viol,
             void(*constr)(int, int, double *, double *, void *),
             void(*obj)(int, int, double *, double *, void *),
             void(*gradob)(int, int, double *, double *,
                           void(*)(int, int, double *, double *, void *), void *),
             void(*gradcn)(int, int, double *, double *,
                           void(*)(int, int, double *, double *, void *), void *),
             void *cd, int feasb)
#else
static void
update_omega(nparam, ncsipl, ncsipn, mesh_pts, nineqn, nob, nobL, nfsip,
             steps, fmax, cs, ob, x, viol, constr, obj, gradob, gradcn, cd, feasb)
int nparam, ncsipl, ncsipn, *mesh_pts, nineqn, nobL, nob, nfsip, feasb;
double *x, steps, fmax;
struct _constraint *cs;
struct _objective *ob;
struct _violation *viol;
void(* constr)();
void(* obj)();
void(* gradob)();
void(* gradcn)();
void    *cd;
#endif
{
    int i, j, i_max, index, offset, nineq, display;
    double epsilon, g_max, fprev, fnow, fnext, fmult;

    epsilon = 1.e0;
    glob_info.tot_actf_sip = glob_info.tot_actg_sip = 0;
    nineq = glob_info.nnineq;
    if (glob_prnt.iter % glob_prnt.iter_mod)
        display = FALSE;
    else
        display = TRUE;
    /* Clear previous constraint sets                   */
    for (i = 1; i <= ncsipl; i++)
        cs[nineq-ncsipl+i].act_sip = FALSE;
    for (i = 1; i <= ncsipn; i++)
        cs[nineqn-ncsipn+i].act_sip = FALSE;
    /* Clear previous objective sets                    */
    for (i = nob - nfsip + 1; i <= nob; i++)
        ob[i].act_sip = FALSE;

    /*--------------------------------------------------*/
    /* Update Constraint Sets Omega_k                */
    /*--------------------------------------------------*/

    if (ncsipn != 0)
    {
        offset = nineqn - ncsipn;
        for (i = 1; i <= glob_info.ncsipn; i++)
        {
            for (j = 1; j <= mesh_pts[glob_info.nfsip+i]; j++)
            {
                offset++;
                if (j == 1)
                {
                    if (cs[offset].val >= -epsilon &&
                        cs[offset].val >= cs[offset+1].val)
                    {
                        cs[offset].act_sip = TRUE;
                        glob_info.tot_actg_sip++;
                        if (cs[offset].mult == 0.e0 && !glob_log.first)
                        {
                            glob_grd.valnom = cs[offset].val;
                            gradcn(nparam, offset, x + 1, cs[offset].grad + 1, constr,
                                   cd);
                        }
                        continue;
                    }
                }
                else if (j == mesh_pts[glob_info.nfsip+i])
                {
                    if (cs[offset].val >= -epsilon &&
                        cs[offset].val > cs[offset-1].val)
                    {
                        cs[offset].act_sip = TRUE;
                        glob_info.tot_actg_sip++;
                        if (cs[offset].mult == 0.e0 && !glob_log.first)
                        {
                            glob_grd.valnom = cs[offset].val;
                            gradcn(nparam, offset, x + 1, cs[offset].grad + 1, constr,
                                   cd);
                        }
                        continue;
                    }
                }
                else
                {
                    if (cs[offset].val >= -epsilon && cs[offset].val >
                        cs[offset-1].val && cs[offset].val >=
                        cs[offset+1].val)
                    {
                        cs[offset].act_sip = TRUE;
                        glob_info.tot_actg_sip++;
                        if (cs[offset].mult == 0.e0 && !glob_log.first)
                        {
                            glob_grd.valnom = cs[offset].val;
                            gradcn(nparam, offset, x + 1, cs[offset].grad + 1, constr,
                                   cd);
                        }
                        continue;
                    }
                }
                if (cs[offset].val >= -glob_grd.epsmac)
                {
                    cs[offset].act_sip = TRUE;
                    glob_info.tot_actg_sip++;
                    if (cs[offset].mult == 0.e0 && !glob_log.first)
                    {
                        glob_grd.valnom = cs[offset].val;
                        gradcn(nparam, offset, x + 1, cs[offset].grad + 1, constr, cd);
                    }
                    continue;
                }
                if (cs[offset].mult > 0.e0)
                {
                    cs[offset].act_sip = TRUE;
                    glob_info.tot_actg_sip++;
                }
                /* Add if binding for d1  */
                if (cs[offset].d1bind)
                {
                    cs[offset].act_sip = TRUE;
                    glob_info.tot_actg_sip++;
                    if (cs[offset].mult == 0.e0 && !glob_log.first)
                    {
                        glob_grd.valnom = cs[offset].val;
                        gradcn(nparam, offset, x + 1, cs[offset].grad + 1, constr, cd);
                    }
                }

            }
        }
    }
    if (ncsipl != 0)
    {
        /* Don't need to get gradients */
        offset = nineq - ncsipl;
        for (i = 1; i <= glob_info.ncsipl; i++)
        {
            if (feasb)
                index = glob_info.nfsip + glob_info.ncsipn + i;
            else
                index = glob_info.ncsipn + i;
            for (j = 1; j <= mesh_pts[index]; j++)
            {
                offset++;
                if (j == 1)
                {
                    if (cs[offset].val >= -epsilon &&
                        cs[offset].val >= cs[offset+1].val)
                    {
                        cs[offset].act_sip = TRUE;
                        glob_info.tot_actg_sip++;
                        continue;
                    }
                }
                else
                    if (j == mesh_pts[index])
                    {
                        if (cs[offset].val >= -epsilon &&
                            cs[offset].val > cs[offset-1].val)
                        {
                            cs[offset].act_sip = TRUE;
                            glob_info.tot_actg_sip++;
                            continue;
                        }
                    }
                    else
                    {
                        if (cs[offset].val >= -epsilon && cs[offset].val >
                            cs[offset-1].val && cs[offset].val >=
                            cs[offset+1].val)
                        {
                            cs[offset].act_sip = TRUE;
                            glob_info.tot_actg_sip++;
                            continue;
                        }
                    }
                if (cs[offset].val >= -glob_grd.epsmac ||
                    cs[offset].mult > 0.e0 || cs[offset].d1bind)
                {
                    cs[offset].act_sip = TRUE;
                    glob_info.tot_actg_sip++;
                }
            }
        }
    }
    /* Include some extra points during 1st iteration        */
    /* (gradients are already evaluated for first iteration) */
    /* Current heuristics: maximizers and end-points.        */
    if (glob_log.first)
    {
        if (feasb)
        {
            offset = nineqn - ncsipn;
            for (i = 1; i <= glob_info.ncsipn; i++)
            {
                i_max = ++offset;
                g_max = cs[i_max].val;
                if (!cs[i_max].act_sip)
                { /* add first point       */
                    cs[i_max].act_sip = TRUE;
                    glob_info.tot_actg_sip++;
                }
                for (j = 2;j <= mesh_pts[glob_info.nfsip+i];j++)
                {
                    offset++;
                    if (cs[offset].val > g_max)
                    {
                        i_max = offset;
                        g_max = cs[i_max].val;
                    }
                }
                if (!cs[i_max].act_sip)
                {
                    cs[i_max].act_sip = TRUE;
                    glob_info.tot_actg_sip++;
                }
                if (!cs[offset].act_sip)
                { /* add last point          */
                    cs[offset].act_sip = TRUE;
                    glob_info.tot_actg_sip++;
                }
            }
        }
        offset = nineq - ncsipl;
        for (i = 1; i <= glob_info.ncsipl; i++)
        {
            i_max = ++offset;
            g_max = cs[i_max].val;
            if (!cs[i_max].act_sip)
            { /* add first point       */
                cs[i_max].act_sip = TRUE;
                glob_info.tot_actg_sip++;
            }
            if (feasb)
                index = glob_info.nfsip + glob_info.ncsipn + i;
            else
                index = glob_info.ncsipn + i;
            for (j = 2;j <= mesh_pts[index]; j++)
            {
                offset++;
                if (cs[offset].val > g_max)
                {
                    i_max = offset;
                    g_max = cs[i_max].val;
                }
            }
            if (!cs[i_max].act_sip)
            {
                cs[i_max].act_sip = TRUE;
                glob_info.tot_actg_sip++;
            }
            if (!cs[offset].act_sip)
            { /* add last point          */
                cs[offset].act_sip = TRUE;
                glob_info.tot_actg_sip++;
            }
        }
    }

    /* If necessary, append xi_bar                 */
    if (steps < 1.e0 && viol->type == CONSTR)
    {
        i = viol->index;
        if (!cs[i].act_sip)
        {
            cs[i].act_sip = TRUE;
            glob_info.tot_actg_sip++;
        }
    }
    if (glob_prnt.iprint >= 2 && display)
        fprintf(glob_prnt.io, " |Xi_k| for g %26d\n",
                glob_info.tot_actg_sip);

    for (i = 1; i <= ncsipl; i++)
        cs[nineq-ncsipl+i].d1bind = FALSE;
    for (i = 1; i <= ncsipn; i++)
        cs[nineqn-ncsipn+i].d1bind = FALSE;

    /*---------------------------------------------------------*/
    /* Update Objective Set Omega_k          */
    /*---------------------------------------------------------*/

    if (nfsip)
    {
        offset = nob - nfsip;
        if (feasb)
            index = glob_info.nfsip;
        else
            index = glob_info.ncsipn;
        for (i = 1; i <= index; i++)
        {
            for (j = 1; j <= mesh_pts[i]; j++)
            {
                offset++;
                if (nobL > nob)
                {
                    fnow = fabs(ob[offset].val);
                    fmult = DMAX1(fabs(ob[offset].mult),
                                  fabs(ob[offset].mult_L));
                }
                else
                {
                    fnow = ob[offset].val;
                    fmult = ob[offset].mult;
                }
                if (j == 1)
                {
                    if (nobL > nob)
                        fnext = fabs(ob[offset+1].val);
                    else
                        fnext = ob[offset+1].val;
                    if ((fnow >= fmax - epsilon) && fnow >= fnext)
                    {
                        ob[offset].act_sip = TRUE;
                        glob_info.tot_actf_sip++;
                        if (fmult == 0.e0 && !glob_log.first)
                        {
                            glob_grd.valnom = ob[offset].val;
                            if (feasb)
                                gradob(nparam, offset, x + 1,
                                       ob[offset].grad + 1, obj, cd);
                            else
                                gradcn(nparam, offset, x + 1, ob[offset].grad + 1,
                                       constr, cd);
                        }
                        continue;
                    }
                }
                else if (j == mesh_pts[i])
                {
                    if (nobL > nob)
                        fprev = fabs(ob[offset-1].val);
                    else
                        fprev = ob[offset-1].val;
                    if ((fnow >= fmax - epsilon) && fnow > fprev)
                    {
                        ob[offset].act_sip = TRUE;
                        glob_info.tot_actf_sip++;
                        if (fmult == 0.e0 && !glob_log.first)
                        {
                            glob_grd.valnom = ob[offset].val;
                            if (feasb)
                                gradob(nparam, offset, x + 1,
                                       ob[offset].grad + 1, obj, cd);
                            else
                                gradcn(nparam, offset, x + 1, ob[offset].grad + 1,
                                       constr, cd);
                        }
                        continue;
                    }
                }
                else
                {
                    if (nobL > nob)
                    {
                        fprev = fabs(ob[offset-1].val);
                        fnext = fabs(ob[offset+1].val);
                    }
                    else
                    {
                        fprev = ob[offset-1].val;
                        fnext = ob[offset+1].val;
                    }
                    if ((fnow >= fmax - epsilon) && fnow > fprev &&
                        fnow >= fnext)
                    {
                        ob[offset].act_sip = TRUE;
                        glob_info.tot_actf_sip++;
                        if (fmult == 0.e0 && !glob_log.first)
                        {
                            glob_grd.valnom = ob[offset].val;
                            if (feasb)
                                gradob(nparam, offset, x + 1,
                                       ob[offset].grad + 1, obj, cd);
                            else
                                gradcn(nparam, offset, x + 1, ob[offset].grad + 1,
                                       constr, cd);
                        }
                        continue;
                    }
                }
                if (fnow >= fmax - glob_grd.epsmac && !ob[offset].act_sip)
                {
                    ob[offset].act_sip = TRUE;
                    glob_info.tot_actf_sip++;
                    if (fmult == 0.e0 && !glob_log.first)
                    {
                        glob_grd.valnom = ob[offset].val;
                        if (feasb)
                            gradob(nparam, offset, x + 1,
                                   ob[offset].grad + 1, obj, cd);
                        else
                            gradcn(nparam, offset, x + 1, ob[offset].grad + 1,
                                   constr, cd);
                    }
                    continue;
                }
                if (fmult != 0.e0 && !ob[offset].act_sip)
                {
                    ob[offset].act_sip = TRUE;
                    glob_info.tot_actf_sip++;
                    continue;
                }
            }
        }
        /* Addition of objectives for first iteration.          */
        /* Current heuristics: maximizers and end-points        */
        if (glob_log.first)
        {
            offset = nob - nfsip;
            if (feasb)
                index = glob_info.nfsip;
            else
                index = glob_info.ncsipn;
            for (i = 1; i <= index; i++)
            {
                i_max = ++offset;
                if (nobL == nob)
                    g_max = ob[i_max].val;
                else
                    g_max = fabs(ob[i_max].val);
                if (!ob[i_max].act_sip)
                { /* add first point       */
                    ob[i_max].act_sip = TRUE;
                    glob_info.tot_actf_sip++;
                }
                for (j = 2;j <= mesh_pts[i];j++)
                {
                    offset++;
                    if (nobL == nob)
                        fnow = ob[offset].val;
                    else
                        fnow = fabs(ob[offset].val);
                    if (fnow > g_max)
                    {
                        i_max = offset;
                        g_max = fnow;
                    }
                }
                if (!ob[i_max].act_sip)
                {
                    ob[i_max].act_sip = TRUE;
                    glob_info.tot_actf_sip++;
                }
                if (!ob[offset].act_sip)
                { /* add last point          */
                    ob[offset].act_sip = TRUE;
                    glob_info.tot_actf_sip++;
                }
            }
        }

        /* If necessary, append omega_bar            */
        if (steps < 1.e0 && viol->type == OBJECT)
        {
            i = viol->index;
            if (!ob[i].act_sip)
            {
                ob[i].act_sip = TRUE;
                glob_info.tot_actf_sip++;
            }
        }
        if (glob_prnt.iprint >= 2 && display)
            fprintf(glob_prnt.io, " |Omega_k| for f %26d\n",
                    glob_info.tot_actf_sip);
    }
    viol->type = NONE;
    viol->index = 0;
    return;
}
/*******************************************************************/
/*   CFSQP : Computation of the search direction                   */
/*******************************************************************/


#ifdef __STDC__
static void
dqp(int, int, int, int, int, int, int, int, int, int, int, int, int,
    int, int, int *, struct _parameter *, double *, int,
    struct _objective *, double, double *, struct _constraint *,
    double **, double *, double *, double *, double *,
    double **, double **, double *, double, int);
static void
di1(int, int, int, int, int, int, int, int, int, int, int, int, int *,
    int, struct _parameter *, double *, struct _objective *,
    double, double *, struct _constraint *, double *,
    double *, double *, double *, double **, double *, double);
#else
static void dqp();
static void di1();
#endif

#ifdef __STDC__
static void
dir(int nparam, int nob, int nobL, int nfsip, int nineqn, int neq, int neqn,
    int nn, int ncsipl, int ncsipn, int ncnstr,
    int feasb, double *steps, double epskt, double epseqn,
    double *sktnom, double *scvneq, double Ck, double *d0nm,
    double *grdftd, int *indxob, int *indxcn, int *iact, int *iskp,
    int *iskip, int *istore, struct _parameter *param, double *di,
    double *d, struct _constraint *cs, struct _objective *ob,
    double *fM, double *fMp, double *fmax, double *psf, double *grdpsf,
    double *penp, double **a, double *bl, double *bu, double *clamda,
    double *cvec, double **hess, double **hess1,
    double *backup, double *signeq, struct _violation *viol,
    void(*obj)(int, int, double *, double *, void *),
    void(*constr)(int, int, double *, double *, void *))
#else
static void
dir(nparam, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn, ncnstr,
    feasb, steps, epskt, epseqn, sktnom, scvneq, Ck, d0nm,
    grdftd, indxob, indxcn, iact, iskp, iskip, istore, param, di, d, cs, ob,
    fM, fMp, fmax, psf, grdpsf, penp, a, bl, bu, clamda, cvec, hess, hess1,
    backup, signeq, viol, obj, constr)
int nparam, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn, ncnstr,
*iskp, feasb;
int *indxob, *indxcn, *iact, *iskip, *istore;
double *steps, epskt, epseqn, *sktnom, Ck, *d0nm, *grdftd, *fM, *fMp,
*fmax, *psf, *scvneq;
double  *di, *d, *grdpsf, *penp, **a, *bl, *bu, *clamda, *cvec, **hess,
**hess1, *backup, *signeq;
struct  _constraint *cs;
struct  _objective  *ob;
struct  _parameter  *param;
struct  _violation  *viol;
void(* obj)(), (* constr)();
#endif
{
    int  i, j, k, kk, ncg, ncf, nqprm0, nclin0, nctot0, infoqp, nqprm1, ncl,
    nclin1, ncc, nff, nrowa0, nrowa1, ninq, nobb, nobbL,
    nncn, ltem1, ltem2, display, need_d1;
    double fmxl, vv, dx, dmx, dnm1, dnm, v0, v1, vk, temp1, temp2, theta,
    rhol, rhog, rho, grdfd0, grdfd1, dummy, grdgd0, grdgd1, thrshd,
    sign, *adummy, dnmtil, *tempv;

    ncg = ncf = *iskp = 0;
    ncl = glob_info.nnineq - nineqn;
    glob_log.local = glob_log.update = FALSE;
    glob_log.rhol_is1 = FALSE;
    thrshd = tolfea;
    adummy = make_dv(1);
    adummy[1] = 0.e0;
    dummy = 0.e0;
    temp1 = temp2 = 0.e0;
    if (glob_prnt.iter % glob_prnt.iter_mod)
        display = FALSE;
    else
        display = TRUE;
    need_d1 = TRUE;

    if (nobL <= 1)
    {
        nqprm0 = nparam;
        nclin0 = ncnstr;
    }
    else
    {
        nqprm0 = nparam + 1;
        nclin0 = ncnstr + nobL;
    }
    nctot0 = nqprm0 + nclin0;
    vv = 0.e0;
    nrowa0 = DMAX1(nclin0, 1);
    for (i = 1; i <= ncnstr; i++)
    {
        if (feasb)
        {
            if (i > nineqn && i <= glob_info.nnineq)
                iskip[glob_info.nnineq+2-i] = i;
            iw[i] = i;
        }
        else
        {
            if (i <= ncl)
                iskip[ncl+2-i] = nineqn + i;
            if (i <= ncl)
                iw[i] = nineqn + i;
            if (i > ncl)
                iw[i] = nineqn + neqn + i;
        }
    }
    for (i = 1; i <= nob; i++)
        iw[ncnstr+i] = i;
    nullvc(nparam, cvec);
    glob_log.d0_is0 = FALSE;
    dqp(nparam, nqprm0, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn,
        ncnstr, nctot0, nrowa0, nineqn, &infoqp, param, di, feasb, ob,
        *fmax, grdpsf, cs, a, cvec, bl, bu, clamda, hess, hess1, di, vv, 0);
    if (infoqp != 0)
    {
        glob_prnt.info = 5;
        if (!feasb)
            glob_prnt.info = 2;
        nstop = 0;
        free_dv(adummy);
        return;
    }
    /*-------------------------------------------------------------*/
    /*    Reorder indexes of constraints & objectives              */
    /*-------------------------------------------------------------*/
    if (nn > 1)
    {
        j = 1;
        k = nn;
        for (i = nn; i >= 1; i--)
        {
            if (fuscmp(cs[indxcn[i]].mult, thrshd))
            {
                iact[j] = indxcn[i];
                j++;
            }
            else
            {
                iact[k] = indxcn[i];
                k--;
            }
        }
    }
    if (nobL > 1)
    {
        j = nn + 1;
        k = nn + nob;
        for (i = nob; i >= 1; i--)
        {
            kk = nqprm0 + ncnstr;
            ltem1 = fuscmp(ob[i].mult, thrshd);
            ltem2 = (nobL != nob) && (fuscmp(ob[i].mult_L, thrshd));
            if (ltem1 || ltem2)
            {
                iact[j] = i;
                j++;
            }
            else
            {
                iact[k] = i;
                k--;
            }
        }
    }
    if (nob > 0)
        vv = ob[iact[nn+1]].val;
    *d0nm = sqrt(scaprd(nparam, di, di));
    if (glob_log.first && nclin0 == 0)
    {
        dx = sqrt(scaprd(nparam, param->x, param->x));
        dmx = DMAX1(dx, 1.e0);
        if (*d0nm > dmx)
        {
            for (i = 1; i <= nparam; i++)
                di[i] = di[i] * dmx / (*d0nm);
            *d0nm = dmx;
        }
    }
    matrvc(nparam, nparam, hess, di, w);
    if (nn == 0)
        *grdftd = -scaprd(nparam, w, di);
    *sktnom = sqrt(scaprd(nparam, w, w));
    if (((*d0nm <= epskt) || ((gLgeps > 0.e0) && (*sktnom <= gLgeps))) &&
        (neqn == 0 || *scvneq <= epseqn))
    {
        /* We are finished! */
        nstop = 0;
        if (feasb && glob_log.first && neqn != 0)
        {
            /*  Finished, but still need to estimate nonlinear equality
                constraint multipliers   */
            glob_log.get_ne_mult = TRUE;
            glob_log.d0_is0 = TRUE;
        }
        if (!feasb)
            glob_prnt.info = 2;
        free_dv(adummy);
        if (glob_prnt.iprint < 3 || !display)
            return;
        if (nobL <= 1)
            nff = 1;
        if (nobL > 1)
            nff = 2;
        sbout1(glob_prnt.io, nparam, "multipliers for x   ", dummy,
               param->mult, 2, 2);
        if (ncnstr != 0)
        {
            fprintf(glob_prnt.io, "\t\t\t %s\t %22.14e\n",
                    "           for g    ", cs[1].mult);
            for (j = 2; j <= ncnstr; j++)
                fprintf(glob_prnt.io, " \t\t\t\t\t\t %22.14e\n", cs[j].mult);
        }
        if (nobL > 1)
        {
            fprintf(glob_prnt.io, "\t\t\t %s\t %22.14e\n",
                    "           for f    ", ob[1].mult);
            for (j = 2; j <= nob; j++)
                fprintf(glob_prnt.io, " \t\t\t\t\t\t %22.14e\n", ob[j].mult);
        }
        return;
    }
    if (glob_prnt.iprint >= 3 && display)
    {
        sbout1(glob_prnt.io, nparam, "d0                  ", dummy, di, 2, 2);
        sbout1(glob_prnt.io, 0, "d0norm              ", *d0nm, adummy, 1, 2);
        sbout1(glob_prnt.io, 0, "ktnorm              ", *sktnom, adummy, 1, 2);
    }
    if (neqn != 0 && *d0nm <= DMIN1(0.5e0*epskt, (0.1e-1)*glob_grd.rteps)
        && *scvneq > epseqn)
    {
        /* d0 is "zero", but equality constraints not satisfied  */
        glob_log.d0_is0 = TRUE;
        return;
    }
    /*--------------------------------------------------------------*/
    /*     Single objective without nonlinear constraints requires  */
    /*     no d1 and dtilde; multi-objectives without nonlinear     */
    /*     constraints requires no d1.       */
    /*--------------------------------------------------------------*/
    if (nn != 0)
        *grdftd = slope(nob, nobL, neqn, nparam, feasb, ob, grdpsf,
                        di, d, *fmax, dummy, 0, adummy, 0);

    if (nn == 0 && nobL <= 1)
    {
        for (i = 1; i <= nparam; i++)
            d[i] = 0.e0;
        dnmtil = 0.e0;
        free_dv(adummy);
        return;
    }
    if (nn == 0)
    {
        dnm = *d0nm;
        rho = 0.e0;
        rhog = 0.e0;
        goto L310;
    }
    /*-------------------------------------------------------------*/
    /*     compute modified first order direction d1    */
    /*-------------------------------------------------------------*/

    /* First check that it is necessary */
    if (glob_info.mode == 1)
    {
        vk = DMIN1(Ck * (*d0nm) * (*d0nm), *d0nm);
        need_d1 = FALSE;
        for (i = 1; i <= nn; i++)
        {
            tempv = cs[indxcn[i]].grad;
            grdgd0 = scaprd(nparam, tempv, di);
            temp1 = vk + cs[indxcn[i]].val + grdgd0;
            if (temp1 > 0.e0)
            {
                need_d1 = TRUE;
                break;
            }
        }
    }
    if (need_d1)
    {
        nqprm1 = nparam + 1;
        if (glob_info.mode == 0)
            nclin1 = ncnstr + DMAX1(nobL, 1);
        if (glob_info.mode == 1)
            nclin1 = ncnstr;
        nrowa1 = DMAX1(nclin1, 1);
        ninq = glob_info.nnineq;
        di1(nparam, nqprm1, nob, nobL, nfsip, nineqn, neq, neqn, ncnstr,
            ncsipl, ncsipn, nrowa1, &infoqp, glob_info.mode,
            param, di, ob, *fmax, grdpsf, cs, cvec, bl, bu, clamda,
            hess1, d, *steps);
        if (infoqp != 0)
        {
            glob_prnt.info = 6;
            if (!feasb)
                glob_prnt.info = 2;
            nstop = 0;
            free_dv(adummy);
            return;
        }
        dnm1 = sqrt(scaprd(nparam, d, d));
        if (glob_prnt.iprint >= 3 && display)
        {
            sbout1(glob_prnt.io, nparam, "d1                  ", dummy, d, 2, 2);
            sbout1(glob_prnt.io, 0, "d1norm              ", dnm1, adummy, 1, 2);
        }
    }
    else
    {
        dnm1 = 0.e0;
        for (i = 1; i <= nparam; i++)
            d[i] = 0.e0;
    }
    if (glob_info.mode != 1)
    {
        v0 = pow(*d0nm, 2.1);
        v1 = DMAX1(0.5e0, pow(dnm1, 2.5));
        rho = v0 / (v0 + v1);
        rhog = rho;
    }
    else
    {
        rhol = 0.e0;
        if (need_d1)
        {
            for (i = 1; i <= nn; i++)
            {
                tempv = cs[indxcn[i]].grad;
                grdgd0 = scaprd(nparam, tempv, di);
                grdgd1 = scaprd(nparam, tempv, d);
                temp1 = vk + cs[indxcn[i]].val + grdgd0;
                temp2 = grdgd1 - grdgd0;
                if (temp1 <= 0.e0)
                    continue;
                if (fabs(temp2) < glob_grd.epsmac)
                {
                    rhol = 1.e0;
                    glob_log.rhol_is1 = TRUE;
                    break;
                }
                rhol = DMAX1(rhol, -temp1 / temp2);
                if (temp2 < 0.e0 && rhol < 1.e0)
                    continue;
                rhol = 1.e0;
                glob_log.rhol_is1 = TRUE;
                break;
            }
        }
        theta = 0.2e0;
        if (rhol == 0.e0)
        {
            rhog = rho = 0.e0;
            dnm = *d0nm;
            goto L310;
        }
        if (nobL > 1)
        {
            rhog = slope(nob, nobL, neqn, nparam, feasb, ob, grdpsf,
                         di, d, *fmax, theta, glob_info.mode, adummy, 0);
            rhog = DMIN1(rhol, rhog);
        }
        else
        {
            grdfd0 = *grdftd;
            if (nob == 1)
                grdfd1 = scaprd(nparam, ob[1].grad, d);
            else
                grdfd1 = 0.e0;
            grdfd1 = grdfd1 - scaprd(nparam, grdpsf, d);
            temp1 = grdfd1 - grdfd0;
            temp2 = (theta - 1.e0) * grdfd0 / temp1;
            if (temp1 <= 0.e0)
                rhog = rhol;
            else
                rhog = DMIN1(rhol, temp2);
        }
        rho = rhog;
        if (*steps == 1.e0 && rhol < 0.5e0)
            rho = rhol;
    }
    for (i = 1; i <= nparam; i++)
    {
        if (rho != rhog)
            cvec[i] = di[i];
        di[i] = (1.e0 - rho) * di[i] + rho * d[i];
    }
    dnm = sqrt(scaprd(nparam, di, di));
    if (!(glob_prnt.iprint < 3 || glob_info.mode == 1 || nn == 0) && display)
    {
        sbout1(glob_prnt.io, 0, "rho                 ", rho, adummy, 1, 2);
        sbout1(glob_prnt.io, nparam, "d                   ", dummy, di, 2, 2);
        sbout1(glob_prnt.io, 0, "dnorm               ", dnm, adummy, 1, 2);
    }
L310:
    for (i = 1; i <= nob; i++)
        bl[i] = ob[i].val;
    if (rho != 1.e0)
    {
        if (!(glob_prnt.iprint != 3 || glob_info.mode == 0 || nn == 0)
            && display)
        {
            sbout1(glob_prnt.io, 0, "Ck                  ", Ck, adummy, 1, 2);
            sbout1(glob_prnt.io, 0, "rhol                ", rho, adummy, 1, 2);
            sbout1(glob_prnt.io, nparam, "dl                  ", dummy, di, 2, 2);
            sbout1(glob_prnt.io, 0, "dlnorm              ", dnm, adummy, 1, 2);
        }
        if (glob_info.mode != 0)
        {
            glob_log.local = TRUE;
            step1(nparam, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn,
                  ncnstr, &ncg, &ncf, indxob, indxcn, iact, iskp, iskip, istore,
                  feasb, *grdftd, ob, fM, fMp, fmax, psf, penp, steps, scvneq, bu,
                  param->x, di, d, cs, backup, signeq, viol, obj, constr, param->cd);
            if (!glob_log.update)
                nstop = 1;
            else
            {
                free_dv(adummy);
                return;
            }
            glob_log.local = FALSE;
            if (rho != rhog && nn != 0)
                for (i = 1; i <= nparam; i++)
                    di[i] = (1 - rhog) * cvec[i] + rhog * d[i];
            dnm = sqrt(scaprd(nparam, di, di));
        }
    }
    if (!(glob_prnt.iprint < 3 || glob_info.mode == 0 || nn == 0) &&
        display)
    {
        sbout1(glob_prnt.io, 0, "rhog                ", rhog, adummy, 1, 2);
        sbout1(glob_prnt.io, nparam, "dg                  ", dummy, di, 2, 2);
        sbout1(glob_prnt.io, 0, "dgnorm              ", dnm, adummy, 1, 2);
    }
    if (rho != 0.e0)
        *grdftd = slope(nob, nobL, neqn, nparam, feasb, ob,
                        grdpsf, di, d, *fmax, theta, 0, bl, 1);
    if (glob_info.mode != 1 || rho != rhog)
        for (i = 1; i <= nparam; i++)
            bu[i] = param->x[i] + di[i];
    x_is_new = TRUE;
    if (rho != rhog)
        ncg = 0;
    ncc = ncg + 1;
    fmxl = -bgbnd;
    ninq = nncn = ncg;
    j = 0;
    /*--------------------------------------------------------------*/
    /*   iskip[1]-iskip[iskp] store the indexes of linear inequality*/
    /*   constraints that are not to be used to compute d~          */
    /*   iskip[nnineq-nineqn+1]-iskip[nnineq-ncn+1-iskp] store      */
    /*   those that are to be used to compute d~                    */
    /*--------------------------------------------------------------*/
    for (i = ncc; i <= ncnstr; i++)
    {
        if (i <= nn)
            kk = iact[i];
        else
            kk = indxcn[i];
        if (kk > nineqn && kk <= glob_info.nnineq)
        {
            iskip[ncl+1-j] = kk;
            j++;
        }
        if (kk <= glob_info.nnineq)
        {
            tempv = cs[kk].grad;
            temp1 = dnm * sqrt(scaprd(nparam, tempv, tempv));
            temp2 = cs[kk].mult;
        }
        if (temp2 != 0.e0 || cs[kk].val >= (-0.2e0*temp1) ||
            kk > glob_info.nnineq)
        {
            ninq++;
            iw[ninq] = kk;
            if (feasb && kk <= nineqn)
                istore[kk] = 1;
            constr(nparam, kk, bu + 1, &(cs[kk].val), param->cd);
            if (!feasb || (feasb && (kk > glob_info.nnineq + neqn)))
                continue;
            if (kk <= nineqn)
                nncn = ninq;
            fmxl = DMAX1(fmxl, cs[kk].val);
            if (feasb && (kk <= nineqn || (kk > glob_info.nnineq
                                           && kk <= (glob_info.nnineq + neqn))))
                glob_info.ncallg++;
            if (fabs(fmxl) > bgbnd)
            {
                for (i = 1; i <= nparam; i++)
                    d[i] = 0.e0;
                dnmtil = 0.e0;
                nstop = 1;
                free_dv(adummy);
                return;
            }
            continue;
        }
        if (kk <= nineqn)
            continue;
        (*iskp)++;
        iskip[*iskp] = kk;
        j--;
    }
    if ((neqn != 0) && (feasb))
        resign(nparam, neqn, psf, grdpsf, penp, cs, signeq, 10, 20);
    ninq -= neq;
    /*  if (!feasb) ninq+=neqn;   BUG???   */
    if (ncg != 0)
        for (i = 1; i <= ncg; i++)
        {
            iw[i] = iact[i];
            istore[iact[i]] = 1;
            fmxl = DMAX1(fmxl, cs[iact[i]].val);
            if (fabs(fmxl) > bgbnd)
            {
                for (i = 1; i <= nparam; i++)
                    d[i] = 0.e0;
                dnmtil = 0.e0;
                nstop = 1;
                free_dv(adummy);
                return;
            }
        }
    if (nobL <= 1)
    {
        iw[1+ninq+neq] = 1;
        nobb = nob;
        goto L1110;
    }
    if (rho != rhog)
        ncf = 0;
    nff = ncf + 1;
    nobb = ncf;
    sign = 1.e0;
    fmxl = -bgbnd;
    if (ob[iact[nn+1]].mult < 0.e0)
        sign = -1.e0;
    for (i = nff; i <= nob; i++)
    {
        kk = iact[nn+i];
        if (!feasb)
            kk = iact[i];
        if (feasb)
            k = nn + 1;
        if (!feasb)
            k = 1;
        for (j = 1; j <= nparam; j++)
            w[nparam+j] = sign * ob[iact[k]].grad[j] - ob[kk].grad[j];
        temp1 = fabs(ob[kk].val - sign * vv);
        temp2 = dnm * sqrt(scaprd(nparam, &w[nparam], &w[nparam]));
        if (temp1 != 0.e0 && temp2 != 0.e0)
        {
            temp1 = temp1 / temp2;
            temp2 = ob[kk].mult;
            if (temp2 == 0.e0 && temp1 > 0.2e0)
                continue;
        }
        nobb++;
        iw[nobb+ninq+neq] = kk;
        if (feasb)
            istore[nineqn+kk] = 1;
        else
            istore[kk] = 1;
        if (!feasb)
        {
            constr(nparam, indxob[kk], bu + 1, &(ob[kk].val), param->cd);
            glob_info.ncallg++;
        }
        else
        {
            obj(nparam, kk, bu + 1, &(ob[kk].val), param->cd);
            glob_info.ncallf++;
            if (nobL != nob)
                fmxl = DMAX1(fmxl, -ob[kk].val);
        }
        fmxl = DMAX1(fmxl, ob[kk].val);
        if (fabs(fmxl) > bgbnd)
        {
            for (i = 1; i <= nparam; i++)
                d[i] = 0.e0;
            dnmtil = 0.e0;
            nstop = 1;
            free_dv(adummy);
            return;
        }
    }
    if (ncf != 0)
    {
        for (i = 1; i <= ncf; i++)
        {
            iw[ninq+neq+i] = iact[i+nn];
            istore[nineqn+iact[i+nn]] = 1;
            fmxl = DMAX1(fmxl, ob[iact[i+nn]].val);
            if (nobL != nob)
                fmxl = DMAX1(fmxl, -ob[iact[i+nn]].val);
            if (fabs(fmxl) > bgbnd)
            {
                for (i = 1; i <= nparam; i++)
                    d[i] = 0.e0;
                dnmtil = 0.e0;
                nstop = 1;
                free_dv(adummy);
                return;
            }
        }
    }
    L1110:
    matrvc(nparam, nparam, hess, di, cvec);
    vv = -DMIN1(0.01e0 * dnm, pow(dnm, 2.5));
    /*--------------------------------------------------------------*/
    /*    compute a correction dtilde to d=(1-rho)d0+rho*d1         */
    /*--------------------------------------------------------------*/
    if (nobL != nob)
        nobbL = 2 * nobb;
    if (nobL == nob)
        nobbL = nobb;
    if (nobbL <= 1)
    {
        nqprm0 = nparam;
        nclin0 = ninq + neq;
    }
    else
    {
        nqprm0 = nparam + 1;
        nclin0 = ninq + neq + nobbL;
    }
    nctot0 = nqprm0 + nclin0;
    nrowa0 = DMAX1(nclin0, 1);
    i = ninq + neq;
    nstop = 1;
    dqp(nparam, nqprm0, nobb, nobbL, nfsip, nncn, neq, neqn, nn, ncsipl, ncsipn, i,
        nctot0, nrowa0, nineqn, &infoqp, param, di, feasb, ob, fmxl,
        grdpsf, cs, a, cvec, bl, bu, clamda, hess, hess1, d, vv, 1);
    dnmtil = sqrt(scaprd(nparam, d, d));
    if (infoqp != 0 || dnmtil > dnm)
    {
        for (i = 1; i <= nparam; i++)
            d[i] = 0.e0;
        dnmtil = 0.e0;
        nstop = 1;
        free_dv(adummy);
        return;
    }
    if (dnmtil != 0.e0)
        for (i = 1; i <= nineqn + nob; i++)
            istore[i] = 0;
    if (glob_prnt.iprint < 3 || !display)
    {
        free_dv(adummy);
        return;
    }
    sbout1(glob_prnt.io, nparam, "dtilde              ", dummy, d, 2, 2);
    sbout1(glob_prnt.io, 0, "dtnorm              ", dnmtil, adummy, 1, 2);
    free_dv(adummy);
    return;
}

/*******************************************************************/
/*     job=0:     compute d0        */
/*     job=1:     compute d~        */
/*******************************************************************/
#ifdef __STDC__
static void
dqp(int nparam, int nqpram, int nob, int nobL, int nfsip, int nineqn,
    int neq, int neqn, int nn, int ncsipl, int ncsipn, int ncnstr,
    int nctotl, int nrowa, int nineqn_tot, int *infoqp,
    struct _parameter *param, double *di, int feasb, struct _objective *ob,
    double fmax, double *grdpsf, struct _constraint *cs,
    double **a, double *cvec, double *bl, double *bu, double *clamda,
    double **hess, double **hess1, double *x,
    double vv, int job)
#else
static void
dqp(nparam, nqpram, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn,
    ncnstr, nctotl, nrowa, nineqn_tot, infoqp, param, di, feasb, ob,
    fmax, grdpsf, cs, a, cvec, bl, bu, clamda, hess, hess1, x, vv, job)
int nparam, nqpram, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn,
ncnstr, nctotl, nrowa, nineqn_tot, *infoqp, job, feasb;
double fmax, vv;
double *di, *grdpsf, **a, *cvec, *bl, *bu, *clamda, **hess, **hess1, *x;
struct  _constraint *cs;
struct  _objective  *ob;
struct  _parameter  *param;
#endif
{
    int i, ii, j, jj, ij, k, iout, mnn, nqnp, zero, temp1, temp2, ncnstr_used,
    numf_used;
    int *iw_hold;
    double x0i, xdi, *bj, *htemp, *atemp;

    iout = 6;
    bj = make_dv(nrowa);
    iw_hold = make_iv(nrowa);
    for (i = 1; i <= nparam; i++)
    {
        x0i = param->x[i];
        if (job == 1)
            xdi = di[i];
        if (job == 0)
            xdi = 0.e0;
        bl[i] = param->bl[i] - x0i - xdi;
        bu[i] = param->bu[i] - x0i - xdi;
        cvec[i] = cvec[i] - grdpsf[i];
    }
    if (nobL > 1)
    {
        bl[nqpram] = -bgbnd;
        bu[nqpram] = bgbnd;
    }
    ii = ncnstr - nn;
    /*---------------------------------------------------------------*/
    /*     constraints are assigned to a in reverse order            */
    /*---------------------------------------------------------------*/
    k = 0;
    for (i = 1; i <= ncnstr; i++)
    {
        jj = iw[ncnstr+1-i];
        if ((jj > glob_info.nnineq) || (jj <= nineqn_tot - ncsipn) ||
            ((jj > nineqn_tot) && (jj <= glob_info.nnineq - ncsipl)) ||
            cs[jj].act_sip)
        {
            k++;
            x0i = vv;
            if (i <= (neq - neqn) || (i > neq && i <= (ncnstr - nineqn)))
                x0i = 0.e0;
            if (!feasb)
                x0i = 0.e0;
            bj[k] = x0i - cs[jj].val;
            for (j = 1; j <= nparam; j++)
                a[k][j] = -cs[jj].grad[j];
            if (nobL > 1)
                a[k][nqpram] = 0.e0;
            iw_hold[k] = jj;
        }
    }
    ncnstr_used = k;
    /*---------------------------------------------------------------*/
    /* Assign objectives for QP         */
    /*---------------------------------------------------------------*/
    if (nobL == 1)
    {
        for (i = 1; i <= nparam; i++)
            cvec[i] = cvec[i] + ob[1].grad[i];
    }
    else if (nobL > 1)
    {
        numf_used = nob - nfsip + glob_info.tot_actf_sip;
        if (job && nfsip)
        {     /* compute # objectives used for dtilde */
            numf_used = 0;
            for (i = 1; i <= nob; i++)
                if (ob[iw[ncnstr+i]].act_sip)
                    numf_used++;
        }
        for (i = 1; i <= nob; i++)
        {
            ij = ncnstr + i;
            if ((i <= nob - nfsip) || ob[iw[ij]].act_sip)
            {
                k++;
                iw_hold[k] = iw[ij];  /* record which are used */
                bj[k] = fmax - ob[iw[ij]].val;
                if (nobL > nob)
                    bj[k+numf_used] = fmax + ob[iw[ij]].val;
                for (j = 1; j <= nparam; j++)
                {
                    a[k][j] = -ob[iw[ij]].grad[j];
                    if (nobL > nob)
                        a[k+numf_used][j] = ob[iw[ij]].grad[j];
                }
                a[k][nqpram] = 1.e0;
                if (nobL > nob)
                    a[k+numf_used][nqpram] = 1.e0;
            }
        }
        cvec[nqpram] = 1.e0;
        if (nobL > nob)
            k = k + numf_used;  /* k=# rows for a         */
    }        /*  =# constraints for QP */
    matrcp(nparam, hess, nparam + 1, hess1);
    nullvc(nqpram, x);

    iw[1] = 1;
    zero = 0;
    temp1 = neq - neqn;
    temp2 = nparam + 1;
    mnn = k + 2 * nqpram;
    htemp = convert(hess1, nparam + 1, nparam + 1);
    atemp = convert(a, nrowa, nqpram);

    ql0001_(&k, &temp1, &nrowa, &nqpram, &temp2, &mnn, (htemp + 1),
            (cvec + 1), (atemp + 1), (bj + 1), (bl + 1), (bu + 1), (x + 1),
            (clamda + 1), &iout, infoqp, &zero, (w + 1), &lenw, (iw + 1), &leniw,
            &glob_grd.epsmac);

    free_dv(htemp);
    free_dv(atemp);
    if (*infoqp != 0 || job == 1)
    {
        free_iv(iw_hold);
        free_dv(bj);
        return;
    }

    /*---------------------------------------------------------------*/
    /*  Save multipliers from the computation of d0                  */
    /*---------------------------------------------------------------*/
    nullvc(nqpram, param->mult);
    if (ncsipl + ncsipn)
        for (i = 1; i <= ncnstr; i++)
            cs[i].mult = 0.e0;
    if (nfsip)
        for (i = 1; i <= nob; i++)
        {
            ob[i].mult = 0.e0;
            ob[i].mult_L = 0.e0;
        }
    for (i = 1; i <= nqpram; i++)
    {
        ii = k + i;
        if (clamda[ii] == 0.e0 && clamda[ii+nqpram] == 0.e0)
            continue;
        else if (clamda[ii] != 0.e0)
            clamda[ii] = -clamda[ii];
        else
            clamda[ii] = clamda[ii+nqpram];
    }
    nqnp = nqpram + ncnstr;
    for (i = 1; i <= nqpram; i++)                /* Simple bounds */
        param->mult[i] = clamda[k+i];
    if (nctotl > nqnp)
    {                         /* Objectives    */
        for (i = 1; i <= numf_used; i++)
        {
            ij = ncnstr_used + i;
            if (nobL != nob)
            {
                ii = k - 2 * numf_used + i;
                ob[iw_hold[ij]].mult = clamda[ii] - clamda[ii+numf_used];
                ob[iw_hold[ij]].mult_L = clamda[ii+numf_used];
            }
            else
            {
                ii = k - numf_used + i;
                ob[iw_hold[ij]].mult = clamda[ii];
            }
        }
    }
    for (i = 1; i <= ncnstr_used; i++)             /* Constraints   */
        cs[iw_hold[i]].mult = clamda[i];
    free_iv(iw_hold);
    free_dv(bj);
    return;
}

/****************************************************************/
/*    Computation of first order direction d1   */
/****************************************************************/
#ifdef __STDC__
static void
di1(int nparam, int nqpram, int nob, int nobL, int nfsip, int nineqn,
    int neq, int neqn, int ncnstr, int ncsipl, int ncsipn,
    int nrowa, int *infoqp, int mode, struct _parameter *param,
    double *d0, struct _objective *ob, double fmax, double
    *grdpsf, struct _constraint *cs, double *cvec, double *bl, double *bu,
    double *clamda, double **hess1, double *x, double steps)
#else
static void
di1(nparam, nqpram, nob, nobL, nfsip, nineqn, neq, neqn, ncnstr, ncsipl,
    ncsipn, nrowa, infoqp, mode, param, d0, ob, fmax, grdpsf, cs,
    cvec, bl, bu, clamda, hess1, x, steps)
int nparam, nqpram, nob, nobL, nfsip, nineqn, neq, neqn, ncnstr,
nrowa, *infoqp, mode, ncsipl, ncsipn;
double fmax, steps;
double *d0, *grdpsf, *cvec, *bl, *bu, *clamda, **hess1, *x;
struct _constraint *cs;
struct _objective  *ob;
struct _parameter  *param;
#endif
{
    int i, k, ii, jj, iout, j, mnn, zero, temp1, temp3, ncnstr_used, numf_used;
    int *iw_hold;
    double x0i, eta, *atemp, *htemp, **a, *bj;

    if ((ncsipl + ncsipn) != 0)
        nrowa = nrowa - (ncsipl + ncsipn) + glob_info.tot_actg_sip;
    if (nfsip)
    {
        if (nobL > nob)
            nrowa = nrowa - 2 * nfsip + 2 * glob_info.tot_actf_sip;
        else
            nrowa = nrowa - nfsip + glob_info.tot_actf_sip;
    }
    nrowa = DMAX1(nrowa, 1);
    a = make_dm(nrowa, nqpram);
    bj = make_dv(nrowa);
    iw_hold = make_iv(nrowa);
    iout = 6;
    if (mode == 0)
        eta = 0.1e0;
    if (mode == 1)
        eta = 3.e0;
    for (i = 1; i <= nparam; i++)
    {
        x0i = param->x[i];
        bl[i] = param->bl[i] - x0i;
        bu[i] = param->bu[i] - x0i;
        if (mode == 0)
            cvec[i] = -eta * d0[i];
        if (mode == 1)
            cvec[i] = 0.e0;
    }
    bl[nqpram] = -bgbnd;
    bu[nqpram] = bgbnd;
    cvec[nqpram] = 1.e0;
    ii = ncnstr - nineqn;
    k = 0;
    for (i = 1; i <= ncnstr; i++)
    {
        jj = ncnstr + 1 - i;
        if ((jj > glob_info.nnineq) || (jj <= nineqn - ncsipn) ||
            ((jj > nineqn) && (jj <= glob_info.nnineq - ncsipl)) ||
            cs[jj].act_sip)
        {
            k++;
            bj[k] = -cs[jj].val;
            for (j = 1; j <= nparam; j++)
                a[k][j] = -cs[jj].grad[j];
            a[k][nqpram] = 0.e0;
            if ((i > (neq - neqn) && i <= neq) || i > ii)
                a[k][nqpram] = 1.e0;
            iw_hold[k] = jj;
        }
    }
    ncnstr_used = k;

    if (mode != 1)
    {
        numf_used = nob - nfsip + glob_info.tot_actf_sip;
        for (i = 1; i <= nob; i++)
        {
            if ((i <= nob - nfsip) || ob[i].act_sip)
            {
                k++;
                bj[k] = fmax - ob[i].val;
                for (j = 1; j <= nparam; j++)
                {
                    a[k][j] = -ob[i].grad[j] + grdpsf[j];
                    if (nobL > nob)
                        a[k+numf_used][j] = ob[i].grad[j] + grdpsf[j];
                }
                a[k][nqpram] = 1.e0;
                if (nobL > nob)
                    a[k+numf_used][nqpram] = 1.e0;
            }
        }
        if (nob == 0)
        {
            k++;
            bj[k] = fmax;
            for (j = 1; j <= nparam; j++)
                a[k][j] = grdpsf[j];
            a[k][nqpram] = 1.e0;
        }
    }
    diagnl(nqpram, eta, hess1);
    nullvc(nqpram, x);
    hess1[nqpram][nqpram] = 0.e0;

    iw[1] = 1;
    zero = 0;
    temp1 = neq - neqn;
    if (nobL > nob)
        temp3 = k + numf_used;
    else
        temp3 = k;
    mnn = temp3 + 2 * nqpram;
    htemp = convert(hess1, nparam + 1, nparam + 1);
    atemp = convert(a, nrowa, nqpram);

    ql0001_(&temp3, &temp1, &nrowa, &nqpram, &nqpram, &mnn, (htemp + 1),
            (cvec + 1), (atemp + 1), (bj + 1), (bl + 1), (bu + 1), (x + 1),
            (clamda + 1), &iout, infoqp, &zero, (w + 1), &lenw, (iw + 1), &leniw,
            &glob_grd.epsmac);

    free_dv(htemp);
    free_dv(atemp);
    free_dm(a, nrowa);
    free_dv(bj);
    /* Determine binding constraints */
    if (ncsipl + ncsipn)
    {
        for (i = 1; i <= ncnstr_used; i++)
            if (clamda[i] > 0.e0)
                cs[iw_hold[i]].d1bind = TRUE;
    }
    free_iv(iw_hold);
    return;
}
/*****************************************************************/
/*     CFSQP : Armijo or nonmonotone line search, with some      */
/*             ad hoc strategies to decrease the number of       */
/*             function evaluations as much as possible.         */
/*****************************************************************/


#ifdef __STDC__
static void
step1(int nparam, int nob, int nobL, int nfsip, int nineqn, int neq, int neqn,
      int nn, int ncsipl, int ncsipn, int ncnstr, int *ncg, int *ncf,
      int *indxob, int *indxcn, int *iact, int *iskp, int *iskip,
      int *istore, int feasb, double grdftd, struct _objective *ob,
      double *fM, double *fMp, double *fmax, double *psf, double *penp,
      double *steps, double *scvneq, double *xnew, double *x, double *di,
      double *d, struct _constraint *cs, double *backup, double *signeq,
      struct _violation *sip_viol,
      void(*obj)(int, int, double *, double *, void *),
      void(*constr)(int, int, double *, double *, void *), void *cd)
#else
static void
step1(nparam, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn, ncnstr,
      ncg, ncf, indxob, indxcn, iact, iskp, iskip, istore, feasb, grdftd, ob,
      fM, fMp, fmax, psf, penp, steps, scvneq, xnew, x, di, d, cs, backup,
      signeq, sip_viol, obj, constr, cd)
int  nparam, nob, nobL, nfsip, nineqn, neq, neqn, nn, ncsipl, ncsipn, ncnstr,
*ncg, *ncf, feasb, *iskp;
int *indxob, *indxcn, *iact, *iskip, *istore;
double  grdftd, *fM, *fMp, *fmax, *steps, *scvneq, *psf;
double  *xnew, *x, *di, *d, *penp, *backup, *signeq;
struct  _constraint *cs;
struct  _objective  *ob;
struct  _violation  *sip_viol;
void(* obj)(), (* constr)();
void    *cd;
#endif
{
    int i, ii, ij, jj, itry, ikeep, j, job, nlin, mnm, ltem1, ltem2, reform,
    fbind, cdone, fdone, eqdone, display, sipldone;
    double prod1, prod, dummy, fmax1, tolfe, ostep, temp, **adummy, fii;

    nlin = glob_info.nnineq - nineqn;
    itry = ii = jj = 1;
    ostep = *steps = 1.e0;
    fbind = cdone = fdone = eqdone = FALSE;
    dummy = 0.e0;
    sipldone = (ncsipl == 0);
    if (glob_log.local)
        glob_log.dlfeas = FALSE;
    ikeep = nlin - *iskp;
    prod1 = (0.1e0) * grdftd;        /* alpha = 0.1e0  */
    tolfe = 0.e0;                  /* feasibility tolerance */
    adummy = make_dm(1, 1);
    adummy[1][1] = 0.e0;
    if (glob_prnt.iter % glob_prnt.iter_mod)
        display = FALSE;
    else
        display = TRUE;
    if (glob_prnt.iprint >= 3 && display)
        sbout1(glob_prnt.io, 0, "directional deriv   ", grdftd, *(adummy + 1),
               1, 2);
    w[1] = *fM;
    for (;;)
    {
        reform = TRUE;
        if (glob_prnt.iprint >= 3 && display)
            fprintf(glob_prnt.io, "\t\t\t trial number            %22d\n",
                    itry);
        prod = prod1 * (*steps);
        if (!feasb || (nobL > 1))
            prod = prod + tolfe;
        for (i = 1; i <= nparam; i++)
        {
            if (glob_log.local)
                xnew[i] = x[i] + (*steps) * di[i];
            else
                xnew[i] = x[i] + (*steps) * di[i] + d[i] * (*steps) * (*steps);
        }
        x_is_new = TRUE;
        if (glob_prnt.iprint >= 3 && display)
        {
            sbout1(glob_prnt.io, 0, "trial step          ", *steps,
                   *(adummy + 1), 1, 2);
            sbout1(glob_prnt.io, nparam, "trial point         ",
                   dummy, xnew, 2, 2);
        }

        /* Generate an upper bound step size using the linear constraints
           not used in the computation of dtilde */
        if (*iskp != 0)
        {
            ostep = *steps;
            for (i = ii; i <= *iskp; i++)
            {
                ij = iskip[i];
                constr(nparam, ij, xnew + 1, &(cs[ij].val), cd);
                if (glob_prnt.iprint >= 3 && display)
                {
                    if (i == 1)
                        fprintf(glob_prnt.io,
                                "\t\t\t trial constraints  %d \t %22.14e\n", ij,
                                cs[ij].val);
                    if (i != 1)
                        fprintf(glob_prnt.io,
                                "\t\t\t\t\t %d \t %22.14e\n", ij, cs[ij].val);
                }
                if (cs[ij].val <= tolfe)
                    continue;
                ii = i;
                if (ncsipl && ii > glob_info.nnineq - ncsipl)
                {
                    sip_viol->type = CONSTR;
                    sip_viol->index = ij;
                }
                else
                {
                    sip_viol->type = NONE; /* non-SIP constraint violated */
                    sip_viol->index = 0;
                }
                goto L1120;
            }
            *iskp = 0;
        }

        /* Refine the upper bound using the linear SI constraints not
           in Omega_k */
        if (!sipldone)
        {
            for (i = jj; i <= ncsipl; i++)
            {
                ij = glob_info.nnineq - ncsipl + i;
                if (cs[ij].act_sip || element(iskip, nlin - ikeep, ij))
                    continue;
                constr(nparam, ij, xnew + 1, &(cs[ij].val), cd);
                if (glob_prnt.iprint >= 3 && display)
                {
                    if (i == 1)
                        fprintf(glob_prnt.io,
                                "\t\t\t trial constraints  %d \t %22.14e\n", ij,
                                cs[ij].val);
                    if (i != 1)
                        fprintf(glob_prnt.io,
                                "\t\t\t\t\t %d \t %22.14e\n", ij, cs[ij].val);
                }
                if (cs[ij].val <= tolfe)
                    continue;
                jj = i;
                sip_viol->type = CONSTR;
                sip_viol->index = ij;
                goto L1120;
            }
            sipldone = TRUE;
        }
        if (nn == 0)
            goto L310;

        /* Check nonlinear constraints                            */
        if (!glob_log.local && fbind)
            goto L315;
        do
        {
            for (i = 1; i <= nn; i++)
            {
                *ncg = i;
                ii = iact[i];
                ij = glob_info.nnineq + neqn;
                if (!((ii <= glob_info.nnineq && istore[ii] == 1) ||
                      (ii > glob_info.nnineq && ii <= ij && eqdone)))
                {
                    temp = 1.e0;
                    if (ii > glob_info.nnineq && ii <= ij)
                        temp = signeq[ii-glob_info.nnineq];
                    constr(nparam, ii, xnew + 1, &(cs[ii].val), cd);
                    cs[ii].val *= temp;
                    glob_info.ncallg++;
                }
                if (glob_prnt.iprint >= 3 && display)
                {
                    if (i == 1 && ikeep == nlin)
                        fprintf(glob_prnt.io,
                                "\t\t\t trial constraints  %d \t %22.14e\n", ii,
                                cs[ii].val);
                    if (i != 1 || ikeep != nlin)
                        fprintf(glob_prnt.io,
                                "\t\t\t\t\t %d \t %22.14e\n", ii, cs[ii].val);
                }
                if (!(glob_log.local || cs[ii].val <= tolfe))
                {
                    shift(nn, ii, iact);
                    if (ncsipn && ii > nineqn - ncsipn)
                    {
                        sip_viol->type = CONSTR;
                        sip_viol->index = ii;
                    }
                    else
                    {
                        sip_viol->type = NONE; /* non-SIP constraint violated */
                        sip_viol->index = 0;
                    }
                    goto L1110;
                }
                if (glob_log.local && cs[ii].val > tolfe)
                {
                    if (ncsipn && ii > nineqn - ncsipn)
                    {
                        sip_viol->type = CONSTR;
                        sip_viol->index = ii;
                    }
                    else
                    {
                        sip_viol->type = NONE; /* non-SIP constraint violated */
                        sip_viol->index = 0;
                    }
                    goto L1500;
                }
            }
L310:
            cdone = eqdone = TRUE;
            if (glob_log.local)
                glob_log.dlfeas = TRUE; /* dl is feasible */
L315:
            if (fdone)
                break;
            if (nob > 0)
                fmax1 = -bgbnd;
            else
                fmax1 = 0.e0;
            for (i = 0; i <= nob; i++)
            {
                if (nob != 0 && i == 0)
                    continue;
                *ncf = i;
                ii = iact[nn+i];
                if (feasb)
                {
                    if (!(eqdone || neqn == 0))
                    {
                        for (j = 1; j <= neqn; j++)
                            constr(nparam, glob_info.nnineq + j, xnew + 1,
                                   &(cs[glob_info.nnineq+j].val), cd);
                        glob_info.ncallg += neqn;
                    }
                    if (neqn != 0)
                    {
                        if (eqdone)
                            job = 20;
                        if (!eqdone)
                            job = 10;
                        resign(nparam, neqn, psf, *(adummy + 1), penp, cs, signeq,
                               job, 10);
                    }
                    if (istore[nineqn+ii] != 1 && i != 0)
                    {
                        obj(nparam, ii, xnew + 1, &(ob[ii].val), cd);
                        glob_info.ncallf++;
                    }
                    if (i == 0)
                        fii = 0.e0;
                    else
                        fii = ob[ii].val;
                    if (i == 0 && glob_prnt.iprint >= 3 && display)
                        fprintf(glob_prnt.io,
                                "\t\t\t trial penalty term \t %22.14e\n", -*psf);
                    if (i == 1 && glob_prnt.iprint >= 3 && display)
                        fprintf(glob_prnt.io,
                                "\t\t\t trial objectives   %d \t %22.14e\n",
                                ii, fii - *psf);
                    if (i > 1 && glob_prnt.iprint >= 3 && display)
                        fprintf(glob_prnt.io,
                                "\t\t\t\t\t %d \t %22.14e\n", ii, fii - *psf);
                }
                else
                {
                    if (istore[ii] != 1)
                    {
                        constr(nparam, indxob[ii], xnew + 1, &(ob[ii].val), cd);
                        glob_info.ncallg++;
                    }
                    if (ob[ii].val > tolfe)
                        reform = FALSE;
                    if (i == 1 && glob_prnt.iprint > 2 && display)
                        fprintf(glob_prnt.io,
                                "\t\t\t trial objectives   %d \t %22.14e\n",
                                indxob[ii], ob[ii].val);
                    if (i != 1 && glob_prnt.iprint > 2 && display)
                        fprintf(glob_prnt.io,
                                "\t\t\t\t\t %d \t %22.14e\n", indxob[ii],
                                ob[ii].val);
                    fii = ob[ii].val;
                }
                fmax1 = DMAX1(fmax1, fii);
                if (nobL != nob)
                    fmax1 = DMAX1(fmax1, -fii);
                if (!feasb && reform)
                    continue;
                if (!glob_log.local)
                {
                    if ((fii - *psf) > (*fMp + prod))
                    {
                        fbind = TRUE;
                        shift(nob, ii, &iact[nn]);
                        if (nfsip && ii > nob - nfsip)
                        {
                            sip_viol->type = OBJECT;
                            sip_viol->index = ii;
                        }
                        else
                        {
                            sip_viol->type = NONE;
                            sip_viol->index = 0;
                        }
                        goto L1110;
                    }
                    if (nobL == nob || (-fii - *psf) <= (*fMp + prod))
                        continue;
                    fbind = TRUE;
                    shift(nob, ii, &iact[nn]);
                    if (nfsip && ii > nob - nfsip)
                    {
                        sip_viol->type = OBJECT;
                        sip_viol->index = ii;
                    }
                    else
                    {
                        sip_viol->type = NONE;
                        sip_viol->index = 0;
                    }
                    goto L1110;
                }
                ltem1 = (fii - *psf) > (*fMp + prod);
                ltem2 = (nobL != nob) && ((-fii - *psf) > (*fMp + prod));
                if (ltem1 || ltem2)
                    goto L1500;
            }
            fbind = FALSE;
            fdone = eqdone = TRUE;
        }
        while (!cdone);
        if (ostep == *steps)
            mnm = ikeep + neq - neqn;
        if (ostep != *steps)
            mnm = ncnstr - nn;
        for (i = 1; i <= mnm; i++)
        {
            ii = indxcn[i+nn];
            if (ikeep != nlin && ostep == *steps)
            {
                if (i <= ikeep)
                    ii = iskip[nlin+2-i];
                else
                    ii = indxcn[nn+i-ikeep+nlin];
            }
            constr(nparam, ii, xnew + 1, &(cs[ii].val), cd);
        }
        *scvneq = 0.e0;
        for (i = 1; i <= ncnstr; i++)
        {
            if (i > glob_info.nnineq && i <= (glob_info.nnineq + neqn))
                *scvneq = *scvneq - cs[i].val;
            backup[i] = cs[i].val;
        }
        for (i = 1; i <= nob; i++)
            backup[i+ncnstr] = ob[i].val;
        if (!feasb && reform)
        {
            for (i = 1; i <= nparam; i++)
                x[i] = xnew[i];
            nstop = 0;
            goto L1500;
        }
        if (glob_log.local)
            *ncg = ncnstr;
        if (glob_log.local)
            glob_log.update = TRUE;
        *fM = fmax1;
        *fMp = fmax1 - *psf;
        *fmax = fmax1;
        for (i = 1; i <= nn; i++)
            iact[i] = indxcn[i];
        for (i = 1; i <= nob; i++)
            iact[nn+i] = i;
        goto L1500;
        L1110:
        cdone = fdone = eqdone = reform = FALSE;
        L1120:
        itry++;
        if (glob_info.modec == 2)
            fbind = FALSE;
        if (*steps >= 1.e0)
            for (i = 1; i <= nob + nineqn; i++)
                istore[i] = 0;
        *steps = *steps * 0.5e0;
        if (*steps < glob_grd.epsmac)
            break;
    }
    glob_prnt.info = 4;
    nstop = 0;
    L1500:
    free_dm(adummy, 1);
    if (*steps < 1.e0)
        return;
    for (i = 1; i <= nob + nineqn; i++)
        istore[i] = 0;
    return;
}
/******************************************************************/
/*   CFSQP : Update the Hessian matrix using BFGS formula with    */
/*           Powell's modification.                               */
/******************************************************************/


#ifdef __STDC__
static void
hessian(int nparam, int nob, int nfsip, int nobL, int nineqn, int neq,
        int neqn, int nn, int ncsipn, int ncnstr, int nfs, int *nstart,
        int feasb, double *xnew, struct _parameter *param,
        struct _objective *ob, double fmax, double *fM, double *fMp,
        double *psf, double *grdpsf, double *penp, struct _constraint *cs,
        double *gm, int *indxob, int *indxcn, double *delta, double *eta,
        double *gamma, double **hess, double *hd, double steps, int *nrst,
        double *signeq, double *span,
        void(*obj)(int, int, double *, double *, void *),
        void(*constr)(int, int, double *, double *, void *),
        void(*gradob)(int, int, double *, double *,
                      void(*)(int, int, double *, double *, void *), void *),
        void(*gradcn)(int, int, double *, double *,
                      void(*)(int, int, double *, double *, void *), void *),
        double **phess, double *psb, double *psmu,
        struct _violation *sip_viol)
#else
static void
hessian(nparam, nob, nfsip, nobL, nineqn, neq, neqn, nn, ncsipn, ncnstr,
        nfs, nstart, feasb, xnew, param, ob, fmax, fM, fMp, psf, grdpsf, penp,
        cs, gm, indxob, indxcn, delta, eta, gamma, hess, hd, steps, nrst, signeq,
        span, obj, constr, gradob, gradcn, phess, psb, psmu, sip_viol)
int nparam, nob, nobL, nineqn, neq, neqn, nn, nfsip, ncsipn, ncnstr,
nfs, *nstart, feasb, *nrst;
int     *indxob, *indxcn;
double steps, *psf, fmax, *fM, *fMp;
double  *xnew, *grdpsf, *penp, *gm, *delta, *eta, *gamma,
**hess, *hd, *signeq, *span, **phess, *psb, *psmu;
struct  _constraint *cs;
struct  _objective  *ob;
struct  _parameter  *param;
struct  _violation  *sip_viol;
void(* obj)(), (* constr)(), (* gradob)(), (* gradcn)();
#endif
{
    int    i, j, k, ifail, np, mnm, done, display;
    double dhd, gammd, etad, dummy, theta, signgj, psfnew, delta_s;
    double *tempv;

    /* Check to see whether user-accessible stopping criterion
       is satisfied. The check of gLgeps is made just after
       computing d0 */

    if (!glob_log.get_ne_mult)
    {
        if (feasb && nstop && !neqn)
            if ((fabs(w[1] - fmax) <= objeps) ||
                (fabs(w[1] - fmax) <= objrep*fabs(w[1])))
                nstop = 0;
        if (!nstop)
        {
            for (i = 1; i <= nparam; i++)
                param->x[i] = xnew[i];
            x_is_new = TRUE;
            return;
        }
    }

    delta_s = glob_grd.rteps;  /* SIP */
    if (glob_prnt.iter % glob_prnt.iter_mod)
        display = FALSE;
    else
        display = TRUE;
    psfnew = 0.e0;
    glob_prnt.ipd = 0;
    done = FALSE;
    dummy = 0.e0;
    nullvc(nparam, delta);
    nullvc(nparam, eta);
    for (j = 1; j <= 2; j++)
    {
        nullvc(nparam, gamma);
        if (nobL > 1)
        {
            for (i = 1; i <= nparam; i++)
            {
                hd[i] = 0.e0;
                for (k = 1; k <= nob; k++)
                    hd[i] = hd[i] + ob[k].grad[i] * ob[k].mult;
            }
        }
        if (feasb)
        {
            if (nineqn != 0)
            {
                for (i = 1; i <= nparam; i++)
                {
                    gamma[i] = 0.e0;
                    for (k = 1; k <= nineqn; k++)
                        gamma[i] = gamma[i] + cs[k].grad[i] * cs[k].mult;
                }
            }
            if (neqn != 0)
            {
                for (i = 1; i <= nparam; i++)
                {
                    eta[i] = 0.e0;
                    for (k = 1; k <= neqn; k++)
                        eta[i] = eta[i] + cs[glob_info.nnineq+k].grad[i] *
                                 cs[glob_info.nnineq+k].mult;
                }
            }
        }
        for (i = 1; i <= nparam; i++)
        {
            if (nobL > 1)
            {
                if (done)
                    psb[i] = hd[i] + param->mult[i] + gamma[i];
                gamma[i] = gamma[i] + hd[i] - grdpsf[i] + eta[i];
            }
            else if (nobL == 1)
            {
                if (done)
                    psb[i] = ob[1].grad[i] + param->mult[i] + gamma[i];
                gamma[i] = gamma[i] + ob[1].grad[i] - grdpsf[i] + eta[i];
            }
            else if (nobL == 0)
            {
                if (done)
                    psb[i] = param->mult[i] + gamma[i];
                gamma[i] = gamma[i] - grdpsf[i] + eta[i];
            }
            if (!done)
                delta[i] = gamma[i];
        }
        if (!done && !glob_log.d0_is0)
        {
            if (nn != 0)
            {
                for (i = 1; i <= nn; i++)
                {
                    if ((feasb) && (i > nineqn))
                        signgj = signeq[i-nineqn];
                    if ((!feasb) || (i <= nineqn))
                        signgj = 1.e0;
                    if ((feasb) && (ncsipn) && (i > nineqn - ncsipn) &&
                        (cs[indxcn[i]].mult == 0.e0))
                        continue;
                    glob_grd.valnom = cs[indxcn[i]].val * signgj;
                    gradcn(nparam, indxcn[i], xnew + 1, cs[indxcn[i]].grad + 1,
                           constr, param->cd);
                }
                resign(nparam, neqn, psf, grdpsf, penp, cs, signeq, 11, 11);
            }
            for (i = 1; i <= nob; i++)
            {
                glob_grd.valnom = ob[i].val;
                if ((i <= nob - nfsip) || (i > nob - nfsip &&
                                           ((ob[i].mult != 0.e0) || (ob[i].mult_L != 0.e0))))
                {
                    if (feasb)
                        gradob(nparam, i, xnew + 1, ob[i].grad + 1, obj, param->cd);
                    else
                        gradcn(nparam, indxob[i], xnew + 1, ob[i].grad + 1,
                               constr, param->cd);
                }
            }
            done = TRUE;
        }
        if (glob_log.d0_is0)
            done = TRUE;
    }
    if (!glob_log.d0_is0)
    {
        if (!(feasb && steps < delta_s && ((sip_viol->type == OBJECT &&
                                            !ob[sip_viol->index].act_sip) || (sip_viol->type == CONSTR &&
                                                                              !cs[sip_viol->index].act_sip))))
        {
            if (*nrst < (5*nparam) || steps > 0.1e0)
            {
                (*nrst)++;
                for (i = 1; i <= nparam; i++)
                {
                    gamma[i] = gamma[i] - delta[i];
                    delta[i] = xnew[i] - param->x[i];
                }
                matrvc(nparam, nparam, hess, delta, hd);
                dhd = scaprd(nparam, delta, hd);
                if (sqrt(scaprd(nparam, delta, delta)) <= glob_grd.epsmac)
                {
                    /* xnew too close to x!! */
                    nstop = 0;
                    glob_prnt.info = 8;
                    return;
                }
                gammd = scaprd(nparam, delta, gamma);
                if (gammd >= (0.2e0*dhd))
                    theta = 1.e0;
                else
                    theta = 0.8e0 * dhd / (dhd - gammd);
                for (i = 1; i <= nparam; i++)
                    eta[i] = hd[i] * (1.e0 - theta) + theta * gamma[i];
                etad = theta * gammd + (1.e0 - theta) * dhd;
                for (i = 1; i <= nparam; i++)
                {
                    for (j = i; j <= nparam; j++)
                    {
                        hess[i][j] = hess[i][j] - hd[i] * hd[j] / dhd +
                                     eta[i] * eta[j] / etad;
                        hess[j][i] = hess[i][j];
                    }
                }
            }
            else
            {
                *nrst = 0;
                diagnl(nparam, 1.e0, hess);
            }
        }
        for (i = 1; i <= nparam; i++)
            param->x[i] = xnew[i];
        x_is_new = TRUE;
    }
    if (neqn != 0 && (feasb))
    {
        i = glob_info.nnineq - nineqn;
        if (i != 0)
        {
            for (j = 1; j <= nparam; j++)
            {
                gamma[j] = 0.e0;
                for (k = 1; k <= i; k++)
                    gamma[j] = gamma[j] + cs[k+nineqn].grad[j] *
                               cs[nineqn+k].mult;
            }
            for (i = 1; i <= nparam; i++)
                psb[i] = psb[i] + gamma[i];
        }
        i = neq - neqn;
        if (i != 0)
        {
            for (j = 1; j <= nparam; j++)
            {
                gamma[j] = 0.e0;
                for (k = 1; k <= i; k++)
                    gamma[j] = gamma[j] + cs[k+neqn+glob_info.nnineq].grad[j] *
                               cs[glob_info.nnineq+neqn+k].mult;
            }
            for (i = 1; i <= nparam; i++)
                psb[i] = psb[i] + gamma[i];
        }
        /* Update penalty parameters for nonlinear equality constraints */
        estlam(nparam, neqn, &ifail, bgbnd, phess, delta, eta,
               gamma, cs, psb, hd, xnew, psmu);
        if (glob_log.get_ne_mult)
            return;
        for (i = 1; i <= neqn; i++)
        {
            if (ifail != 0 || glob_log.d0_is0)
                penp[i] = 2.e0 * penp[i];
            else
            {
                etad = psmu[i] + penp[i];
                if (etad < 1.e0)
                    penp[i] = DMAX1((1.e0 - psmu[i]), (2.e0 * penp[i]));
            }
            if (penp[i] > bgbnd)
            {
                nstop = 0;
                glob_prnt.info = 9;
                return;
            }
        }
        resign(nparam, neqn, psf, grdpsf, penp, cs, signeq, 20, 12);
        *fMp = *fM - *psf;
    }
    if (nfs != 0)
    {
        (*nstart)++;
        np = indexs(*nstart, nfs);
        span[np] = fmax;
        for (i = 1; i <= neqn; i++)
            gm[(np-1)*neqn+i] = cs[glob_info.nnineq+i].val;
        if (neqn != 0)
        {
            psfnew = 0.e0;
            for (i = 1; i <= neqn; i++)
                psfnew = psfnew + gm[i]*penp[i];
        }
        *fM = span[1];
        *fMp = span[1] - psfnew;
        mnm = DMIN1(*nstart, nfs);
        for (i = 2; i <= mnm; i++)
        {
            if (neqn != 0)
            {
                psfnew = 0.e0;
                for (j = 1; j <= neqn; j++)
                    psfnew = psfnew + gm[(i-1)*neqn +j]*penp[j];
            }
            *fM = DMAX1(*fM, span[i]);
            *fMp = DMAX1(*fMp, span[i] - psfnew);
        }
    }
    if (glob_prnt.iprint < 3 || !display)
        return;
    for (i = 1; i <= nob; i++)
    {
        if (!feasb)
        {
            sbout2(glob_prnt.io, nparam, indxob[i], "gradg(j,",
                   ")", ob[i].grad);
            continue;
        }
        if (nob > 1)
            sbout2(glob_prnt.io, nparam, i, "gradf(j,", ")",
                   ob[i].grad);
        if (nob == 1)
            sbout1(glob_prnt.io, nparam, "gradf(j)            ",
                   dummy, ob[1].grad, 2, 2);
    }
    if (ncnstr != 0)
    {
        for (i = 1; i <= ncnstr; i++)
        {
            tempv = cs[indxcn[i]].grad;
            sbout2(glob_prnt.io, nparam, indxcn[i], "gradg(j,", ")", tempv);
        }
        if (neqn != 0)
        {
            sbout1(glob_prnt.io, nparam, "grdpsf(j)           ",
                   dummy, grdpsf, 2, 2);
            sbout1(glob_prnt.io, neqn, "P                   ", dummy,
                   penp, 2, 2);
            sbout1(glob_prnt.io, neqn, "psmu                ", dummy,
                   psmu, 2, 2);
        }
    }
    sbout1(glob_prnt.io, nparam, "multipliers for x   ", dummy,
           param->mult, 2, 2);
    if (ncnstr != 0)
    {
        fprintf(glob_prnt.io, "\t\t\t %s\t %22.14e\n",
                "            for g   ", cs[1].mult);
        for (j = 2; j <= ncnstr; j++)
            fprintf(glob_prnt.io, " \t\t\t\t\t\t %22.14e\n", cs[j].mult);
    }
    if (nobL > 1)
    {
        fprintf(glob_prnt.io, "\t\t\t %s\t %22.14e\n",
                "            for f   ", ob[1].mult);
        for (j = 2; j <= nob; j++)
            fprintf(glob_prnt.io, " \t\t\t\t\t\t %22.14e\n", ob[j].mult);
    }
    for (i = 1; i <= nparam; i++)
    {
        tempv = colvec(hess, i, nparam);
        sbout2(glob_prnt.io, nparam, i, "hess (j,", ")", tempv);
        free_dv(tempv);
    }
    return;
}
/**************************************************************/
/*   CFSQP : Output                                           */
/**************************************************************/


#ifdef __STDC__
static void
out(int miter, int nparam, int nob, int nobL, int nfsip, int ncn,
    int nn, int nineqn, int ncnstr, int ncsipl, int ncsipn,
    int *mesh_pts, double *x, struct _constraint *cs,
    struct _objective *ob, double fM, double fmax,
    double steps, double sktnom, double d0norm, int feasb)
#else
static void
out(miter, nparam, nob, nobL, nfsip, ncn, nn, nineqn, ncnstr, ncsipl, ncsipn,
    mesh_pts, x, cs, ob, fM, fmax, steps, sktnom, d0norm, feasb)
int miter, nparam, nob, nobL, nfsip, ncn, nn, ncnstr, feasb,
ncsipl, ncsipn, nineqn, *mesh_pts;
double  fM, fmax, steps, sktnom, d0norm;
double  *x;
struct  _constraint *cs;
struct  _objective  *ob;
#endif
{
    int i, j, index, display, offset;
    double SNECV, dummy, *adummy, gmax;

    adummy = make_dv(1);
    adummy[1] = 0.e0;
    dummy = 0.e0;
    if (glob_prnt.iter >= miter && nstop != 0)
    {
        glob_prnt.info = 3;
        nstop = 0;
        if (glob_prnt.iprint == 0)
            goto L9000;
    }
    if (glob_prnt.iprint == 0 && glob_prnt.iter < miter)
    {
        glob_prnt.iter++;
        goto L9000;
    }
    if ((glob_prnt.info > 0 && glob_prnt.info < 3) || glob_prnt.info == 7)
        goto L120;
    if (glob_prnt.iprint == 1 && nstop != 0)
    {
        glob_prnt.iter++;
        if (glob_prnt.initvl == 0)
            goto L9000;
        if (feasb && nob > 0)
        {
            fprintf(glob_prnt.io, " objectives\n");
            for (i = 1; i <= nob - nfsip; i++)
            {
                if (nob == nobL)
                    fprintf(glob_prnt.io, " \t\t\t %22.14e\n", ob[i].val);
                else
                    fprintf(glob_prnt.io, " \t\t\t %22.14e\n", fabs(ob[i].val));
            }
            if (nfsip)
            {
                offset = nob - nfsip;
                for (i = 1; i <= glob_info.nfsip; i++)
                {
                    if (nob == nobL)
                        gmax = ob[++offset].val;
                    else
                        gmax = fabs(ob[++offset].val);
                    for (j = 2; j <= mesh_pts[i]; j++)
                    {
                        offset++;
                        if (nob == nobL && ob[offset].val > gmax)
                            gmax = ob[offset].val;
                        else if (nob != nobL && fabs(ob[offset].val) > gmax)
                            gmax = fabs(ob[offset].val);
                    }
                    fprintf(glob_prnt.io, " \t\t\t %22.14e\n", gmax);
                }
            }
        }
        if (glob_info.mode == 1 && glob_prnt.iter > 1 && feasb)
            sbout1(glob_prnt.io, 0, "objective max4      ", fM, adummy, 1, 1);
        if (nob > 1)
            sbout1(glob_prnt.io, 0, "objmax              ", fmax, adummy, 1, 1);
        if (ncnstr == 0)
            fprintf(glob_prnt.io, "\n");
        else
        {
            fprintf(glob_prnt.io, " constraints\n");
            for (i = 1; i <= nineqn - ncsipn; i++)
                fprintf(glob_prnt.io, " \t\t\t %22.14e\n", cs[i].val);
            if (ncsipn)
            {
                offset = nineqn - ncsipn;
                for (i = 1; i <= glob_info.ncsipn; i++)
                {
                    gmax = cs[++offset].val;
                    for (j = 2; j <= mesh_pts[glob_info.nfsip+i]; j++)
                    {
                        offset++;
                        if (cs[offset].val > gmax)
                            gmax = cs[offset].val;
                    }
                    fprintf(glob_prnt.io, " \t\t\t %22.14e\n", gmax);
                }
            }
            for (i = nineqn + 1; i <= glob_info.nnineq - ncsipl; i++)
                fprintf(glob_prnt.io, " \t\t\t %22.14e\n", cs[i].val);
            if (ncsipl)
            {
                offset = glob_info.nnineq - ncsipl;
                for (i = 1; i <= glob_info.ncsipl; i++)
                {
                    gmax = cs[++offset].val;
                    if (feasb)
                        index = glob_info.nfsip + glob_info.ncsipn + i;
                    else
                        index = glob_info.ncsipn + i;
                    for (j = 2; j <= mesh_pts[index]; j++)
                    {
                        offset++;
                        if (cs[offset].val > gmax)
                            gmax = cs[offset].val;
                    }
                    fprintf(glob_prnt.io, " \t\t\t %22.14e\n", gmax);
                }
            }
            for (i = glob_info.nnineq + 1; i <= ncnstr; i++)
                fprintf(glob_prnt.io, " \t\t\t %22.14e\n", cs[i].val);
        }
        if (ncnstr != 0)
            fprintf(glob_prnt.io, "\n");
        goto L9000;
    }
    if (glob_prnt.iprint == 1 && nstop == 0)
        fprintf(glob_prnt.io, " iteration           %26d\n",
                glob_prnt.iter);
    if (glob_prnt.iprint <= 2 && nstop == 0)
        fprintf(glob_prnt.io, " inform              %26d\n",
                glob_prnt.info);
    if (glob_prnt.iprint == 1 && nstop == 0 && (ncsipl + ncsipn) != 0)
        fprintf(glob_prnt.io, " |Xi_k|              %26d\n",
                glob_info.tot_actg_sip);
    if (glob_prnt.iprint == 1 && nstop == 0 && nfsip != 0)
        fprintf(glob_prnt.io, " |Omega_k|           %26d\n",
                glob_info.tot_actf_sip);
    glob_prnt.iter++;
    if (!((glob_prnt.iter) % glob_prnt.iter_mod))
        display = TRUE;
    else
        display = (nstop == 0);
    if (glob_prnt.iter_mod != 1 && display)
        fprintf(glob_prnt.io, "\n iteration           %26d\n",
                glob_prnt.iter - 1);
    if (glob_prnt.initvl == 0 && display)
        sbout1(glob_prnt.io, nparam, "x                   ", dummy, x, 2, 1);
    if (display)
    {
        if (nob > 0)
        {
            fprintf(glob_prnt.io, " objectives\n");
            for (i = 1; i <= nob - nfsip; i++)
            {
                if (nob == nobL)
                    fprintf(glob_prnt.io, " \t\t\t %22.14e\n", ob[i].val);
                else
                    fprintf(glob_prnt.io, " \t\t\t %22.14e\n",
                            fabs(ob[i].val));
            }
        }
        if (nfsip)
        {
            offset = nob - nfsip;
            if (feasb)
                index = glob_info.nfsip;
            else
                index = glob_info.ncsipn;
            for (i = 1; i <= index; i++)
            {
                if (nob == nobL)
                    gmax = ob[++offset].val;
                else
                    gmax = fabs(ob[++offset].val);
                for (j = 2; j <= mesh_pts[i]; j++)
                {
                    offset++;
                    if (nob == nobL && ob[offset].val > gmax)
                        gmax = ob[offset].val;
                    else if (nob != nobL && fabs(ob[offset].val) > gmax)
                        gmax = fabs(ob[offset].val);
                }
                fprintf(glob_prnt.io, " \t\t\t %22.14e\n", gmax);
            }
        }
    }
    if (glob_info.mode == 1 && glob_prnt.iter > 1 && display)
        sbout1(glob_prnt.io, 0, "objective max4      ", fM, adummy, 1, 1);
    if (nob > 1 && display)
        sbout1(glob_prnt.io, 0, "objmax              ", fmax, adummy, 1, 1);
    if (ncnstr != 0 && display)
    {
        fprintf(glob_prnt.io, " constraints\n");
        for (i = 1; i <= nineqn - ncsipn; i++)
            fprintf(glob_prnt.io, " \t\t\t %22.14e\n", cs[i].val);
        if (ncsipn)
        {
            offset = nineqn - ncsipn;
            for (i = 1; i <= glob_info.ncsipn; i++)
            {
                gmax = cs[++offset].val;
                for (j = 2; j <= mesh_pts[glob_info.nfsip+i]; j++)
                {
                    offset++;
                    if (cs[offset].val > gmax)
                        gmax = cs[offset].val;
                }
                fprintf(glob_prnt.io, " \t\t\t %22.14e\n", gmax);
            }
        }
        for (i = nineqn + 1; i <= glob_info.nnineq - ncsipl; i++)
            fprintf(glob_prnt.io, " \t\t\t %22.14e\n", cs[i].val);
        if (ncsipl)
        {
            offset = glob_info.nnineq - ncsipl;
            for (i = 1; i <= glob_info.ncsipl; i++)
            {
                gmax = cs[++offset].val;
                if (feasb)
                    index = glob_info.nfsip + glob_info.ncsipn + i;
                else
                    index = glob_info.ncsipn + i;
                for (j = 2; j <= mesh_pts[index];
                     j++)
                {
                    offset++;
                    if (cs[offset].val > gmax)
                        gmax = cs[offset].val;
                }
                fprintf(glob_prnt.io, " \t\t\t %22.14e\n", gmax);
            }
        }
        for (i = glob_info.nnineq + 1; i <= ncnstr; i++)
            fprintf(glob_prnt.io, " \t\t\t %22.14e\n", cs[i].val);
        if (feasb)
        {
            SNECV = 0.e0;
            for (i = glob_info.nnineq + 1; i <= glob_info.nnineq + nn - nineqn; i++)
                SNECV = SNECV + fabs(cs[i].val);
            if (glob_prnt.initvl == 0 && (nn - nineqn) != 0)
                sbout1(glob_prnt.io, 0, "SNECV               ",
                       SNECV, adummy, 1, 1);
        }
    }
    if (glob_prnt.iter <= 1 && display)
    {
        fprintf(glob_prnt.io, " \n");
        fprintf(glob_prnt.io, " iteration           %26d\n",
                glob_prnt.iter);
        goto L9000;
    }
    if (glob_prnt.iprint >= 2 && glob_prnt.initvl == 0 && display)
        sbout1(glob_prnt.io, 0, "step                ", steps, adummy, 1, 1);
    if (glob_prnt.initvl == 0 && display &&
        (nstop == 0 || glob_prnt.info != 0 || glob_prnt.iprint == 2))
    {
        sbout1(glob_prnt.io, 0, "d0norm              ", d0norm, adummy, 1, 1);
        sbout1(glob_prnt.io, 0, "ktnorm              ", sktnom, adummy, 1, 1);
    }
    if (glob_prnt.initvl == 0 && feasb && display)
        fprintf(glob_prnt.io, " ncallf              %26d\n",
                glob_info.ncallf);
    if (glob_prnt.initvl == 0 && (nn != 0 || !feasb) && display)
        fprintf(glob_prnt.io, " ncallg              %26d\n",
                glob_info.ncallg);
    if (glob_prnt.iprint >= 3 && glob_prnt.iter_mod != 1 && nstop != 0
        && !(glob_prnt.iter % glob_prnt.iter_mod))
        fprintf(glob_prnt.io,
                "\n The following was calculated during iteration %5d:\n",
                glob_prnt.iter);
    if (nstop != 0 && (glob_prnt.iter_mod == 1))
        fprintf(glob_prnt.io, "\n iteration           %26d\n",
                glob_prnt.iter);
L120:
    if (nstop != 0 || glob_prnt.iprint == 0)
        goto L9000;
    fprintf(glob_prnt.io, "\n");
    if (glob_prnt.iprint >= 3)
        fprintf(glob_prnt.io, " inform              %26d\n",
                glob_prnt.info);
    if (glob_prnt.info == 0)
        fprintf(glob_prnt.io,
                "\nNormal termination: You have obtained a solution !!\n");
    if (glob_prnt.info == 0 && sktnom > 0.1e0)
        fprintf(glob_prnt.io,
                "Warning: Norm of Kuhn-Tucker vector is large !!\n");
    if (glob_prnt.info == 3)
    {
        fprintf(glob_prnt.io,
                "\nWarning: Maximum iterations have been reached before\n");
        fprintf(glob_prnt.io, "obtaining a solution !!\n\n");
    }
    if (glob_prnt.info == 4)
    {
        fprintf(glob_prnt.io,
                "\nError : Step size has been smaller than the computed\n");
        fprintf(glob_prnt.io, "machine precision !!\n\n");
    }
    if (glob_prnt.info == 5)
        fprintf(glob_prnt.io,
                "\nError: Failure in constructing d0 !!\n\n");
    if (glob_prnt.info == 6)
        fprintf(glob_prnt.io,
                "\nError: Failure in constructing d1 !!\n\n");
    if (glob_prnt.info == 8)
    {
        fprintf(glob_prnt.io,
                "\nError: The new iterate is numerically equivalent to the\n");
        fprintf(glob_prnt.io,
                "previous iterate, though the stopping criterion is not \n");
        fprintf(glob_prnt.io, "satisfied\n");
    }
    if (glob_prnt.info == 9)
    {
        fprintf(glob_prnt.io,
                "\nError: Could not satisfy nonlinear equality constraints -\n");
        fprintf(glob_prnt.io, "       Penalty parameter too large\n");
    }
    fprintf(glob_prnt.io, "\n");
L9000:
    free_dv(adummy);
    glob_prnt.initvl = 0;
    return;
}
/*************************************************************/
/*   CFSQP : Computation of gradients of objective           */
/*           functions by forward finite differences         */
/*************************************************************/


#ifdef __STDC__
void grobfd(int nparam, int j, double *x, double *gradf,
            void(*obj)(int, int, double *, double *, void *), void *cd)
#else
void grobfd(nparam, j, x, gradf, obj, cd)
int nparam, j;
double *x, *gradf;
void(*obj)();
void   *cd;
#endif
{
    int i;
    double xi, delta;

    for (i = 0; i <= nparam - 1; i++)
    {
        xi = x[i];
        delta = DMAX1(glob_grd.udelta,
                      glob_grd.rteps * DMAX1(1.e0, fabs(xi)));
        if (xi < 0.e0)
            delta = -delta;
        if (!(glob_prnt.ipd == 1 || j != 1 || glob_prnt.iprint < 3))
        {
            /*  formats are not set yet...  */
            if (i == 0)
                fprintf(glob_prnt.io, "\tdelta(i)\t %22.14f\n", delta);
            if (i != 0)
                fprintf(glob_prnt.io, "\t\t\t %22.14f\n", delta);
        }
        x[i] = xi + delta;
        x_is_new = TRUE;
        (*obj)(nparam, j, x, &gradf[i], cd);
        gradf[i] = (gradf[i] - glob_grd.valnom) / delta;
        x[i] = xi;
        x_is_new = TRUE;
    }
    return;
}

/***********************************************************/
/*   CFSQP : Computation of gradients of constraint        */
/*           functions by forward finite differences       */
/***********************************************************/

#ifdef __STDC__
void grcnfd(int nparam, int j, double *x, double *gradg,
            void(*constr)(int, int, double *, double *, void *), void *cd)
#else
void grcnfd(nparam, j, x, gradg, constr, cd)
int nparam, j;
double *x, *gradg;
void(*constr)();
void   *cd;
#endif
{
    int i;
    double xi, delta;

    for (i = 0; i <= nparam - 1; i++)
    {
        xi = x[i];
        delta = DMAX1(glob_grd.udelta,
                      glob_grd.rteps * DMAX1(1.e0, fabs(xi)));
        if (xi < 0.e0)
            delta = -delta;
        if (!(j != 1 || glob_prnt.iprint < 3))
        {
            /*  formats are not set yet...  */
            if (i == 0)
                fprintf(glob_prnt.io, "\tdelta(i)\t %22.14f\n", delta);
            if (i != 0)
                fprintf(glob_prnt.io, "\t\t\t %22.14f\n", delta);
            glob_prnt.ipd = 1;
        }
        x[i] = xi + delta;
        x_is_new = TRUE;
        (*constr)(nparam, j, x, &gradg[i], cd);
        gradg[i] = (gradg[i] - glob_grd.valnom) / delta;
        x[i] = xi;
        x_is_new = TRUE;
    }
    return;
}
/************************************************************/
/*    Utility functions used by CFSQP -                     */
/*    Available functions:                                  */
/*      diagnl        error         estlam                  */
/*      colvec        scaprd        small                   */
/*      fool          matrvc        matrcp                  */
/*      nullvc        resign        sbout1                  */
/*      sbout2        shift         slope                   */
/*      fuscmp        indexs        element                 */
/************************************************************/


#ifdef __STDC__
static void fool(double, double, double *);
#else
static void fool();
#endif

/************************************************************/
/*    Set a=diag*I, a diagonal matrix                       */
/************************************************************/

#ifdef __STDC__
static void diagnl(int nrowa, double diag, double **a)
#else
static void diagnl(nrowa, diag, a)
int nrowa;
double **a, diag;
#endif
{
    int i, j;

    for (i = 1; i <= nrowa; i++)
    {
        for (j = i; j <= nrowa; j++)
        {
            a[i][j] = 0.e0;
            a[j][i] = 0.e0;
        }
        a[i][i] = diag;
    }
    return;
}

/***********************************************************/
/*    Display error messages                               */
/***********************************************************/

#ifdef __STDC__
static void error(const char string[], int *inform)
#else
static void error(string, inform)
const char string[];
int *inform;
#endif
{
    if (glob_prnt.iprint > 0)
        fprintf(stderr, "%s\n", string);
    *inform = 7;
    return;
}

/***********************************************************/
/*    Compute an estimate of multipliers for updating      */
/*    penalty parameter (nonlinear equality constraints)   */
/***********************************************************/

#ifdef __STDC__
static void
estlam(int nparam, int neqn, int *ifail, double bigbnd, double **hess,
       double *cvec, double *a, double *b, struct _constraint *cs,
       double *psb, double *bl, double *bu, double *x)
#else
static void
estlam(nparam, neqn, ifail, bigbnd, hess, cvec, a, b, cs, psb, bl, bu, x)
int nparam, neqn, *ifail;
double bigbnd, **hess, *cvec, *a, *b, *psb, *bl, *bu, *x;
struct _constraint *cs;
#endif
{
    int i, j, zero, one, lwar2, mnn, iout;
    double *ctemp;

    for (i = 1; i <= neqn; i++)
    {
        bl[i] = (-bigbnd);
        bu[i] = bigbnd;
        cvec[i] = scaprd(nparam, cs[i+glob_info.nnineq].grad, psb);
        x[i] = 0.e0;
        for (j = i; j <= neqn; j++)
        {
            hess[i][j] = scaprd(nparam, cs[i+glob_info.nnineq].grad,
                                cs[j+glob_info.nnineq].grad);
            hess[j][i] = hess[i][j];
        }
    }
    zero = 0;
    one = 1;
    iw[1] = 1;
    mnn = 2 * neqn;
    ctemp = convert(hess, neqn, neqn);
    lwar2 = lenw - 1;
    iout = 6;

    ql0001_(&zero, &zero, &one, &neqn, &neqn, &mnn, (ctemp + 1), (cvec + 1),
            (a + 1), (b + 1), (bl + 1), (bu + 1), (x + 1), (w + 1), &iout, ifail,
            &zero, (w + 3), &lwar2, (iw + 1), &leniw, &glob_grd.epsmac);

    free_dv(ctemp);
    return;
}

/**************************************************************/
/*   Extract a column vector from a matrix                    */
/**************************************************************/

#ifdef __STDC__
static double *colvec(double **a, int col, int nrows)
#else
static double *colvec(a, col, nrows)
double **a;
int col, nrows;
#endif
{
    double *temp;
    int i;

    temp = make_dv(nrows);
    for (i = 1;i <= nrows;i++)
        temp[i] = a[i][col];
    return temp;
}

/************************************************************/
/*    Compute the scalar product z=x'y                      */
/************************************************************/

#ifdef __STDC__
static double scaprd(int n, double *x, double *y)
#else
static double scaprd(n, x, y)
double *x, *y;
int n;
#endif
{
    int i;
    double z;

    z = 0.e0;
    for (i = 1;i <= n;i++)
        z = z + x[i] * y[i];
    return z;
}

/***********************************************************/
/*    Used by smallNumber()                                      */
/***********************************************************/

#ifdef __STDC__
static void fool(double x, double y, double *z)
#else
static void fool(x, y, z)
double x, y, *z;
#endif
{
    *z = x * y + y;
    return;
}

/**********************************************************/
/*    Computes the machine precision                      */
/**********************************************************/

static double smallNumber()
{
    double one, two, z, tsmall;

    one = 1.e0;
    two = 2.e0;
    tsmall = one;
    do
    {
        tsmall = tsmall / two;
        fool(tsmall, one, &z);
    }
    while (z > 1.e0);
    return tsmall*two*two;
}

/**********************************************************/
/*     Compares value with threshold to see if exceeds    */
/**********************************************************/

#ifdef __STDC__
static int fuscmp(double val, double thrshd)
#else
static int fuscmp(val, thrshd)
double val, thrshd;
#endif
{
    int temp;

    if (fabs(val) <= thrshd)
        temp = FALSE;
    else
        temp = TRUE;
    return temp;
}

/**********************************************************/
/*     Find the residue of i with respect to nfs          */
/**********************************************************/

#ifdef __STDC__
static int indexs(int i, int nfs)
#else
static int indexs(i, nfs)
int i, nfs;
#endif
{
    int mm = i;

    while (mm > nfs)
        mm -= nfs;
    return mm;
}

/*********************************************************/
/*     Copies matrix a to matrix b                       */
/*********************************************************/

#ifdef __STDC__
static void matrcp(int ndima, double **a, int ndimb, double **b)
#else
static void matrcp(ndima, a, ndimb, b)
double **a, **b;
int ndima, ndimb;
#endif
{
    int i, j;

    for (i = 1; i <= ndima; i++)
        for (j = 1; j <= ndima; j++)
            b[i][j] = a[i][j];
    if (ndimb <= ndima)
        return;
    for (i = 1; i <= ndimb; i++)
    {
        b[ndimb][i] = 0.e0;
        b[i][ndimb] = 0.e0;
    }
    return;
}

/*******************************************************/
/*     Computes y=ax                                   */
/*******************************************************/

#ifdef __STDC__
static void matrvc(int la, int na, double **a, double *x, double *y)
#else
static void matrvc(la, na, a, x, y)
double **a, *x, *y;
int la, na;
#endif
{
    int i, j;
    double yi;

    for (i = 1; i <= la; i++)
    {
        yi = 0.e0;
        for (j = 1; j <= na; j++)
            yi = yi + a[i][j] * x[j];
        y[i] = yi;
    }
    return;
}

/******************************************************/
/*      Set x=0                                       */
/******************************************************/

#ifdef __STDC__
static void nullvc(int nparam, double *x)
#else
static void nullvc(nparam, x)
int nparam;
double *x;
#endif
{
    int i;

    for (i = 1; i <= nparam; i++)
        x[i] = 0.e0;
    return;
}

/*********************************************************/
/*   job1=10: g*signeq,   job1=11: gradg*signeq,         */
/*                        job1=12: job1=10&11            */
/*   job1=20: do not change sign                         */
/*   job2=10: psf,        job2=11: grdpsf,               */
/*     job2=12: job2=10&11            */
/*   job2=20: do not change sign                         */
/*********************************************************/

#ifdef __STDC__
static void
resign(int n, int neqn, double *psf, double *grdpsf, double *penp,
       struct _constraint *cs, double *signeq, int job1, int job2)
#else
static void
resign(n, neqn, psf, grdpsf, penp, cs, signeq, job1, job2)
int job1, job2, n, neqn;
double *psf, *grdpsf, *penp, *signeq;
struct _constraint *cs;
#endif
{
    int i, j, nineq;

    nineq = glob_info.nnineq;
    if (job2 == 10 || job2 == 12)
        *psf = 0.e0;
    for (i = 1; i <= neqn; i++)
    {
        if (job1 == 10 || job1 == 12)
            cs[i+nineq].val =
                signeq[i] * cs[i+nineq].val;
        if (job2 == 10 || job2 == 12)
            *psf = *psf + cs[i+nineq].val * penp[i];
        if (job1 == 10 || job1 == 20)
            continue;
        for (j = 1; j <= n; j++)
            cs[i+nineq].grad[j] = cs[i+nineq].grad[j] * signeq[i];
    }
    if (job2 == 10 || job2 == 20)
        return;
    nullvc(n, grdpsf);
    for (i = 1; i <= n; i++)
        for (j = 1; j <= neqn; j++)
            grdpsf[i] = grdpsf[i] + cs[j+nineq].grad[i] * penp[j];
    return;
}

/**********************************************************/
/*      Write output to file                              */
/**********************************************************/

#ifdef __STDC__
static void
sbout1(FILE *io, int n, const char *s1, double z, double *z1, int job, int level)
#else
static void sbout1(io, n, s1, z, z1, job, level)
FILE *io;
int n, job, level;
double z, *z1;
const char *s1;
#endif
{
    int j;

    if (job != 2)
    {
        if (level == 1)
            fprintf(io, " %s\t %22.14e\n", s1, z);
        if (level == 2)
            fprintf(io, "\t\t\t %s\t %22.14e\n", s1, z);
        return;
    }
    if (n == 0)
        return;
    if (level == 1)
        fprintf(io, " %s\t %22.14e\n", s1, z1[1]);
    if (level == 2)
        fprintf(io, "\t\t\t %s\t %22.14e\n", s1, z1[1]);
    for (j = 2; j <= n; j++)
    {
        if (level == 1)
            fprintf(io, " \t\t\t %22.14e\n", z1[j]);
        if (level == 2)
            fprintf(io, " \t\t\t\t\t\t %22.14e\n", z1[j]);
    }
    return;
}

/*********************************************************/
/*      Write output to file                             */
/*********************************************************/

#ifdef __STDC__
static void
sbout2(FILE *io, int n, int i, const char *s1, const char *s2, double *z)
#else
static void sbout2(io, n, i, s1, s2, z)
FILE *io;
int n, i;
double *z;
const char *s1, *s2;
#endif
{
    int j;

    fprintf(io, "\t\t\t %8s %5d %1s\t %22.14e\n", s1, i, s2, z[1]);
    for (j = 2; j <= n; j++)
        fprintf(io, "\t\t\t\t\t\t %22.14e\n", z[j]);
    return;
}

/*********************************************************/
/*      Extract ii from iact and push in front           */
/*********************************************************/

#ifdef __STDC__
static void shift(int n, int ii, int *iact)
#else
static void shift(n, ii, iact)
int n, ii, *iact;
#endif
{
    int j, k;

    if (ii == iact[1])
        return;
    for (j = 1; j <= n; j++)
    {
        if (ii != iact[j])
            continue;
        for (k = j; k >= 2; k--)
            iact[k] = iact[k-1];
        break;
    }
    if (n != 0)
        iact[1] = ii;
    return;
}

/****************************************************************/
/*      job=0 : Compute the generalized gradient of the minimax */
/*      job=1 : Compute rhog in mode = 1                        */
/****************************************************************/

#ifdef __STDC__
static double
slope(int nob, int nobL, int neqn, int nparam, int feasb,
      struct _objective *ob, double *grdpsf, double *x, double *y,
      double fmax, double theta, int job, double *prev, int old)
#else
static double
slope(nob, nobL, neqn, nparam, feasb, ob, grdpsf, x, y, fmax, theta, job,
      prev, old)
int nob, nobL, neqn, nparam, job, feasb, old;
double fmax, theta;
double *grdpsf, *x, *y, * prev;
struct _objective *ob;
#endif
{
    int i;
    double slope1, rhs, rhog, grdftx, grdfty, diff, grpstx, grpsty;
    double tslope;

    tslope = -bgbnd;
    if (feasb && nob == 0)
        tslope = 0.e0;
    if (neqn == 0 || !feasb)
    {
        grpstx = 0.e0;
        grpsty = 0.e0;
    }
    else
    {
        grpstx = scaprd(nparam, grdpsf, x);
        grpsty = scaprd(nparam, grdpsf, y);
    }
    for (i = 1; i <= nob; i++)
    {
        if (old)
            slope1 = prev[i] + scaprd(nparam, ob[i].grad, x);
        else
            slope1 = ob[i].val + scaprd(nparam, ob[i].grad, x);
        tslope = DMAX1(tslope, slope1);
        if (nobL != nob)
            tslope = DMAX1(tslope, -slope1);
    }
    tslope = tslope - fmax - grpstx;
    if (job == 0)
        return tslope;
    rhs = theta * tslope + fmax;
    rhog = 1.e0;
    for (i = 1; i <= nob; i++)
    {
        grdftx = scaprd(nparam, ob[i].grad, x) - grpstx;
        grdfty = scaprd(nparam, ob[i].grad, y) - grpsty;
        diff = grdfty - grdftx;
        if (diff <= 0.e0)
            continue;
        rhog = DMIN1(rhog, (rhs - ob[i].val - grdftx) / diff);
        if (nobL != nob)
            rhog = DMIN1(rhog, -(rhs + ob[i].val + grdftx) / diff);
    }
    tslope = rhog;
    return tslope;
}

/************************************************************/
/*  Determine whether index is in set                       */
/************************************************************/

#ifdef __STDC__
static int element(int *set, int length, int index)
#else
static int element(set, length, index)
int *set;
int length, index;
#endif
{
    int i, temp;

    temp = 0;
    for (i = 1; i <= length; i++)
    {
        if (set[i] == 0)
            break;
        if (set[i] == index)
        {
            temp = 1;
            return temp;
        }
    }
    return temp;
}
/*************************************************************/
/*     Memory allocation utilities for CFSQP                 */
/*                   */
/*     All vectors and matrices are intended to              */
/*     be subscripted from 1 to n, NOT 0 to n-1.             */
/*     The addreses returned assume this convention.         */
/*************************************************************/


/*************************************************************/
/*     Create double precision vector                        */
/*************************************************************/

#ifdef __STDC__
static double *
make_dv(int len)
#else
static double *
make_dv(len)
int len;
#endif
{
    double *v;

    if (!len)
        len = 1;
    v = (double *)calloc(len, sizeof(double));
    if (!v)
    {
        fprintf(stderr, "Run-time error in make_dv");
        exit(1);
    }
    return --v;
}

/*************************************************************/
/*     Create integer vector                                 */
/*************************************************************/

#ifdef __STDC__
static int *
make_iv(int len)
#else
static int *
make_iv(len)
int len;
#endif
{
    int *v;

    if (!len)
        len = 1;
    v = (int *)calloc(len, sizeof(int));
    if (!v)
    {
        fprintf(stderr, "Run-time error in make_iv");
        exit(1);
    }
    return --v;
}

/*************************************************************/
/*     Create a double precision matrix                      */
/*************************************************************/

#ifdef __STDC__
static double **
make_dm(int rows, int cols)
#else
static double **
make_dm(rows, cols)
int rows, cols;
#endif
{
    double **temp;
    int i;

    if (rows == 0)
        rows = 1;
    if (cols == 0)
        cols = 1;
    temp = (double **)calloc(rows, sizeof(double *));
    if (!temp)
    {
        fprintf(stderr, "Run-time error in make_dm");
        exit(1);
    }
    temp--;
    for (i = 1; i <= rows; i++)
    {
        temp[i] = (double *)calloc(cols, sizeof(double));
        if (!temp[i])
        {
            fprintf(stderr, "Run-time error in make_dm");
            exit(1);
        }
        temp[i]--;
    }
    return temp;
}

/*************************************************************/
/*     Free a double precision vector                        */
/*************************************************************/

#ifdef __STDC__
static void
free_dv(double *v)
#else
static void
free_dv(v)
double *v;
#endif
{
    free((char *)(v + 1));
}

/*************************************************************/
/*     Free an integer vector                                */
/*************************************************************/

#ifdef __STDC__
static void
free_iv(int *v)
#else
static void
free_iv(v)
int *v;
#endif
{
    free((char *)(v + 1));
}

/*************************************************************/
/*     Free a double precision matrix                        */
/*************************************************************/

#ifdef __STDC__
static void
free_dm(double **m, int rows)
#else
static void
free_dm(m, rows)
double **m;
int rows;
#endif
{
    int i;

    if (!rows)
        rows = 1;
    for (i = 1; i <= rows; i++)
        free((char *)(m[i] + 1));
    free((char *)(m + 1));
}

/*************************************************************/
/*     Converts matrix a into a form that can easily be      */
/*     passed to a FORTRAN subroutine.                       */
/*************************************************************/

#ifdef __STDC__
static double *
convert(double **a, int m, int n)
#else
static double *
convert(a, m, n)
double **a;
int m, n;
#endif
{
    double *temp;
    int i, j;

    temp = make_dv(m * n);

    for (i = 1; i <= n; i++)     /* loop thru columns */
        for (j = 1; j <= m; j++)  /* loop thru row     */
            temp[(m*(i-1)+j)] = a[j][i];

    return temp;
}


// Wavelets ----------------------------------------------------------------
void wt1(double a[], unsigned long n, int isign,
         void(*wtstep)(double [], unsigned long, int))
{
    unsigned long nn;

    if (n < 4)
        return;
    if (isign >= 0)
    {
        for (nn = n;nn >= 4;nn >>= 1)
            (*wtstep)(a, nn, isign);
    }
    else
    {
        for (nn = 4;nn <= n;nn <<= 1)
            (*wtstep)(a, nn, isign);
    }
}

void wtn(double a[], unsigned long nn[], int ndim, int isign,
         void(*wtstep)(double [], unsigned long, int))
{
    unsigned long i1, i2, i3, k, n, nnew, nprev = 1, nt, ntot = 1;
    int idim;
    double *wksp;

    for (idim = 1;idim <= ndim;idim++)
        ntot *= nn[idim];
    ask_Tvector(wksp, 1, ntot);
    for (idim = 1;idim <= ndim;idim++)
    {
        n = nn[idim];
        nnew = n * nprev;
        if (n > 4)
        {
            for (i2 = 0;i2 < ntot;i2 += nnew)
            {
                for (i1 = 1;i1 <= nprev;i1++)
                {
                    for (i3 = i1 + i2, k = 1;k <= n;++k, i3 += nprev)
                        wksp[k] = a[i3];
                    if (isign >= 0)
                    {
                        for (nt = n;nt >= 4;nt >>= 1)
                            (*wtstep)(wksp, nt, isign);
                    }
                    else
                    {
                        for (nt = 4;nt <= n;nt <<= 1)
                            (*wtstep)(wksp, nt, isign);
                    }

                    for (i3 = i1 + i2, k = 1;k <= n;++k, i3 += nprev)
                        a[i3] = wksp[k];
                }
            }
        }
        nprev = nnew;
    }
    free_Tvector(wksp, 1, ntot);
}

typedef struct
{
    unsigned int ncof, ioff, joff;
    double *cc, *cr;
}
wavefilt;

wavefilt wfilt;

void pwtset(int n)
{
    int k;
    float sig = -1.0;
    static double c2[3] =
        {
            0.0, 0.707106781186547, 0.707106781186547
        };
    static double c4[5] =
        {
            0.0, 0.4829629131445341, 0.8365163037378079,
            0.2241438680420134, -0.1294095225512604
        };
    static double c12[13] =
        {
            0.0, 0.111540743350, 0.494623890398, 0.751133908021,
            0.315250351709, -0.226264693965, -0.129766867567,
            0.097501605587, 0.027522865530, -0.031582039318,
            0.000553842201, 0.004777257511, -0.001077301085
        };
    static double c20[21] =
        {
            0.0, 0.026670057901, 0.188176800078, 0.527201188932,
            0.688459039454, 0.281172343661, -0.249846424327,
            -0.195946274377, 0.127369340336, 0.093057364604,
            -0.071394147166, -0.029457536822, 0.033212674059,
            0.003606553567, -0.010733175483, 0.001395351747,
            0.001992405295, -0.000685856695, -0.000116466855,
            0.000093588670, -0.000013264203
        };
    static double c2r[2], c4r[5], c12r[13], c20r[21];

    wfilt.ncof = n;
    if (n == 2)
    {
        wfilt.cc = c2;
        wfilt.cr = c2r;
    }
    else if (n == 4)
    {
        wfilt.cc = c4;
        wfilt.cr = c4r;
    }
    else if (n == 12)
    {
        wfilt.cc = c12;
        wfilt.cr = c12r;
    }
    else if (n == 20)
    {
        wfilt.cc = c20;
        wfilt.cr = c20r;
    }
    else
        nrerror("unimplemented value n in pwtset");
    for (k = 1;k <= n;k++)
    {
        wfilt.cr[wfilt.ncof+1-k] = sig * wfilt.cc[k];
        sig = -sig;
    }
    wfilt.ioff = wfilt.joff = -(n >> 1);
}

void pwt(double a[], unsigned long n, int isign)
{
    double ai, ai1, *wksp;
    unsigned long i, ii, jf, jr, k, n1, ni, nj, nh, nmod;

    if (n < 4)
        return;
    ask_Tvector(wksp, 1, n);
    nmod = wfilt.ncof * n;
    n1 = n - 1;
    nh = n >> 1;
    memset(wksp+1,0,n*sizeof(double));
    if (isign >= 0)
    {
        for (ii = 1, i = 1;i <= n;i += 2, ii++)
        {
            ni = i + nmod + wfilt.ioff;
            nj = i + nmod + wfilt.joff;
            double &aux1=wksp[ii];
            double &aux2=wksp[ii+nh];
            unsigned long kmax=4*(wfilt.ncof/4);
            // Loop unrolling (every 4 coefficients)
            for (k = 1;k <= kmax;k+=4)
            {
                unsigned long k_1=k+1;
                unsigned long k_2=k+2;
                unsigned long k_3=k+3;
                jf = n1 & (ni + k);
                jr = n1 & (nj + k);
                aux1 += wfilt.cc[k] * a[jf+1];
                aux2 += wfilt.cr[k] * a[jr+1];
                unsigned long jf_1 = n1 & (ni + k_1);
                unsigned long jr_1 = n1 & (nj + k_1);
                aux1 += wfilt.cc[k_1] * a[jf_1+1];
                aux2 += wfilt.cr[k_1] * a[jr_1+1];
                unsigned long jf_2 = n1 & (ni + k_2);
                unsigned long jr_2 = n1 & (nj + k_2);
                aux1 += wfilt.cc[k_2] * a[jf_2+1];
                aux2 += wfilt.cr[k_2] * a[jr_2+1];
                unsigned long jf_3 = n1 & (ni + k_3);
                unsigned long jr_3 = n1 & (nj + k_3);
                aux1 += wfilt.cc[k_3] * a[jf_3+1];
                aux2 += wfilt.cr[k_3] * a[jr_3+1];
            }
            // The rest of coefficients
            for (k = kmax+1;k <= wfilt.ncof;++k)
            {
                jf = n1 & (ni + k);
                jr = n1 & (nj + k);
                aux1 += wfilt.cc[k] * a[jf+1];
                aux2 += wfilt.cr[k] * a[jr+1];
            }
        }
    }
    else
    {
        for (ii = 1, i = 1;i <= n;i += 2, ii++)
        {
            ai = a[ii];
            ai1 = a[ii+nh];
            ni = i + nmod + wfilt.ioff;
            nj = i + nmod + wfilt.joff;
            for (k = 1;k <= wfilt.ncof;k++)
            {
                jf = (n1 & (ni + k)) + 1;
                jr = (n1 & (nj + k)) + 1;
                wksp[jf] += wfilt.cc[k] * ai;
                wksp[jr] += wfilt.cr[k] * ai1;
            }
        }
    }
    memcpy(&a[1],&wksp[1],n*sizeof(double));
    free_Tvector(wksp, 1, n);
}

/* Gamma function ---------------------------------------------------------- */
#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
    int n;
    double sum, del, ap;

    *gln = gammln(a);
    if (x <= 0.0)
    {
        if (x < 0.0)
            nrerror("x less than 0 in routine gser");
        *gamser = 0.0;
        return;
    }
    else
    {
        ap = a;
        del = sum = 1.0 / a;
        for (n = 1;n <= ITMAX;n++)
        {
            ++ap;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS)
            {
                *gamser = sum * exp(-x + a * log(x) - (*gln));
                return;
            }
        }
        nrerror("a too large, ITMAX too small in routine gser");
        return;
    }
}
#undef ITMAX
#undef EPS

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
    int i;
    double an, b, c, d, del, h;

    *gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1;i <= ITMAX;i++)
    {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = b + an / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS)
            break;
    }
    if (i > ITMAX)
        nrerror("a too large, ITMAX too small in gcf");
    *gammcf = exp(-x + a * log(x) - (*gln)) * h;
}
#undef ITMAX
#undef EPS
#undef FPMIN

double gammp(double a, double x)
{
    double gamser, gammcf, gln;

    if (x < 0.0 || a <= 0.0)
        nrerror("Invalid arguments in routine gammp");
    if (x < (a + 1.0))
    {
        gser(&gamser, a, x, &gln);
        return gamser;
    }
    else
    {
        gcf(&gammcf, a, x, &gln);
        return 1.0 -gammcf;
    }
}

/* Solving linear equation systems via Cholesky ---------------------------- */
void choldc(double *a, int n, double *p)
{
    int i, j, k;
    double sum;

    for (i = 1;i <= n;i++)
    {
        for (j = i;j <= n;j++)
        {
            for (sum = a[i*n+j], k = i - 1;k >= 1;k--)
                sum -= a[i*n+k] * a[j*n+k];
            if (i == j)
            {
                if (sum <= 0.0)
                    nrerror("choldc failed");
                p[i] = sqrt(sum);
            }
            else
                a[j*n+i] = sum / p[i];
        }
    }
}

void cholsl(double *a, int n, double *p, double *b, double *x)
{
    int i, k;
    double sum;

    for (i = 1;i <= n;i++)
    {
        for (sum = b[i], k = i - 1;k >= 1;k--)
            sum -= a[i*n+k] * x[k];
        x[i] = sum / p[i];
    }
    for (i = n;i >= 1;i--)
    {
        for (sum = x[i], k = i + 1;k <= n;k++)
            sum -= a[k*n+i] * x[k];
        x[i] = sum / p[i];
    }
}

/* Polynomial interpolation ------------------------------------------------ */
void polint(double *xa, double *ya, int n, double x, double &y, double &dy)
{
    int i, m, ns = 1;
    double den, dif, dift, ho, hp, w;
    double *c, *d;
    dif = fabs(x - xa[1]);
    ask_Tvector(c, 1, n);
    ask_Tvector(d, 1, n);
    for (i = 1;i <= n;i++)
    {
        if ((dift = fabs(x - xa[i])) < dif)
        {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns--];
    for (m = 1;m < n;m++)
    {
        for (i = 1;i <= n - m;i++)
        {
            ho = xa[i] - x;
            hp = xa[i+m] - x;
            w = c[i+1] - d[i];
            if ((den = ho - hp) == 0.0)
            {
                nrerror("error in routine polint\n");
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        y += (dy = (2 * ns < (n - m) ? c[ns+1] : d[ns--]));
    }
    free_Tvector(d, 1, n);
    free_Tvector(c, 1, n);
}

/* Simulated annealing ----------------------------------------------------- */
double amotsa(double **p, double y[], double psum[], int ndim, double pb[],
              double *yb, double (*funk)(double []), int ihi, double *yhi, double fac,
              double tt, int &idum)
{
    int j;
    double fac1,fac2,yflu,ytry,*ptry;

    ask_Tvector(ptry,1,ndim);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=1;j<=ndim;j++)
        ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=(*funk)(ptry);
    if (ytry <= *yb)
    {
        for (j=1;j<=ndim;j++)
            pb[j]=ptry[j];
        *yb=ytry;
    }
    yflu=ytry-tt*log(ran1(&idum));
    if (yflu < *yhi)
    {
        y[ihi]=ytry;
        *yhi=yflu;
        for (j=1;j<=ndim;j++)
        {
            psum[j] += ptry[j]-p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }
    free_Tvector(ptry,1,ndim);
    return yflu;
}

void amebsa(double **p, double y[], int ndim, double pb[], double *yb,
            double ftol, double (*funk)(double []), int *iter, double temptr)
{
    int i,ihi,ilo,j,m,n,mpts=ndim+1;
    double rtol,sum,swap,yhi,ylo,ynhi,ysave,yt,ytry,*psum;
    int idum=-1;

    ask_Tvector(psum,1,ndim);
    double tt = -temptr;
    for (n=1;n<=ndim;n++)
    {
        for (sum=0.0,m=1;m<=mpts;m++)
            sum += p[m][n];
        psum[n]=sum;
    }
    for (;;)
    {
        ilo=1;
        ihi=2;
        ynhi=ylo=y[1]+tt*log(ran1(&idum));
        yhi=y[2]+tt*log(ran1(&idum));
        if (ylo > yhi)
        {
            ihi=1;
            ilo=2;
            ynhi=yhi;
            yhi=ylo;
            ylo=ynhi;
        }
        for (i=3;i<=mpts;i++)
        {
            yt=y[i]+tt*log(ran1(&idum));
            if (yt <= ylo)
            {
                ilo=i;
                ylo=yt;
            }
            if (yt > yhi)
            {
                ynhi=yhi;
                ihi=i;
                yhi=yt;
            }
            else if (yt > ynhi)
            {
                ynhi=yt;
            }
        }
        rtol=2.0*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo));
        if (rtol < ftol || *iter < 0)
        {
            swap=y[1];
            y[1]=y[ilo];
            y[ilo]=swap;
            for (n=1;n<=ndim;n++)
            {
                swap=p[1][n];
                p[1][n]=p[ilo][n];
                p[ilo][n]=swap;
            }
            break;
        }
        *iter -= 2;
        ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,-1.0,tt,idum);
        if (ytry <= ylo)
        {
            ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,2.0,tt,idum);
        }
        else if (ytry >= ynhi)
        {
            ysave=yhi;
            ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,0.5,tt,idum);
            if (ytry >= ysave)
            {
                for (i=1;i<=mpts;i++)
                {
                    if (i != ilo)
                    {
                        for (j=1;j<=ndim;j++)
                        {
                            psum[j]=0.5*(p[i][j]+p[ilo][j]);
                            p[i][j]=psum[j];
                        }
                        y[i]=(*funk)(psum);
                    }
                }
                *iter -= ndim;
                for (n=1;n<=ndim;n++)
                {
                    for (sum=0.0,m=1;m<=mpts;m++)
                        sum += p[m][n];
                    psum[n]=sum;
                }
            }
        }
        else
            ++(*iter);
    }
    free_Tvector(psum,1,ndim);
}
