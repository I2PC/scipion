/***************************************************************************
 *
 * Authors:    Monica Chagoyen
 *             Carlos Oscar            coss@cnb.csic.es (2011)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
/****************************************************************************/
/* Program for finding the center of an image. The original in Fortran was  */
/* written by A. Santisteban. For any information about the program, contact*/
/* him. Translated into C by Juan P. Secilla (MSC)  Jun/86      */
/****************************************************************************/

#include <data/xmipp_program.h>

// Old code -----------------------------------------------------------------
#define NATURAL       1                /* natural, 1 byte/pixel format */
#define INTFMT        2                /* integer, 2 byte/pixel format */
#define LONGFMT       3                /* long, 4 bytes/pixel format   */
#define FLOATFMT      4                /* float, 4 bytes/pixel format  */
#define FOURIER       5                /* Fourier transform format     */
#define SPIDER        6                /* Spider (header) format       */
#define DEF_IT     5  /*Default iterations number */
#define DEF_DEL    2
#define DEF_IN     2
#define ALLOC_SIZE  65000              /* Allocation unit size         */
typedef unsigned char  BYTE;       /*** Only this and float are used ***/
typedef unsigned short UWORD;
typedef unsigned long  ULONG;

float coseno [1281];
BYTE **imagen;
float r1, r2, r3, cxp1, cxp2, cxp3, cxm1, cxm2, cxm3, rh = 1000., xc0, yc0, del,
        rbajo, ralto, zzz0;
int largo, lancho, indmul, in, ni, mu, ncic;
int ir, m, mu1, mu4, ntot,  mt, idz, ncic2, ntot4;

// Routines needed --------------------------------------------------------
/**************************************************************************/
/* Image allocation routine. Returns a pointer to pointers to individually*/
/* allocated rows of the image. (row, col) are the image dimensions.      */
/* Format indicates the format of the resultant image.                    */
/* Returns: the pointer to the pointers if everything OK, NULL if there   */
/* is not enough memory. In this latter case, it leaves everything as was */
/* before the call.                                                       */
/* Version 2.0: allocates blocks of rows to improve speed                 */
/* Version 3.0: modified for OS/2 1.0 ==> halloc used instead of malloc   */
/*    just because malloc does not give more then 4Mb, and halloc does.   */
/*    should work in any environment just changing halloc for malloc      */
/**************************************************************************/

void **imalloc(int row, int col, int format)
{
    size_t i, j, k;                 /* Counters                             */
    unsigned element_size;          /* Size of element to allocate          */
    unsigned pointer_size;          /* Id. of pointers                      */
    char **temp;                    /* Temporal value to work with          */
    unsigned no_blocks;             /* No. of ALLOC_SIZE blocks to allocate */
    long tot_size;                  /* Total allocation size                */
    unsigned row_alloc;             /* No. of rows to alloc at the same time*/
    unsigned all_size;              /* No. of bits to alloc at one time     */
    unsigned rest_size;             /* Rest of bytes to alloc               */
    unsigned rest_rows;             /* Rest of rows to alloc.               */
    char *aux;                      /* Aux. pointer                         */


    /******************* Assign appropriate value to size flag ******************/

    if (format == FOURIER)
        col += 2;      /* Special treatment for FFT format (see foutrans.c) */

    switch (format)
    {
    case NATURAL:
        element_size = sizeof(BYTE);
        pointer_size = sizeof(BYTE *);
        break;
    case INTFMT:
        element_size = sizeof(UWORD);
        pointer_size = sizeof(UWORD *);
        break;
    case LONGFMT:
        element_size = sizeof(ULONG);
        pointer_size = sizeof(ULONG *);
        break;
    case FLOATFMT:
    case FOURIER:
        element_size = sizeof(float);
        pointer_size = sizeof(float *);
        break;
    default:
        return NULL;
    }

    row_alloc = ALLOC_SIZE / (col * element_size);  /* No. of rows to alloc */
    all_size = row_alloc * col * element_size;
    tot_size = ((long) element_size) * ((long) row) * ((long) col);
    no_blocks = tot_size / all_size;
    rest_size = tot_size - ((long) no_blocks) * ((long) all_size);
    rest_rows = rest_size / (col * element_size);

    /********************* Allocate base pointer ********************************/

    if ((temp = (char **) malloc(row * pointer_size)) == NULL)
        return NULL;               /* Not even this little bit of memory */

    /*********************** Allocate most blocks *******************************/

    j = 0;
    for (i = 0; i < no_blocks; i++)
    {
        if ((aux = (char *)malloc((long)all_size)) == NULL)
        {
            for (j = 0; j < i; j++)
                free(temp[j*row_alloc]);
            free((char *) temp);
            return NULL;
        }
        for (k = 0; k < row_alloc; k++, j++)
            temp [j] = aux + k * col * element_size;
    }

    /*************************** Alloc the last block ***************************/

    if (rest_size != 0)
    {
        if ((aux = (char *)malloc((long)rest_size)) == NULL)
        {
            for (j = 0; j < no_blocks; j++)
                free(temp[j*row_alloc]);
            free(temp);
            return NULL;
        }
        for (k = 0; k < rest_rows; k++, j++)
            temp [j] = aux + k * col * element_size;
    }

    /************************* return OK pointer value  **********************/

    return (void **)temp;
}

/**************************************************************************/
/* This function frees an image previously allocated with image_alloc.    */
/* hfree used instead of free (see imalloc header). Change hfree to free  */
/* for portability                                                        */
/**************************************************************************/

void imfree(char **image, int row, int  col, int format)

{
    size_t i;                       /* Counters                             */
    unsigned element_size;          /* Size of element to allocate          */
    unsigned no_blocks;             /* No. of ALLOC_SIZE blocks to allocate */
    long tot_size;                  /* Total allocation size                */
    unsigned row_alloc;             /* No. of rows to alloc at the same time*/
    unsigned all_size;              /* No. of bits to alloc at one time     */

    if (image == NULL)                  /* No allocation at the moment */
        return;

    if (format == FOURIER)
        col += 2;      /* Special treatment for FFT format (see foutrans.c) */

    switch (format)
    {
    case NATURAL:
        element_size = sizeof(BYTE);
        break;
    case INTFMT:
        element_size = sizeof(UWORD);
        break;
    case LONGFMT:
        element_size = sizeof(ULONG);
        break;
    case FLOATFMT:
    case FOURIER:
        element_size = sizeof(float);
        break;
    default:
        return;
    }

    row_alloc = ALLOC_SIZE / (col * element_size);  /* No. of rows to free  */
    all_size = row_alloc * col * element_size;
    tot_size = ((long) element_size) * ((long) row) * ((long) col);
    no_blocks = tot_size / all_size;

    if (image == NULL)  /* No allocation at the moment */
        return;

    /*************************** Free most blocks *******************************/

    for (i = 0; i < no_blocks; i++)
        if (image [i*row_alloc] != NULL)
            free(image [i*row_alloc]);

    /*************************** Free the last block ****************************/

    if (image [i*row_alloc] != NULL)
        free(image [i*row_alloc]);

    free(image);
}

float conv1x(double y, double x)
/**************************************************************************/
/* Returns the value of the image at the point (y, x) using linear interp-*/
/* olation. For higher accuracy, use "conv3x" (available in 370 assembly  */
/* languaje.          */
/**************************************************************************/
{   short i, j;   /* Row and column */
    float intfila1, intfila2; /* Partial row interpolations */
    float escala;   /* Scale factor */

    j = (short int) y;     /* Trunc the x, y coordinates to short */
    i = (short int) x;

    escala =  y - j;
    /* 1st row interpolation */
    intfila1 = imagen[i][j] + escala * ((short)imagen[i][j+1] - (short)imagen[i][j]);
    /* 2nd row interpolation */
    intfila2 = imagen[i+1][j] + escala * ((short)imagen[i+1][j+1] - (short)imagen [i+1][j]);
    /* Column interpolation */
    return intfila1 + (x - i)*(intfila2 - intfila1);
}

void ergrot(double xc0, double yc0, float* zzz)
{
    //double xc0, yc0;
    //float *zzz; /* It hides global variable, handle with care */

    static float a[266], b[257];
    double za, zb;
    double xp1, xp2, xp3, xm1, xm2, xm3;
    double axp1, axp2, axp3, axm1, axm2, axm3;
    short i7, iz2, iz1, kv, l1, l2, i1, i, j;
    float r, x, y, zz, ai, bi;
    float bj, aj, z;

    r = r1 - r3;
    axp1 = 0.;
    axm1 = 0.;
    axp2 = 0.;
    axm2 = 0.;
    axp3 = 0.;
    axm3 = 0.;
    for (i7 = 1; i7 <= ir; i7++) /* do 31 */
    { r += r3;
        iz2 = ntot - m + mu4;
        iz1 = iz2 - 1;
        for (i = 1; i <= mu; i++) /* do 13 */
        {   iz1 += idz;
            iz2 += idz;
            za = 0.;
            zb = 0.;
            for (kv = 1; kv <= ncic; kv++) /* do 11 */
            { iz1 += m;
                x = xc0 + r * coseno[iz1];
                y = yc0 + r * coseno[iz1-mu4];
                z = conv1x(y, x);
                iz1 += m;
                x = xc0 + r * coseno[iz1];
                y = yc0 + r * coseno[iz1-mu4];
                zz = conv1x(y, x);
                zb += (z - zz);
                iz2 += m;
                x = xc0 + r * coseno[iz2];
                y = yc0 + r * coseno[iz2-mu4];
                z = conv1x(y, x);
                iz2 += m;
                x = xc0 + r * coseno[iz2];
                y = yc0 + r * coseno[iz2-mu4];
                zz = conv1x(y, x);
                za += (z - zz);
            }
            b[i] = zb;
            a[i] = za;
        }
        xp1 = a[mu];
        xp2 = b[mu];
        xm2 = xp2 * xp2;
        xp3 = xp1 * xp2 * coseno[mu4-ncic];
        xm3 = xp3;
        xp2 = xm2 * coseno[mu4-ncic2];
        xp1 = xp1 * xp1;
        xm1 = xp1;
        for (i = 1; i <= mu1; i++) /* do 14 */
        {   l1 = 4 * i * ncic + mu4;
            l2 = mu4;
            ai = a[i];
            bi = b[i];
            xp1 += (ai * ai * coseno[l1]);
            xp2 += (bi * bi * coseno[l1-ncic2]);
            xp3 += (ai * bi * coseno[l1-ncic]);
            xm1 += (ai * ai);
            xm2 += (bi * bi);
            xm3 += (ai * bi * coseno[mu4+ncic]);
            i1 = i + 1;
            ai = a[i];
            bi = b[i];
            for (j = i1; j <= mu; j++) /* do 15 */
            { l1 += ncic2;
                l2 += ncic2;
                aj = a[j];
                bj = b[j];
                double ajai2=2.0*aj*ai;
                double bjbi2=2.0*bj*bi;
                double aibj=ai*bj;
                double ajbi=aj*bi;
                xp1 += (ajai2 * coseno[l1]);
                xm1 += (ajai2 * coseno[l2]);
                xp2 += (bjbi2 * coseno[l1-ncic2]);
                xm2 += (bjbi2 * coseno[l2]);
                xp3 += ((aibj + ajbi) * coseno[l1-ncic]);
                xm3 += (aibj * coseno[l2-ncic] + ajbi * coseno[l2+ncic]);
            }
        }
        axp1 += xp1 * r;
        axm1 += xm1 * r;
        axp2 += xp2 * r;
        axm2 += xm2 * r;
        axp3 += xp3 * r;
        axm3 += xm3 * r;
    }
    (*zzz) = axp1 * cxp1 + axp2 * cxp2 + axp3 * cxp3 + axm1 * cxm1 + axm2 * cxm2 + axm3 * cxm3;
    (*zzz) /= (zzz0 * ir);
}

void suaviza()
{
    unsigned char pix;
    short i, j, k, n;
    long isuma;
    float racua, rbcua, dr, r;

    isuma = 0;
    n = 0;
    i = 0;
    racua = ralto * ralto;
    rbcua = rbajo * rbajo;
    dr = 3.141592653 / (ralto - rbajo);
    for (k = 1; k <= lancho; k++)
    {
        if (k % 16 == 0)
            for (j = 1; j <= largo; j++)
            {
                r = (xc0 - k) * (xc0 - k) + (yc0 - j) * (yc0 - j);
                if (r > racua)
                {
                    isuma += imagen [k][j];
                    n++;
                }
            }
    }
    if (n != 0)
        i = (short int)(((float) isuma) / n + .5);
    m = i;
    pix = i;
    for (k = 1; k <= lancho; k++)
    {
        if (k % 16 == 0)
            for (j = 1; j <= largo; j++)
            {
                r = (xc0 - k) * (xc0 - k) + (yc0 - j) * (yc0 - j);
                if (r < rbcua)
                    continue;
                if (r > racua)
                {
                    imagen[k][j] = pix;
                    continue;
                }
                r = sqrt(r);
                r -= rbajo;
                imagen[k][j] = (unsigned char)((0.5 + 0.5 * cos(dr * r)) * (imagen[k][j] - m) + m);
            }
    }
    printf("\nImage smoothed. Outer mean = %d\n", m);
}

void busca()
{
    static float f[22][22];
    float pi = 3.141592653;
    short in1, mu5, i, mu3, ind, ind1, indice, iflag, j, ix, iy, k;
    float h, th, costh, sinth, fi, anc, xk0, yk0, hpi, x, y, z, xx, yy;
    float b, g, d1, e1, ext;
    int count = 0;

    suaviza();
    in = XMIPP_MIN(in, 10);
    mu1 = mu - 1;
    m = 2 * mu;
    mu4 = mu * ncic;
    mt = 2 * m;
    ntot = mt * ncic;
    ntot4 = ntot + mu4 + ncic;
    zzz0 = ntot;
    ncic2 = 2 * ncic;
    in1 = in + 1;
    h = 2. * pi / ntot;
    idz = 2 - ntot;
    th = ncic * h;
    costh = cos(th);
    sinth = sin(th);
    b = 2. / (th * th) * (1. + costh * costh - sin(2. * th) / th);
    g = 4. / (th * th) * (sinth / th - costh);
    d1 = 2. * th / 45.;
    e1 = d1 * sinth * 2.;
    d1 *= sin(2. * th);
    mu5 = mu4 - 1;
    for (i = 1; i <= mu5; i++)
    {
        fi = i * h;
        coseno[i] = sin(fi);
    }
    coseno[mu4] = 1.;
    mu3 = 2 * mu4;
    for (i = 1; i <= mu5; i++)
    {
        coseno[mu3-i] = coseno[i];
        coseno[mu3+i] = -coseno[i];
        coseno[ntot-i] = -coseno[i];
        coseno[ntot+i] = coseno[i];
    }
    coseno[mu3] = 0.;
    coseno[mu3+mu4] = -1.;
    coseno[ntot] = 0.;
    coseno[ntot+mu4] = 1.;
    ind = 2 * in + 1;
    ind1 = ind - 1;
    indice = 0;
    iflag = 0;
e9:
    if (indice >= ni)
        goto e10;
    if (iflag == 2)
        goto e22;
    anc = (int)(3. + in * del);
    xk0 = xc0;
    yk0 = yc0;
    rh = XMIPP_MIN(rh, r2);
    rh = XMIPP_MIN(rh, xk0 - anc - 1.);
    rh = XMIPP_MIN(rh, yk0 - anc - 1.);
    rh = XMIPP_MIN(rh, largo - xk0 - anc);
    rh = XMIPP_MIN(rh, lancho - yk0 - anc);
    ir = (int)((rh - r1) / r3 + 1);
    hpi = h / 2. / pi / ir;
    cxp1 = 2. * b * e1 * hpi;
    cxp2 = -2. * g * d1 * hpi;
    cxp3 = 2. * (g * e1 - b * d1) * hpi;
    cxm1 = (b * b + e1 * e1) * hpi;
    cxm2 = (g * g + d1 * d1) * hpi;
    cxm3 = 2. * (b * g - d1 * e1) * hpi;
    /*    if(iflag == 1) goto 21            */
    if (iflag == 2)
        goto e22;
    x = xc0 - in * del;
    for (i = 1; i <= ind; i++) /* do 5 */
    { y = yc0 - in * del;
        for (j = 1; j <= ind; j++) /*do 4 */
        {   ergrot(x, y, &z);
            /*     printf ("%10.3f%10.3f%10.5f%10.3f%10.3f AA\n", x,y,z,del,rh);*/
            printf(".");
            fflush(stdout);
            z = 100. + indmul * z;
            f[i][j] = z;
            y += del;
        }
        x += del;
    }
    // Introduced by Sjors dd 28.9.2004 to avoid infinite loops
    count++;
    if (count > 1000)
        goto s1;
e23:
    ext = -1000000.;
    for (i = 1; i <= ind; i++) /* do 7 */
        for (j = 1; j <= ind; j++) /* do 7 */
            {   if (f[i][j] > ext)
            {
                ix = i;
                iy = j;
                ext = f[i][j];
            }
        }
    xc0 = xc0 + (ix - in - 1) * del;
    yc0 = yc0 + (iy - in - 1) * del;
    if (ix == 1 || ix == ind || iy == 1 || iy == ind)
        goto e8;
    del /= in;
    iflag = 2;
    indice++;
    goto e9;
e8:
    iflag = 1;
    goto e9;
e22:
    f[1][1] = f[ix-1][iy-1];
    f[ind][1] = f[ix+1][iy-1];
    f[1][ind] = f[ix-1][iy+1];
    f[ind][ind] = f[ix+1][iy+1];
    f[1][in1] = f[ix-1][iy];
    f[in1][1] = f[ix][iy-1];
    f[ind][in1] = f[ix+1][iy];
    f[in1][ind] = f[ix][iy+1];
    f[in1][in1] = f[ix][iy];
    x = xc0 - (in - 1) * del;
    for (i = 2; i < ind1; i++) /* do 11 */
    { y = yc0 - (in - 1) * del;
        for (j = 2;j <= ind1; j++) /* do 12 */
            {   if (i == in1 && j == in1)
                goto e13;
            ergrot(x, y, &z);
            /*     printf ("%10.3f%10.3f%10.5f%10.3f BB \n", x,y,z,del);*/
            printf(".");
            fflush(stdout);
            z = 100. + indmul * z;
            f[i][j] = z;
e13:
            y += del;
        }
        x += del;
    }
    x = xc0 - (in - 1) * del;
    y = yc0 - (in - 1) * del;

    for (k = 2; k <= ind1; k++) /* do 16 */
        { if (k == in1)
            goto e17;
        xx = xc0 - in * del;
        ergrot(xx, y, &z);
        /* printf ("%10.3f%10.3f%10.5f%10.3f CC \n", xx,y,z,del);*/
        printf(".");
        fflush(stdout);
        z = 100. + indmul * z;
        f[1][k] = z;
        yy = yc0 - in * del;
        ergrot(x, yy, &z);
        /* printf ("%10.3f%10.3f%10.5f%10.3f DD \n", xx,y,z,del);*/
        printf(".");
        fflush(stdout);
        z = 100. + indmul * z;
        f[k][1] = z;
        xx = xc0 + in * del;
        ergrot(xx, y, &z);
        /* printf ("%10.3f%10.3f%10.5f%10.3f EE \n", xx,y,z,del);*/
        printf(".");
        fflush(stdout);
        z = 100. + indmul * z;
        f[ind][k] = z;
        yy = yc0 + in * del;
        ergrot(x, yy, &z);
        /* printf ("%10.3f%10.3f%10.5f%10.3f FF\n", x,yy,z,del);*/
        printf(".");
        fflush(stdout);
        z = 100. + indmul * z;
        f[k][ind] = z;
e17:
        x += del;
        y += del;
    }
    goto e23;
e10:
    std::cout << "\nOptimal center coordinates: x= " << yc0-1 << " ,y= " << xc0-1 << " " << std::endl;
    return;
s1:
    yc0=xc0=-1;
    printf("\nNot-converged\n");
    return;
}

// Wrapper to old code ----------------------------------------------------
class ProgFindCenter2D: public XmippProgram
{
public:
    /// Filenames
    FileName fnIn, fnOroot;

    /// Radii
    double _r1, _r2, _r3, _r4;

    /// Starting center
    double x0, y0;

    /// Harmonic
    int _ncic;

    /// Optimization type
    int _indmul;

    /// Define parameters
    void defineParams()
    {
        addUsageLine("Find the best symmetry of rotation of an image or collection of images");
        addUsageLine("+It is very useful when you want to calculate the rotational symmetry of ");
        addUsageLine("+a set of images (in fact it computes the center of the average image). ");
        addUsageLine("+Image dimensions must be less than 512x512.");
        addSeeAlsoLine("xmipp_image_rotational_spectrum");
        addParamsLine(" -i <file>          : Image, image stack or image selfile");
        addParamsLine(" --oroot <rootname> : Rootname for output files");
        addParamsLine("                    :+ <rootname>_center.xmd contains the image center");
        addParamsLine("                    :+ <rootname>_analyzed_image.xmp (if verbose>=2) contains the image that was actually analyzed");
        addParamsLine("[--r1+ <radius=15>]  : Lowest integration radius (% of the image radius)");
        addParamsLine("[--r2+ <radius=80>]  : Highest integration radius (% of the image radius)");
        addParamsLine("[--r3+ <radius=90>]  : Lowest smoothing radius (% of the image radius)");
        addParamsLine("[--r4+ <radius=100>] : Highest smoothing radius (% of the image radius)");
        addParamsLine("[--x0 <x>]          : Initial center of rotation");
        addParamsLine("[--y0 <y>]          : Initial center of rotation");
        addParamsLine("[--harm+ <n=1>]      : Harmonic to optimize");
        addParamsLine("[--opt+ <opt=-1>]    : Type of optimization (-1=minimize, +1=maximize)");
        addExampleLine("xmipp_image_find_center -i image.xmp");
    }

    void readParams()
    {
        fnIn=getParam("-i");
        fnOroot=getParam("--oroot");
        _r1=getDoubleParam("--r1");
        _r2=getDoubleParam("--r2");
        _r3=getDoubleParam("--r3");
        _r4=getDoubleParam("--r4");
        if (checkParam("--x0"))
            x0=getDoubleParam("--x0");
        else
            x0=-1;
        if (checkParam("--y0"))
            y0=getDoubleParam("--y0");
        else
            y0=-1;
        _ncic=getIntParam("--harm");
        _indmul=getIntParam("--opt");
    }

    void show()
    {
        if (verbose==0)
            return;
        std::cout << "Input:          " << fnIn    << std::endl
        << "Output root:    " << fnOroot << std::endl
        << "R1:             " << _r1     << std::endl
        << "R2:             " << _r2     << std::endl
        << "R3:             " << _r3     << std::endl
        << "R4:             " << _r4     << std::endl
        << "Harmonic:       " << _ncic   << std::endl
        << "Opt:            " << _indmul << std::endl
        << "Initial center: (" << x0 << "," << y0 << ")\n"
        ;
    }

    void run()
    {
        show();
        size_t id;
        // Get the input image or the average of the input images
        Image<double> I, Iaux;
        if (fnIn.isMetaData())
        {
            MetaData MD(fnIn);
            int N=0;
            FOR_ALL_OBJECTS_IN_METADATA(MD)
            {
                Iaux.readApplyGeo(MD,__iter.objId);
                if (N==0)
                    I()=Iaux();
                else
                    I()+=Iaux();
                ++N;
            }
            I()/=N;
        }
        else
        {
            I.read(fnIn, HEADER);
            size_t Xdim, Ydim, Zdim, Ndim;
            I.getDimensions(Xdim, Ydim, Zdim, Ndim);
            I.clear();
            if (Ndim>1)
            {
                for (size_t  n = FIRST_IMAGE; n <= Ndim; n++)
                {
                    Iaux.read(fnIn,DATA,n);
                    if (n == FIRST_IMAGE)
                        I=Iaux;
                    else
                        I()+=Iaux();
                }
                I()/=Ndim;
            }
            else
                I.read(fnIn);
        }
        I().rangeAdjust(0,255);
        if (verbose==2)
            I.write(fnOroot+"_analyzed_image.xmp");

        // Adapt to old code
        if ((imagen = (unsigned char **)imalloc(YSIZE(I()) + 1, XSIZE(I()) + 1, NATURAL)) == NULL)
            REPORT_ERROR(ERR_MEM_NOTENOUGH,"");
        FOR_ALL_ELEMENTS_IN_ARRAY2D(I())
        imagen[i+1][j+1]=(unsigned char)IMGPIXEL(I,i,j);
        if (x0>=0)
            xc0=(float)x0+1; //+1 because of Fortran indexing
        else
            xc0=XSIZE(I())/2+1;
        if (y0>=0)
            yc0=(float)y0+1;
        else
            yc0=YSIZE(I())/2+1;
        r1=(float)(_r1/100.0*XSIZE(I())/2.0);
        r2=(float)(_r2/100.0*XSIZE(I())/2.0);
        r3=1;
        rbajo=(float)(_r3/100.0*XSIZE(I())/2.0);
        ralto=(float)(_r4/100.0*XSIZE(I())/2.0);
        ncic=_ncic;
        indmul=_indmul;
        del = DEF_DEL;
        in = DEF_IN;
        ni = DEF_IT;
        mu = (int)(PI / 2 * r2 / ncic);
        if (mu < 3)
            REPORT_ERROR(ERR_ARG_INCORRECT,"A higher integration radius is needed (r2>6*harm/pi)");
        largo=YSIZE(I());
        lancho=XSIZE(I());
        busca();
        MetaData MD;
        id=MD.addObject();
        if (yc0>0)
        {
            MD.setValue(MDL_X,(double)(yc0-1),id);
            MD.setValue(MDL_Y,(double)(xc0-1),id);
        }
        else
        {
            MD.setValue(MDL_X,(double)(XSIZE(I())/2),id);
            MD.setValue(MDL_Y,(double)(YSIZE(I())/2),id);
        }
        MD.write(fnOroot+"_center.xmd");
        imfree((char**)imagen, YSIZE(I()) + 1, XSIZE(I()) + 1, NATURAL);
    }
};

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgFindCenter2D program;
    program.read(argc, argv);
    program.tryRun();
    return 0;
}
