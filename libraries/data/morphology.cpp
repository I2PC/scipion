/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Pedro A. de Alarc�n     (pedro@cnb.csic.es)
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
#include "morphology.h"
#include "filters.h"

#include <vector>

/* Dilate/Erode 2D steps --------------------------------------------------- */
void dilate2D_step(const MultidimArray<double> &in, MultidimArray<double> &out, int neig,
                   int count)
{
    double sum = 0;
    double dcount=count;
    for (int i = STARTINGY(in) + 1;i < FINISHINGY(in); i++)
        for (int j = STARTINGX(in) + 1;j < FINISHINGX(in); j++)
        {
            if (A2D_ELEM(in,i, j) == 0)
            {
                // 4-environment
                A2D_ELEM(out, i, j) = 0;
                sum = A2D_ELEM(in,i - 1, j) + A2D_ELEM(in,i + 1, j) +
                	  A2D_ELEM(in,i, j - 1) + A2D_ELEM(in,i, j + 1);
                if (sum > dcount)
                { //change the value to foreground
                    A2D_ELEM(out, i, j) = 1;
                }
                else if (neig == 8)
                { //8-environment
                    sum +=A2D_ELEM(in,i - 1, j - 1) + A2D_ELEM(in,i - 1, j + 1) +
                          A2D_ELEM(in,i + 1, j - 1) + A2D_ELEM(in,i + 1, j + 1);
                    if (sum > dcount)
                    { //change the value to foreground
                        A2D_ELEM(out, i, j) = 1;
                    }
                }
            }
            else
            {
                A2D_ELEM(out, i, j) = A2D_ELEM(in,i, j);
            }
            sum = 0;
        }
}

void erode2D_step(const MultidimArray<double> &in, MultidimArray<double> &out, int neig,
                  int count)
{
    double dcount=count;
    double sum;
    for (int i = STARTINGY(in) + 1;i < FINISHINGY(in); i++)
        for (int j = STARTINGX(in) + 1;j < FINISHINGX(in); j++)
        {
            if (A2D_ELEM(in, i, j) == 1)
            {
                // 4-environment
                A2D_ELEM(out,i, j) = 1;
                sum = A2D_ELEM(in,i - 1, j) + A2D_ELEM(in,i + 1, j) +
                      A2D_ELEM(in, i, j - 1) + A2D_ELEM(in, i, j + 1);
                if ((4 - sum) > dcount)
                { //change the value to background
                    A2D_ELEM(out, i, j) = 0;
                }
                else if (neig == 8)
                { //8-environment
                    sum += A2D_ELEM(in, i - 1, j - 1) + A2D_ELEM(in, i - 1, j + 1) +
                    	   A2D_ELEM(in, i + 1, j - 1) + A2D_ELEM(in, i + 1, j + 1);
                    if ((neig - sum) > dcount)
                    { //change the value to background
                        A2D_ELEM(out,i, j) = 0;
                    }
                }
            }
            else
            {
                A2D_ELEM(out, i, j) = A2D_ELEM(in, i, j);
            }
            sum = 0;
        }
}

/* Dilate/Erode 2D --------------------------------------------------------- */
void dilate2D(const MultidimArray<double> &in, MultidimArray<double> &out, int neig, int count,
              int size)
{
    MultidimArray<double> tmp;
    int i;
    tmp = in;
    for (i = 0;i < size;i++)
    {
        dilate2D_step(tmp, out, neig, count);
        tmp = out;
    }
}

void erode2D(const MultidimArray<double> &in, MultidimArray<double> &out, int neig, int count,
             int size)
{
    MultidimArray<double> tmp;
    int i;
    tmp = in;
    for (i = 0;i < size;i++)
    {
        erode2D_step(tmp, out, neig, count);
        tmp = out;
    }
}

/* Opening and closing 2D -------------------------------------------------- */
void closing2D(const MultidimArray<double> &in, MultidimArray<double> &out, int neig,
               int count, int size)
{
    MultidimArray<double> tmp;
    int i;
    tmp = in;
    for (i = 0;i < size;i++)
    { //dilate
        dilate2D_step(tmp, out, neig, count);
        tmp = out;
    }
    for (i = 0;i < size;i++)
    { // erode
        erode2D_step(tmp, out, neig, count);
        tmp = out;
    }
}

void opening2D(const MultidimArray<double> &in, MultidimArray<double> &out, int neig,
               int count, int size)
{
    MultidimArray<double> tmp;
    int i;
    tmp = in;
    for (i = 0;i < size;i++)
    { // erode
        erode2D_step(tmp, out, neig, count);
        tmp = out;
    }
    for (i = 0;i < size;i++)
    { //dilate
        dilate2D_step(tmp, out, neig, count);
        tmp = out;
    }
}

/* Border ------------------------------------------------------------------ */
void border(const MultidimArray<double> &img, MultidimArray<double> &border)
{
    /*
    border.initZeros(img);
    erode2D(img,border,8,0,1);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(border)
       border(i,j)=img(i,j)-border(i,j);
    */
    border.initZeros(img);
    for (int i = STARTINGY(border) + 1; i <= FINISHINGY(border) - 1; i++)
        for (int j = STARTINGX(border) + 1; j <= FINISHINGX(border) - 1; j++)
            if (img(i, j))
            {
                // Number of 4 neighbours
                int N4 = 0;
                if (img(i + 1, j))
                    N4++;
                if (img(i - 1, j))
                    N4++;
                if (img(i, j + 1))
                    N4++;
                if (img(i, j - 1))
                    N4++;
                if (N4 != 4)
                    border(i, j) = 1;
            }
}

/* Simplified Border ------------------------------------------------------- */
void simplify_border(const MultidimArray<double> &border,
                     MultidimArray<double> &simplified_border)
{
    MultidimArray<double> aux;
    aux.initZeros(border);
    for (int i = STARTINGY(border) + 1; i <= FINISHINGY(border) - 1; i++)
        for (int j = STARTINGX(border) + 1; j <= FINISHINGX(border) - 1; j++)
            if (border(i, j))
            {
                // Number of 4 neighbours
                int N4 = 0;
                if (border(i + 1, j))
                    N4++;
                if (border(i - 1, j))
                    N4++;
                if (border(i, j + 1))
                    N4++;
                if (border(i, j - 1))
                    N4++;
                if (N4 <= 2)
                    aux(i, j) = 1;
            }

    // Again removing all those without any 4 neighbour
    simplified_border.initZeros(aux);
    for (int i = STARTINGY(border) + 1; i <= FINISHINGY(border) - 1; i++)
        for (int j = STARTINGX(border) + 1; j <= FINISHINGX(border) - 1; j++)
            if (aux(i, j))
            {
                // Number of 4 neighbours
                int N4 = 0;
                if (aux(i + 1, j))
                    N4++;
                if (aux(i - 1, j))
                    N4++;
                if (aux(i, j + 1))
                    N4++;
                if (aux(i, j - 1))
                    N4++;
                if (N4 >= 1)
                    simplified_border(i, j) = 1;
            }
}

/* Random convex hull ------------------------------------------------------ */
void random_convex_hull(const MultidimArray<double> &img, MultidimArray<double> &hull,
                        long N)
{
    hull = img;

    std::vector<int> full_tx, full_ty;
    // Build the list of points
    FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
    if (img(i, j) > 0)
    {
        full_tx.push_back(j);
        full_ty.push_back(i);
    }

    long n = 0;
    int idx_max = full_tx.size() - 1;
    if (idx_max < 2)
        return;
    while (n++ < N)
    {
        // Generate 3 random points
        int i = ROUND(rnd_unif(0, idx_max));
        int j = ROUND(rnd_unif(0, idx_max));
        int k = ROUND(rnd_unif(0, idx_max));

        // Get their positions
        int tx[3], ty[3];
        tx[0] = full_tx[i];
        ty[0] = full_ty[i];
        tx[1] = full_tx[j];
        ty[1] = full_ty[j];
        tx[2] = full_tx[k];
        ty[2] = full_ty[k];

        // Fill the triangle
        fillTriangle(hull, tx, ty, 1);
    }
    /*
    MultidimArray<double> aux;
    aux.initZeros(hull);
    closing2D(hull,aux,8,0,1);
    hull=aux;
    */
    fillBinaryObject(hull);
}

/* Dilate/erode 3D steps --------------------------------------------------- */
void dilate3D_step(const MultidimArray<double> &in, MultidimArray<double> &out, int neig,
                   int count)
{
    int sum = 0;
    for (int k = STARTINGZ(in) + 1;k < FINISHINGZ(in); k++)
        for (int i = STARTINGY(in) + 1;i < FINISHINGY(in); i++)
            for (int j = STARTINGX(in) + 1;j < FINISHINGX(in); j++)
            {
                if (A3D_ELEM(in,k, i, j) == 0)
                {
                    // 6-environment
                    A3D_ELEM(out,k, i, j) = 0;
                    sum = (int)(A3D_ELEM(in,k - 1, i, j) + A3D_ELEM(in,k + 1, i, j) +
                    		A3D_ELEM(in,k, i - 1, j) + A3D_ELEM(in,k, i + 1, j)
                                + A3D_ELEM(in,k, i, j - 1) + A3D_ELEM(in,k, i, j + 1));
                    if (sum > count)
                    { //change the value to foreground
                        A3D_ELEM(out,k, i, j) = 1;
                    }
                    else if (neig == 18)
                    { //18-environment
                        sum = (int)(sum + A3D_ELEM(in,k - 1, i, j - 1) + A3D_ELEM(in,k - 1, i, j + 1) + A3D_ELEM(in,k + 1, i, j - 1) + A3D_ELEM(in,k + 1, i, j + 1) +
                                    A3D_ELEM(in,k, i + 1, j + 1) + A3D_ELEM(in,k, i + 1, j - 1) + A3D_ELEM(in,k, i - 1, j + 1) + A3D_ELEM(in,k, i - 1, j - 1) +
                                    A3D_ELEM(in,k - 1, i + 1, j) + A3D_ELEM(in,k - 1, i - 1, j) + A3D_ELEM(in,k + 1, i + 1, j) + A3D_ELEM(in,k + 1, i - 1, j));
                        if (sum > count)
                        { //change the value to foreground
                            A3D_ELEM(out,k, i, j) = 1;
                        }
                    }
                    else if (neig == 26)
                    { //26-environment
                        sum = (int)(sum + A3D_ELEM(in,k - 1, i + 1, j + 1) + A3D_ELEM(in,k - 1, i + 1, j - 1) +
                                    A3D_ELEM(in,k - 1, i - 1, j + 1) + A3D_ELEM(in,k - 1, i - 1, j - 1) +
                                    A3D_ELEM(in,k + 1, i + 1, j + 1) + A3D_ELEM(in,k + 1, i + 1, j - 1) +
                                    A3D_ELEM(in,k + 1, i - 1, j + 1) + A3D_ELEM(in,k + 1, i - 1, j - 1));
                        if (sum > count)
                        { //change the value to foreground
                            A3D_ELEM(out,k, i, j) = 1;
                        }
                    }

                }
                else
                {
                    A3D_ELEM(out,k, i, j) = A3D_ELEM(in,k, i, j);
                }
                sum = 0;
            }
}

void erode3D_step(const MultidimArray<double> &in, MultidimArray<double> &out, int neig,
                  int count)
{
    int sum = 0;
    for (int k = STARTINGZ(in) + 1;k < FINISHINGZ(in); k++)
        for (int i = STARTINGY(in) + 1;i < FINISHINGY(in); i++)
            for (int j = STARTINGX(in) + 1;j < FINISHINGX(in); j++)
            {
                if (A3D_ELEM(in,k, i, j) == 1)
                {
                    // 6-environment
                    A3D_ELEM(out,k, i, j) = 1;

                    sum = (int)(A3D_ELEM(in,k - 1, i, j) + A3D_ELEM(in,k + 1, i, j) + A3D_ELEM(in,k, i - 1, j) + A3D_ELEM(in,k, i + 1, j)
                                + A3D_ELEM(in,k, i, j - 1) + A3D_ELEM(in,k, i, j + 1));
                    if ((6 - sum) > count)
                    { //change the value to background
                        A3D_ELEM(out,k, i, j) = 0;
                    }
                    else if (neig == 18)
                    { //18-environment
                        sum = (int)(sum + A3D_ELEM(in,k - 1, i, j - 1) + A3D_ELEM(in,k - 1, i, j + 1) + A3D_ELEM(in,k + 1, i, j - 1) + A3D_ELEM(in,k + 1, i, j + 1) +
                                    A3D_ELEM(in,k, i + 1, j + 1) + A3D_ELEM(in,k, i + 1, j - 1) + A3D_ELEM(in,k, i - 1, j + 1) + A3D_ELEM(in,k, i - 1, j - 1) +
                                    A3D_ELEM(in,k - 1, i + 1, j) + A3D_ELEM(in,k - 1, i - 1, j) + A3D_ELEM(in,k + 1, i + 1, j) + A3D_ELEM(in,k + 1, i - 1, j));
                        if ((neig - sum) > count)
                        { //change the value to background
                            A3D_ELEM(out,k, i, j) = 0;
                        }
                    }
                    else if (neig == 26)
                    { //26-environment
                        sum = (int)(sum + A3D_ELEM(in,k - 1, i, j - 1) + A3D_ELEM(in,k - 1, i, j + 1) + A3D_ELEM(in,k + 1, i, j - 1) + A3D_ELEM(in,k + 1, i, j + 1) +
                                    A3D_ELEM(in,k, i + 1, j + 1) + A3D_ELEM(in,k, i + 1, j - 1) + A3D_ELEM(in,k, i - 1, j + 1) + A3D_ELEM(in,k, i - 1, j - 1) +
                                    A3D_ELEM(in,k - 1, i + 1, j) + A3D_ELEM(in,k - 1, i - 1, j) + A3D_ELEM(in,k + 1, i + 1, j) + A3D_ELEM(in,k + 1, i - 1, j) +
                                    A3D_ELEM(in,k - 1, i + 1, j + 1) + A3D_ELEM(in,k - 1, i + 1, j - 1) +
                                    A3D_ELEM(in,k - 1, i - 1, j + 1) + A3D_ELEM(in,k - 1, i - 1, j - 1) +
                                    A3D_ELEM(in,k + 1, i + 1, j + 1) + A3D_ELEM(in,k + 1, i + 1, j - 1) +
                                    A3D_ELEM(in,k + 1, i - 1, j + 1) + A3D_ELEM(in,k + 1, i - 1, j - 1));

                        if ((neig - sum) > count)
                        { //change the value to background
                            A3D_ELEM(out,k, i, j) = 0;
                        }
                    }

                }
                else
                {
                    A3D_ELEM(out,k, i, j) = A3D_ELEM(in,k, i, j);
                }
                sum = 0;
            }

}

/* Dilate/Erode 3D --------------------------------------------------------- */
void dilate3D(const MultidimArray<double> &in, MultidimArray<double> &out, int neig, int count,
              int size)
{
    MultidimArray<double> tmp;
    int i;
    tmp = in;
    for (i = 0;i < size;i++)
    {
        dilate3D_step(tmp, out, neig, count);
        tmp = out;
    }
}

void erode3D(const MultidimArray<double> &in, MultidimArray<double> &out, int neig, int count,
             int size)
{
    MultidimArray<double> tmp;
    int i;
    tmp = in;
    for (i = 0;i < size;i++)
    {
        erode3D_step(tmp, out, neig, count);
        tmp = out;
    }

}

/* Opening/Closing 3D ------------------------------------------------------ */
void closing3D(const MultidimArray<double> &in, MultidimArray<double> &out, int neig,
               int count, int size)
{
    MultidimArray<double> tmp;
    int i;
    tmp = in;
    for (i = 0;i < size;i++)
    { //dilate
        dilate3D_step(tmp, out, neig, count);
        tmp = out;
    }
    for (i = 0;i < size;i++)
    { // erode
        erode3D_step(tmp, out, neig, count);
        tmp = out;
    }
}

void opening3D(const MultidimArray<double> &in, MultidimArray<double> &out, int neig,
               int count, int size)
{
    MultidimArray<double> tmp;
    int i;
    tmp=in;
    for (i = 0;i < size;i++)
    { // erode
        erode3D_step(tmp, out, neig, count);
        tmp = out;
    }
    for (i = 0;i < size;i++)
    { //dilate
        dilate3D_step(tmp, out, neig, count);
        tmp = out;
    }
}

// Grey operations ---------------------------------------------------------
void dilate3D(const MultidimArray<double> &in,
              const MultidimArray<double> &structuringElement,
              MultidimArray<double> &out)
{
    out.initZeros(in);
    double maxval=in.computeMax();
    for (size_t kk=0; kk<ZSIZE(out); kk++)
        for (size_t ii=0; ii<YSIZE(out); ii++)
            for (size_t jj=0; jj<XSIZE(out); jj++)
            {
                double maxLocal=DIRECT_A3D_ELEM(in,kk,ii,jj)+
                    A3D_ELEM(structuringElement,0,0,0);
                int k0=XMIPP_MAX(0,kk+STARTINGZ(structuringElement))-kk;
                int kF=XMIPP_MIN(ZSIZE(out)-1,kk+FINISHINGZ(structuringElement))-kk;
                int i0=XMIPP_MAX(0,ii+STARTINGY(structuringElement))-ii;
                int iF=XMIPP_MIN(YSIZE(out)-1,ii+FINISHINGY(structuringElement))-ii;
                int j0=XMIPP_MAX(0,jj+STARTINGX(structuringElement))-jj;
                int jF=XMIPP_MIN(XSIZE(out)-1,jj+FINISHINGX(structuringElement))-jj;
                for (int k=k0; k<=kF; k++)
                    for (int i=i0; i<=iF; i++)
                        for (int j=j0; j<=jF; j++)
                        {
                            double val=DIRECT_A3D_ELEM(in,kk+k,ii+i,jj+j)+
                                A3D_ELEM(structuringElement,k,i,j);
                            maxLocal=XMIPP_MAX(maxLocal,val);
                        }
                maxLocal=XMIPP_MIN(maxLocal,maxval);
                DIRECT_A3D_ELEM(out,kk,ii,jj)=maxLocal;
            }
}

void erode3D(const MultidimArray<double> &in,
              const MultidimArray<double> &structuringElement,
              MultidimArray<double> &out)
{
    out.initZeros(in);
    double minval=in.computeMin();
    for (size_t kk=0; kk<ZSIZE(out); kk++)
        for (size_t ii=0; ii<YSIZE(out); ii++)
            for (size_t jj=0; jj<XSIZE(out); jj++)
            {
                double minLocal=DIRECT_A3D_ELEM(in,kk,ii,jj)-
                    A3D_ELEM(structuringElement,0,0,0);
                int k0=XMIPP_MAX(0,kk+STARTINGZ(structuringElement))-kk;
                int kF=XMIPP_MIN(ZSIZE(out)-1,kk+FINISHINGZ(structuringElement))-kk;
                int i0=XMIPP_MAX(0,ii+STARTINGY(structuringElement))-ii;
                int iF=XMIPP_MIN(YSIZE(out)-1,ii+FINISHINGY(structuringElement))-ii;
                int j0=XMIPP_MAX(0,jj+STARTINGX(structuringElement))-jj;
                int jF=XMIPP_MIN(XSIZE(out)-1,jj+FINISHINGX(structuringElement))-jj;
                for (int k=k0; k<=kF; k++)
                    for (int i=i0; i<=iF; i++)
                        for (int j=j0; j<=jF; j++)
                        {
                            double val=DIRECT_A3D_ELEM(in,kk+k,ii+i,jj+j)-
                                A3D_ELEM(structuringElement,k,i,j);
                            minLocal=XMIPP_MIN(minLocal,val);
                        }
                minLocal=XMIPP_MAX(minLocal,minval);
                DIRECT_A3D_ELEM(out,kk,ii,jj)=minLocal;
            }
}

/* Sharpening -------------------------------------------------------------- */
void sharpening(const MultidimArray<double> &in, double width, double strength,
    MultidimArray<double> &out)
{
    // Build the quadratic kernel
    int diameter=(int)(2*width+1);
    MultidimArray<double> kernel(diameter,diameter,diameter);
    kernel.setXmippOrigin();
    
    double width2=width*width;
    double minval, maxval;
    in.computeDoubleMinMax(minval,maxval);
    double c=minval+(maxval-minval)*strength/100;
    double a=(minval-c)/width2;

    FOR_ALL_ELEMENTS_IN_ARRAY3D(kernel)
    {
        double r2=k*k+i*i+j*j;
        A3D_ELEM(kernel,k,i,j)=a*r2+c;
    }
    
    // Create the dilated and eroded versions
    MultidimArray<double> dilated, eroded;
    dilate3D(in,kernel,dilated);
    erode3D(in,kernel,eroded);
#ifdef DEBUG
    Image<double> save;
    save()=dilated; save.write("PPPdilated.vol");
    save()=eroded; save.write("PPPeroded.vol");
#endif
    
    // Sharpen
    out=in;
    double eps=(maxval-minval)*0.01;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(in)
    {
        double threshold=0.5*(A3D_ELEM(dilated,k,i,j)+A3D_ELEM(eroded,k,i,j));
        if      (A3D_ELEM(in,k,i,j)>threshold+eps)  A3D_ELEM(out,k,i,j)=A3D_ELEM(dilated,k,i,j);
        else if (A3D_ELEM(in,k,i,j)<threshold-eps)  A3D_ELEM(out,k,i,j)=A3D_ELEM(eroded,k,i,j);
    }
}
