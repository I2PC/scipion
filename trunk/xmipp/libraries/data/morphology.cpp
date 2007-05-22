/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Pedro A. de Alarc�n     (pedro@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
#include "morphology.h"
#include "filters.h"

#include <vector>

/* Dilate/Erode 2D steps --------------------------------------------------- */
void dilate2D_step(const matrix2D<double> &in, matrix2D<double> &out, int neig,
                   int count)
{
    int sum = 0;
    for (int i = STARTINGX(in) + 1;i < FINISHINGX(in); i++)
        for (int j = STARTINGY(in) + 1;j < FINISHINGY(in); j++)
        {
            if (in(i, j) == 0)
            {
                // 4-environment
                out(i, j) = 0;
                sum = (int)(in(i - 1, j) + in(i + 1, j) + in(i, j - 1) + in(i, j + 1));
                if (sum > count)
                { //change the value to foreground
                    out(i, j) = 1;
                }
                else if (neig == 8)
                { //8-environment
                    sum = (int)(sum + in(i - 1, j - 1) + in(i - 1, j + 1) + in(i + 1, j - 1) + in(i + 1, j + 1));
                    if (sum > count)
                    { //change the value to foreground
                        out(i, j) = 1;
                    }
                }
            }
            else
            {
                out(i, j) = in(i, j);
            }
            sum = 0;
        }
}

void erode2D_step(const matrix2D<double> &in, matrix2D<double> &out, int neig,
                  int count)
{
    int sum = 0;
    for (int i = STARTINGX(in) + 1;i < FINISHINGX(in); i++)
        for (int j = STARTINGY(in) + 1;j < FINISHINGY(in); j++)
        {
            if (in(i, j) == 1)
            {
                // 4-environment
                out(i, j) = 1;
                sum = (int)(in(i - 1, j) + in(i + 1, j) + in(i, j - 1) + in(i, j + 1));
                if ((4 - sum) > count)
                { //change the value to background
                    out(i, j) = 0;
                }
                else if (neig == 8)
                { //8-environment
                    sum = (int)(sum + in(i - 1, j - 1) + in(i - 1, j + 1) + in(i + 1, j - 1) + in(i + 1, j + 1));
                    if ((neig - sum) > count)
                    { //change the value to background
                        out(i, j) = 0;
                    }
                }
            }
            else
            {
                out(i, j) = in(i, j);
            }
            sum = 0;
        }
}

/* Dilate/Erode 2D --------------------------------------------------------- */
void dilate2D(const matrix2D<double> &in, matrix2D<double> &out, int neig, int count,
              int size)
{
    matrix2D<double> tmp;
    int i;
    tmp.resize(in);
    tmp = in;
    for (i = 0;i < size;i++)
    {
        dilate2D_step(tmp, out, neig, count);
        tmp = out;
    }
}

void erode2D(const matrix2D<double> &in, matrix2D<double> &out, int neig, int count,
             int size)
{
    matrix2D<double> tmp;
    int i;
    tmp.resize(in);
    tmp = in;
    for (i = 0;i < size;i++)
    {
        erode2D_step(tmp, out, neig, count);
        tmp = out;
    }
}

/* Opening and closing 2D -------------------------------------------------- */
void closing2D(const matrix2D<double> &in, matrix2D<double> &out, int neig,
               int count, int size)
{
    matrix2D<double> tmp;
    int i;
    tmp.resize(in);
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

void opening2D(const matrix2D<double> &in, matrix2D<double> &out, int neig,
               int count, int size)
{
    matrix2D<double> tmp;
    int i;
    tmp.resize(in);
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
void border(const matrix2D<double> &img, matrix2D<double> &border)
{
    /*
    border.init_zeros(img);
    erode2D(img,border,8,0,1);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(border)
       border(i,j)=img(i,j)-border(i,j);
    */
    border.init_zeros(img);
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
void simplify_border(const matrix2D<double> &border,
                     matrix2D<double> &simplified_border)
{
    matrix2D<double> aux;
    aux.init_zeros(border);
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
    simplified_border.init_zeros(aux);
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
void random_convex_hull(const matrix2D<double> &img, matrix2D<double> &hull,
                        long N)
{
    hull = img;

    vector<int> full_tx, full_ty;
    // Build the list of points
    FOR_ALL_ELEMENTS_IN_MATRIX2D(img)
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
        fill_triangle(hull, tx, ty, 1);
    }
    /*
    matrix2D<double> aux;
    aux.init_zeros(hull);
    closing2D(hull,aux,8,0,1);
    hull=aux;
    */
    fill_binary_object(hull);
}

/* Dilate/erode 3D steps --------------------------------------------------- */
void dilate3D_step(const matrix3D<double> &in, matrix3D<double> &out, int neig,
                   int count)
{
    int sum = 0;
    for (int k = STARTINGZ(in) + 1;k < FINISHINGZ(in); k++)
        for (int i = STARTINGY(in) + 1;i < FINISHINGY(in); i++)
            for (int j = STARTINGX(in) + 1;j < FINISHINGX(in); j++)
            {
                if (in(k, i, j) == 0)
                {
                    // 6-environment
                    out(k, i, j) = 0;
                    sum = (int)(in(k - 1, i, j) + in(k + 1, i, j) + in(k, i - 1, j) + in(k, i + 1, j)
                                + in(k, i, j - 1) + in(k, i, j + 1));
                    if (sum > count)
                    { //change the value to foreground
                        out(k, i, j) = 1;
                    }
                    else if (neig == 18)
                    { //18-environment
                        sum = (int)(sum + in(k - 1, i, j - 1) + in(k - 1, i, j + 1) + in(k + 1, i, j - 1) + in(k + 1, i, j + 1) +
                                    in(k, i + 1, j + 1) + in(k, i + 1, j - 1) + in(k, i - 1, j + 1) + in(k, i - 1, j - 1) +
                                    in(k - 1, i + 1, j) + in(k - 1, i - 1, j) + in(k + 1, i + 1, j) + in(k + 1, i - 1, j));
                        if (sum > count)
                        { //change the value to foreground
                            out(k, i, j) = 1;
                        }
                    }
                    else if (neig == 26)
                    { //26-environment
                        sum = (int)(sum + in(k - 1, i + 1, j + 1) + in(k - 1, i + 1, j - 1) +
                                    in(k - 1, i - 1, j + 1) + in(k - 1, i - 1, j - 1) +
                                    in(k + 1, i + 1, j + 1) + in(k + 1, i + 1, j - 1) +
                                    in(k + 1, i - 1, j + 1) + in(k + 1, i - 1, j - 1));
                        if (sum > count)
                        { //change the value to foreground
                            out(k, i, j) = 1;
                        }
                    }

                }
                else
                {
                    out(k, i, j) = in(k, i, j);
                }
                sum = 0;
            }
}

void erode3D_step(const matrix3D<double> &in, matrix3D<double> &out, int neig,
                  int count)
{
    int sum = 0;
    for (int k = STARTINGZ(in) + 1;k < FINISHINGZ(in); k++)
        for (int i = STARTINGY(in) + 1;i < FINISHINGY(in); i++)
            for (int j = STARTINGX(in) + 1;j < FINISHINGX(in); j++)
            {
                if (in(k, i, j) == 1)
                {
                    // 6-environment
                    out(k, i, j) = 1;

                    sum = (int)(in(k - 1, i, j) + in(k + 1, i, j) + in(k, i - 1, j) + in(k, i + 1, j)
                                + in(k, i, j - 1) + in(k, i, j + 1));
                    if ((6 - sum) > count)
                    { //change the value to background
                        out(k, i, j) = 0;
                    }
                    else if (neig == 18)
                    { //18-environment
                        sum = (int)(sum + in(k - 1, i, j - 1) + in(k - 1, i, j + 1) + in(k + 1, i, j - 1) + in(k + 1, i, j + 1) +
                                    in(k, i + 1, j + 1) + in(k, i + 1, j - 1) + in(k, i - 1, j + 1) + in(k, i - 1, j - 1) +
                                    in(k - 1, i + 1, j) + in(k - 1, i - 1, j) + in(k + 1, i + 1, j) + in(k + 1, i - 1, j));
                        if ((neig - sum) > count)
                        { //change the value to background
                            out(k, i, j) = 0;
                        }
                    }
                    else if (neig == 26)
                    { //26-environment
                        sum = (int)(sum + in(k - 1, i, j - 1) + in(k - 1, i, j + 1) + in(k + 1, i, j - 1) + in(k + 1, i, j + 1) +
                                    in(k, i + 1, j + 1) + in(k, i + 1, j - 1) + in(k, i - 1, j + 1) + in(k, i - 1, j - 1) +
                                    in(k - 1, i + 1, j) + in(k - 1, i - 1, j) + in(k + 1, i + 1, j) + in(k + 1, i - 1, j) +
                                    in(k - 1, i + 1, j + 1) + in(k - 1, i + 1, j - 1) +
                                    in(k - 1, i - 1, j + 1) + in(k - 1, i - 1, j - 1) +
                                    in(k + 1, i + 1, j + 1) + in(k + 1, i + 1, j - 1) +
                                    in(k + 1, i - 1, j + 1) + in(k + 1, i - 1, j - 1));

                        if ((neig - sum) > count)
                        { //change the value to background
                            out(k, i, j) = 0;
                        }
                    }

                }
                else
                {
                    out(k, i, j) = in(k, i, j);
                }
                sum = 0;
            }

}

/* Dilate/Erode 3D --------------------------------------------------------- */
void dilate3D(const matrix3D<double> &in, matrix3D<double> &out, int neig, int count,
              int size)
{
    matrix3D<double> tmp;
    int i;
    tmp.resize(in);
    tmp = in;
    for (i = 0;i < size;i++)
    {
        dilate3D_step(tmp, out, neig, count);
        tmp = out;
    }
}

void erode3D(const matrix3D<double> &in, matrix3D<double> &out, int neig, int count,
             int size)
{
    matrix3D<double> tmp;
    int i;
    tmp.resize(in);
    tmp = in;
    for (i = 0;i < size;i++)
    {
        erode3D_step(tmp, out, neig, count);
        tmp = out;
    }

}

/* Opening/Closing 3D ------------------------------------------------------ */
void closing3D(const matrix3D<double> &in, matrix3D<double> &out, int neig,
               int count, int size)
{
    matrix3D<double> tmp;
    int i;
    tmp.resize(in);
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

void opening3D(const matrix3D<double> &in, matrix3D<double> &out, int neig,
               int count, int size)
{
    matrix3D<double> tmp;
    int i;
    tmp.resize(in);
    tmp = in;
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
