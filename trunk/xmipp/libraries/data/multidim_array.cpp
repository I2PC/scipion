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

#include "multidim_array.h"


// Show a complex array ---------------------------------------------------
std::ostream& operator<<(std::ostream& ostrm,
                         const MultidimArray< std::complex<double> >& v)
{
    if (v.xdim == 0)
        ostrm << "NULL MultidimArray\n";
    else
        ostrm << std::endl;

    for (int l = 0; l < NSIZE(v); l++)
    {
        if (NSIZE(v)>1)
            ostrm << "Image No. " << l << std::endl;
        for (int k = STARTINGZ(v); k <= FINISHINGZ(v); k++)
        {
            if (ZSIZE(v)>1)
                ostrm << "Slice No. " << k << std::endl;
            for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
            {
                for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
                    ostrm << A3D_ELEM(v, k, i, j) << ' ';
                ostrm << std::endl;
            }
        }
    }

    return ostrm;
}

/** Force positive -------------------------------------------------------- */
void forcePositive(MultidimArray<double> &V)
{
    bool negativeRemaining;

    if (V.getDim()==2) // IMAGE
    {
        do
        {
            negativeRemaining=false;

            FOR_ALL_ELEMENTS_IN_ARRAY2D(V)
            if (V(i, j)<=0)
            {
                std::vector<double> neighbours;
                for (int ii=-2; ii<=2; ii++)
                {
                    int iii=i+ii;
                    if (iii<0 || iii>=YSIZE(V))
                        continue;
                    for (int jj=-2; jj<=2; jj++)
                    {
                        int jjj=j+jj;
                        if (jjj<0 || jjj>=XSIZE(V))
                            continue;
                        double val=V(iii,jjj);
                        if (val>0)
                            neighbours.push_back(val);
                    }
                }
                int N=neighbours.size();
                if (N==0)
                    negativeRemaining=true;
                else
                {
                    std::sort(neighbours.begin(),neighbours.end());
                    if (N%2==0)
                        V(i,j)=0.5*(neighbours[N/2-1]+neighbours[N/2]);
                    else
                        V(i,j)=neighbours[N/2];
                }
            }
        }
        while (negativeRemaining);
    }
    else if (V.getDim()==3) // VOLUME
    {
        do
        {
            negativeRemaining=false;

            FOR_ALL_ELEMENTS_IN_ARRAY3D(V)
            if (V(k, i, j)<=0)
            {
                std::vector<double> neighbours;
                for (int kk=-2; kk<=2; kk++)
                {
                    int kkk=k+kk;
                    if (kkk<0 || kkk>=ZSIZE(V))
                        continue;
                    for (int ii=-2; ii<=2; ii++)
                    {
                        int iii=i+ii;
                        if (iii<0 || iii>=YSIZE(V))
                            continue;
                        for (int jj=-2; jj<=2; jj++)
                        {
                            int jjj=j+jj;
                            if (jjj<0 || jjj>=XSIZE(V))
                                continue;
                            double val=V(kkk,iii,jjj);
                            if (val>0)
                                neighbours.push_back(val);
                        }
                    }
                    int N=neighbours.size();
                    if (N==0)
                        negativeRemaining=true;
                    else
                    {
                        std::sort(neighbours.begin(),neighbours.end());
                        if (N%2==0)
                            V(k,i,j)=0.5*(neighbours[N/2-1]+
                                          neighbours[N/2]);
                        else
                            V(k,i,j)=neighbours[N/2];
                    }
                }
            }
        }
        while (negativeRemaining);
    }
    else
    {
        REPORT_ERROR(ERR_NOT_IMPLEMENTED,"");
    }
}
