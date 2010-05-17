/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

//-----------------------------------------------------------------------------
// xmippUmatrix.cc: Unified Distance Matrix (Umatrix)
//-----------------------------------------------------------------------------

#include "umatrix.h"

/**
*   Same as () operator.
*   It takes a SOM as input and returns the calculated Umatrix
*/
void xmippUmatrix::getUmatrix(const In& in, Out& out) const
{
    operator()(in, out);
};

/**
*   Gets the Umatrix but removing a variable from the analysis.
*   It takes a SOM as input and returns the calculated Umatrix
*/

void xmippUmatrix::getUmatrix(const In& in, Out& out, const unsigned& _varOut) const
{
    In tmpSOM = in;
    for (unsigned i = 0 ; i < tmpSOM.height() ; i++)
        for (unsigned j = 0 ; j < tmpSOM.width() ; j++)
            tmpSOM.itemAtPos(SomPos(j, i))[_varOut] = 0.;
    operator()(tmpSOM, out);
}

/**
*   Gets the Umatrix but removing a list variables from the analysis.
*   It takes a SOM as input and returns the calculated Umatrix
*/

void xmippUmatrix::getUmatrix(const In& in, Out& out, const std::vector<unsigned>& _varsOut) const
{
    In tmpSOM = in;
    for (unsigned i = 0 ; i < tmpSOM.height() ; i++)
        for (unsigned j = 0 ; j < tmpSOM.width() ; j++)
            for (unsigned v = 0 ; v < _varsOut.size(); v++)
                tmpSOM.itemAtPos(SomPos(j, i))[_varsOut[v]] = 0.;
    operator()(tmpSOM, out);
}

/**
*   () operator.
*    It takes a SOM as input and returns the calculated Umatrix
*/

void xmippUmatrix::operator()(const In& in, Out& out) const
{
    bool bx, by, bz;
    xmippFeature dx, dy, dz1, dz2, dz;
    unsigned i, j;


    // Check Maps dimensions

    if ((out.height() != 2*in.height() - 1) || (out.width() != 2*in.width() - 1))
        throw std::invalid_argument("xmippUmatrix: Invalid Umatrix dimensions. (OutMap should be 2*InMap - 1)");


    // Set calibrated tag.

    if (in.calibrated())
        out.calibrated(true);

    // Calculate UMatrix

    if (in.layout() == "RECT")

        // rectangular topology
    {
        for (j = 0; j < in.height(); j++)
            for (i = 0; i < in.width(); i++)
            {

                bx = false;
                by = false;
                bz = false;
                if (i < (in.width() - 1))
                {
                    dx = eDist(in.itemAtPos(SomPos(i, j)), in.itemAtPos(SomPos(i + 1, j)));
                    bx = true;
                }
                if (j < (in.height() - 1))
                {
                    dy = eDist(in.itemAtPos(SomPos(i, j)), in.itemAtPos(SomPos(i, j + 1)));
                    by = true;
                }
                if (j < (in.height() - 1) && i < (in.width() - 1))
                {
                    dz1 = eDist(in.itemAtPos(SomPos(i, j)), in.itemAtPos(SomPos(i + 1, j + 1)));
                    dz2 = eDist(in.itemAtPos(SomPos(i, j + 1)), in.itemAtPos(SomPos(i + 1, j)));
                    bz = true;
                }

                dz = (dz1 / sqrt((xmippFeature) 2.0) + dz2 / sqrt((xmippFeature) 2.0)) / 2;

                if (bx)
                    out.itemAtPos(SomPos(2*i + 1, 2*j))[0] = dx;
                if (by)
                    out.itemAtPos(SomPos(2*i, 2*j + 1))[0] = dy;
                if (bz)
                    out.itemAtPos(SomPos(2*i + 1, 2*j + 1))[0] = dz;
            } // for i
    } // if Layout
    else
        /* hexagonal topology */
    {
        for (j = 0; j < in.height(); j++)
            for (i = 0; i < in.width(); i++)
            {

                bx = false;
                by = false;
                bz = false;
                if (i < (in.width() - 1))
                {
                    dx = eDist(in.itemAtPos(SomPos(i, j)), in.itemAtPos(SomPos(i + 1, j)));
                    bx = true;
                }

                if (j < (in.height() - 1))
                {
                    if (j % 2)
                    {
                        dy = eDist(in.itemAtPos(SomPos(i, j)), in.itemAtPos(SomPos(i, j + 1)));
                        by = true;
                    }
                    else
                    {
                        if (i > 0)
                        {
                            dy = eDist(in.itemAtPos(SomPos(i, j)), in.itemAtPos(SomPos(i - 1, j + 1)));
                            by = true;
                        }
                    }
                }


                if (j < (in.height() - 1))
                {
                    if (!(j % 2))
                    {
                        dz = eDist(in.itemAtPos(SomPos(i, j)), in.itemAtPos(SomPos(i, j + 1)));
                        bz = true;
                    }
                    else
                    {
                        if (i < (in.width() - 1))
                        {
                            dz = eDist(in.itemAtPos(SomPos(i, j)), in.itemAtPos(SomPos(i + 1, j + 1)));
                            bz = true;
                        }
                    }
                }


                if (bx)
                    out.itemAtPos(SomPos(2*i + 1, 2*j))[0] = dx;

                if (by)
                {
                    if (j % 2)
                        out.itemAtPos(SomPos(2*i, 2*j + 1))[0] = dy;
                    else
                        out.itemAtPos(SomPos(2*i - 1, 2*j + 1))[0] = dy;
                }

                if (bz)
                {
                    if (j % 2)
                        out.itemAtPos(SomPos(2*i + 1, 2*j + 1))[0] = dz;
                    else
                        out.itemAtPos(SomPos(2*i, 2*j + 1))[0] = dz;
                }


            } // for i

    } // Else


    /* Set the values corresponding to the model vectors themselves
       to medians of the surrounding values */

    if (in.layout() == "RECT")

        // Rectangular Topology
    {

        /* medians of the 4-neighborhood */
        for (j = 0; j < out.height(); j += 2)
            for (i = 0; i < out.width(); i += 2)
            {
                if (i > 0 && j > 0 && i < out.width() - 1 && j < out.height() - 1)
                {
                    // In the Middle of the Map
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = (v[1] + v[2]) / 2.0;
                }
                else if (j == 0 && i > 0 && i < out.width() - 1)
                {
                    /* in the upper edge */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[1];
                }
                else if (j == out.height() - 1 && i > 0 && i < out.width() - 1)
                {
                    /* in the lower edge */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[1];
                }
                else if (i == 0 && j > 0 && j < out.height() - 1)
                {
                    /* in the left edge */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[1];
                }
                else if (i == out.width() - 1 && j > 0 && j < out.height() - 1)
                {
                    /* in the right edge */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[1];
                }
                else if (i == 0 && j == 0)
                    /* the upper left-hand corner */
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 2.0;
                else if (i == out.width() - 1 && j == 0)
                {
                    /* the upper right-hand corner */
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 2.0;
                }
                else if (i == 0 && j == out.height() - 1)
                {
                    /* the lower left-hand corner */
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j - 1))[0]) / 2.0;
                }
                else if (i == out.width() - 1 && j == out.height() - 1)
                {
                    /* the lower right-hand corner */
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j - 1))[0]) / 2.0;
                }

                if (in.calibrated())
                    out.targetAtPos(SomPos(i, j)) = in.targetAtPos(SomPos(i / 2, j / 2));
            } // for i
    } else
        // Hexagonal Topology

    {

        for (j = 0; j < out.height(); j += 2)
            for (i = 0; i < out.width(); i += 2)
            {

                if (i > 0 && j > 0 && i < out.width() - 1 && j < out.height() - 1)
                {
                    /* in the middle of the map */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    if (!(j % 4))
                    {
                        v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    }
                    else
                    {
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i + 1, j + 1))[0]);
                    }
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = (v[2] + v[3]) / 2.0;
                }
                else if (j == 0 && i > 0 && i < out.width() - 1)
                {
                    /* in the upper edge */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = (v[1] + v[2]) / 2.0;
                }
                else if (j == out.height() - 1 && i > 0 && i < out.width() - 1)
                {
                    /* in the lower edge */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    if (!(j % 4))
                    {
                        v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    }
                    else
                    {
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                    }
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = (v[1] + v[2]) / 2.0;
                }
                else if (i == 0 && j > 0 && j < out.height() - 1)
                {
                    /* in the left edge */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    if (!(j % 4))
                    {
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                        sort(v.begin(), v.end());
                        out.itemAtPos(SomPos(i, j))[0] = v[1];
                    }
                    else
                    {
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i + 1, j + 1))[0]);
                        sort(v.begin(), v.end());
                        out.itemAtPos(SomPos(i, j))[0] = v[2];
                    }
                }
                else if (i == out.width() - 1 && j > 0 && j < out.height() - 1)
                {
                    /* in the right edge */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    if (j % 4)
                    {
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                        sort(v.begin(), v.end());
                        out.itemAtPos(SomPos(i, j))[0] = v[1];
                    }
                    else
                    {
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
                        sort(v.begin(), v.end());
                        out.itemAtPos(SomPos(i, j))[0] = v[2];
                    }
                }
                else if (i == 0 && j == 0)
                    /* the upper left-hand corner */
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 2.0;
                else if (i == out.width() - 1 && j == 0)
                {
                    /* the upper right-hand corner */
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[1];
                }
                else if (i == 0 && j == out.height() - 1)
                {
                    /* the lower left-hand corner */
                    if (!(j % 4))
                        out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i + 1, j))[0] +
                                                          out.itemAtPos(SomPos(i, j - 1))[0]) / 2.0;
                    else
                    {
                        std::vector<xmippFeature> v;
                        v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                        sort(v.begin(), v.end());
                        out.itemAtPos(SomPos(i, j))[0] = v[1];
                    }
                }
                else if (i == out.width() - 1 && j == out.height() - 1)
                {
                    /* the lower right-hand corner */
                    if (j % 4)
                        out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j))[0] +
                                                          out.itemAtPos(SomPos(i, j - 1))[0]) / 2.0;
                    else
                    {
                        std::vector<xmippFeature> v;
                        v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                        v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                        v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                        sort(v.begin(), v.end());
                        out.itemAtPos(SomPos(i, j))[0] = v[1];
                    }
                }

                if (in.calibrated())
                    out.targetAtPos(SomPos(i, j)) = in.targetAtPos(SomPos(i / 2, j / 2));
            } // for i

    } // else


    /* find the minimum and maximum values */

    xmippFeature Max = -MAXFLOAT;
    xmippFeature Min = MAXFLOAT;

    for (i = 0;i < out.width();i++)
        for (j = 0;j < out.height();j++)
        {
            if (out.itemAtPos(SomPos(i, j))[0] > Max)
                Max = out.itemAtPos(SomPos(i, j))[0];
            if (out.itemAtPos(SomPos(i, j))[0] < Min)
                Min = out.itemAtPos(SomPos(i, j))[0];
        }

    xmippFeature bw = Max - Min;

    /* scale values to [0,1] */
    for (i = 0;i < out.width();i++)
        for (j = 0;j < out.height();j++)
            out.itemAtPos(SomPos(i, j))[0] = 1.0 - (out.itemAtPos(SomPos(i, j))[0] - Min) / bw;


    // Smooth map depending of the type of smoothing


    if (smooth == MEDIAN)
        medianSmoothing(out);
    else if (smooth == AVERAGE)
        averageSmoothing(out);


} // operator ()


//-----------------------------------------------------------------------------

void xmippUmatrix::medianSmoothing(Out& out) const
{
    unsigned i, j;

    /* rectangular topology */
    if (out.layout() == "RECT")
    {
        for (j = 0;j < out.height();j++)
            for (i = 0;i < out.width();i++)
                if (i && j && (j < out.height() - 1) && (i < out.width() - 1))
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[2];
                }
                else if (i && (i < out.width() - 1) && !j)
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[2];
                }
                else if (!i && j && (j < out.height() - 1))
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[2];
                }
                else if (i && (i < out.width() - 1) && (j == out.height() - 1))
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[2];
                }
                else if (j && (j < out.height() - 1) && (i == out.width() - 1))
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[2];
                }

        std::vector<xmippFeature> v;
        v.push_back(out.itemAtPos(SomPos(1, out.height() - 1))[0]);
        v.push_back(out.itemAtPos(SomPos(0, out.height() - 1))[0]);
        v.push_back(out.itemAtPos(SomPos(0, out.height() - 2))[0]);
        sort(v.begin(), v.end());
        out.itemAtPos(SomPos(0, out.height() - 1))[0] = v[1];

        v.clear();
        v.push_back(out.itemAtPos(SomPos(out.width() - 2, out.height() - 1))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 1, out.height() - 1))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 1, out.height() - 2))[0]);
        sort(v.begin(), v.end());
        out.itemAtPos(SomPos(out.width() - 1, out.height() - 1))[0] = v[1];

        v.clear();
        v.push_back(out.itemAtPos(SomPos(out.width() - 2, 0))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 1, 0))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 1, 1))[0]);
        sort(v.begin(), v.end());
        out.itemAtPos(SomPos(out.width() - 1, 0))[0] = v[1];

        v.clear();
        v.push_back(out.itemAtPos(SomPos(1, 0))[0]);
        v.push_back(out.itemAtPos(SomPos(0, 1))[0]);
        v.push_back(out.itemAtPos(SomPos(0, 0))[0]);
        sort(v.begin(), v.end());
        out.itemAtPos(SomPos(0, 0))[0] = v[1];

    }
    else
        /* hexagonal topology */
    {
        /*else*/
        for (j = 1;j < out.height() - 1;j++)
            for (i = 1;i < out.width() - 1;i++)
                /* Non-borders */
            {
                if ((j % 4) == 1)
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[3];
                }
                else if ((j % 4) == 2)
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[3];
                }
                else if ((j % 4) == 3)
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[3];
                }
                else if ((j % 4) == 0)
                {
                    std::vector<xmippFeature> v;
                    v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                    v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
                    v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                    sort(v.begin(), v.end());
                    out.itemAtPos(SomPos(i, j))[0] = v[3];
                }
            }
        /* north border */
        j = 0;
        for (i = 1;i < out.width() - 1;i++)
        {
            std::vector<xmippFeature> v;
            v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
            v.push_back(out.itemAtPos(SomPos(i, j))[0]);
            v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
            v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
            v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
            sort(v.begin(), v.end());
            out.itemAtPos(SomPos(i, j))[0] = v[2];
        }
        /*south border*/
        j = out.height() - 1;
        for (i = 1;i < out.width() - 1;i++)
        {
            if ((j % 4) == 1)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
            else if ((j % 4) == 2)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
            else if ((j % 4) == 3)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
            else if ((j % 4) == 0)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
        }
        /*east border*/
        i = out.width() - 1;
        for (j = 1;j < out.height() - 1;j++)
        {
            if ((j % 4) == 1)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
            else if ((j % 4) == 2)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
            else if ((j % 4) == 3)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
            else if ((j % 4) == 0)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i - 1, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i - 1, j + 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[3];
            }
        }
        i = 0;
        for (j = 1;j < out.height() - 1;j++)

            /*west border*/
        {
            if ((j % 4) == 1)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
            else if ((j % 4) == 2)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j + 1))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[3];
            }
            else if ((j % 4) == 3)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j + 1))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
            else if ((j % 4) == 0)
            {
                std::vector<xmippFeature> v;
                v.push_back(out.itemAtPos(SomPos(i, j - 1))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i + 1, j))[0]);
                v.push_back(out.itemAtPos(SomPos(i, j + 1))[0]);
                sort(v.begin(), v.end());
                out.itemAtPos(SomPos(i, j))[0] = v[2];
            }
        }


        /*Corners*/

        std::vector<xmippFeature> v;
        v.push_back(out.itemAtPos(SomPos(1, 0))[0]);
        v.push_back(out.itemAtPos(SomPos(0, 0))[0]);
        v.push_back(out.itemAtPos(SomPos(0, 1))[0]);
        sort(v.begin(), v.end());
        out.itemAtPos(SomPos(0, 0))[0] = v[1];

        v.clear();
        v.push_back(out.itemAtPos(SomPos(out.width() - 1, 0))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 1, 1))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 2, 0))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 2, 1))[0]);
        sort(v.begin(), v.end());
        out.itemAtPos(SomPos(out.width() - 1, 0))[0] = v[2];

        v.clear();
        v.push_back(out.itemAtPos(SomPos(out.width() - 1, out.height() - 1))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 1, out.height() - 2))[0]);
        v.push_back(out.itemAtPos(SomPos(out.width() - 2, out.height() - 1))[0]);
        sort(v.begin(), v.end());
        out.itemAtPos(SomPos(out.width() - 1, out.height() - 1))[0] = v[1];

        v.clear();
        v.push_back(out.itemAtPos(SomPos(0, out.height() - 1))[0]);
        v.push_back(out.itemAtPos(SomPos(1, out.height() - 1))[0]);
        v.push_back(out.itemAtPos(SomPos(0, out.height() - 2))[0]);
        sort(v.begin(), v.end());
        out.itemAtPos(SomPos(0, out.height() - 1))[0] = v[1];
    }


} // median Smoothing


//-----------------------------------------------------------------------------

void xmippUmatrix::averageSmoothing(Out& out) const
{

    unsigned i, j;

    /* rectangular topology */
    if (out.layout() == "RECT")
    {
        for (j = 0;j < out.height();j++)
            for (i = 0;i < out.width();i++)
                if (i && j && (j < out.height() - 1) && (i < out.width() - 1))
                {
                    /* Non borders */
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 5.0;
                }
                else if (i && (i < out.width() - 1) && !j)
                {
                    /* West brdr*/
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 4.0;
                }
                else if (!i && j && (j < out.height() - 1))
                {
                    /*north*/
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 4.0;
                }
                else if (i && (i < out.width() - 1) && (j == out.height() - 1))
                {
                    /* south */
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j))[0]) / 4.0;
                }
                else if (j && (j < out.height() - 1) && (i == out.width() - 1))
                {
                    /* east*/
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 4.0;
                }
        /*corners*/

        out.itemAtPos(SomPos(0, out.height() - 1))[0] = (out.itemAtPos(SomPos(1, out.height() - 1))[0] +
                out.itemAtPos(SomPos(0, out.height() - 1))[0] +
                out.itemAtPos(SomPos(0, out.height() - 2))[0]) / 3.0;

        out.itemAtPos(SomPos(out.width() - 1, out.height() - 1))[0] = (out.itemAtPos(SomPos(out.width() - 2, out.height() - 1))[0] +
                out.itemAtPos(SomPos(out.width() - 1, out.height() - 1))[0] +
                out.itemAtPos(SomPos(out.width() - 1, out.height() - 2))[0]) / 3.0;

        out.itemAtPos(SomPos(out.width() - 1, 0))[0] = (out.itemAtPos(SomPos(out.width() - 2, 0))[0] +
                out.itemAtPos(SomPos(out.width() - 1, 0))[0] +
                out.itemAtPos(SomPos(out.width() - 1, 1))[0]) / 3.0;


        out.itemAtPos(SomPos(0, 0))[0] = (out.itemAtPos(SomPos(1, 0))[0] +
                                          out.itemAtPos(SomPos(0, 1))[0] +
                                          out.itemAtPos(SomPos(0, 0))[0]) / 3.0;

    }
    else
        /* hexagonal topology */
    {
        /*else*/
        for (j = 1;j < out.height() - 1;j++)
            for (i = 1;i < out.width() - 1;i++)
                /* Non-borders */
            {
                if ((j % 4) == 1)
                {
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j + 1))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 7.0;
                }
                else if ((j % 4) == 2)
                {
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j + 1))[0]) / 7.0;
                }
                else if ((j % 4) == 3)
                {
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j + 1))[0]) / 7.0;
                }
                else if ((j % 4) == 0)
                {
                    out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i, j - 1))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j))[0] +
                                                      out.itemAtPos(SomPos(i, j))[0] +
                                                      out.itemAtPos(SomPos(i + 1, j))[0] +
                                                      out.itemAtPos(SomPos(i - 1, j + 1))[0] +
                                                      out.itemAtPos(SomPos(i, j + 1))[0]) / 7.0;
                }
            }
        /* north border */
        j = 0;
        for (i = 1;i < out.width() - 1;i++)
            out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j))[0] +
                                              out.itemAtPos(SomPos(i, j))[0] +
                                              out.itemAtPos(SomPos(i + 1, j))[0] +
                                              out.itemAtPos(SomPos(i - 1, j + 1))[0] +
                                              out.itemAtPos(SomPos(i, j + 1))[0]) / 5.0;
        /*south border*/
        j = out.height() - 1;
        for (i = 1;i < out.width() - 1;i++)
        {
            if ((j % 4) == 1)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j))[0]) / 5.0;
            }
            else if ((j % 4) == 2)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j))[0]) / 5.0;
            }
            else if ((j % 4) == 3)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j))[0]) / 5.0;
            }
            else if ((j % 4) == 0)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j))[0]) / 5.0;
            }
        }
        /*east border*/
        i = out.width() - 1;
        for (j = 1;j < out.height() - 1;j++)
        {
            if ((j % 4) == 1)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j + 1))[0] +
                                                  out.itemAtPos(SomPos(i, j + 1))[0]) / 5.0;
            }
            else if ((j % 4) == 2)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i, j + 1))[0]) / 4.0;
            }
            else if ((j % 4) == 3)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i, j + 1))[0]) / 5.0;
            }
            else if ((j % 4) == 0)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i - 1, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i - 1, j + 1))[0] +
                                                  out.itemAtPos(SomPos(i, j + 1))[0]) / 6.0;
            }
        }
        i = 0;
        for (j = 1;j < out.height() - 1;j++)

            /*west border*/
        {
            if ((j % 4) == 1)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j + 1))[0]) / 5.0;
            }
            else if ((j % 4) == 2)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j + 1))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j + 1))[0]) / 6.0;
            }
            else if ((j % 4) == 3)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j + 1))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j + 1))[0]) / 5.0;
            }
            else if ((j % 4) == 0)
            {
                out.itemAtPos(SomPos(i, j))[0] = (out.itemAtPos(SomPos(i, j - 1))[0] +
                                                  out.itemAtPos(SomPos(i, j))[0] +
                                                  out.itemAtPos(SomPos(i + 1, j))[0] +
                                                  out.itemAtPos(SomPos(i, j + 1))[0]) / 4.0;
            }
        }


        /*Corners*/

        out.itemAtPos(SomPos(0, 0))[0] = (out.itemAtPos(SomPos(1, 0))[0] +
                                          out.itemAtPos(SomPos(0, 0))[0] +
                                          out.itemAtPos(SomPos(0, 1))[0]) / 3.0;

        out.itemAtPos(SomPos(out.width() - 1, 0))[0] = (out.itemAtPos(SomPos(out.width() - 1, 0))[0] +
                out.itemAtPos(SomPos(out.width() - 1, 1))[0] +
                out.itemAtPos(SomPos(out.width() - 2, 0))[0] +
                out.itemAtPos(SomPos(out.width() - 2, 1))[0]) / 4.0;


        out.itemAtPos(SomPos(out.width() - 1, out.height() - 1))[0] = (out.itemAtPos(SomPos(out.width() - 1, out.height() - 1))[0] +
                out.itemAtPos(SomPos(out.width() - 1, out.height() - 2))[0] +
                out.itemAtPos(SomPos(out.width() - 2, out.height() - 1))[0]) / 3.0;


        out.itemAtPos(SomPos(0, out.height() - 1))[0] = (out.itemAtPos(SomPos(0, out.height() - 1))[0] +
                out.itemAtPos(SomPos(1, out.height() - 1))[0] +
                out.itemAtPos(SomPos(0, out.height() - 2))[0]) / 3.0;
    }


} // average Smoothing
