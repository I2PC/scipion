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

#include "filters.h"
#include "fftw.h"
#include "polar.h"
#include <list>

/* Substract background ---------------------------------------------------- */
void substract_background_plane(Matrix2D<double> &I)
{

    I.checkDimension(2);

    Matrix2D<double> A(3, 3);
    Matrix1D<double> x(3), b(3);

    // Solve the plane 'x'
    A.initZeros();
    b.initZeros();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
    {
        A(0, 0) += j * j;
        A(0, 1) += j * i;
        A(0, 2) += j;
        A(1, 1) += i * i;
        A(1, 2) += i;
        A(2, 2) += 1;
        b(0)    += j * MAT_ELEM(I, i, j);
        b(1)    += i * MAT_ELEM(I, i, j);
        b(2)    += MAT_ELEM(I, i, j);
    }
    A(1, 0)  = A(0, 1);
    A(2, 0)  = A(0, 2);
    A(2, 1)  = A(1, 2);
    solve(A, b, x);

    // Now substract the plane
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I)
        MAT_ELEM(I, i, j) -= x(0) * i + x(1) * j + x(2);
}

/* Substract background ---------------------------------------------------- */
void substract_background_rolling_ball(Matrix2D<double> &I, int radius)
{

    I.checkDimension(2);

    // Build the ball
    int arcTrimPer;
    int shrinkFactor;
    if (radius<=10) {
        shrinkFactor = 1;
        arcTrimPer = 24; // trim 24% in x and y
    } else if (radius<=30) {
        shrinkFactor = 2;
        arcTrimPer = 24; // trim 24% in x and y
    } else if (radius<=100) {
        shrinkFactor = 4;
        arcTrimPer = 32; // trim 32% in x and y
    } else {
        shrinkFactor = 8;
        arcTrimPer = 40; // trim 40% in x and y
    }
    
    double smallballradius = radius/shrinkFactor;
    if (smallballradius<1) smallballradius = 1;
    double r2 = smallballradius*smallballradius;
    int xtrim = (int)(arcTrimPer*smallballradius)/100; // only use a patch of the rolling ball
    int halfWidth = ROUND(smallballradius - xtrim);
    int ballWidth = 2*halfWidth+1;
    Matrix2D<double> ball(ballWidth,ballWidth);
    ball.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(ball)
    {
        double temp=r2-i*i-j*j;
        ball(i,j)=temp>0. ? sqrt(temp):0.;
    }

    // Shrink the image: each point in the shrinked image is the
    // minimum of its neighbourhood
    int sXdim = (XSIZE(I)+shrinkFactor-1)/shrinkFactor;
    int sYdim = (YSIZE(I)+shrinkFactor-1)/shrinkFactor;
    Matrix2D<double> shrinkI(sYdim,sXdim);
    shrinkI.setXmippOrigin();
    for (int ySmall=0; ySmall<sYdim; ySmall++) {
        for (int xSmall=0; xSmall<sXdim; xSmall++) {
            double minVal = 1e38;
            for (int j=0, y=shrinkFactor*ySmall; j<shrinkFactor&&y<YSIZE(I); j++, y++)
                for (int k=0, x=shrinkFactor*xSmall; k<shrinkFactor&&x<XSIZE(I); k++, x++) {
                    double thispixel = DIRECT_MAT_ELEM(I,y,x);
                    if (thispixel<minVal)
                        minVal = thispixel;
                }
            DIRECT_MAT_ELEM(shrinkI,ySmall,xSmall) = minVal;
        }
    }
    
    // Now roll the ball
    radius=ballWidth/2;
    Matrix2D<double> Irolled;
    Irolled.resize(shrinkI);
    Irolled.initConstant(-500);
    for (int yb=-radius; yb<YSIZE(shrinkI)+radius; yb++) {
        // Limits of the ball
        int y0 = yb-radius;
        if (y0 < 0) y0 = 0;
        int y0b = y0-yb+radius; //y coordinate in the ball corresponding to y0
        int yF = yb+radius;
        if (yF>=YSIZE(shrinkI)) yF = YSIZE(shrinkI)-1;

        for (int xb=-radius; xb<XSIZE(shrinkI)+radius; xb++) {
            // Limits of the ball
            int x0 = xb-radius;
            if (x0 < 0) x0 = 0;
            int x0b = x0-xb+radius;
            int xF = xb+radius;
            if (xF>=XSIZE(shrinkI)) xF = XSIZE(shrinkI)-1;

            double z = 1e38;
            for (int yp=y0, ybp=y0b; yp<=yF; yp++,ybp++)
                for (int xp=x0, xbp=x0b; xp<=xF; xp++, xbp++) {
                    double zReduced=DIRECT_MAT_ELEM(shrinkI,yp,xp) -
                        DIRECT_MAT_ELEM(ball,ybp,xbp);
                    if (z > zReduced) z = zReduced;
                }
            for (int yp=y0, ybp=y0b; yp<=yF; yp++,ybp++)
                for (int xp=x0, xbp=x0b; xp<=xF; xp++, xbp++) {
                    double zMin = z + DIRECT_MAT_ELEM(ball,ybp,xbp);
                    if (DIRECT_MAT_ELEM(Irolled,yp,xp) < zMin)
                        DIRECT_MAT_ELEM(Irolled,yp,xp) = zMin;
                }
        }
    }
    
    // Now rescale the background
    Matrix2D<double> bgEnlarged;
    scaleToSize(1, Irolled, bgEnlarged, YSIZE(I),XSIZE(I));
    bgEnlarged.copyShape(I);
    I-=bgEnlarged;
}

/* Contranst enhancement --------------------------------------------------- */
void contrast_enhancement(Image<double> *I)
{
    (*I)().rangeAdjust(0, 255);
}

/* Region growing for images ----------------------------------------------- */
void region_growing2D(const Matrix2D<double> &I_in, Matrix2D<double> &I_out,
                      int i, int j,
                      float stop_colour, float filling_colour, bool less, int neighbourhood)
{

    I_in.checkDimension(2);

    std::list<int> iNeighbours;   /* A list for neighbour pixels */
    int iCurrenti, iCurrentj;     /* Coordinates of the current pixel considered */

    /* First task is copying the input image into the output one */
    I_out = I_in;

    /**** Then the region growing is done ****/
    /* Insert at the beginning of the list the seed coordinates */
    iNeighbours.push_front(j);
    iNeighbours.push_front(i);

    /* Fill the seed coordinates */
    MAT_ELEM(I_out, i, j) = filling_colour;

    while (!iNeighbours.empty())
    {
        Matrix1D<double> r(2);

        /* Take the current pixel to explore */
        iCurrenti = iNeighbours.front();
        iNeighbours.pop_front();
        iCurrentj = iNeighbours.front();
        iNeighbours.pop_front();

#define CHECK_POINT(i,j) \
    XX(r)=j; YY(r)=i; \
    if (!I_out.outside(r))  { \
        if (MAT_ELEM(I_out,i,j)!=filling_colour) \
            if ((less && MAT_ELEM (I_out,i,j) < stop_colour) || \
                (!less && MAT_ELEM (I_out,i,j) > stop_colour)) { \
                MAT_ELEM (I_out,i,j)=filling_colour; \
                iNeighbours.push_front(j); \
                iNeighbours.push_front(i); \
            } \
    }

        /* Make the exploration of the pixel's neighbours */
        CHECK_POINT(iCurrenti  , iCurrentj - 1);
        CHECK_POINT(iCurrenti  , iCurrentj + 1);
        CHECK_POINT(iCurrenti - 1, iCurrentj);
        CHECK_POINT(iCurrenti + 1, iCurrentj);
        if (neighbourhood == 8)
        {
            CHECK_POINT(iCurrenti - 1, iCurrentj - 1);
            CHECK_POINT(iCurrenti - 1, iCurrentj + 1);
            CHECK_POINT(iCurrenti + 1, iCurrentj - 1);
            CHECK_POINT(iCurrenti + 1, iCurrentj + 1);
        }
    }
}

/* Region growing for volumes ----------------------------------------------- */
void region_growing3D(const Matrix3D<double> &V_in, Matrix3D<double> &V_out,
                      int k, int i, int j,
                      float stop_colour, float filling_colour, bool less)
{
    V_in.checkDimension(3);

    std::list<int> iNeighbours;       /* A list for neighbour voxels */
    int iCurrentk, iCurrenti, iCurrentj;     /* Coordinates of the current voxel considered */

    /* First task is copying the input volume into the output one */
    V_out = V_in;

    /**** Then the region growing is done in output volume ****/
    /* Insert at the beginning of the list the seed coordinates */
    iNeighbours.push_front(j);
    iNeighbours.push_front(i);
    iNeighbours.push_front(k);

    /* Fill the seed coordinates */
    VOL_ELEM(V_out, k, i, j) = filling_colour;

    while (!iNeighbours.empty())
    {
        Matrix1D<double> r(3);

        /* Take the current pixel to explore */
        iCurrentk = iNeighbours.front();
        iNeighbours.pop_front();
        iCurrenti = iNeighbours.front();
        iNeighbours.pop_front();
        iCurrentj = iNeighbours.front();
        iNeighbours.pop_front();

        /* a macro for doing exploration of a voxel. If the voxel has a value
        lower than stop_colour, its filled with filling colour and added to the
        list for exploring its neighbours */
#define CHECK_POINT_3D(k,i,j) \
    XX(r)=j; YY(r)=i; ZZ(r)=k; \
    if (!V_out.outside(r))  { \
        if (VOL_ELEM(V_out,k,i,j)!=filling_colour) \
            if ((less && VOL_ELEM (V_out,k,i,j) < stop_colour)|| \
                (!less &&VOL_ELEM (V_out,k,i,j) > stop_colour)) { \
                VOL_ELEM (V_out,k,i,j)=filling_colour; \
                iNeighbours.push_front(j); \
                iNeighbours.push_front(i); \
                iNeighbours.push_front(k); \
            } \
    }

        /* Make the exploration of the pixel�s neighbours */
        CHECK_POINT_3D(iCurrentk  , iCurrenti  , iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk  , iCurrenti  , iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk  , iCurrenti - 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk  , iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk  , iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk  , iCurrenti + 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk  , iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk  , iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti  , iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti  , iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti  , iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk - 1, iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti  , iCurrentj);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti  , iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti  , iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti - 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj + 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj - 1);
        CHECK_POINT_3D(iCurrentk + 1, iCurrenti + 1, iCurrentj + 1);
    }
}

void distance_transform(const Matrix2D<int> &in,
    Matrix2D<int> &out, bool wrap)
{
    std::list<int> toExplore;   /* A list of points to explore */
    
    in.checkDimension(2);

    out.resize(in);
    out.initConstant(XSIZE(in)+YSIZE(in));

#define CHECK_POINT_DIST(i,j,proposedDistance) \
    { \
        int ip=i; \
        int jp=j; \
        if (wrap) \
        { \
            ip=intWRAP(ip,STARTINGY(out),FINISHINGY(out));\
            jp=intWRAP(jp,STARTINGX(out),FINISHINGX(out));\
        } \
        if (ip>=STARTINGY(out) && ip<=FINISHINGY(out) && \
            jp>=STARTINGX(out) && jp<=FINISHINGX(out)) \
            if (out(ip,jp)>proposedDistance) \
            { \
                out(ip,jp)=proposedDistance; \
                toExplore.push_back(ip); \
                toExplore.push_back(jp); \
                toExplore.push_back(proposedDistance); \
            } \
    }

    // Look for all elements in the binary mask and set the corresponding
    // distance to 0
    FOR_ALL_ELEMENTS_IN_MATRIX2D(in)
        if (in(i,j))
        {
            out(i,j)=0;
            CHECK_POINT_DIST(i-1,j,1);
            CHECK_POINT_DIST(i+1,j,1);
            CHECK_POINT_DIST(i,j-1,1);
            CHECK_POINT_DIST(i,j+1,1);
        }

    while (!toExplore.empty())
    {
        int i=toExplore.front(); toExplore.pop_front();
        int j=toExplore.front(); toExplore.pop_front();
        int proposedDistance=toExplore.front(); toExplore.pop_front();
        
        if (proposedDistance==out(i,j))
        {
            // If this is the current distance (i.e., it has not
            // been supersceded by a later value
            CHECK_POINT_DIST(i-1,j,proposedDistance+1);
            CHECK_POINT_DIST(i+1,j,proposedDistance+1);
            CHECK_POINT_DIST(i,j-1,proposedDistance+1);
            CHECK_POINT_DIST(i,j+1,proposedDistance+1);
        }
    }
}

/* Label image ------------------------------------------------------------ */
int label_image2D(const Matrix2D<double> &I, Matrix2D<double> &label,
                  int neighbourhood)
{
    I.checkDimension(2);

    label = I;
    int colour = 32000;
    bool found;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
    {
        if (label(i, j) != 1)
            continue;
        region_growing2D(label, label, i, j, 0, colour, false, neighbourhood);
        colour++;
    }
    FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
    if (label(i, j) != 0)
        label(i, j) = label(i, j) - 31999;
    return colour -32000;
}

/* Label volume ------------------------------------------------------------ */
int label_image3D(const Matrix3D<double> &V, Matrix3D<double> &label)
{
    V.checkDimension(2);

    label = V;
    int colour = 32000;
    bool found;
    FOR_ALL_ELEMENTS_IN_MATRIX3D(label)
    {
        if (label(k, i, j) != 1)
            continue;
        region_growing3D(label, label, k, i, j, 0, colour, false);
        colour++;
    }
    FOR_ALL_ELEMENTS_IN_MATRIX3D(label)
    if (label(k, i, j) != 0)
        label(k, i, j) = label(k, i, j) - 31999;
    return colour -32000;
}

/* Remove small components ------------------------------------------------- */
void remove_small_components(Matrix2D<double> &I, int size,
                             int neighbourhood)
{
    I.checkDimension(2);

    Matrix2D<double> label;
    int imax = label_image2D(I, label, neighbourhood);
    Matrix1D<int> nlabel(imax + 1);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(label) nlabel((int)(label(i, j)))++;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
    if (nlabel((int)(label(i, j))) < size)
        I(i, j) = 0;
}

/* Keep biggest component -------------------------------------------------- */
void keep_biggest_component(Matrix2D<double> &I, double percentage,
                            int neighbourhood)
{
    I.checkDimension(2);

    Matrix2D<double> label;
    int imax = label_image2D(I, label, neighbourhood);
    Matrix1D<int> nlabel(imax + 1);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
    {
        int idx = (int)(label(i, j));
        if (idx == 0)
            continue;
        nlabel(idx)++;
    }
    Matrix1D<int> best = nlabel.indexSort();
    best -= 1;
    int nbest = XSIZE(best) - 1;
    double total = nlabel.sum();
    double explained = nlabel(best(nbest));
    while (explained < percentage*total)
    {
        nbest--;
        explained += nlabel(best(nbest));
    }

    FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
    {
        bool among_the_best = false;
        for (int n = nbest; n < imax + 1; n++)
            if (label(i, j) == best(n))
            {
                among_the_best = true;
                break;
            }
        if (!among_the_best)
            I(i, j) = 0;
    }
}

/* Fill object ------------------------------------------------------------- */
void fill_binary_object(Matrix2D<double> &I, int neighbourhood)
{
    I.checkDimension(2);

    Matrix2D<double> label;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(I) I(i, j) = 1 - I(i, j);
    int imax = label_image2D(I, label, neighbourhood);
    double l0 = label(STARTINGY(I), STARTINGX(I));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(label)
    if (label(i, j) == l0)
        I(i, j) = 0;
    else
        I(i, j) = 1;
}

/* Otsu Segmentation ------------------------------------------------------- */
void OtsuSegmentation(Matrix3D<double> &V)
{
    V.checkDimension(3);

    // Compute the probability density function
    histogram1D hist;
    hist.clear();
    compute_hist(V,hist,200);
    hist/=hist.sum();

    // Compute the cumulative 0th and 1st order moments
    double x;
    Matrix1D<double> mom0, mom1;
    mom0.initZeros(XSIZE(hist));
    mom1.initZeros(XSIZE(hist));
    mom0(0)=hist(0);
    hist.index2val(0,x); mom1(0)=hist(0)*x;
    for (int i=1; i<XSIZE(mom0); i++)
    {
        mom0(i)=mom0(i-1)+hist(i);
        hist.index2val(i,x); mom1(i)=mom1(i-1)+hist(i)*x;
    }

    // Maximize sigma2B
    double bestSigma2B=-1;
    int ibestSigma2B=-1;
    for (int i=0; i<XSIZE(hist)-1; i++)
    {
        double w1=mom0(i);
        double w2=1-mom0(i);
        double mu1=mom1(i);
        double mu2=mom1(XSIZE(mom1)-1)-mom1(i);
        double sigma2B=w1*w2*(mu1-mu2)*(mu1-mu2);
        if (sigma2B>bestSigma2B)
        {
             bestSigma2B=sigma2B;
             ibestSigma2B=i;
        }
    }

    hist.index2val(ibestSigma2B,x);
    V.binarize(x); 
}

/* Entropy Segmentation ---------------------------------------------------- */
void EntropySegmentation(Matrix3D<double> &V)
{
    V.checkDimension(3);

    // Compute the probability density function
    histogram1D hist;
    hist.clear();
    compute_hist(V,hist,200);
    hist/=hist.sum();

    // Compute the cumulative 0th and 1st order moments
    double x;
    Matrix1D<double> mom0;
    mom0.initZeros(XSIZE(hist));
    mom0(0)=hist(0);
    for (int i=1; i<XSIZE(mom0); i++)
        mom0(i)=mom0(i-1)+hist(i);

    // Entropy for black and white parts of the histogram
    const double epsilon = 1e-15;
    Matrix1D<double> h1, h2;
    h1.initZeros(XSIZE(hist));
    h2.initZeros(XSIZE(hist));
    for (int i=0; i<XSIZE(hist); i++) {
        // Entropy h1
        double w1=mom0(i);
        if (w1>epsilon)
            for (int ii=0; ii<=i; ii++)
                if (hist(ii)>epsilon)
                {
                    double aux=hist(ii)/w1;
                    h1(i) -= aux*log10(aux);
                }

        // Entropy h2
        double w2=1-mom0(i);
        if (w2>epsilon)
            for (int ii=i+1; ii<XSIZE(hist); ii++)
                if (hist(ii)>epsilon)
                {
                    double aux=hist(ii)/w2;
                    h2(i) -= aux*log10(aux);
                }
    }

    // Find histogram index with maximum entropy
    double Hmax=h1(0)+h2(0);
    int iHmax=0;
    for (int i=1; i<XSIZE(hist)-1; i++) {
        double H = h1(i)+h2(i);
        if (H > Hmax) {
            Hmax = H;
            iHmax = i;
        }
    }

    hist.index2val(iHmax,x);
    V.binarize(x); 
}

/* Otsu+Entropy Segmentation ----------------------------------------------- */
void EntropyOtsuSegmentation(Matrix3D<double> &V, double percentil)
{
    V.checkDimension(3);

    // Compute the probability density function
    histogram1D hist;
    hist.clear();
    compute_hist(V,hist,200);
    hist/=hist.sum();

    // Compute the cumulative 0th and 1st order moments
    double x;
    Matrix1D<double> mom0,mom1;
    mom0.initZeros(XSIZE(hist));
    mom1.initZeros(XSIZE(hist));
    mom0(0)=hist(0);
    hist.index2val(0,x); mom1(0)=hist(0)*x;
    for (int i=1; i<XSIZE(mom0); i++)
    {
        mom0(i)=mom0(i-1)+hist(i);
        hist.index2val(i,x); mom1(i)=mom1(i-1)+hist(i)*x;
    }

    // Entropy for black and white parts of the histogram
    const double epsilon = 1e-15;
    Matrix1D<double> h1, h2;
    h1.initZeros(XSIZE(hist));
    h2.initZeros(XSIZE(hist));
    for (int i=0; i<XSIZE(hist); i++) {
        // Entropy h1
        double w1=mom0(i);
        if (w1>epsilon)
            for (int ii=0; ii<=i; ii++)
                if (hist(ii)>epsilon)
                {
                    double aux=hist(ii)/w1;
                    h1(i) -= aux*log10(aux);
                }

        // Entropy h2
        double w2=1-mom0(i);
        if (w2>epsilon)
            for (int ii=i+1; ii<XSIZE(hist); ii++)
                if (hist(ii)>epsilon)
                {
                    double aux=hist(ii)/w2;
                    h2(i) -= aux*log10(aux);
                }
    }

    // Compute sigma2B and H
    Matrix1D<double> sigma2B, H, HSigma2B;
    sigma2B.initZeros(XSIZE(hist)-1);
    H.initZeros(XSIZE(hist)-1);
    HSigma2B.initZeros(XSIZE(hist)-1);
    for (int i=0; i<XSIZE(hist)-1; i++)
    {
        double w1=mom0(i);
        double w2=1-mom0(i);
        double mu1=mom1(i);
        double mu2=mom1(XSIZE(mom1)-1)-mom1(i);
        sigma2B(i)=w1*w2*(mu1-mu2)*(mu1-mu2);
        H(i) = h1(i)+h2(i);
        HSigma2B(i)=-log10(sigma2B(i))/H(i);
           // The logic behind this expression is
           // Otsu:    max sigma2B -> max log10(sigma2B) -> min -log10(sigma2B)
           // Entropy: max H       -> max H              -> min 1/H
    }
    
    // Sort HSigma2B and take a given percentage of it
    Matrix1D<double> HSigma2Bsorted=HSigma2B.sort();
    int iTh=ROUND(XSIZE(HSigma2B)*percentil);
    double threshold=HSigma2Bsorted(iTh);
    
    // Find the first value within HSigma2B falling below this threshold
    iTh=0;
    while (HSigma2B(iTh)>threshold) iTh++;
    iTh--;
    
    hist.index2val(iTh,x);
    V.binarize(x); 
}

/* Best shift -------------------------------------------------------------- */
void best_shift(const Matrix2D<double> &I1, const Matrix2D<double> &I2,
                double &shiftX, double &shiftY, const Matrix2D<int> *mask)
{
    I1.checkDimension(2);
    I2.checkDimension(2);

    int              imax, jmax, i_actual, j_actual;
    double           max, xmax, ymax, sumcorr, avecorr, stdcorr, dummy;
    float            xshift, yshift, shift;
    int              n_max = -1;
    bool             neighbourhood = true;
    Matrix2D<double> Mcorr;

    correlation_matrix(I1, I2, Mcorr);

    /*
      Warning: for masks with a small number of non-zero pixels, this routine is NOT reliable...
      Anyway, maybe using a mask is not a good idea at al...
     */

    // Adjust statistics within shiftmask to average 0 and stddev 1
    if (mask != NULL)
    {
        if ((*mask).sum() < 2)
        {
            shiftX = shiftY = 0.;
            return;
        }
        else
        {
            computeStats_within_binary_mask(*mask, Mcorr, dummy, dummy, avecorr, stdcorr);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(Mcorr)
            if (MAT_ELEM(*mask, i, j))
                MAT_ELEM(Mcorr, i, j) = (MAT_ELEM(Mcorr, i, j) - avecorr) / stdcorr;
            else
                MAT_ELEM(Mcorr, i, j) = 0.;
        }
    }
    else
        Mcorr.statisticsAdjust(0, 1);
    Mcorr.maxIndex(imax, jmax);
    max = MAT_ELEM(Mcorr, imax, jmax);

    while (neighbourhood)
    {
        n_max ++;
        for (int i = -n_max; i <= n_max; i++)
            for (int j = -n_max; j <= n_max; j++)
            {
                i_actual = i + imax;
                j_actual = j + jmax;
                if (i_actual < STARTINGY(Mcorr)  || j_actual < STARTINGX(Mcorr) ||
                    i_actual > FINISHINGY(Mcorr) || j_actual > FINISHINGX(Mcorr))
                    neighbourhood = false;
                else if (max / 1.414 > MAT_ELEM(Mcorr, i_actual, j_actual))
                    neighbourhood = false;
            }
    }

    // We have the neighbourhood => looking for the gravity centre
    xmax = ymax = sumcorr = 0.;
    for (int i = -n_max; i <= n_max; i++)
        for (int j = -n_max; j <= n_max; j++)
        {
            i_actual = i + imax;
            j_actual = j + jmax;
            if (i_actual >= STARTINGY(Mcorr)  && j_actual >= STARTINGX(Mcorr) &&
                i_actual <= FINISHINGY(Mcorr) && j_actual <= FINISHINGX(Mcorr))
            {
                ymax += i_actual * MAT_ELEM(Mcorr, i_actual, j_actual);
                xmax += j_actual * MAT_ELEM(Mcorr, i_actual, j_actual);
                sumcorr += MAT_ELEM(Mcorr, i_actual, j_actual);
            }
        }
    shiftX = xmax / sumcorr;
    shiftY = ymax / sumcorr;
}

/* Best non-wrapping shift ------------------------------------------------- */
//#define DEBUG
void best_nonwrapping_shift(const Matrix2D<double> &I1,
    const Matrix2D<double> &I2, double &shiftX, double &shiftY)
{
    I1.checkDimension(2);
    I2.checkDimension(2);

    best_shift(I1, I2, shiftX, shiftY);
    double bestCorr, corr;
    Matrix2D<double> Iaux;
    
    translate(1, Iaux, I1, vectorR2(-shiftX,-shiftY), DONT_WRAP);
    //I1.translate(vectorR2(-shiftX,-shiftY),Iaux, DONT_WRAP);
    bestCorr=corr=correlation_index(I2,Iaux);
    double finalX=shiftX;
    double finalY=shiftY;
    #ifdef DEBUG
        std::cout << "shiftX=" << shiftX << " shiftY=" << shiftY
                  << " corr=" << corr << std::endl;
        ImageXmipp save;
        save()=I1;   save.write("PPPI1.xmp");
        save()=I2;   save.write("PPPI2.xmp");
        save()=Iaux; save.write("PPPpp.xmp");
    #endif        
    
    Iaux.initZeros();
    double testX=(shiftX>0) ? (shiftX-XSIZE(I1)):(shiftX+XSIZE(I1));
    double testY=shiftY;
    translate(1, Iaux, I1, vectorR2(-testX,-testY), DONT_WRAP);
    //I1.translate(vectorR2(-testX,-testY),Iaux,DONT_WRAP);
    corr=correlation_index(I2,Iaux);
    if (corr>bestCorr) finalX=testX;
    #ifdef DEBUG
        std::cout << "shiftX=" << testX << " shiftY=" << testY
                  << " corr=" << corr << std::endl;
        save()=Iaux; save.write("PPPmp.xmp");
    #endif        

    Iaux.initZeros();
    testX=shiftX;
    testY=(shiftY>0) ? (shiftY-YSIZE(I1)):(shiftY+YSIZE(I1));
    translate(1, Iaux, I1, vectorR2(-testX,-testY), DONT_WRAP);
    //I1.translate(vectorR2(-testX,-testY),Iaux,DONT_WRAP);
    corr=correlation_index(I2,Iaux);
    if (corr>bestCorr)
        finalY=testY;
    #ifdef DEBUG
        std::cout << "shiftX=" << testX << " shiftY=" << testY
                  << " corr=" << corr << std::endl;
        save()=Iaux; save.write("PPPpm.xmp");
    #endif        

    Iaux.initZeros();
    testX=(shiftX>0) ? (shiftX-XSIZE(I1)):(shiftX+XSIZE(I1));
    testY=(shiftY>0) ? (shiftY-YSIZE(I1)):(shiftY+YSIZE(I1));
    translate(1, Iaux, I1, vectorR2(-testX,-testY), DONT_WRAP);
    //I1.translate(vectorR2(-testX,-testY),Iaux,DONT_WRAP);
    corr=correlation_index(I2,Iaux);
    if (corr>bestCorr) {
        finalX=testX;
        finalY=testY;
    }
    #ifdef DEBUG
        std::cout << "shiftX=" << testX << " shiftY=" << testY
                  << " corr=" << corr << std::endl;
        save()=Iaux; save.write("PPPmm.xmp");
    #endif        
    
    shiftX=finalX;
    shiftY=finalY;
}
#undef DEBUG

/* Align two images -------------------------------------------------------- */
void alignImages(const Matrix2D< double >& Iref, Matrix2D< double >& I,
    Matrix2D< double >&M)
{
    Iref.checkDimension(2);
    I.checkDimension(2);

    Matrix2D<double> ARS, ASR, R;
    ARS.initIdentity(3);
    ASR.initIdentity(3);
    Matrix2D<double> IauxSR=I, IauxRS=I;

    Polar_fftw_plans *plans=NULL;
    Polar< std::complex<double> > polarFourierIref;
    normalizedPolarFourierTransform(
	    Iref,
	    polarFourierIref,
	    false,
            XSIZE(Iref)/5,
	    XSIZE(Iref)/2,
            plans,
            1);

    XmippFftw local_transformer;
    Matrix1D<double> rotationalCorr;
    rotationalCorr.resize(2*polarFourierIref.getSampleNoOuterRing()-1);
    local_transformer.setReal(rotationalCorr);

    // Align the image with the reference
    for (int i=0; i<2; i++)
    {
        double shiftX, shiftY;
	
        // Shift then rotate	
        best_nonwrapping_shift(I,IauxSR,shiftX,shiftY);
        ASR(0,2)+=shiftX;
        ASR(1,2)+=shiftY;
        applyGeometry(1, IauxSR, I, ASR, IS_NOT_INV, WRAP);
        
        Polar< std::complex<double> > polarFourierI;
	normalizedPolarFourierTransform(
		IauxSR,
		polarFourierI,
		true,
                XSIZE(Iref)/5,
		XSIZE(Iref)/2,
                plans,
                1);
        
        double bestRot = best_rotation(polarFourierIref,polarFourierI,
            local_transformer);
	R=rotation2DMatrix(-bestRot);
        ASR=R*ASR;
        applyGeometry(1, IauxSR, I, ASR, IS_NOT_INV, WRAP);

        // Rotate then shift
	normalizedPolarFourierTransform(
		IauxRS,
		polarFourierI,
		true,
                XSIZE(Iref)/5,
		XSIZE(Iref)/2,
                plans,
                1);
        bestRot = best_rotation(polarFourierIref,polarFourierI,
            local_transformer);
	R=rotation2DMatrix(-bestRot);
        ARS=R*ARS;
        applyGeometry(1, IauxRS, I, ARS, IS_NOT_INV, WRAP);

        best_nonwrapping_shift(Iref,IauxRS,shiftX,shiftY);
        ARS(0,2)+=shiftX;
        ARS(1,2)+=shiftY;
        applyGeometry(1, IauxRS, I, ARS, IS_NOT_INV, WRAP);
    }
    
    double corrRS=correlation_index(IauxRS,Iref);
    double corrSR=correlation_index(IauxSR,Iref);
    if (corrRS>corrSR)
    {
        I=IauxRS;
        M=ARS;
    }
    else
    {
        I=IauxSR;
        M=ASR;
    }
}

/* Estimate 2D Gaussian ---------------------------------------------------- */
/* See Brandle, Chen, Bischof, Lapp. Robust parametric and semi-parametric
   spot fitting for spot array images. 2000 */
double unnormalizedGaussian2D(const Matrix1D<double> &r,
                  const Matrix1D<double> &mu,
                  const Matrix2D<double> &sigmainv)
{
    double x=XX(r)-XX(mu);
    double y=YY(r)-YY(mu);
    return exp(-0.5*(DIRECT_MAT_ELEM(sigmainv,0,0)*x*x+
                   2*DIRECT_MAT_ELEM(sigmainv,0,1)*x*y+
                     DIRECT_MAT_ELEM(sigmainv,1,1)*y*y));
}

void estimateGaussian2D(const Matrix2D<double> &I,
    double &a, double &b, Matrix1D<double> &mu, Matrix2D<double> &sigma,
    bool estimateMu, int iterations)
{
    I.checkDimension(2);

    Matrix2D<double> z(I);

    // Estimate b
    histogram1D hist;
    compute_hist(z,hist,100);
    b=hist.percentil(5);

    // Iteratively estimate all parameters
    for (int n=0; n<iterations; n++)
    {
        // Reestimate z
        FOR_ALL_ELEMENTS_IN_MATRIX2D(z)
           z(i,j)=XMIPP_MAX(I(i,j)-b,0);

        // Sum of z
        double T=z.sum();
        
        // Estimate center
        mu.initZeros(2);
        if (estimateMu)
        {
            FOR_ALL_ELEMENTS_IN_MATRIX2D(z)
            {
                double val=z(i,j);
                XX(mu)+=val*j;
                YY(mu)+=val*i;
            }
            mu/=T;
        }
        
        // Estimate sigma
        sigma.initZeros(2,2);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(z)
        {
            double val=z(i,j);
            double j_mu=j-XX(mu);
            double i_mu=i-YY(mu);
            DIRECT_MAT_ELEM(sigma,0,0)+=val*j_mu*j_mu;
            DIRECT_MAT_ELEM(sigma,0,1)+=val*i_mu*j_mu;
            DIRECT_MAT_ELEM(sigma,1,1)+=val*i_mu*i_mu;
        }
        DIRECT_MAT_ELEM(sigma,1,0)=DIRECT_MAT_ELEM(sigma,0,1);
        sigma/=T;
        
        // Estimate amplitude
        Matrix2D<double> sigmainv=sigma.inv();
        Matrix1D<double> r(2);
        double G2=0;
        a=0;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(z)
        {
            XX(r)=j;
            YY(r)=i;
            double G=unnormalizedGaussian2D(r,mu,sigmainv);
            a+=z(i,j)*G;
            G2+=G*G;
        }
        a/=G2;

        // Reestimate b
        FOR_ALL_ELEMENTS_IN_MATRIX2D(z)
        {
            XX(r)=j;
            YY(r)=i;
            double G=unnormalizedGaussian2D(r,mu,sigmainv);
            z(i,j)=I(i,j)-a*G;
        }
        compute_hist(z,hist,100);
        b=hist.percentil(5);
    }
}

/* Fourier-Bessel decomposition. ------------------------------------------- */
void Fourier_Bessel_decomposition(const Matrix2D<double> &img_in,
                                  Matrix2D<double> &m_out, double r1, double r2, int k1, int k2)
{
    img_in.checkDimension(2);

    for (int k = k1; k <= k2; k++)
    {
        int k_1 = k - 1;

        // Compute h and a,b coefficients
        double coefca = 0, coefcb = 0, coefsa = 0, coefsb = 0;
        double h = 0, my5 = 0;
        if (k_1 != 0)
        {
            double my = 1 + PI * r2 / 2 / k_1;
            double my2 = 2 * my;
            double my4 = my * k_1;
            double my5 = my4 - 1;
            double ntot = 4 * my4;
            h = 2 * PI / ntot;
            double hdpi = h / PI;
            double th = k_1 * h;
            double ys = sin(th);
            double zs = cos(th);
            double ys2 = sin(2 * th);
            double b1 = 2 / (th * th) * (1 + zs * zs - ys2 / th);
            double g1 = 4 / (th * th) * (ys / th - zs);
            double d1 = 2 * th / 45 * ys2;
            double e1 = d1 * ys * 2;
            coefca = (b1 + e1) * hdpi;
            coefcb = (g1 - d1) * hdpi;
            coefsa = (b1 - e1) * hdpi;
            coefsb = (g1 + d1) * hdpi;
        }
        else
        {
            double my = 1 + PI * r2 / 2;
            double my2 = 2 * my;
            double my4 = my;
            double my5 = my4 - 1;
            double ntot = 4 * my4;
            h = 2 * PI / ntot;
            coefca = h / PI / 2.;
        }

        Matrix1D<double> sine(CEIL(my5));
        FOR_ALL_ELEMENTS_IN_MATRIX1D(sine) sine(i) = sin((i + 1) * h);

    }
}

/* Harmonic decomposition. ------------------------------------------------- */
void harmonic_decomposition(const Matrix2D<double> &img_in,
                            Matrix1D<double> &v_out)
{}

/* Shah energy ------------------------------------------------------------- */
/* This function computes the current functional energy */
double Shah_energy(const Matrix2D<double> &img,
                   const Matrix2D<double> &surface_strength,
                   const Matrix2D<double> &edge_strength,
                   double K, const Matrix1D<double> &W)
{
    img.checkDimension(2);

    int Ydim1 = YSIZE(img) - 1;
    int Xdim1 = XSIZE(img) - 1;

    double Kinv = 1.0 / K;

    /* Calculate surface energy */
    double E1 = 0.0, E2 = 0.0, E3 = 0.0, E4 = 0.0;
    for (int i = 1; i < Ydim1; i++)
        for (int j = 1; j < Xdim1; j++)
        {
            /* Calculate data matching terms */
            double D = dMij(img, i, j);
            double F = dMij(surface_strength, i, j);
            double S = dMij(edge_strength, i, j);
            E1 += W(0) * (D - F) * (D - F);
            E3 += W(2) * K * S * S;

            /* Calculate first derivative terms */
            double Fx = (dMij(surface_strength, i, j + 1) - dMij(surface_strength, i, j - 1)) / 2;
            double Fy = (dMij(surface_strength, i + 1, j) - dMij(surface_strength, i - 1, j)) / 2;
            double Sx = (dMij(edge_strength, i, j + 1)    - dMij(edge_strength, i, j - 1)) / 2;
            double Sy = (dMij(edge_strength, i + 1, j)    - dMij(edge_strength, i - 1, j)) / 2;
            E2 += W(1) * (1 - S) * (1 - S) * (Fx * Fx + Fy * Fy);
            E4 += W(3) * Kinv * (Sx * Sx + Sy * Sy);
        }

    return E1 + E2 + E3 + E4; // Total energy
}

/* Update Surface Shah ----------------------------------------------------- */
/* This routine performs one update to the edge estimate based
   on a finite differences solution to the following equation:
       0 = dE/df = dF/df - d(dF/dfx)/dx - d(dF/dfy)/dy
           + dd(dF/dfxx)/dxx + dd(dF/dfxy)/dxy + dd(dF/dfyy)/dyy */
double Update_surface_Shah(Matrix2D<double> &img,
                           Matrix2D<double> &surface_strength,
                           Matrix2D<double> &edge_strength,
                           const Matrix1D<double> &W)
{
    img.checkDimension(2);

    double Diff = 0.0, Norm = 0.0;
    int Ydim1 = YSIZE(img) - 1;
    int Xdim1 = XSIZE(img) - 1;

    /* Update surface estimate */
    for (int i = 1; i < Ydim1; i++)
        for (int j = 1; j < Xdim1; j++)
        {
            /* Calculate edge partial derivative terms */
            double S  =  dMij(edge_strength, i, j);
            double Sx = (dMij(edge_strength, i, j + 1)    - dMij(edge_strength, i, j - 1)) / 2;
            double Sy = (dMij(edge_strength, i + 1, j)    - dMij(edge_strength, i - 1, j)) / 2;

            double nS  = 1 - S;
            double nS2 = nS * nS;

            /* Calculate surface partial derivative terms (excluding central pixel) */
            double F, D;
            F = D = dMij(img, i, j);
            double Fx = (dMij(surface_strength, i, j + 1) - dMij(surface_strength, i, j - 1)) / 2;
            double Fy = (dMij(surface_strength, i + 1, j) - dMij(surface_strength, i - 1, j)) / 2;
            double Fxx =  dMij(surface_strength, i, j + 1) + dMij(surface_strength, i, j - 1);
            double Fyy =  dMij(surface_strength, i + 1, j) + dMij(surface_strength, i - 1, j);

            /* Calculate surface partial derivative weights */
            double wFx = 4 * W(1) * nS * Sx;
            double wFy = 4 * W(1) * nS * Sy;
            double wFxx = -2 * W(1) * nS2;
            double wFyy = -2 * W(1) * nS2;

            /* Calculate new surface value */
            double Constant = -2 * W(0) * D;
            double Central  = -2 * W(0) + 2 * wFxx + 2 * wFyy;
            double Neighbors = wFx * Fx + wFy * Fy + wFxx * Fxx + wFyy * Fyy;

            if (ABS(Central) > XMIPP_EQUAL_ACCURACY)
                F = (Constant + Neighbors) / Central;
            F = CLIP(F, 0, 1);

            // Compute the difference.
            Diff += ABS(dMij(surface_strength, i, j) - F);
            Norm += ABS(dMij(surface_strength, i, j));

            // Update the new value.
            dMij(surface_strength, i, j) = F;
        }
    return Diff / Norm; // Return the relative difference.
}

/* Update Edge Shah -------------------------------------------------------- */
/* This routine performs one update to the edge estimate based
   on a finite differences solution to the following equation:
   0 = dE/ds = dF/ds - d(dF/dsx)/dx - d(dF/dsy)/dy */
double Update_edge_Shah(Matrix2D<double> &img,
                        Matrix2D<double> &surface_strength,
                        Matrix2D<double> &edge_strength,
                        double K,
                        const Matrix1D<double> &W)
{
    img.checkDimension(2);

    double Diff = 0.0, Norm = 0.0;
    int Ydim1 = YSIZE(img) - 1;
    int Xdim1 = XSIZE(img) - 1;
    double Kinv = 1.0 / K;

    /* Update edge estimate */
    for (int i = 1; i < Ydim1; i++)
        for (int j = 1; j < Xdim1; j++)
        {
            /* Calculate first and second derivative terms */
            double Fx = (dMij(surface_strength, i, j + 1) - dMij(surface_strength, i, j - 1)) / 2;
            double Fy = (dMij(surface_strength, i + 1, j) - dMij(surface_strength, i - 1, j)) / 2;
            double Constant = W(1) * (Fx * Fx + Fy * Fy);

            /* Calculate weights for central pixel and neighbors */
            double Central   = W(2) * K + W(3) * Kinv * 4;
            double Neighbors = W(3) * Kinv * (
                                   dMij(edge_strength, i - 1, j) + dMij(edge_strength, i + 1, j)
                                   + dMij(edge_strength, i, j - 1) + dMij(edge_strength, i, j + 1));

            /* Calculate new S value */
            double Old_edge_strength = dMij(edge_strength, i, j);
            double S = (Constant + Neighbors) / (Constant + Central);
            if (S < 0)
                dMij(edge_strength, i, j) /= 2;
            else if (S > 1)
                dMij(edge_strength, i, j) = (dMij(edge_strength, i, j) + 1) / 2;
            else
                dMij(edge_strength, i, j) = S;

            // Compute the difference.
            Diff += ABS(dMij(edge_strength, i, j) - Old_edge_strength);
            Norm += ABS(Old_edge_strength);
        }
    return Diff / Norm; // Return the relative difference.
}

/* Smoothing Shah ---------------------------------------------------------- */
#define SHAH_CONVERGENCE_THRESHOLD  0.0001
void Smoothing_Shah(Matrix2D<double> &img,
                    Matrix2D<double> &surface_strength,
                    Matrix2D<double> &edge_strength,
                    const Matrix1D<double> &W,
                    int OuterLoops,
                    int InnerLoops,
                    int RefinementLoops,
                    bool adjust_range)
{

    img.checkDimension(2);

    typeCast(img, surface_strength);
    if (adjust_range)
        surface_strength.rangeAdjust(0, 1);
    edge_strength.resize(img);

    for (int k = 1; k <= RefinementLoops; k++)
    {
        // Initialize Edge Image.
        edge_strength.initZeros();

        double diffsurface = MAXFLOAT; // Reset surface difference
        for (int i = 0; ((i < OuterLoops) && OuterLoops) ||
             ((diffsurface > SHAH_CONVERGENCE_THRESHOLD) && !OuterLoops); i++)
        {

            /* std::cout << "Iteration ..." << i+1;*/
            /* Iteratively update surface estimate */
            for (int j = 0; j < InnerLoops; j++)
                diffsurface =
                    Update_surface_Shah(img, surface_strength, edge_strength, W);

            /* Iteratively update edge estimate */
            double diffedge;
            for (int j = 0; j < InnerLoops; j++)
                diffedge =
                    Update_edge_Shah(img, surface_strength, edge_strength, k, W);

            /* Calculate new functional energy */
            double energy = Shah_energy(img, surface_strength, edge_strength, k, W);
            /* std::cout << " Energy " << energy
                 << " ... Relative Diff " << diffsurface
                 << " ... Edge diff " << diffedge
                 << std::endl; */
        }
    }
}

/* Tomographic diffusion --------------------------------------------------- */
//#define DEBUG
double tomographicDiffusion(Matrix3D< double >& V,
    const Matrix1D< double >& alpha, double lambda)
{
    V.checkDimension(3);

    double alphax=XX(alpha);
    double alphay=YY(alpha);
    double alphaz=ZZ(alpha);
    double diffx, diffy, diffz;

    // Compute regularization error
    double regError=0;
    for (int z=1; z<ZSIZE(V)-1; z++)
        for (int y=1; y<YSIZE(V)-1; y++)
            for (int x=1; x<XSIZE(V)-1; x++)
            {
                diffx=DIRECT_VOL_ELEM(V,z,y,x+1)-DIRECT_VOL_ELEM(V,z,y,x-1);
                diffy=DIRECT_VOL_ELEM(V,z,y+1,x)-DIRECT_VOL_ELEM(V,z,y-1,x);
                diffz=DIRECT_VOL_ELEM(V,z+1,y,x)-DIRECT_VOL_ELEM(V,z-1,y,x);
                regError+=sqrt(alphax*diffx*diffx+
                               alphay*diffy*diffy+
                               alphaz*diffz*diffz);
            }
    regError*=0.5;
    
    // Compute the gradient of the regularization error
    Matrix3D<double> gradient;
    gradient.initZeros(V);
    for (int z=2; z<ZSIZE(V)-2; z++)
        for (int y=2; y<YSIZE(V)-2; y++)
            for (int x=2; x<XSIZE(V)-2; x++)
            {
                // First term
                double V000=DIRECT_VOL_ELEM(V,z,y,x);
                double V_200=DIRECT_VOL_ELEM(V,z,y,x-2);
                double V_110=DIRECT_VOL_ELEM(V,z,y+1,x-1);
                double V_1_10=DIRECT_VOL_ELEM(V,z,y-1,x-1);
                double V_101=DIRECT_VOL_ELEM(V,z+1,y,x-1);
                double V_10_1=DIRECT_VOL_ELEM(V,z-1,y,x-1);
                diffx=V000-V_200;
                diffy=V_110-V_1_10;
                diffz=V_101-V_10_1;
                double t1=diffx/sqrt(alphax*diffx*diffx+
                    alphay*diffy*diffy+alphaz*diffz*diffz);
                
                // Second term
                double V200=DIRECT_VOL_ELEM(V,z,y,x+2);
                double V110=DIRECT_VOL_ELEM(V,z,y+1,x+1);
                double V1_10=DIRECT_VOL_ELEM(V,z,y-1,x+1);
                double V101=DIRECT_VOL_ELEM(V,z+1,y,x+1);
                double V10_1=DIRECT_VOL_ELEM(V,z-1,y,x+1);
                diffx=V200-V000;
                diffy=V110-V1_10;
                diffz=V101-V10_1;
                double t2=diffx/sqrt(alphax*diffx*diffx+
                    alphay*diffy*diffy+alphaz*diffz*diffz);
                
                // Third term
                double V0_20=DIRECT_VOL_ELEM(V,z,y-2,x);
                double V0_11=DIRECT_VOL_ELEM(V,z+1,y-1,x);
                double V0_1_1=DIRECT_VOL_ELEM(V,z-1,y-1,x);
                diffx=V1_10-V_1_10;
                diffy=V000-V0_20;
                diffz=V0_11-V0_1_1;
                double t3=diffy/sqrt(alphax*diffx*diffx+
                    alphay*diffy*diffy+alphaz*diffz*diffz);
                
                // Fourth term
                double V020=DIRECT_VOL_ELEM(V,z,y+2,x);
                double V011=DIRECT_VOL_ELEM(V,z+1,y+1,x);
                double V01_1=DIRECT_VOL_ELEM(V,z-1,y+1,x);
                diffx=V110-V_110;
                diffy=V020-V000;
                diffz=V011-V01_1;
                double t4=diffy/sqrt(alphax*diffx*diffx+
                    alphay*diffy*diffy+alphaz*diffz*diffz);

                // Fifth term
                double V00_2=DIRECT_VOL_ELEM(V,z-2,y,x);
                diffx=V10_1-V_10_1;
                diffy=V01_1-V0_1_1;
                diffz=V000-V00_2;
                double t5=diffz/sqrt(alphax*diffx*diffx+
                    alphay*diffy*diffy+alphaz*diffz*diffz);

                // Sixth term
                double V002=DIRECT_VOL_ELEM(V,z+2,y,x);
                diffx=V101-V_101;
                diffy=V011-V0_11;
                diffz=V002-V000;
                double t6=diffz/sqrt(alphax*diffx*diffx+
                    alphay*diffy*diffy+alphaz*diffz*diffz);
                
                // Compute gradient
                DIRECT_VOL_ELEM(gradient,z,y,x)=
                    0.5*(alphax*(t1-t2)+alphay*(t3-t4)+alphaz*(t5-t6));
            }
    #ifdef DEBUG
        VolumeXmipp save;
        save()=V;
        save.write("PPPvolume.vol");
        save()=gradient;
        save.write("PPPgradient.vol");
        std::cout << "Press any key\n";
        char c; std::cin >> c;
    #endif
    
    // Update volume
    for (int z=2; z<ZSIZE(V)-2; z++)
        for (int y=2; y<YSIZE(V)-2; y++)
            for (int x=2; x<XSIZE(V)-2; x++)
                DIRECT_VOL_ELEM(V,z,y,x)-=
                    lambda*DIRECT_VOL_ELEM(gradient,z,y,x);

    // Finish
    return regError;
}
#undef DEBUG

/* Rotational invariant moments -------------------------------------------- */
void rotational_invariant_moments(const Matrix2D<double> &img,
                                  const Matrix2D<int> *mask,
                                  Matrix1D<double> &v_out)
{
    img.checkDimension(2);

    // Prepare some variables
    double m_11 = 0, m_02 = 0, m_20 = 0, m_12 = 0, m_21 = 0, m_03 = 0, m_30 = 0; //, m_00=0;
    double normalize_x = 2.0 / XSIZE(img);
    double normalize_y = 2.0 / YSIZE(img);

    // Compute low-level moments
    FOR_ALL_ELEMENTS_IN_MATRIX2D(img)
    {
        if (mask != NULL)
            if (!(*mask)(i, j))
                continue;
        double val = img(i, j);
        double dx = j * normalize_x;
        double dy = i * normalize_y;
        double dx2 = dx * dx, dx3 = dx2 * dx;
        double dy2 = dy * dy, dy3 = dy2 * dy;
        // m_00+=val;
        m_11 += val * dx * dy;
        m_02 += val * dy2;
        m_20 += val * dx2;
        m_12 += val * dx * dy2;
        m_21 += val * dx2 * dy;
        m_03 += val * dy3;
        m_30 += val * dx3;
    }
    //m_11/=m_00; m_02/=m_00; m_20/=m_00;
    //m_12/=m_00; m_21/=m_00; m_03/=m_00; m_30/=m_00;

    // Compute high-level rotational invariant moments
    v_out.resize(5);
    v_out(0) = m_20 + m_02;
    v_out(1) = (m_20 - m_02) * (m_20 - m_02) + 4 * m_11 * m_11;
    v_out(2) = (m_30 - 3 * m_12) * (m_30 - 3 * m_12) +
               (3 * m_21 - m_03) * (3 * m_21 - m_03);
    v_out(3) = (m_30 + m_12) * (m_30 + m_12) +
               (m_21 + m_03) * (m_21 + m_03);
    v_out(4) = (m_30 - 3 * m_12) * (m_30 + m_12) * (
                   (m_30 + m_12) * (m_30 + m_12)
                   - 3 * (m_21 + m_03) * (m_21 + m_03)) +
               (3 * m_21 - m_03) * (m_21 + m_03) * (
                   3 * (m_30 + m_12) * (m_30 + m_12)
                   - (m_21 + m_03) * (m_21 + m_03));
    /*
       v_out( 5)=(m_20+m_02)*(
                    (m_30+m_12)*(m_30+m_12)
            -3*(m_21+m_03)*(m_21+m_03))
          +4*m_11*(m_30+m_12)*(m_03+m_21);
       v_out( 6)=(3*m_21-m_03)*(m_12+m_30)*(
                      (m_30+m_12)*(m_30+m_12)
            -3*(m_21+m_03)*(m_21+m_03))
          -(m_30-3*m_12)*(m_12+m_03)*(
             3*(m_30+m_12)*(m_30+m_12)
         -(m_21+m_03)*(m_21+m_03));
    */
}

/* Inertia moments --------------------------------------------------------- */
void inertia_moments(const Matrix2D<double> &img,
                     const Matrix2D<int> *mask,
                     Matrix1D<double> &v_out,
                     Matrix2D<double> &u)
{
    img.checkDimension(2);

    // Prepare some variables
    double m_11 = 0, m_02 = 0, m_20 = 0;
    double normalize_x = 2.0 / XSIZE(img);
    double normalize_y = 2.0 / YSIZE(img);

    // Compute low-level moments
    FOR_ALL_ELEMENTS_IN_MATRIX2D(img)
    {
        if (mask != NULL)
            if (!(*mask)(i, j))
                continue;
        double val = img(i, j);
        double dx = j * normalize_x;
        double dy = i * normalize_y;
        double dx2 = dx * dx;
        double dy2 = dy * dy;
        m_11 += val * dx * dy;
        m_02 += val * dy2;
        m_20 += val * dx2;
    }

    // Compute the eigen values of the inertia matrix
    // [m_02 m_11
    //  m_11 m_20]
    Matrix2D<double> A(2, 2);
    A(0, 0) = m_02;
    A(0, 1) = A(1, 0) = m_11;
    A(1, 1) = m_20;
    Matrix2D<double> v;
    svdcmp(A, u, v_out, v);
    v_out = v_out.sort();
}

/* Fill triangle ----------------------------------------------------------- */
void fill_triangle(Matrix2D<double> &img, int *tx, int *ty, double color)
{
    img.checkDimension(2);

    /*
     * Order in y values
     */
    int y1 = 0;
    while (!y1)
    {
        y1 = 1;
        for (int y2 = 0; y2 < 2; y2++)
        {
            if (ty[y2] > ty[y2 + 1] ||
                (ty[y2] == ty[y2 + 1] && tx[y2] < tx[y2 + 1]))
            {
                int x1 = ty[y2];
                ty[y2] = ty[y2 + 1];
                ty[y2 + 1] = x1;
                x1 = tx[y2];
                tx[y2] = tx[y2 + 1];
                tx[y2 + 1] = x1;
                y1 = 0;
            }
        }
    }

    int dx1 = tx[1] - tx[0];
    int dx2 = tx[2] - tx[0];
    int dy1 = ty[1] - ty[0];
    int dy2 = ty[2] - ty[0];

    int sx1 = SGN0(dx1);
    int sx2 = SGN0(dx2);
    int sy1 = SGN0(dy1);

    int ix1 = ABS(dx1);
    int ix2 = ABS(dx2);
    int iy1 = ABS(dy1);
    int iy2 = ABS(dy2);

    int inc1 = XMIPP_MAX(ix1, iy1);
    int inc2 = XMIPP_MAX(ix2, iy2);

    int x1, x2, y2, xl, xr;
    x1 = x2 = y1 = y2 = 0;
    xl = xr = tx[0];
    int y = ty[0];

    while (y != ty[1])
    {
        int doit1 = 0;
        int doit2 = 0;

        while (!doit1)
        {   /* Wait until y changes */
            x1 += ix1;
            y1 += iy1;
            if (x1 > inc1)
            {
                x1 -= inc1;
                xl += sx1;
            }
            if (y1 > inc1)
            {
                y1 -= inc1;
                y += sy1;
                doit1 = 1;
            }
        }

        while (!doit2)
        {   /* Wait until y changes */
            x2 += ix2;
            y2 += iy2;
            if (x2 > inc2)
            {
                x2 -= inc2;
                xr += sx2;
            }
            if (y2 > inc2)
            {
                y2 -= inc2;
                doit2 = 1;
            }
        }

        for (int myx = xl; myx <= xr; myx++)
            img(y, myx) = color;
    }

    dx1 = tx[2] - tx[1];
    dy1 = ty[2] - ty[1];

    sx1 = SGN0(dx1);
    sy1 = SGN0(dy1);

    ix1 = ABS(dx1);
    iy1 = ABS(dy1);

    inc1 = XMIPP_MAX(ix1, iy1);
    xl = tx[1];
    x1 = 0;

    while (y != ty[2])
    {
        int doit1 = 0;
        int doit2 = 0;

        while (!doit1)
        {   /* Wait until y changes */
            x1 += ix1;
            y1 += iy1;
            if (x1 > inc1)
            {
                x1 -= inc1;
                xl += sx1;
            }
            if (y1 > inc1)
            {
                y1 -= inc1;
                y += sy1;
                doit1 = 1;
            }
        }

        while (!doit2)
        {   /* Wait until y changes */
            x2 += ix2;
            y2 += iy2;
            if (x2 > inc2)
            {
                x2 -= inc2;
                xr += sx2;
            }
            if (y2 > inc2)
            {
                y2 -= inc2;
                doit2 = 1;
            }
        }

        for (int myx = xl; myx <= xr; myx++)
            img(y, myx) = color;
    }
}

/* Local thresholding ------------------------------------------------------ */
void local_thresholding(Matrix2D<double> &img,
                        double C,
                        double dimLocal,
                        Matrix2D<int> &result,
                        Matrix2D<int> *mask)
{
    img.checkDimension(2);

    // Convolve the input image with the kernel
    Matrix2D<double> convolved;
    convolved.initZeros(img);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(convolved)
    {
        if (mask != NULL)
            if (!(*mask)(i, j))
                continue;
        int ii0 = XMIPP_MAX(STARTINGY(convolved), FLOOR(i - dimLocal));
        int jj0 = XMIPP_MAX(STARTINGX(convolved), FLOOR(j - dimLocal));
        int iiF = XMIPP_MIN(FINISHINGY(convolved), CEIL(i + dimLocal));
        int jjF = XMIPP_MIN(FINISHINGX(convolved), CEIL(j + dimLocal));
        double N = 0;
        for (int ii = ii0; ii <= iiF; ii++)
            for (int jj = jj0; jj <= jjF; jj++)
            {
                if (mask == NULL)
                {
                    convolved(i, j) += img(ii, jj);
                    ++N;
                }
                else if ((*mask)(i, j))
                {
                    convolved(i, j) += img(ii, jj);
                    ++N;
                }
            }
        if (N != 0)
            convolved(i, j) /= N;
    }

    // Substract the original from the convolved image and threshold
    result.initZeros(img);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(convolved)
    {
        if (mask != NULL)
            if (!(*mask)(i, j))
                continue;
        if (img(i, j) - convolved(i, j) > C)
            result(i, j) = 1;
    }
}

/* Center translationally -------------------------------------------------- */
void centerImageTranslationally(Matrix2D<double> &I)
{
    I.checkDimension(2);

    Matrix2D<double> Ix  = I; Ix.selfReverseX();   Ix.setXmippOrigin();
    Matrix2D<double> Iy  = I; Iy.selfReverseY();   Iy.setXmippOrigin();
    Matrix2D<double> Ixy = Ix; Ixy.selfReverseY(); Ixy.setXmippOrigin();
    
    double meanShiftX=0, meanShiftY=0, shiftX, shiftY;
    best_nonwrapping_shift(I,Ix, meanShiftX,meanShiftY);
    best_nonwrapping_shift(I,Iy, shiftX,shiftY);
    meanShiftX+=shiftX; meanShiftY+=shiftY;
    best_nonwrapping_shift(I,Ixy,shiftX,shiftY);
    meanShiftX+=shiftX; meanShiftY+=shiftY;
    meanShiftX/=3; meanShiftY/=3;
    
    Matrix1D<double> shift(2);
    VECTOR_R2(shift,-meanShiftX,-meanShiftY);
    Matrix2D<double> aux = I;
    translate(3, I, aux, shift);
    //I.selfTranslateBSpline(3,shift);
}

/* Center rotationally ----------------------------------------------------- */
void centerImageRotationally(Matrix2D<double> &I)
{
    I.checkDimension(2);

    Matrix2D<double> Ix  = I; Ix.selfReverseX();
    Ix.setXmippOrigin();

    Polar_fftw_plans *plans=NULL;
    Polar< std::complex<double> > polarFourierI, polarFourierIx;
    normalizedPolarFourierTransform(Ix,polarFourierIx,false,XSIZE(Ix)/5,
        XSIZE(Ix)/2,plans);
    normalizedPolarFourierTransform(I, polarFourierI, true, XSIZE(I)/5,
        XSIZE(I)/2,plans);

    XmippFftw local_transformer;
    Matrix1D<double> rotationalCorr;
    rotationalCorr.resize(2*polarFourierI.getSampleNoOuterRing()-1);
    local_transformer.setReal(rotationalCorr);
    double bestRot = best_rotation(polarFourierIx,polarFourierI,
        local_transformer);

    Matrix2D<double> aux = I;
    rotate(3, I, aux, -bestRot/2,WRAP);
    //I.selfRotateBSpline(3,-bestRot/2,WRAP);
}

/* Center both rotationally and translationally ---------------------------- */
//#define DEBUG
void centerImage(Matrix2D<double> &I, int Niter, bool limitShift)
{
    I.checkDimension(2);

    I.setXmippOrigin();
    double avg=I.computeAvg();
    I-=avg;

    Matrix2D<double> Ix, Iy, Ixy, Iaux, A;
    A.initIdentity(3);
    Iaux=I;
    
    Matrix2D<int> mask;
    mask.initZeros(I);
    BinaryCircularMask(mask,XSIZE(I)/2);
    
    Matrix1D<double> lineY, lineX;
    lineY.initZeros(YSIZE(I)); STARTINGX(lineY)=STARTINGY(I);
    lineX.initZeros(XSIZE(I)); STARTINGX(lineX)=STARTINGX(I);

    Polar_fftw_plans *plans=NULL;
    for (int i=0; i<Niter; i++)
    {
        // Mask Iaux
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (!mask(i,j)) Iaux(i,j)=0;

        // Center translationally
        Ix  = Iaux;  Ix.selfReverseX();  Ix.setXmippOrigin();
        Iy  = Iaux;  Iy.selfReverseY();  Iy.setXmippOrigin();
        Ixy = Ix;    Ixy.selfReverseY(); Ixy.setXmippOrigin();
        
        double meanShiftX=0, meanShiftY=0, shiftX, shiftY, Nx=0, Ny=0;
        best_nonwrapping_shift(Iaux,Ix, shiftX,shiftY);
        #ifdef DEBUG
            ImageXmipp save;
            save()=Ix; save.write("PPPx.xmp");
            std::cout << "con Ix: " << shiftX << " " << shiftY << std::endl;
        #endif
        if (ABS(shiftX)<XSIZE(I)/3 || !limitShift) {meanShiftX+=shiftX; Nx++;}
        if (ABS(shiftY)<YSIZE(I)/3 || !limitShift) {meanShiftY+=shiftY; Ny++;}
        best_nonwrapping_shift(Iaux,Iy, shiftX,shiftY);
        #ifdef DEBUG
            save()=Iy; save.write("PPPy.xmp");
            std::cout << "con Iy: " << shiftX << " " << shiftY << std::endl;
        #endif
        if (ABS(shiftX)<XSIZE(I)/3 || !limitShift) {meanShiftX+=shiftX; Nx++;}
        if (ABS(shiftY)<YSIZE(I)/3 || !limitShift) {meanShiftY+=shiftY; Ny++;}
        best_nonwrapping_shift(Iaux,Ixy,shiftX,shiftY);
        #ifdef DEBUG
            save()=Ixy; save.write("PPPxy.xmp");
            std::cout << "con Ixy: " << shiftX << " " << shiftY << std::endl;
        #endif
        if (ABS(shiftX)<XSIZE(I)/3 || !limitShift) {meanShiftX+=shiftX; Nx++;}
        if (ABS(shiftY)<YSIZE(I)/3 || !limitShift) {meanShiftY+=shiftY; Ny++;}
        if (Nx>0) meanShiftX/=Nx;
        if (Ny>0) meanShiftY/=Ny;

        A(0,2)+=-meanShiftX/2;
        A(1,2)+=-meanShiftY/2;
        Iaux.initZeros();
        applyGeometry(1, Iaux, I, A, IS_NOT_INV, WRAP);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (!mask(i,j)) Iaux(i,j)=0;
        
        #ifdef DEBUG
            std::cout << "Iter " << i << std::endl;
            std::cout << "shift=" << -meanShiftX << "," << -meanShiftY << std::endl;
            save()=I; save.write("PPP.xmp");
            save()=Iaux; save.write("PPPshift.xmp");
        #endif

        // Center rotationally
        Ix  = Iaux;  Ix.selfReverseX();  Ix.setXmippOrigin();
        
        Polar< std::complex<double> > polarFourierI;
        normalizedPolarFourierTransform(Iaux, polarFourierI, true, XSIZE(I)/5,
            XSIZE(I)/2,plans);
        XmippFftw local_transformer;
        Matrix1D<double> rotationalCorr;
        rotationalCorr.resize(2*polarFourierI.getSampleNoOuterRing()-1);
        local_transformer.setReal(rotationalCorr);

        Polar< std::complex<double> > polarFourierIx;
        normalizedPolarFourierTransform(Ix,polarFourierIx,false,XSIZE(Ix)/5,
            XSIZE(Ix)/2,plans);
        double bestRot = best_rotation(polarFourierIx,polarFourierI,
            local_transformer);
        bestRot = realWRAP(bestRot,0,180);
        if (bestRot>90) bestRot=bestRot-180;
        
        A=rotation2DMatrix(-bestRot/2)*A;
        Iaux.initZeros();
        applyGeometry(1, Iaux, I, A, IS_NOT_INV, WRAP);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (!mask(i,j)) Iaux(i,j)=0;

        #ifdef DEBUG
            std::cout << "rot=" << -bestRot/2 << std::endl;
            save()=Iaux; save.write("PPProt.xmp");
        #endif
        
        // Remove horizontal/vertical ambiguity
        lineX.initZeros();
        lineY.initZeros();
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Iaux)
        {
            double val=Iaux(i,j);
            if (j==0) lineY(i)=val;
            else if (lineY(i)<val) lineY(i)=val;
            if (i==0) lineX(j)=val;
            else if (lineX(j)<val) lineX(j)=val;
        }
        
        double thX=lineX.computeMin()+0.75*(lineX.computeMax()-lineX.computeMin());
        double thY=lineY.computeMin()+0.75*(lineY.computeMax()-lineY.computeMin());
        int x0=STARTINGX(lineX);  while (lineX(x0)<thX) x0++;
        int y0=STARTINGX(lineY);  while (lineY(y0)<thY) y0++;
        int xF=FINISHINGX(lineX); while (lineX(xF)<thX) xF--;
        int yF=FINISHINGX(lineY); while (lineY(yF)<thY) yF--;
        if ((xF-x0)>(yF-y0))
            A=rotation2DMatrix(90)*A;
        applyGeometry(1, Iaux, I, A, IS_NOT_INV, WRAP);
        #ifdef DEBUG
            lineX.write("PPPlineX.txt");
            lineY.write("PPPlineY.txt");
            std::cout << "dev X=" << xF-x0 << std::endl;
            std::cout << "dev Y=" << yF-y0 << std::endl;
            save()=Iaux; save.write("PPPhorver.xmp");
            std::cout << "Press any key\n";
            char c; std::cin >> c;
        #endif
    }
    applyGeometry(3, Iaux, I, A, IS_NOT_INV,WRAP);
    I=Iaux;
    I+=avg;
    if (plans!=NULL) delete plans;
}
#undef DEBUG
