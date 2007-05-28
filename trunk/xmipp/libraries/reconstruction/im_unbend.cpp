/***************************************************************************
 *
 * Authors:     Debora Gil
                Roberto Marabini (roberto@mipg.upenn.edu)
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

#include "im_unbend.h"

#include <sys/types.h>
#include <sys/times.h>

//////////////////////////// AUXILIAR FUNCTIONS

// Comparison function for qsort of Lat. Pts
// in increasing x values.
// Needed for Computation of DELAUNAY Triang.
/*int XYZCompare(const void *v1,const void *v2)
{
      XYZ *p1,*p2;
      p1 = (XYZ *) v1;
      p2 = (XYZ *) v2;

      return (p1->x - p2->x);

}*/

int LatPtCompare(const void *v1, const void *v2)
{
    LatPoint *p1, *p2;
    p1 = (LatPoint *) v1;
    p2 = (LatPoint *) v2;

    return (int)(p1->x - p2->x);
}

//////////////////////////////////////////////////////////////
//////////////////////////////  I/O FUNCTIONS


void ImUmbend::ReadMRCCord()
{
    ExpLat.read(FN_Correlation);
    return;
}

/////////////////////////////////////////////////////
//////////////////////////// PEAKS CORRESPONDANCE


////Peaks Correspondance.
void ImUmbend::PeaksCorresp()
{
    double auxX, auxY;
    LatPoint valCoord;
    //Threshold over CCvalue
    double Th = cc_peak_factor * ExpLat.cc_max;


    for (int kk = 0; kk < ExpLat.MRC_Xcoord.size(); kk++)
    {

        //Compute ideal position in regular grid
        valCoord.i = ExpLat.MRC_Xindex[kk];
        valCoord.j = ExpLat.MRC_Yindex[kk];
        valCoord.x = ExpLat.O[0] + ExpLat.MRC_Xindex[kk] * ExpLat.a(0) + ExpLat.MRC_Yindex[kk] * ExpLat.b(0);
        valCoord.y = ExpLat.O[1] + ExpLat.MRC_Xindex[kk] * ExpLat.a(1) + ExpLat.MRC_Yindex[kk] * ExpLat.b(1);
        //Deviation of experimental position
        valCoord.Incrx = ExpLat.MRC_Xcoord[kk] - valCoord.x;
        valCoord.Incry = ExpLat.MRC_Ycoord[kk] - valCoord.y;
        valCoord.Interp = (ExpLat.MRC_CCcoheficiente[kk] < Th);
        INCR_coord.push_back(valCoord);
    }

#undef DEBUG
//  #define DEBUG
#ifdef DEBUG
    for (int kk = 0; kk < INCR_coord.size(); kk++)
        cout <<  INCR_coord[kk].x << " " << INCR_coord[kk].y << " "
        <<  INCR_coord[kk].Incrx << " " << INCR_coord[kk].Incry << endl;

#endif
#undef DEBUG
}





//////////////////////////////////////////////////////////////
///////////////////////////// IMAGE UNBENDING
// OBS: Interpolation based on Delanauy is not efficient by the point on triangle search (FindNearest)
// We implement instead interpolation of shifts to the whole regular grid (using Delanauy) and then linear
// intrepolation to the whole image.

void ImUmbend::UnBending()
{

    //Del. Triang.
    vector <ITRIANGLE> LatTri;
    //Interpolation Parameters
    LatPoint TargetPt;
    //Auxiliar variables
    float Tx, Ty;
    int Ti, Tj;
    float A, A1, A2, A3, A4;
    float O[2];

    //Read Input Image
    try
    {
        inIm.read(inImfile);
    }
    catch (Xmipp_error XE)
    {
        cerr << XE;
        exit(1);
    }
    //xmipp coordinates in image
    //inIm().setXmippOrigin();
    O[1] = inIm().ColNo() * 0.5;
    O[2] = inIm().RowNo() * 0.5;
    outIm().resize(inIm());

    //Lattice Triangulation (MRC coordinates)
    LatTriang(LatTri);

#undef DEBUG


    // Extend Shifts to Regular Grid (MRC Coordinates)
    // (Remove it when you find an efficient triangle search)
    int Nh, Nk;
    Nk = (int)ceil(min(abs(ExpLat.dim[0] / ExpLat.a(0)), abs(ExpLat.dim[0] / ExpLat.b(0))));
    Nh = (int)ceil(min(abs(ExpLat.dim[1] / ExpLat.a(1)), abs(ExpLat.dim[1] / ExpLat.b(1))));
    matrix2D <double> MatIncrX, MatIncrY;
    MatIncrX.initZeros(Nh, Nk);
    MatIncrY.initZeros(Nh, Nk);

#undef DEBUG
    //   #define DEBUG
#ifdef DEBUG
    cout <<  MatIncrX.RowNo() << " " << MatIncrX.ColNo() << endl;
#endif
#undef DEBUG

    Scattered2Regular(MatIncrX, MatIncrY, LatTri);




    //Interpolate Values for Transformation inIm(i+shift_Y,j+shift_X)
    FOR_ALL_ELEMENTS_IN_MATRIX2D(outIm())
    {

        //Interpolate Experimental Shifts


        // Change coordinates from xmipp to MRC (Only for interpolations based on Delanauy triangulation of MRC mesh)
        //TargetPt.x=j+O[1];
        //TargetPt.y=i+O[2];
        TargetPt.x = j;
        TargetPt.y = i;
        TargetPt.Incrx = 0;
        TargetPt.Incry = 0;



        ShiftsInterpReg(MatIncrX, MatIncrY, TargetPt);
        // Not efficient
        // ShiftsInterp(TargetPt,LatTri);

#ifdef TIMES
        times(&after);
        cout << "Total Shift time " << after.tms_utime - before.tms_utime << endl;
#endif

#undef DEBUG
        //    #define DEBUG
#ifdef DEBUG
        cout <<  TargetPt.x << " " << TargetPt.y << " " << TargetPt.Incrx << " " << TargetPt.Incry << endl;
#endif
#undef DEBUG

        //New pixel position (in xmipp coordinates)
        Ty = i + TargetPt.Incry;
        Tx = j + TargetPt.Incrx;
        Ti = (int)floor(Ty);
        Tj = (int)floor(Tx);
        //Check not outside input image
        if ((!inIm().outside(Ti, Tj)) && (!inIm().outside(Ti + 1, Tj + 1)))
        {
            //BiLinear Interpolation for image transformation
            A1 = fabs(Ty - Ti - 1) * fabs(Tx - Tj - 1);
            A2 = fabs(Ty - Ti - 1) * fabs(Tx - Tj);
            A3 = fabs(Ty - Ti) * fabs(Tx - Tj - 1);
            A4 = fabs(Ty - Ti) * fabs(Tx - Tj);
            A = A1 + A2 + A3 + A4;
            //A=1;
            outIm(i, j) = A1 * inIm(Ti, Tj) + A2 * inIm(Ti, Tj + 1) + A3 * inIm(Ti + 1, Tj) + A4 * inIm(Ti + 1, Tj + 1);
            outIm(i, j) = outIm(i, j) / A;

        }
        else
        {
            outIm(i, j) = 0;
        }
    }

    /* save image  */
    try
    {
        outIm.write(outImfile);
    }
    catch (Xmipp_error XE)
    {
        cerr << XE;
        exit(1);
    }

}

/////////////////////////////////////////////////////
////////////////////////// INTERPOLATION

//Shifts Interpolation from Regular grid
void ImUmbend::ShiftsInterpReg(matrix2D <double> & MatIncrX, matrix2D <double> & MatIncrY, LatPoint &TargetPt)
{
    int i, j;
    float det;
    float Tx, Ty, Ti, Tj, TiM, TjM;
    float A, A1, A2, A3, A4;
    float ACoeff[4];



    //Indexes of nearest point in grid
    det = ExpLat.a(0) * ExpLat.b(1) - ExpLat.a(1) * ExpLat.b(0);
    Tx = TargetPt.x;
    Ty = TargetPt.y;
    j = (int)((ExpLat.b(1) * Tx - ExpLat.b(0) * Ty) / det); //Coordenada x
    i = (int)((-ExpLat.a(1) * Tx + ExpLat.a(0) * Ty) / det); //Coordenada y

    //Nearest neighbors
    Tj = j * ExpLat.a(0) + i * ExpLat.b(0);
    Ti = j * ExpLat.a(1) + i * ExpLat.b(1);
    TjM = (j + 1) * ExpLat.a(0) + (i + 1) * ExpLat.b(0);
    TiM = (j + 1) * ExpLat.a(1) + (i + 1) * ExpLat.b(1);

    //Area computation (CHOSE INTERPOLATION method here)

    Interp2D(Tx, Ty, Ti, Tj, TiM, TjM, ACoeff);
    A1 = ACoeff[0];
    A2 = ACoeff[1];
    A3 = ACoeff[2];
    A4 = ACoeff[3];
    A = A1 + A2 + A3 + A4;

    //Shift Interpolation
    if ((!MatIncrX.outside(i, j)) && (!MatIncrX.outside(i + 1, j + 1)))
    {
        TargetPt.Incrx = A1 * MatIncrX(i, j) + A2 * MatIncrX(i, j + 1) + A3 * MatIncrX(i + 1, j) + A4 * MatIncrX(i + 1, j + 1);
        TargetPt.Incrx = TargetPt.Incrx / A;
        TargetPt.Incry = A1 * MatIncrY(i, j) + A2 * MatIncrY(i, j + 1) + A3 * MatIncrY(i + 1, j) + A4 * MatIncrY(i + 1, j + 1);
        TargetPt.Incry = TargetPt.Incry / A;

#undef DEBUG
        //  #define DEBUG
#ifdef DEBUG
        cout <<  Tx << " " << Ty << " " << i << " " << j << endl;
        cout << TargetPt.Incrx << " " << TargetPt.Incry << "  " << MatIncrX(i, j) << " " << MatIncrY(i, j) << endl;
#endif
#undef DEBUG
    }
    else if ((MatIncrX.outside(i + 1, j)) &&
             (MatIncrX.outside(i, j + 1)) &&
             (MatIncrX.outside(i + 1, j + 1))
            )
    {
        TargetPt.Incrx = MatIncrX(i, j);
        TargetPt.Incry = MatIncrY(i, j);
    }
    else if ((MatIncrX.outside(i + 1, j)) &&
             (MatIncrX.outside(i + 1, j + 1))
            )
    {
        A3 = A4 = 0;
        A = A1 + A2 + A3 + A4;
        TargetPt.Incrx = A1 * MatIncrX(i, j) + A2 * MatIncrX(i, j + 1);
        TargetPt.Incrx = TargetPt.Incrx / A;
        TargetPt.Incry = A1 * MatIncrY(i, j) + A2 * MatIncrY(i, j + 1);
        TargetPt.Incry = TargetPt.Incry / A;
    }
    else if ((MatIncrX.outside(i, j + 1)) &&
             (MatIncrX.outside(i + 1, j + 1))
            )
    {
        A2 = A4 = 0;
        A = A1 + A2 + A3 + A4;
        TargetPt.Incrx = A1 * MatIncrX(i, j) + A3 * MatIncrX(i + 1, j);
        TargetPt.Incrx = TargetPt.Incrx / A;
        TargetPt.Incry = A1 * MatIncrY(i, j) + A3 * MatIncrY(i + 1, j);
        TargetPt.Incry = TargetPt.Incry / A;
    }

}
//2D Interpolation on Square acoording to InterpModel
void ImUmbend::Interp2D(float Tx, float Ty, float Ti, float Tj, float TiM, float TjM, float * ACoeff)
{
    double R0, h, x0, y0;
    double Na, Nb;
    int indR;
    int Bdim = 300;

    Na = sqrt(ExpLat.a(0) * ExpLat.a(0) + ExpLat.a(1) * ExpLat.a(1));
    Nb = sqrt(ExpLat.b(0) * ExpLat.b(0) + ExpLat.b(1) * ExpLat.b(1));

    h = 0.01;
    if (strcmp(InterpModel.c_str(), "Bessel") == 0)
    {

        x0 = SIG * (Tx - TjM) / Na;
        y0 = SIG * (Ty - TiM) / Nb;
        R0 = sqrt(x0 * x0 + y0 * y0);
        indR = ROUND(R0 / h);
        if (indR < Bdim)
            ACoeff[3] = B[indR];
        else
            ACoeff[3] = 0;

        x0 = SIG * (Tx - Tj) / Na;
        y0 = SIG * (Ty - TiM) / Nb;
        R0 = sqrt(x0 * x0 + y0 * y0);
        indR = ROUND(R0 / h);
        if (indR < Bdim)
            ACoeff[2] = B[indR];
        else
            ACoeff[2] = 0;

        x0 = SIG * (Tx - TjM) / Na;
        y0 = SIG * (Ty - Ti) / Nb;
        R0 = sqrt(x0 * x0 + y0 * y0);
        indR = ROUND(R0 / h);
        if (indR < Bdim)
            ACoeff[1] = B[indR];
        else
            ACoeff[1] = 0;

        x0 = SIG * (Tx - Tj) / Na;
        y0 = SIG * (Ty - Ti) / Nb;
        R0 = sqrt(x0 * x0 + y0 * y0);
        indR = ROUND(R0 / h);
        if (indR < Bdim)
            ACoeff[0] = B[indR];
        else
            ACoeff[0] = 0;
    }
    else
    {
        ACoeff[0] = fabs(Ty - TiM) * fabs(Tx - TjM);
        ACoeff[1] = fabs(Ty - TiM) * fabs(Tx - Tj);
        ACoeff[2] = fabs(Ty - Ti) * fabs(Tx - TjM);
        ACoeff[3] = fabs(Ty - Ti) * fabs(Tx - Tj);
    }
}


///////////////////////////////////////////////////////////////////////////
//Linear Interpolation from scattered data set to regular grid
void ImUmbend::Scattered2Regular(matrix2D <double> & MatIncrX, matrix2D <double> & MatIncrY, vector <ITRIANGLE> &LatTri)
{

    int k;
    int N = INCR_coord.size();
    //Interpolation Parameters
    LatPoint TargetPt;

    //Interpolation to the whole grid
    for (k = 0;k < N;k++)
    {
        if (INCR_coord[k].Interp)
        {
            TargetPt.x = INCR_coord[k].x;
            TargetPt.y = INCR_coord[k].y;
            TargetPt.Incrx = 0;
            TargetPt.Incry = 0;
            ShiftsInterp(TargetPt, LatTri);
            INCR_coord[k] = TargetPt;
        }
    }

    //Matrix Form
    int i, j;
    double Tx, Ty, det;
    for (k = 0;k < N;k++)
    {

        det = ExpLat.a(0) * ExpLat.b(1) - ExpLat.a(1) * ExpLat.b(0);
        Tx = INCR_coord[k].x;
        Ty = INCR_coord[k].y;
        //index coordinates in lattice
        j =  ROUND((ExpLat.b(1) * Tx - ExpLat.b(0) * Ty) / det); //Coordenada x
        i =  ROUND((-ExpLat.a(1) * Tx + ExpLat.a(0) * Ty) / det); //Coordenada y



        MatIncrX(i, j) = INCR_coord[k].Incrx;
        MatIncrY(i, j) = INCR_coord[k].Incry;
#undef DEBUG
        //  #define DEBUG
#ifdef DEBUG
        cout <<  i << " " << j << " " << INCR_coord[k].x << " " << INCR_coord[k].y << "  " << MatIncrX(i, j) << " " << MatIncrY(i, j) << endl;
#endif
#undef DEBUG


    }

}


// Displacement Interpolation from Triangulation of Irregular Grid
// Border values are set to nearest neighbour
void ImUmbend::ShiftsInterp(LatPoint &TargetPt, vector <ITRIANGLE> &  LatTri)
{
    //Interpolation Parameters
    int Tind, Ptind;
    int n;
    float w[3], del;
    float x[3], y[3], Incrx[3], Incry[3];
    float xi, yi;



    //#undef TIMES
    //#define TIMES
#ifdef TIMES
    struct tms before, after;
    times(&before);
#endif


    //find nearest Triangle
    Tind = FindNearestTri(TargetPt, LatTri);
#ifdef TIMES
    times(&after);
    cout << "Triang. search time " << after.tms_utime - before.tms_utime << endl;
#endif

    if (Tind > 0)
    {
        //Lattice Position of Target
        xi = TargetPt.x;
        yi = TargetPt.y;
        //Lattice Position of Interpolant nearest Pts
        x[0] = INCR_coord[LatTri[Tind].p1].x;
        y[0] = INCR_coord[LatTri[Tind].p1].y;
        Incrx[0] = INCR_coord[LatTri[Tind].p1].Incrx;
        Incry[0] = INCR_coord[LatTri[Tind].p1].Incry;
        x[1] = INCR_coord[LatTri[Tind].p2].x;
        y[1] = INCR_coord[LatTri[Tind].p2].y;
        Incrx[1] = INCR_coord[LatTri[Tind].p2].Incrx;
        Incry[1] = INCR_coord[LatTri[Tind].p2].Incry;
        x[2] = INCR_coord[LatTri[Tind].p3].x;
        y[2] = INCR_coord[LatTri[Tind].p3].y;
        Incrx[2] = INCR_coord[LatTri[Tind].p3].Incrx;
        Incry[2] = INCR_coord[LatTri[Tind].p3].Incry;


        //Linear Interpolation (from MatLab griddata)
        //Barycentric Coord.
        //IMPORTANT: Need that (x[k],y[k]) define a triangle (del!=0)
        del = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
        w[2] = ((x[0] - xi) * (y[1] - yi) - (x[1] - xi) * (y[0] - yi)) / del;
        w[1] = ((x[2] - xi) * (y[0] - yi) - (x[0] - xi) * (y[2] - yi)) / del;
        w[0] = ((x[1] - xi) * (y[2] - yi) - (x[2] - xi) * (y[1] - yi)) / del;

        //Displacement Interpolation
        for (n = 0;n < 3;n++)
        {
            TargetPt.Incrx = TargetPt.Incrx + Incrx[n] * w[n];
            TargetPt.Incry = TargetPt.Incry + Incry[n] * w[n];
        }
    }
    else
    {
        //Use nearest neighbour interpolation
        Ptind = FindNearestPt(TargetPt);
        TargetPt.Incrx = INCR_coord[Ptind].Incrx;
        TargetPt.Incry = INCR_coord[Ptind].Incry;
    }
}

//////////////////////////////////////////////// TRIANGULATION HANDLINNG ////////////////////////////////////

/////////////Delaunay Triangulation of Lat. Pts
void ImUmbend::LatTriang(vector <ITRIANGLE> & LatTri)
{
    int i;
    int N = 0;
    int ntri = 0;
    int incr = 1;
    ITRIANGLE *v = NULL;
    XYZ *p = NULL;
    LatPoint * pLat = NULL;
    N = INCR_coord.size();
    p = (XYZ *) malloc((N + 3) * sizeof(XYZ));
    v = (ITRIANGLE *) malloc(3 * N * sizeof(ITRIANGLE));
    pLat = (LatPoint *)malloc(N * sizeof(LatPoint));



    //Sort vertixes in x increasing order
    for (i = 0;i < N;i++)
    {

        pLat[i] = INCR_coord[i];
    }
    qsort(pLat, N, sizeof(LatPoint), LatPtCompare);
    //Copy sorted values to INCR_coord vector
    for (i = 0;i < N;i++)
    {
        INCR_coord[i] = pLat[i];
    }


    //Set Lat. Pts Coord.
    for (i = 0;i < N;i++)
    {
        p[i].x = INCR_coord[i].x;
        p[i].y = INCR_coord[i].y;
        p[i].z = 0;
    }

    ///// Triangulation routine in DelTriang.h
    Triangulate(N, p, v, &ntri);


    //Copy Triangulation to vector format
    for (i = 0;i < ntri;i++)
    {
        LatTri.push_back(v[i]);
#undef DEBUG
        //  #define DEBUG
#ifdef DEBUG
        cout <<  i << " " << v[i].p1 << " " << v[i].p2 << " " << v[i].p3 << endl;
#endif
#undef DEBUG

    }



    //free memory
    free(pLat);
    free(v);
    free(p);

}

/////////////////// Computation of Nearest Triangl. to TargetPt
// Returns index of triangle in vector LatTri
int ImUmbend::FindNearestTri(LatPoint &TargetPt, vector <ITRIANGLE> &  LatTri)
{

    int k, Vind, t;
    XYZ Pi, P1, P2, P3;
    float x1, x2, x3, y1, y2, y3;
    int N = LatTri.size();

    //Target point position
    Pi.x = TargetPt.x;
    Pi.y = TargetPt.y;


    //Nearest Triangle Computation
    for (k = 0;k < N;k++)
    {

        //Triangle Vertex position
        Vind = LatTri[k].p1;
        P1.x = INCR_coord[Vind].x;
        P1.y = INCR_coord[Vind].y;
        Vind = LatTri[k].p2;
        P2.x = INCR_coord[Vind].x;
        P2.y = INCR_coord[Vind].y;
        Vind = LatTri[k].p3;
        P3.x = INCR_coord[Vind].x;
        P3.y = INCR_coord[Vind].y;

        //Inside Condition
        x1 = P1.x - Pi.x;
        y1 = P1.y - Pi.y;
        x2 = P2.x - Pi.x;
        y2 = P2.y - Pi.y;
        x3 = P3.x - Pi.x;
        y3 = P3.y - Pi.y;
        t = (x1 * y2 > x2 * y1) + (x2 * y3 > x3 * y2) + (x3 * y1 > x1 * y3);

        // printf("Inside Condition of Triangle %d is %d:\n",k,t);

        if ((t == 3) || (t == 0)) return k;
    }

    return -1;
}

/////////////////// Computation of Nearest Point to TargetPt
// Returns index of nearest Vertex in vector INCR_coord
int ImUmbend::FindNearestPt(LatPoint &TargetPt)
{
    int k, indMin;
    int N = INCR_coord.size();
    float Dist, MinDist;

    MinDist = -1;

    //Nearest Pt is the one achieving minimum distance
    for (k = 0;k < N;k++)
    {
        //Distance Computation
        Dist = (TargetPt.x - INCR_coord[k].x) * (TargetPt.x - INCR_coord[k].x) + (TargetPt.y - INCR_coord[k].y) * (TargetPt.y - INCR_coord[k].y);
        if (Dist < MinDist)
        {
            MinDist = Dist;
            indMin = k;
        }

    }


    return indMin;

}
