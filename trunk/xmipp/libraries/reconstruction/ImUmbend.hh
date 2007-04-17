/***************************************************************************
 *
 * Authors:     Debora Gil
                Roberto Marabini
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

#ifndef _IMUMBEND_HH
#  define _IMUMBEND_HH

///////////////////////////// COMMON LIBRARIES /////////////////////////
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippMatrices2D.hh>
#include <XmippData/xmippImages.hh>
#include <XmippData/xmippProjection.hh>
#include <XmippInterface/xmippAPHorigmerg.hh>
#include <XmippData/xmippMacros.hh>

//////////////////////////// SPECIFIC LIBRARIES ////////////////////////

#include <XmippInterface/xmippCCLattice_IO.hh>
//#include <../Tools.h>
#include <XmippData/xmippDelTriang.hh>


/////////////////////////// CONSTANTS //////////////////////////////////

#define BORDER_POINT 0
//Bessel Interpolation
#define SIG 3
static double B[] = {0.0046,
    0.0046,
    0.0046,
    0.0046,
    0.0046,
    0.0045,
    0.0045,
    0.0045,
    0.0045,
    0.0045,
    0.0044,
    0.0044,
    0.0044,
    0.0043,
    0.0043,
    0.0043,
    0.0042,
    0.0042,
    0.0041,
    0.0041,
    0.0041,
    0.0040,
    0.0040,
    0.0039,
    0.0039,
    0.0039,
    0.0038,
    0.0038,
    0.0037,
    0.0037,
    0.0037,
    0.0036,
    0.0036,
    0.0035,
    0.0035,
    0.0035,
    0.0034,
    0.0034,
    0.0034,
    0.0033,
    0.0033,
    0.0032,
    0.0032,
    0.0032,
    0.0031,
    0.0031,
    0.0031,
    0.0030,
    0.0030,
    0.0030,
    0.0029,
    0.0029,
    0.0029,
    0.0028,
    0.0028,
    0.0028,
    0.0028,
    0.0027,
    0.0027,
    0.0027,
    0.0026,
    0.0026,
    0.0026,
    0.0026,
    0.0025,
    0.0025,
    0.0025,
    0.0024,
    0.0024,
    0.0024,
    0.0024,
    0.0023,
    0.0023,
    0.0023,
    0.0023,
    0.0022,
    0.0022,
    0.0022,
    0.0022,
    0.0021,
    0.0021,
    0.0021,
    0.0021,
    0.0020,
    0.0020,
    0.0020,
    0.0020,
    0.0020,
    0.0019,
    0.0019,
    0.0019,
    0.0019,
    0.0018,
    0.0018,
    0.0018,
    0.0018,
    0.0018,
    0.0017,
    0.0017,
    0.0017,
    0.0017,
    0.0017,
    0.0016,
    0.0016,
    0.0016,
    0.0016,
    0.0016,
    0.0016,
    0.0015,
    0.0015,
    0.0015,
    0.0015,
    0.0015,
    0.0014,
    0.0014,
    0.0014,
    0.0014,
    0.0014,
    0.0014,
    0.0013,
    0.0013,
    0.0013,
    0.0013,
    0.0013,
    0.0013,
    0.0013,
    0.0012,
    0.0012,
    0.0012,
    0.0012,
    0.0012,
    0.0012,
    0.0011,
    0.0011,
    0.0011,
    0.0011,
    0.0011,
    0.0011,
    0.0011,
    0.0011,
    0.0010,
    0.0010,
    0.0010,
    0.0010,
    0.0010,
    0.0010,
    0.0010,
    0.0010,
    0.0009,
    0.0009,
    0.0009,
    0.0009,
    0.0009,
    0.0009,
    0.0009,
    0.0009,
    0.0008,
    0.0008,
    0.0008,
    0.0008,
    0.0008,
    0.0008,
    0.0008,
    0.0008,
    0.0008,
    0.0008,
    0.0007,
    0.0007,
    0.0007,
    0.0007,
    0.0007,
    0.0007,
    0.0007,
    0.0007,
    0.0007,
    0.0007,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0006,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0005,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0004,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0003,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0002,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0001,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
         0};

/////////////////////////////// DATA TYPES //////////////////////////////

/**@name ImUmbend */

// Struct to store corelation peak related information
typedef struct{
	//Lattice indexes
	int i,j;
	//Point Position
	float x,y;
	//Deviation from Theoric Position
	float Incrx,Incry;
	//Interpolation flag
	bool Interp;
}LatPoint;


//@{
// ImUmbend structure ----------------------------------------------------------
/** ImUmbend class.
    ImUmbend unbends a 2D crystal image.
*/

class ImUmbend
{

 public:

   //////////////////////////////// Umbend Input Parameters
   /** Input file (MRC cor) */
   FileName FN_Correlation;
   /** Input Image ( SPI file) */
   char * inImfile;
   /** Output Image **/
   char * outImfile;

   /** THRESHOLD OF CROSS-CORRELATION PEAK HEIGHT CALCULATED AS;
       max of croscorrelation * FACTOR (READ FROM UNIT 5)*/
   double cc_peak_factor;

   /** Interpolation Model for Extension of Experimental Shifts to whole crystal image*/
   string  InterpModel;

  //////////////////////////////// Umbend Variables

  /** Experimental Lattice and Peaks Coordinates */
   CCLattice_IO ExpLat;


  /**  Experimental Displacements */
  vector <LatPoint> INCR_coord;


  /** Crystal Image  */
  ImageXmipp    inIm,outIm;

////////////////////////////////  FUNCTIONS //////////////////////////

public:

////////////////////////////  Constructors




//////////////////////////// I/O functions

/** Read *.cor file */
void ReadMRCCord();

//////////////////////////  Peaks Correspondance
/**Compute diference between ideal red spots and experimental ones */
void PeaksCorresp();

/////////////////////////// Unbending

/**Image Unbending*/
void UnBending( );

///////////////////////////// Interpolation

/** Shifts Interpolation from Regular grid*/
void ShiftsInterpReg(matrix2D <double> & ,matrix2D <double> & ,LatPoint & );
/** 2D Interpolation on Square acoording to InterpModel*/
void Interp2D(float Tx,float Ty,float Ti,float Tj,float TiM,float TjM, float * ACoeff);

/**Linear Interpolation from scattered data set to regular grid*/
void Scattered2Regular(matrix2D <double> & ,matrix2D <double> & ,vector <ITRIANGLE> & );
/**Displacement Interpolation from Triangulation of Irregular Grid*/
void ShiftsInterp( LatPoint &, vector <ITRIANGLE> & );

////////////////////////////  Triangulation

/**Lattice Triangulation (using DelTriang.hh)*/
void LatTriang(vector <ITRIANGLE> & );

/**Returns index of triangle in vector LatTri*/
int FindNearestTri(LatPoint &,  vector <ITRIANGLE> & );

/** Returns index of nearest Vertex in vector INCR_coord*/
int FindNearestPt(LatPoint &);
///////////////////////////  I/O
/** Show parameters */
   friend ostream & operator << (ostream &o, const ImUmbend &prm) {
      o << "Input Correlation File      : " << prm.FN_Correlation << endl
        << "Correlation_Peak_Threshold               : " << prm.cc_peak_factor << endl
	<< "Input Image           :"<< prm.inImfile<< endl
        ;
      return o;
   };

//@}
};


#endif
