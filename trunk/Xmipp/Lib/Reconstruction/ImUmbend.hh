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


//////////////////////////// SPECIFIC LIBRARIES ////////////////////////

#include <XmippInterface/xmippCCLattice_IO.hh>
//#include <../Tools.h>
#include <XmippData/xmippDelTriang.hh>


/////////////////////////// CONSTANTS //////////////////////////////////

#define BORDER_POINT 0


/////////////////////////////// DATA TYPES //////////////////////////////


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

// Read *.cor file
void ReadMRCCord();

//////////////////////////  Peaks Correspondance
/**Compute diference between ideal red spots and experimental ones */
void PeaksCorresp();

/////////////////////////// Unbending 

//Image Unbending
void UnBending( );

///////////////////////////// Interpolation

//Bilinear Interpolation from Regular grid
void ShiftsInterpReg(matrix2D <double> & ,matrix2D <double> & ,LatPoint & );

//Linear Interpolation from scattered data set to regular grid
void Scattered2Regular(matrix2D <double> & ,matrix2D <double> & ,vector <ITRIANGLE> & );

//Displacement Interpolation
void ShiftsInterp( LatPoint &, vector <ITRIANGLE> & );

////////////////////////////  Triangulation

//Lattice Triangulation (using DelTriang.hh)
void LatTriang(vector <ITRIANGLE> & );

// Returns index of triangle in vector LatTri
int FindNearest(LatPoint &,  vector <ITRIANGLE> & );

};
   /* Show parameters --------------------------------------------------------- */
   ostream & operator << (ostream &o, const ImUmbend &prm) {
      o << "Input Correlation File      : " << prm.FN_Correlation << endl
        << "Correlation_Peak_Threshold               : " << prm.cc_peak_factor << endl
	<< "Input Image           :"<< prm.inImfile<< endl
        ;
      return o;
   };

  
//@}

#endif
