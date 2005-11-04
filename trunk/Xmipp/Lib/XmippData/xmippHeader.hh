/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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
/* ------------------------------------------------------------------------- */
/* XMIPP HEADER                                                                    */
/* ------------------------------------------------------------------------- */
/*

   Valid operations:
   
   headerXMipp A(IMG_XMIPP);              // Creates an XMipp Image header
   headerXMipp A(VOL_XMIPP);              // Creates an XMipp volume header

   A.rename("newName");		          // Renames the header title

   A.write(fp);			          // Writes current header to disk
   A.read(fp);			  	  // Loads a header from disk	

   A.clear();				  // Clears the header

   cout << A;                             // Prints information in the header   
   
   set_eulerAngles(25, 40, 50);           // Sets Euler Angles
   get_eulerAngles(Phi, Theta, Psi);      // Gets Euler Angles
   get_eulerAngles1(Phi1, Theta1, Psi1);      // Gets Euler Angles stored in alternate location 1
   get_eulerAngles2(Phi2, Theta2, Psi2);      // Gets Euler Angles stored in alternate location 2
   clear_fFlag_flag();                    // Clears fFlag flag
   
   rotAngle() = 30;   			  // Sets Rotation Angle to 30
   get_title();                           // Gets the title
   get_time();				  // Gets the time of creation 
   get_date();				  // Gets the date of creation
   
  
*/

#ifndef _XMIPPHEADER_H
#define _XMIPPHEADER_H

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */
#include <iostream>
#include <stdio.h>
#include "xmippFuncs.hh"
#include "xmippMatrices2D.hh"

/* ************************************************************************* */
/* CLASS DEFINITION AND PROTOTYPES                                           */
/* ************************************************************************* */
// The header class is a general class which only contains general information
// for the Xmipp image. 

class headerXmipp {

public:

  //	Types of initialization of the Header
  
  typedef enum { IMG_BYTE = 0, IMG_XMIPP  = 1,
                 IMG_INT  = 9, /*IMG_SHORT  = , IMG_UCHAR = ,*/
                 VOL_BYTE = 2,  VOL_XMIPP  = 3, 
                 VOL_INT  = 10 , /*VOL_SHORT  = , VOL_UCHAR =, */
  		 IMG_FOURIER = -1, VOL_FOURIER = -3} initMode;

   // Constructors ..............................................

   headerXmipp(initMode _im = IMG_BYTE){clear(); im =_im;};   // constructor
                                                 
   // Some operators streams .......................................................
   friend ostream& operator << (ostream& o, const headerXmipp &I);  // cout << I
   void print_hard(ostream &o) const;

   // Input/Output
   
   // reversed is only used in case that the type_check is skipped
   int read(FILE *fp, bool skip_type_check=false, bool reversed=false,
      bool skip_extra_checkings=false);
   int read(const FileName &fn, bool skip_type_check=false,
      bool reversed=false, bool skip_extra_checkings=false);
   void write(FILE *fp, bool force_reversed=false);
   void write(const FileName &fn, bool force_reversed=false);

   // Other header operators
   
   initMode& headerType() {return im;};        // Set type of header
   initMode headerType() const {return im;};   // Get type of header
   
   void clear();
   void set_header();                // Sets a consistent header 

   // Interaction with data
   bool  reversed() const {return __reversed;}
   bool& reversed() {return __reversed;}

  // Dimension of each slice (image)

   void set_dimension(float Ydim, float Xdim);
   void get_dimension(float &Ydim, float &Xdim) const;
   
   int get_header_size() const
      {return (int)header.fNcol*(int)header.fLabrec*sizeof(float);}

   float& Ydim() {return header.fNrow;};      		// Set Ydim
   int iYdim() const 
    {return (int) header.fNrow;};   		// Get Ydim
   float& Xdim() {return header.fNcol;};      		// Set Xdim
   int  iXdim() const 
    {return (int) header.fNcol;};   		// Get Xdim

  // Number of Slices  
   float& Slices() {return header.fNslice;};          // Set Slices
   float  Slices()  const {return header.fNslice;};   // Get Slices
   int    iSlices() const {return (int) header.fNslice;};// Get Slices
   int    iZdim()   const {return (int) header.fNslice;}

  // Rotation Angle 
   float  old_rot() const {return header.fAngle1;};   // Get Rot. Angle
   float& old_rot(){return header.fAngle1;};          // Set Rot. Angle

  // Scale  
   float  Scale() const {return header.fScale;};       // Get Scale
   float& Scale()       {return header.fScale;}   

   // For Xmipp Maximum-Likelihood refinement
   float  Flip()     const {return header.Flip;}
   float& Flip()           {return header.Flip;}
   float  Weight()     const {return header.Weight;}
   float& Weight()           {return header.Weight;}

   // Other ugly fields
   float  fNrec()   const {return header.fNrec;}
   float& fNrec()         {return header.fNrec;}
   float  fNlabel() const {return header.fNlabel;}
   float& fNlabel()       {return header.fNlabel;}
   float  fIform()  const {return header.fIform;}
   float& fIform()        {return header.fIform;}
   float  fImami()  const {return header.fImami;}
   float& fImami()        {return header.fImami;}
   float  fFmax()   const {return header.fFmax;}
   float& fFmax()         {return header.fFmax;}
   float  fFmin()   const {return header.fFmin;}
   float& fFmin()         {return header.fFmin;}
   float  fAv()     const {return header.fAv;}
   float& fAv()           {return header.fAv;}
   float  fSig()    const {return header.fSig;}
   float& fSig()          {return header.fSig;}
   float  fIhist()  const {return header.fIhist;}
   float& fIhist()        {return header.fIhist;}
   float  fLabrec() const {return header.fLabrec;}
   float& fLabrec()       {return header.fLabrec;}
   float  fIangle() const {return header.fIangle;}
   float& fIangle()       {return header.fIangle;}
   float  fXoff()   const {return header.fXoff;}
   float& fXoff()         {return header.fXoff;}
   float  fYoff()   const {return header.fYoff;}
   float& fYoff()         {return header.fYoff;}
   float  fZoff()   const {return header.fZoff;}
   float& fZoff()         {return header.fZoff;}
   float  fLabbyt() const {return header.fLabbyt;}
   float& fLabbyt()       {return header.fLabbyt;}
   float  fLenbyt() const {return header.fLenbyt;}
   float& fLenbyt()       {return header.fLenbyt;}
   matrix2D<double> fGeo_matrix();

  // Origin offsets
  void set_originOffsets(float Xoff, float Yoff);
  void get_originOffsets(float &Xoff, float &Yoff) const;

  // Euler angles
   void set_eulerAngles(float Phi, float Theta, float Psi);
   void set_eulerAngles1(float Phi1, float Theta1, float Psi1);
   void set_eulerAngles2(float Phi2, float Theta2, float Psi2);
/** Clears fFlag flag. The number of triads of Euler angles stored in the header (up to three) is stored here. set_eulerAngles2 makes fFlag=2, set_eulerAngles1 makes fFlag=max(fFlag,1), set_eulerAngles does not change fFlag
        \\Ex: header.clear_fFlag_flag(); */
   void clear_fFlag_flag() {header.fFlag=0.f;};
   template <class T>
      void get_eulerAngles(T &Phi, T &Theta, T &Psi) const {
        Phi = (T) header.fPhi;
        Theta = (T) header.fTheta;
        Psi = (T) header.fPsi;
      }

   template <class T>
      void get_eulerAngles1(T &Phi1, T &Theta1, T &Psi1) const {
        Phi1 = (T) header.fPhi1;
        Theta1 = (T) header.fTheta1;
        Psi1 = (T) header.fPsi1;
      }
      
   template <class T>
      void get_eulerAngles2(T &Phi2, T &Theta2, T &Psi2) const {
        Phi2 = (T) header.fPhi2;
        Theta2 = (T) header.fTheta2;
        Psi2 = (T) header.fPsi2;
      }
   
   float& Phi() {header.fIangle = 1; return header.fPhi;};          // Set Phi
   float  Phi() const {return header.fPhi;};      // Get Phi
   float& Theta() {header.fIangle = 1; return header.fTheta;};      // Set Theta
   float  Theta() const {return header.fTheta;};  // Get Theta
   float& Psi() {header.fIangle = 1; return header.fPsi;};          // Set Psi
   float  Psi() const {return header.fPsi;};      // Get Psi

   float& Phi1() {header.fFlag = 1.f; return header.fPhi1;};      // Set Phi
   float  Phi1() const {return header.fPhi1;};      // Get Phi
   float& Theta1() {header.fFlag = 1.f; return header.fTheta1;};  // Set Theta
   float  Theta1() const {return header.fTheta1;};  // Get Theta
   float& Psi1() {header.fFlag = 1.f; return header.fPsi1;};      // Set Psi
   float  Psi1() const {return header.fPsi1;};      // Get Psi

   float& Phi2() {header.fFlag = 2.f; return header.fPhi2;};      // Set Phi
   float  Phi2() const {return header.fPhi2;};      // Get Phi
   float& Theta2() {header.fFlag = 2.f; return header.fTheta2;};  // Set Theta
   float  Theta2() const {return header.fTheta2;};  // Get Theta
   float& Psi2() {header.fFlag = 2.f; return header.fPsi2;};      // Set Psi
   float  Psi2() const {return header.fPsi2;};      // Get Psi
   
   float   Is_flag_set(void){ return(header.fFlag);}
  // Date and Time

   char* get_date() const;
   char* get_time() const;
   void set_date();
   void set_time();
   
  // Title 
   char* get_title() const;
   void  set_title(FileName newName);

private:

   /* Spider header defined as Private because it must be encapsulated 
      to avoid missuse of it*/
      
   typedef struct {
    /* 1 */  float fNslice;  // NUMBER OF SLICES (PLANES) IN VOLUME        
	                     // (=1 FOR AN IMAGE)  FOR NEW LONG LABEL    
      		             // FORMAT THE VALUE OF NSLICE STORED IN     
		             // THE FILE IS NEGATIVE.                       
    /* 2 */  float fNrow;    // NUMBER OF ROWS PER SLICE (Y)                   
    /* 3 */  float fNrec;    // TOTAL NUMBER OF RECORDS (SEE NOTE #3).   
    /* 4 */  float fNlabel;  // AUXILIARY NUMBER TO COMPUTE TOTAL NUMBER OF RECS
    /* 5 */  float fIform;   // FILE TYPE SPECIFIER.                    
			     // +3 FOR A 3-D FILE  (FLOAT)                   
			     // +1 FOR A 2-D IMAGE (FLOAT)                      
			     // -1 FOR A 2-D FOURIER TRANSFORM 
			     // -3 FOR A 3-D FOURIER TRANSFORM 
			     // -5 FOR A NEW 2-D FOURIER TRANSFORM   
			     // -7 FOR A NEW 3-D FOURIER TRANSFORM
			     // +8 FOR A 2-D EIGHT BIT IMAGE FILE
       			     // +9 FOR A 2-D INT IMAGE FILE
                             // 10 FOR A 3-D INT IMAGE FILE
			     // 11 FOR A 2-D EIGHT BIT COLOR IMAGE FILE                                  
    /* 6 */  float fImami;   // MAXIMUM/MINIMUM FLAG. IS SET AT 0 WHEN THE
    			     // FILE IS CREATED, AND AT 1 WHEN THE MAXIMUM AND 
			     // MINIMUM HAVE BEEN COMPUTED, AND HAVE BEEN STORED 
			     // INTO THIS LABEL RECORD (SEE FOLLOWING WORDS)
    /* 7 */  float fFmax;    // MAXIMUM VALUE
    /* 8 */  float fFmin;    // MINIMUM VALUE
    /* 9 */  float fAv;      // AVERAGE VALUE
    /* 10*/  float fSig;     // STANDARD DEVIATION. A VALUE OF -1. INDICATES 
    			     // THAT SIG HAS NOT BEEN COMPUTED PREVIOUSLY.
    /* 11*/  float fIhist;   // FLAG INDICATING IF THE HISTOGRAM HAS BE 
    			     // COMPUTED. NOT USED IN 3D FILES!
    /* 12*/  float fNcol;    // NUMBER OF PIXELS PER LINE (Columns X)
    /* 13*/  float fLabrec;  // NUMBER OF LABEL RECORDS IN FILE HEADER
    /* 14*/  float fIangle;  // FLAG THAT TILT ANGLES HAVE BEEN FILLED
    /* 15*/  float fPhi;     // EULER: ROTATIONAL ANGLE
    /* 16*/  float fTheta;   // EULER: TILT ANGLE
    /* 17*/  float fPsi;     // EULER: PSI  = TILT ANGLE
    /* 18*/  float fXoff;    // X TRANSLATION
    /* 19*/  float fYoff;    // Y TRANSLATION
    /* 20*/  float fZoff;    // Z TRANSLATION
    /* 21*/  float fScale;   // SCALE
    /* 22*/  float fLabbyt;  // TOTAL NUMBER OF BYTES IN LABEL
    /* 23*/  float fLenbyt;  // RECORD LENGTH IN BYTES
             char  fNada[24];// this is a spider incongruence
    /* 30*/  float fFlag;    // THAT ANGLES ARE SET. 1 = ONE ADDITIONAL 
    			     //	ROTATION IS PRESENT, 2 = ADDITIONAL ROTATION 
			     // THAT PRECEEDS THE ROTATION THAT WAS STORED IN
			     // 15 FOR DETAILS SEE MANUAL CHAPTER VOCEUL.MAN
    /* 31*/  float fPhi1;
    /* 32*/  float fTheta1;
    /* 33*/  float fPsi1;
    /* 34*/  float fPhi2;
    /* 35*/  float fTheta2;
    /* 36*/  float fPsi2;

             double fGeo_matrix[3][3];  // x9 = 72 bytes: Geometric info
             float fAngle1;             // angle info

             float fr1;
             float fr2;                 // lift up cosine mask parameters

	    /** Fraga 23/05/97  For Radon transforms **/

             float RTflag;              // 1=RT, 2=FFT(RT)
             float Astart;
             float Aend;
             float Ainc;
             float Rsigma;  		// 4*7 = 28 bytes
             float Tstart;
             float Tend;
             float Tinc;   		// 4*3 = 12, 12+28 = 40B

             /** Sjors Scheres 17/12/04 **/

             float Weight;          // For Maximum-Likelihood refinement
             float Flip;            // 0=no flipping operation (false), 1=flipping (true)

             char  fNada2[576];         // empty     700-76-40=624-40-8= 576 bytes

 /*212-214*/ char szIDat[12];   // LOGICAL * 1 ARRAY DIMENSIONED 10, CONTAINING 
 			        // THE DATE OF CREATION (10 CHARS)
 /*215-216*/ char szITim[8];    // LOGICAL * 1 ARRAY DIMENSIONED 8, CONTAINING 
 			        // THE TIME OF CREATION (8 CHARS)
 /*217-256*/ char szITit[160];  // LOGICAL * 1 ARRAY DIMENSIONED 160
     
   } SpiderHeader;

   SpiderHeader header;
   initMode im;
   bool     __reversed;
};
#endif
