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
#include "../xmippHeader.hh"
#include <sys/stat.h>
#include <time.h>


/* Some operators *******************************************************/

ostream& operator << (ostream& o, const headerXmipp &I) {
   o << "Image type   : ";
    switch ((int) I.header.fIform) {
      case headerXmipp::IMG_BYTE:          // for a 2D Byte Image
	o << "2D Byte image";
        break;
      case headerXmipp::IMG_XMIPP:         // for a 2D Xmipp image. 
	o << "2D Xmipp image";
        break;
      case headerXmipp::IMG_INT:         // for a 2D Xmipp image. 
	o << "2D INT_Xmipp image";
        break;
//      case headerXmipp::IMG_SHORT:         // for a 2D Xmipp image. 
//	o << "2D SHORT_Xmipp image";
//        break;
//      case headerXmipp::IMG_UCHAR:         // for a 2D Xmipp image. 
//	o << "2D UCHAR_Xmipp image";
//        break;
      case headerXmipp::VOL_BYTE:	   // for a 3D Byte volume.
	o << "3D Byte volume";
        break;
      case headerXmipp::VOL_XMIPP:	   // for a 3D Xmipp volume.	
	o << "3D Xmipp volume";
        break;
      case headerXmipp::VOL_INT:	   // for a 3D Xmipp volume.	
	o << "3D INT_Xmipp volume";
        break;
//      case headerXmipp::VOL_SHORT:	   // for a 3D Xmipp volume.	
//	o << "3D SHORT_Xmipp volume";
//        break;
//      case headerXmipp::VOL_UCHAR:	   // for a 3D Xmipp volume.	
//	o << "3D UCHAR_Xmipp volume";
//       break;
      case headerXmipp::IMG_FOURIER:	   // for a 2D Fourier transform.	
	o << "2D Fourier image";
        break;
      case headerXmipp::VOL_FOURIER:       // for a 3D Fourier transform.
	o << "3D Fourier volume";
        break;
    }
    
   o << endl;
   o << "Reversed     : ";
     if (I.__reversed) o << "TRUE"  << endl;
     else              o << "FALSE" << endl;
   o << "dimensions   : " << I.header.fNslice << " x " << I.header.fNrow << " x " << I.header.fNcol;
   o << "  (slices x rows x columns)" << endl;
   o << "Euler angles : " << endl;
   o << "  Phi   (rotation around Z axis) = " << I.header.fPhi << endl;
   o << "  theta (tilt, second rotation around new Y axis) = " << I.header.fTheta << endl;
   o << "  Psi   (third rotation around new Z axis) = " << I.header.fPsi << endl;
  if(I.header.fFlag == 1.0f || I.header.fFlag == 2.0f ){ 
     o << "  Phi1   = " << I.header.fPhi1 ;
     o << "  theta1 = " << I.header.fTheta1 ;
     o << "  Psi1   = " << I.header.fPsi1 << endl;}
  if(I.header.fFlag == 2.0f) {
     o << "  Phi2   = " << I.header.fPhi2 ;
     o << "  theta2 = " << I.header.fTheta2 ;
     o << "  Psi2   = " << I.header.fPsi2 << endl;}

   o << "Date         : " << I.get_date() << endl; 
   o << "Time         : " << I.get_time() << endl; 
   o << "Title        : " << I.get_title() << endl;
   o << "Header size  : " << I.get_header_size() << endl;
   return o;
}  

void headerXmipp::print_hard(ostream &o) const {
   o << "fNslice=" << header.fNslice << endl;
   o << "fNrow=" <<   header.fNrow << endl;
   o << "fNrec=" <<   header.fNrec << endl;
   o << "fNlabel=" << header.fNlabel << endl;
   o << "fIform=" <<  header.fIform << endl;
   o << "fImami=" <<  header.fImami << endl;
   o << "fFmax=" <<   header.fFmax << endl;
   o << "fFmin=" <<   header.fFmin << endl;
   o << "fAv=" <<     header.fAv << endl;
   o << "fSig=" <<    header.fSig << endl;
   o << "fIhist=" <<  header.fIhist << endl;
   o << "fNcol=" <<   header.fNcol << endl;
   o << "fLabrec=" << header.fLabrec << endl;
   o << "fIangle=" << header.fIangle << endl;
   o << "fPhi=" <<    header.fPhi << endl;
   o << "fTheta=" <<  header.fTheta << endl;
   o << "fPsi=" <<    header.fPsi << endl;
   o << "fXoff=" <<   header.fXoff << endl;
   o << "fYoff=" <<   header.fYoff << endl;
   o << "fZoff=" <<   header.fZoff << endl;
   o << "fScale=" <<  header.fScale << endl;
   o << "fLabbyt=" << header.fLabbyt << endl;
   o << "fLenbyt=" << header.fLenbyt << endl;
}
   
/* Input (read) *******************************************************/
int headerXmipp::read(FILE *fp, bool skip_type_check, bool force_reversed,
   bool skip_extra_checkings) {
  float tmp;
  unsigned tmpSize;       
  struct stat info;
  unsigned long size;
  
   /* Determine reverse status according to the next table .................
      (computed in Linux) */
  #define TYPE_TABLE_SIZE 10
  int type_table[TYPE_TABLE_SIZE][5]={
      0, 0,  48,  65, 11,
      0, 0,  32,  65, 10,
      0, 0,  16,  65,  9,
      0, 0,   0,  65,  8,
      0, 0,  64,  64,  3,
      0, 0, 128,  63,  1,
      0, 0, 128, 191, -1,
      0, 0,  64, 192, -3,
      0, 0, 160, 192, -5,
      0, 0, 224, 192, -7};
  union {
     unsigned char c[4];
     float         f;
  } file_type;
  
   if (!skip_type_check) {
     // Read file type
     fseek(fp, 16, SEEK_SET);
     for (int i=0; i<4; i++)
         fread(&(file_type.c[i]), sizeof(unsigned char), 1, fp);
     fseek(fp,  0, SEEK_SET);

     // Select type
     int correct_type=0;
     __reversed=FALSE;
     #define IS_TYPE(n) \
      (type_table[n][0]==file_type.c[0] && type_table[n][1]==file_type.c[1] && \
         type_table[n][2]==file_type.c[2] && type_table[n][3]==file_type.c[3])
     #define IS_REVERSE_TYPE(n) \
      (type_table[n][0]==file_type.c[3] && type_table[n][1]==file_type.c[2] && \
         type_table[n][2]==file_type.c[1] && type_table[n][3]==file_type.c[0])
   
     for (int i=0; i<TYPE_TABLE_SIZE; i++)
         if      (IS_TYPE(i))
            {correct_type=type_table[i][4]; break;}
         else if (IS_REVERSE_TYPE(i))
            {correct_type=type_table[i][4]; __reversed=TRUE; break;}
     if (correct_type==0) return FALSE;
     
     // Now check this machine type
     file_type.f=1; if (IS_REVERSE_TYPE(5)) __reversed=!__reversed;
  } else __reversed=force_reversed;

  // Read header
  if (!__reversed)
       fread(&header, sizeof(headerXmipp::SpiderHeader), 1, fp);
  else {
       FREAD(&header,             sizeof(float),  36, fp, TRUE);
       FREAD(&header.fGeo_matrix, sizeof(double),  9, fp, TRUE);
       FREAD(&header.fAngle1,     sizeof(float),  11, fp, TRUE);
       FREAD(&header.fNada2,      sizeof(char),  764, fp, TRUE); 
  }

  if (!skip_extra_checkings)
     if (!fstat(fileno(fp), &info)) {
       switch (im) {
       case IMG_XMIPP: 
         size = (unsigned long) (get_header_size() +
             		         header.fNcol*header.fNrow*sizeof(float));
         if ((size != info.st_size) || (header.fIform != 1)) return FALSE;
         else if (skip_type_check) header.fIform=1; // This is done to recover
                                                    // files which are not
                                                    // properly converted from
                                                    // other packages
         break;
       case IMG_INT: 
         size = (unsigned long) (get_header_size() +
             		         header.fNcol*header.fNrow*sizeof(float));
         if ((size != info.st_size) || (header.fIform != 9)) return FALSE;
         else if (skip_type_check) header.fIform=9; // This is done to recover
                                                    // files which are not
                                                    // properly converted from
                                                    // other packages
         break;
       case VOL_BYTE:
         size = (unsigned long) (header.fNslice*header.fNcol*header.fNrow*sizeof(float));
         if ((size != info.st_size)) return FALSE;
         break;
       case VOL_XMIPP: 
         size = (unsigned long) (get_header_size() +
             		         header.fNslice*header.fNcol*header.fNrow*sizeof(float));
         if ((size != info.st_size) || (header.fIform != 3)) return FALSE;
         else if (skip_type_check) header.fIform=3;
         break;
       case VOL_INT: 
         size = (unsigned long) (get_header_size() +
             		         header.fNslice*header.fNcol*header.fNrow*sizeof(float));
         if ((size != info.st_size) || (header.fIform != 10)) return FALSE;
         else if (skip_type_check) header.fIform=10;
         break;
       case IMG_FOURIER: 
         size = (unsigned long) (get_header_size() +
             		         2*header.fNcol*header.fNrow*sizeof(float)); 
				 // The term 2 is to take into account that IMG_FOURIER 
				 // stores complex numbers with 2 floats for each one.
         if ((size != info.st_size) ||
             (header.fIform != -5 && header.fIform != -1)) return FALSE;
         else if (skip_type_check) header.fIform=-1;
         break;
       case VOL_FOURIER: 
         size = (unsigned long) (get_header_size() +
             		         2*header.fNslice*header.fNcol*header.fNrow*sizeof(float));
				 // The term 2 is to take into account that VOL_FOURIER 
				 // stores complex numbers with 2 floats for each one.
         if ((size != info.st_size) ||
             (header.fIform != -7 && header.fIform != -3)) return FALSE;
         else if (skip_type_check) header.fIform=-3;
         break;
       }    
     } else return FALSE;

  /* This is Spider stuff, it's an acient trick.
     The only thing we need to know is that Spider images
     contains a header with the size of the struct SpiderHeader, 
     a "filling" empty space and the pixel data of size
     cols*rows*sizeof(float).  
  */
  
  // Now read and throw empty filling space
//  header.fLabrec = (float) ceil((float) 256/header.fNcol); 
  tmpSize = (int)(header.fNcol * header.fLabrec * 4); //Size of whole header
  tmpSize -= sizeof(headerXmipp::SpiderHeader);      	      //Decrease the real header  
  for (unsigned i = 0; i < tmpSize/4; i++)
      FREAD(&tmp, sizeof(float), 1, fp, __reversed);
  return TRUE;
}

/* Input (read) *******************************************************/

int headerXmipp::read(const FileName &name, bool skip_type_check,
  bool force_reversed, bool skip_extra_checkings) _THROW {
  FILE *fp; 

   // Clear Header first
   clear(); 

   // Open file
   if ((fp = fopen(name.c_str(), "rb")) == NULL)
     REPORT_ERROR(1501,"HeaderXmipp::read: File " + name + " not found");

   int retval=read(fp, skip_type_check, force_reversed, skip_extra_checkings);
   fclose(fp);
   return retval;
}
  
/* Output (write) *******************************************************/
void headerXmipp::write(FILE *fp, bool force_reversed) {
  float tmp;
  unsigned tmpSize;
       
  if (Ydim()==0 || Xdim()==0 || Slices()==0) return;
  
  // Set consistent header before saving
  set_header();
    
  // Since we are writing a new header set time and date
  set_time();
  set_date();

  // Write header
  if (XOR(__reversed,force_reversed)) {
     __reversed=TRUE;
     FWRITE(&header,             sizeof(float),  36, fp, TRUE);
     FWRITE(&header.fGeo_matrix, sizeof(double),  9, fp, TRUE);
     FWRITE(&header.fAngle1,     sizeof(float),  11, fp, TRUE);
     fwrite(&header.fNada2,      sizeof(char),  764, fp); 
  } else {
     __reversed=FALSE;
     fwrite(&header, sizeof(headerXmipp::SpiderHeader), 1, fp);
  }

  /* This is Spider stuff, it's an acient trick.
     The only thing we need to know is that Spider images
     contains a header with the size of the struct SpiderHeader, 
     a "filling" empty space and the pixel data of size
     cols*rows*sizeof(float).  
  */
  
  // Now write empty filling space (filled with zeros)
  tmpSize = (int)(header.fNcol * header.fLabrec * 4); //Size of whole header
  tmpSize -= sizeof(headerXmipp::SpiderHeader);             //Decrease the real header  
  tmp = (float) 0.0;
  for (unsigned i = 0; i < tmpSize/4; i++) {
    fwrite(&tmp, sizeof(float), 1, fp);
  }  
}

/* Clear header *******************************************************/
void headerXmipp::clear() {
  /* Set all header to zero */
  memset(&header, 0, sizeof(headerXmipp::SpiderHeader));    
  set_header(); // Set a consistent header necessary for XMipp
  __reversed=FALSE;
}

/* Set header *******************************************************/
void headerXmipp::set_header() {
int Rows, Cols;
int headrec;
  /* Compute "dark" stuff needed in the header.
     Things were defined in acient time for ancient computers
     and operating systems so never try to understand it. ;-) 
  */
  
  if ((header.fNcol != 0) && (header.fNrow != 0)) {
    header.fNlabel =  (float) ((int) (256/header.fNcol+1)); 
    header.fLabrec = (float) ceil((float) (256/(float)header.fNcol)); 
    
    headrec = (int) 1024 / ((int)header.fNcol * 4); // temporary variable
  
    if(  (1024%(int)header.fNcol !=0)) {
      header.fNrec = header.fNrow+1;
      headrec++;
    } else 
      header.fNrec=header.fNrow; 
 
    header.fLabbyt = header.fNcol*header.fLabrec*4;
    header.fLenbyt = (float) header.fNcol * 4;
  }
  
 // file type
 
    switch (im) {
      case IMG_XMIPP:
    	header.fIform = 1;    // for a 2D Xmipp image. 
        break;
      case IMG_INT:
    	header.fIform = 9;    // for a 2D int Xmipp image. 
        break;
      case VOL_XMIPP:
    	header.fIform = 3;    // for a 3D Xmipp volume. 
        break;
      case VOL_INT:
    	header.fIform = 10;    // for a 3D int Xmipp volume. 
        break;
      case IMG_FOURIER:
    	header.fIform = -1;    // for a 2D Fourier transform. 
        break;
      case VOL_FOURIER:
    	header.fIform = -3;    // for a 3D Fourier transform. 
        break;
    }

    
 // Set scale to 1 (never used by XMipp)
 
    header.fScale = 1;

 // Set Angle of rotation used by Xmipp = 0 (won't be used again) 
 // header.fAngle1 = 0;
 // CO: it's used sometimes to interact with old Xmipp programs

 // Set Geometrical transformation Matrix to Identity (won't be used again)
 
   for (unsigned i=0; i<3; i++)
      for (unsigned j=0; j<3; j++)
          if (i==j) header.fGeo_matrix[i][j]=(double)1.0;
          else      header.fGeo_matrix[i][j]=(double)0.0;

 // Sets statistical information flags to be compatible with Spider
    
    header.fSig = -1;
    header.fImami = 0;
}


// Interaction with data
void headerXmipp::set_dimension(float Ydim, float Xdim){
  header.fNrow = Ydim;
  header.fNcol = Xdim;
}

void headerXmipp::get_dimension(float &Ydim, float &Xdim) const{
  Ydim = header.fNrow;
  Xdim = header.fNcol;
}

char* headerXmipp::get_date() const{
  return (char*) header.szIDat;
}

void headerXmipp::set_date() {
  time_t lTime;
  struct tm *tmTime;
  
  time(&lTime);                // Get seconds since January 1 1970 
  tmTime = localtime(&lTime);  // Get current time
  tmTime->tm_mon++;            // Months are from 0..11 so put ti from 1..12
  sprintf(header.szIDat,"%d%s%d%s%d",tmTime->tm_mday,"-",tmTime->tm_mon,"-",tmTime->tm_year);
}

char* headerXmipp::get_time() const{
  return (char*) header.szITim;
}

void headerXmipp::set_time() {
  time_t lTime;
  struct tm *tmTime;

  time(&lTime);                // Get seconds since January 1 1970 
  tmTime = localtime(&lTime);  // Get current time
  sprintf(header.szITim,"%d%s%d",tmTime->tm_hour,":",tmTime->tm_min);
}

char* headerXmipp::get_title() const{
  return (char*) header.szITit;
}

void headerXmipp::set_title(FileName newName) {
  strcpy(header.szITit, newName.c_str()); // Set title of image in the header
}

void headerXmipp::set_eulerAngles(float Phi, float Theta, float Psi){
  header.fIangle = 1; // sets flag
  header.fPhi = Phi;
  header.fTheta = Theta;
  header.fPsi = Psi;
}

void headerXmipp::set_eulerAngles1(float Phi1, float Theta1, float Psi1){
  if(header.fFlag==2.f)  ;
  else
    header.fFlag = 1.f; // sets flag
  header.fPhi1 = Phi1;
  header.fTheta1 = Theta1;
  header.fPsi1 = Psi1;
}

void headerXmipp::set_eulerAngles2(float Phi2, float Theta2, float Psi2){
  header.fFlag = 2; // sets flag
  header.fPhi2 = Phi2;
  header.fTheta2 = Theta2;
  header.fPsi2 = Psi2;
}

template <class T>
   void headerXmipp::get_eulerAngles(T &Phi, T &Theta, T &Psi) const{
     Phi = (T) header.fPhi;
     Theta = (T) header.fTheta;
     Psi = (T) header.fPsi;
   }

template <class T>
   void headerXmipp::get_eulerAngles1(T &Phi1, T &Theta1, T &Psi1) const{
     Phi1 = (T) header.fPhi1;
     Theta1 = (T) header.fTheta1;
     Psi1 = (T) header.fPsi1;
   }

template <class T>
   void headerXmipp::get_eulerAngles2(T &Phi2, T &Theta2, T &Psi2) const{
     Phi2 = (T) header.fPhi2;
     Theta2 = (T) header.fTheta2;
     Psi2 = (T) header.fPsi2;
   }

matrix2D<double> headerXmipp::fGeo_matrix() {
   matrix2D<double> retval(3,3);
   retval(0,0)=header.fGeo_matrix[0][0];
   retval(0,1)=header.fGeo_matrix[0][1];
   retval(0,2)=header.fGeo_matrix[0][2];
   retval(1,0)=header.fGeo_matrix[1][0];
   retval(1,1)=header.fGeo_matrix[1][1];
   retval(1,2)=header.fGeo_matrix[1][2];
   retval(2,0)=header.fGeo_matrix[2][0];
   retval(2,1)=header.fGeo_matrix[2][1];
   retval(2,2)=header.fGeo_matrix[2][2];
   return retval;
}

void instantiate_header() {
   headerXmipp header;
   float  f;
      header.get_eulerAngles(f,f,f);
      header.get_eulerAngles1(f,f,f);
      header.get_eulerAngles2(f,f,f);
   double d;
      header.get_eulerAngles(d,d,d);
      header.get_eulerAngles1(d,d,d);
      header.get_eulerAngles2(d,d,d);
}
