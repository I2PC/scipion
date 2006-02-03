/***************************************************************************
 *
 * Authors:     Debora gil 
 *              Roberto Marabini
 *
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
 *                                      <                               
 * You should have received a copy of the GNU General Public License   
 * along with this program; if not, write to the Free Software         
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA            
 * 02111-1307  USA                                                     
 *                                                                     
 *  All comments concerning this program package may be sent to the    
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

#include "../xmippCCP4.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippImages.hh>
#include <fstream>
#include <iomanip>
#define GCC_VERSION (__GNUC__ * 10000 \
   + __GNUC_MINOR__ * 100 \
   + __GNUC_PATCHLEVEL__)
/* Test for GCC > 3.3.0 */
#if GCC_VERSION >= 30300
   #include <sstream>
#else
   #include <strstream.h>
#endif

#define VERBOSE
//#define DEBUG

/* ------------------------------------------------------------------------- */
void CCP4::write(const FileName &fn_out, const ImageXmipp &I, bool reversed) {
   FILE *fp;

   //fill mrc header and reverse if needed
   fill_header_from_xmippimage(I, reversed);
   
   //open file
   if ((fp = fopen(fn_out.c_str(), "wb")) == NULL)
     REPORT_ERROR(1503,"CCP4::write: File " + fn_out + " cannot be saved");

   //write header. note that FWRITE can not be used because 
   //floats and longs are invoved
   if(fwrite(&my_mrc_header, sizeof(char), SIZEOF_MRC_HEADER, fp) !=
                            SIZEOF_MRC_HEADER)
     REPORT_ERROR(1503,"CCP4::write: Header of file " + fn_out + " cannot be saved");

   //data, 
   float f; //only float are suported
   FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(I()){
      f=(float) MULTIDIM_ELEM(I(),i);
      FWRITE (&f, sizeof(float), 1, fp, reversed);
      }
   
   fclose(fp);
}

/* ------------------------------------------------------------------------- */
void CCP4::read(const FileName &fn_in, 
                      ImageXmipp &I, bool reversed) {

   FILE *fp;
   read_header_from_file(fn_in, reversed);
   //
   //open file
   if ((fp = fopen(fn_in.c_str(), "rb")) == NULL)
     REPORT_ERROR(1503,"CCP4::read: File " + fn_in + " cannot be read");

   I().resize(my_mrc_header.ny, my_mrc_header.nx);
  //MRC data is made by
   int mode_size;
   switch ( my_mrc_header.mode ) {
   case MODE_BYTE:
     mode_size = sizeof(unsigned char);
     break;
   case MODE_SHORT:
     mode_size = sizeof(short int);
     break;
   case MODE_FLOAT:
     mode_size = sizeof(float);
     break;
   default:
     REPORT_ERROR(1503,"CCP4::read: I do not know how to read this mrc file \
     format, try the reverse_endian flag");
     break;
   }
      

   // Get header size
   struct stat info;
   if (fstat(fileno(fp), &info)) 
      EXIT_ERROR(1,(string)"CCP4: Cannot get size of "+fn_in);
   int header_size=info.st_size-my_mrc_header.nx*
                                my_mrc_header.ny*
				my_mrc_header.nz*mode_size;

   // Skip header
   fseek(fp,header_size,SEEK_SET);

   switch ( my_mrc_header.mode ) {
   case MODE_BYTE:
      unsigned char c;
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(I()) {
	 fread(&c, sizeof(unsigned char), 1, fp);
	 MULTIDIM_ELEM(I(),i)=(double)c;
      }
      break;
   case MODE_SHORT:
      short int si;
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(I()) {
	 FREAD(&si, sizeof(short int), 1, fp,reversed);
	 MULTIDIM_ELEM(I(),i)=(double)si;
      }
      break;
   case MODE_FLOAT:
      float f;
      FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(I()) {
	 FREAD(&f, sizeof(float), 1, fp, reversed);
	 MULTIDIM_ELEM(I(),i)=(double)f;
      }
     break;
   default:
     REPORT_ERROR(1503,"CCP4::read: I do not know how to read this mrc file format");
     break;
   }
   fclose(fp);
}

/* ------------------------------------------------------------------------- */
void CCP4::clear(){
   memset(&my_mrc_header, '\0', SIZEOF_MRC_HEADER);
} /*clear*/


/* ------------------------------------------------------------------------- */
/** Fill mrc header from xmipp image. */
   void CCP4::fill_header_from_xmippimage(ImageXmipp I, bool reversed){
    clear();
   if(reversed==false){
    my_mrc_header.nx	= my_mrc_header.mx = I().ColNo();
    my_mrc_header.ny	= my_mrc_header.my = I().RowNo();
    my_mrc_header.nz	= my_mrc_header.mz = 1;
    my_mrc_header.mode  = MODE_FLOAT;
    my_mrc_header.mapc  = X_AXIS; 
    my_mrc_header.mapr  = Y_AXIS; 
    my_mrc_header.maps  = Z_AXIS; 
    my_mrc_header.amin  = (float)(I().compute_min()); 
    my_mrc_header.amax  = (float)(I().compute_max()); 
    my_mrc_header.amean = (float)(I().compute_avg()); 
    }
   else{
    my_mrc_header.nx    = I().ColNo();
    ByteSwap5(my_mrc_header.nx);
    my_mrc_header.mx	= my_mrc_header.nx;

    my_mrc_header.ny    = I().RowNo();
    ByteSwap5(my_mrc_header.ny);
    my_mrc_header.my	= my_mrc_header.ny;
    
    my_mrc_header.nz	= 1;
    ByteSwap5(my_mrc_header.nz);
    my_mrc_header.mz	= my_mrc_header.nz;
   
    my_mrc_header.mode  = MODE_FLOAT;
    ByteSwap5((my_mrc_header.mode));
    
    my_mrc_header.mapc  = X_AXIS;
    ByteSwap5(my_mrc_header.mapc); 
    my_mrc_header.mapr  = Y_AXIS;
    ByteSwap5(my_mrc_header.mapr); 
    my_mrc_header.maps  = Z_AXIS;
    ByteSwap5(my_mrc_header.maps); 
    
    my_mrc_header.amin = I().compute_min();
    ByteSwap5(my_mrc_header.amin); 
    my_mrc_header.amax = I().compute_max();
    ByteSwap5(my_mrc_header.amax); 
    my_mrc_header.amean = I().compute_avg();
    ByteSwap5(my_mrc_header.amean); 

    } 
}
/* ------------------------------------------------------------------------- */
/** Fill mrc header from mrc file. */
   void CCP4::read_header_from_file(const FileName &fn_in, bool reversed){
    clear();
    FILE *fp;
    //open file
    if ((fp = fopen(fn_in.c_str(), "rb")) == NULL)
      REPORT_ERROR(1503,"CCP4::read_header_from_file: File " + fn_in + " cannot be saved");

   // Get  size
   FREAD(&(my_mrc_header.nx), sizeof(int), 1, fp, reversed);
   FREAD(&(my_mrc_header.ny), sizeof(int), 1, fp, reversed);
   FREAD(&(my_mrc_header.nz), sizeof(int), 1, fp, reversed);
   FREAD(&(my_mrc_header.mode), sizeof(int), 1, fp, reversed); 

   fclose(fp);   
}

