/***************************************************************************
 *
 * Authors: Pedro Antonio de Alarcón (pedro@cnb.uam.es)
 *          Carlos Oscar S. Sorzano
 *          Alberto Pascual Montano
 *          Roberto Marabini
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

#include "../xmippVolumes.hh"

/*
// Specialization for complex numbers --------------------------------------
template <>
void VolumeT<complex<double> >::read(FILE *fh,
  int Zdim, int Ydim, int Xdim, bool reversed, Volume_Type volume_type)
{
  img.resize(Zdim, Ydim,Xdim);
  FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(img)
  {
            float a,b;
			// read real part of a complex number
            FREAD (&a, sizeof(float), 1, fh, reversed);
			// read imaginary part of a complex number
			FREAD (&b, sizeof(float), 1, fh, reversed);			
			// Assign the number
			complex<double> c(a,b);
            MULTIDIM_ELEM(img,i)=c;
  }
}
*/




/* Is Xmipp image? --------------------------------------------------------- */
int Is_VolumeXmipp(const FileName &fn, bool skip_type_check, bool force_reversed) {
   FILE *fp; 
   int result=0;
   headerXmipp header(headerXmipp::VOL_XMIPP);
   headerXmipp header1(headerXmipp::VOL_INT);

   // Open file
   if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501,"Is_VolumeXmipp: File " + fn + " not found");

   // Read header
   result = header.read(fp, skip_type_check, force_reversed);
   if(result==1) {fclose(fp); return (headerXmipp::VOL_XMIPP);}


   result = header1.read(fp, skip_type_check, force_reversed);
   fclose(fp);
   return (headerXmipp::VOL_INT * result);
}

/* Is Xmipp image? --------------------------------------------------------- */
int Is_FourierVolumeXmipp(const FileName &fn, bool skip_type_check, bool force_reversed) {
   FILE *fp; 
   int result=0;
   headerXmipp header(headerXmipp::VOL_FOURIER);

   // Open file
   if ((fp = fopen(fn.c_str(), "rb")) == NULL)
        REPORT_ERROR(1501,"Is_FourierVolumeXmipp: File " + fn + " not found");

   // Read header
   result = header.read(fp, skip_type_check, force_reversed);
   fclose(fp);
   return result;
}

/* Get Volume size ---------------------------------------------------------- */
void GetXmippVolumeSize(const FileName &fn, int &Zdim, int &Ydim, int &Xdim) {
   FILE *fp; 
   int result;
   headerXmipp header(headerXmipp::VOL_XMIPP);

   // Open file
   if ((fp = fopen(fn.c_str(), "rb")) == NULL)
     REPORT_ERROR(1501,"Is_ImageXmipp: File " + fn + " not found");

   // Read header
   result = header.read(fp, false, false);

   fclose(fp);
   if (result) {
      Zdim=header.iZdim();
      Ydim=header.iYdim();
      Xdim=header.iXdim();
   } else {
      Zdim=Ydim=Xdim=-1;
   }
}
