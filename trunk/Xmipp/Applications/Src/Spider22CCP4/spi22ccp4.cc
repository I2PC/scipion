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

#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippVolumes.hh>
#include <XmippData/xmippImages.hh>
#include <sys/stat.h>
#include <XmippInterface/xmippCCP4.hh>


/* Prototypes -============================================================= */
void Usage ();

int main (int argc, char *argv[]) {
/* Input Parameters ======================================================== */
FileName       fn_in;    // input file
FileName       fn_out;   // output file
bool           reverse_endian;
CCP4 mrcimage;

/* Parameters ============================================================== */
   try {
       fn_in  = get_param(argc, argv, "-i");
       fn_out = get_param(argc, argv, "-o");
       reverse_endian=check_param(argc,argv,"-reverse_endian");
   }
   catch (Xmipp_error XE) {cout << XE; Usage();}
  
   ImageXmipp  I;
   try {
      if(Is_ImageXmipp(fn_in)){
	 I.read(fn_in);
	 mrcimage.write(fn_out,I,reverse_endian);
	 }
      else{
	 mrcimage.read(fn_in,I,reverse_endian);
	 I.write(fn_out);
	 }
   } catch (Xmipp_error XE) {cout << XE;}
   
   exit(0);
} //main

/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage () {
  printf (
     "Usage: ccp42spi [Purpose and Parameters]"
     "\nPurpose: Convert from a 3D CCP4 images to Xmipp ones"
     "\n    -i    file_in        input CCP4 file"
     "\n    -o    file_out       output Xmipp file"
     "\n   [-reverse_endian]     by default, output has the same endiannes as input"
     "\n                         use this option to change endianness\n");
}
