/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

void usage() {
   cerr << "Spi2CCP4\n"
        << "   -i <Xmipp volume>           : Volume to convert\n"
        << "   -o <CCP4 volume>            : Converted volume\n"
        << "  [-sam=<Ts=1>]                : Sampling rate (angstrom/pixel)\n"
   ;
}

int main (int argc, char **argv) {
   FileName fn_in, fn_out;
   double   Ts;
   try {
      fn_in =get_param(argc,argv,"-i");
      fn_out=get_param(argc,argv,"-o");
   } catch (Xmipp_error XE) {cerr << XE; usage(); return 1;}
   
   try {
      VolumeXmipp V(fn_in);
      write_as_CCP4(&V,fn_out,Ts);
   } catch (Xmipp_error XE) {cerr << XE;}
   return 0;
}
