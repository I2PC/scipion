/***************************************************************************
 *
 * Authors:    Sjors Scheres
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
 * MERCHANTABILITY or FITNESS FO A PARTICULAR PURPOSE.  See the
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

#include <data/args.h>
#include <data/image.h>
#include <data/selfile.h>

void Usage();

/* MAIN -------------------------------------------------------------------- */
int main (int argc,char *argv[]) {
   SelFile         SF;
   ImageXmipp      img;

   try 
   {
       SF.read(get_param(argc,argv,"-i"));
       SF.go_beginning();
       cerr <<" Resetting all angles, origin offsets, weights and mirror flags to zero ... "<<endl;
       while (!SF.eof()) 
       {
	   img.read(SF.NextImg());
	   img.clear_header();
	   img.write(img.name());
       }
       cerr <<" done!"<<endl;
   }
   } catch (Xmipp_error XE) {cout << XE; Usage();}
}

/* Usage ------------------------------------------------------------------- */
void Usage() {
    printf("Purpose:\n");
    printf(" Reset the geometric transformation (angles & shifts) in the header of 2D-images.\n");
    printf("Usage:\n");
    printf("   header_reset \n");
    printf("        -i  <selfile>      : selfile with input images \n");
    exit(1);
}
