/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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

/* Includes ---------------------------------------------------------------- */
#include <XmippData/xmippArgs.hh>
#include "qapplication.h"
#include "QtMainWidgetMark.hh"
#include "QtWidgetPSD.hh"

/* Prototypes -------------------------------------------------------------- */
void Usage();

/* Main -------------------------------------------------------------------- */
int main( int argc, char **argv ) {
   FileName fnRaw;
   FileName fnRawTilted;
   bool     reversed;
   FileName fn_assign_CTF;
   bool     ctf_mode=false;

   // Get input parameters .................................................
   try {
       fnRaw         = get_param( argc, argv, "-i" );
       fnRawTilted   = get_param( argc, argv, "-tilted", "" );
       reversed      = check_param( argc, argv, "-reverse_endian");
       fn_assign_CTF = get_param( argc, argv, "-psd", "");
       if (check_param(argc, argv, "-ctf")) {
          ctf_mode=true;
          fn_assign_CTF = get_param( argc, argv, "-ctf");
       }
   } catch ( Xmipp_error XE ) { cout << XE; Usage(); exit( 1 ); }

   try {
      Micrograph m, mTilted;
      
      m.open_micrograph( fnRaw, reversed );
      m.compute_8_bit_scaling();
      if ( fnRawTilted != "" ) {
         mTilted.open_micrograph( fnRawTilted, reversed );
	 mTilted.compute_8_bit_scaling();
      }
      
      // Configure application .............................................
      QApplication app( argc, argv );
      QtMainWidgetMark *mainWidget;
      
      if ( fnRawTilted == "" ) mainWidget = new QtMainWidgetMark( &m );
      else mainWidget = new QtMainWidgetMark( &m, &mTilted );
      
      // Check if the PSDs must be shown ...................................
      if (fn_assign_CTF!="") {
         QtWidgetPSD PSDshow;
         if (ctf_mode) PSDshow.set_CTF_mode();
         PSDshow.set_assign_CTF_file(m,fn_assign_CTF);
         PSDshow.show();
      }

      // Run application ...................................................
      app.setMainWidget( mainWidget );
      mainWidget->show();
      
      app.exec();

      // Finish ............................................................
      m.close_micrograph();
      if ( fnRawTilted != "" ) mTilted.close_micrograph();
      delete mainWidget;
   } catch ( Xmipp_error XE ) { cout << XE; }
   return 0;
}

/* Usage ------------------------------------------------------------------- */
void Usage() {
   cerr << "Purpose: Mark particles in a Raw image\n"
        << "         There must exist the image and the corresponding .inf file\n"
        << "\n"
        << "Usage: mark [options]\n"
        << "   -i <input raw file>                : File with the image\n"
        << "  [-tilted <tilted raw file>]         : Image with the tilted pair\n"
	<< "  [-reverse_endian]                   : Raw 16-bit file with reversed endian\n"
        << "  [-psd <assign_CTF_prm_file>]        : Show the PSDs\n"
        << "  [-ctf <assign_CTF_prm_file>]        : Show the CTF models\n"
   ;
}

/* Colimate menu =========================================================== */
/*Colimate:
   PROGRAM Mark {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Mark/Help/mark.html";
      help="Mark particles in raw micrographs";
      OPEN MENU Mark;
      COMMAND LINES {
         + usual: xmipp_mark -i $UNTILTED
                  [-tilted $TILTED]
                  [-reverse_endian]
      }
      PARAMETER DEFINITIONS {
         $UNTILTED {
            label="Untilted or single micrograph";
            help="Raw file";
            type=file existing;
         }
         $TILTED {
            label="Tilted micrograph";
            help="Raw file";
            type=file existing;
         }
         OPT(-reverse_endian) {
            label="Reverse endianness";
            help="The endiannes is reversed in both micrographs";
         }
      }
   }
   MENU Mark {
      "I/O parameters"
      $UNTILTED
      OPT($TILTED)
      OPT(-reverse_endian)
   }
*/
