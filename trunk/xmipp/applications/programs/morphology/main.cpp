/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Pedro A. de Alarc�n (pedro@cnb.uam.es)
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

#include <data/progs.h>
#include <data/args.h>
#include <data/morphology.h>

class Morphology_parameters: public Prog_parameters {
public:
   #define DILATION  	  1
   #define EROSION   	  2
   #define OPENING   	  3
   #define CLOSING        4
   int operation;

   int size;
   int count;
   int neig;
public:
   void read(int argc, char **argv) {
      Prog_parameters::read(argc,argv);
      if (check_param(argc,argv,"-dil"))     	operation=DILATION;
      if (check_param(argc,argv,"-ero"))     	operation=EROSION;
      if (check_param(argc,argv,"-clo"))     	operation=CLOSING;
      if (check_param(argc,argv,"-ope"))     	operation=OPENING;

      size=AtoI(get_param(argc,argv,"-size","1"));
      neig=AtoI(get_param(argc,argv,"-neig","-1"));
      count=AtoI(get_param(argc,argv,"-count","0"));
   }

   void show() {
      Prog_parameters::show();
      cout << "Performing a ";
      switch (operation) {
         case DILATION	      : cout << "Dilation\n"; break;
	 case EROSION	      : cout << "Erosion\n"; break;
	 case OPENING	      : cout << "Opening\n"; break;
	 case CLOSING	      : cout << "Closing\n"; break;
      }
      cout << "Size=" << size << endl
           << "Neighbourhood=" << neig << endl
	   << "Count=" << count << endl;
   }

   void usage() {
      Prog_parameters::usage();
      cerr << "  [-dil]       	     : Apply dilation\n"
	   << "  [-ero]       	     : Apply erosion\n"
	   << "  [-clo]       	     : Apply closing\n"
	   << "  [-ope]       	     : Apply opening\n"
	   << "  [-neig <n=8 | 18>]  : Neighborhood considered \n"
	   << "                        (2D:4,8 3D:6,18,26)\n"
	   << "  [-size <s=1>]       : Size of the Strutural element\n"
	   << "  [-count <c=0>]      : Minimum required neighbors with \n"
	   << "                        distinct value\n"
      ;
   }
};


bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Morphology_parameters *eprm=(Morphology_parameters *) prm;
   if (eprm->neig==-1) eprm->neig=8;
   ImageXmipp retval; retval()=img();

   cout << "Initially the image has " << img().sum() << " pixels set to 1\n";
   switch (eprm->operation) {
      case DILATION:
            dilate2D(img(),retval(),eprm->neig,eprm->count,eprm->size);
         break;
      case EROSION:
            erode2D(img(),retval(),eprm->neig,eprm->count,eprm->size);
         break;
      case OPENING:
            opening2D(img(),retval(),eprm->neig,eprm->count,eprm->size);
         break;
      case CLOSING:
            closing2D(img(),retval(),eprm->neig,eprm->count,eprm->size);
         break;
   }

   img()=retval();
   cout << "Finally the image has " << img().sum() << " pixels set to 1\n";
   return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Morphology_parameters *eprm=(Morphology_parameters *) prm;
   if (eprm->neig==-1) eprm->neig=18;
   VolumeXmipp retval; retval()=vol();

   cout << "Initially the volume has " << vol().sum() << " voxels set to 1\n";
   switch (eprm->operation) {
      case DILATION:
            dilate3D(vol(),retval(),eprm->neig,eprm->count,eprm->size);
         break;
      case EROSION:
            erode3D(vol(),retval(),eprm->neig,eprm->count,eprm->size);
         break;
      case OPENING:
            opening3D(vol(),retval(),eprm->neig,eprm->count,eprm->size);
         break;
      case CLOSING:
            closing3D(vol(),retval(),eprm->neig,eprm->count,eprm->size);
         break;
   }

   vol()=retval();
   cout << "Finally the volume has " << vol().sum() << " voxels set to 1\n";
   return true;
}

int main (int argc, char **argv) {
   Morphology_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus =================================================================== */
/*Colimate:
   PROGRAM Morphology {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Morphology/Help/morphology.html";
      help="Apply morphology filters to binary images and volumes";
      OPEN MENU menu_morphology;
      COMMAND LINES {
	+ usual: xmipp_morphology
               #include "prog_line.mnu"
	         $OP [-neig $NEIG $NEIG_LIST] [-size $SIZE] [-count $COUNT]
      }
      PARAMETER DEFINITIONS {
         #include "prog_vars.mnu"
         $OP {
	    label="Filter type";
	    type=exclusion {
	       "Dilation" {-dil}
	       "Erosion"  {-ero}
	       "Opening"  {-ope}
	       "Closing"  {-clo}
	    };
	 }
	 $NEIG {shown=no; type=natural;}
	 $NEIG_LIST {
	    label="Neighbourhood";
	    type=list {
	        "4" {$NEIG=4;}
	        "8" {$NEIG=8;}
	        "6" {$NEIG=6;}
	       "18" {$NEIG=18;}
	       "26" {$NEIG=26;}
	    };
	
	 }
	 $SIZE {type=float; label="Structuring box size"; by default=1;}
	 $COUNT {type=natural;
	    label="Minimum no. of neighbours to be considered different";
	    by default=0;}
      }
   }

   MENU menu_morphology {
      #include "prog_menu.mnu"
      $OP
      OPT($NEIG)
      OPT($SIZE)
      OPT($COUNT)
   }
*/
