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

#include <XmippData/xmippProgs.hh>
#include <XmippData/xmippArgs.hh>

class Window_parameters: public Prog_parameters {
public:
   bool size_mode;
   bool physical_coords;
   int sizeX, sizeY, sizeZ;
   int x0, y0, z0;
   int xF, yF, zF;
   double init_value;
   int wrong_parameters;
   int average_pad;
   int corner_pad;

   
   void read(int argc, char **argv) _THROW {
      init_value=0.;
      wrong_parameters=0;
      average_pad=FALSE;
      corner_pad=FALSE;
      
      Prog_parameters::read(argc,argv);
      
      if (check_param(argc,argv,"-pad_value"))
         {
	 init_value = AtoF(get_param(argc, argv, "-pad_value"));
	 wrong_parameters=1;
	 }
      if (check_param(argc,argv,"-corner_pad_value"))
         {
         corner_pad=TRUE;
	 if (wrong_parameters>0)  wrong_parameters=-1;
	 else  wrong_parameters=1;
	 }
      if (check_param(argc,argv,"-average_pad_value"))
         {
         average_pad=TRUE;
	 if (wrong_parameters>0)  wrong_parameters=-1;
	 else  wrong_parameters=1;
	 }
      if(wrong_parameters==-1)	 
         REPORT_ERROR(1,"incompatible options");
	 
      if (check_param(argc,argv,"-size")) {
         size_mode=TRUE;
         int i=position_param(argc,argv,"-size");
         if (i+2>=argc) {sizeZ=sizeY=sizeX=AtoI(argv[i+1]);}
         else if (i+3>=argc) {sizeZ=sizeY=AtoI(argv[i+2]); sizeX=AtoI(argv[i+1]);}
         else if (i+4>=argc) {sizeZ=AtoI(argv[i+3]); sizeY=AtoI(argv[i+2]);
	                sizeX=AtoI(argv[i+1]);}
         else REPORT_ERROR(1,"Not enough parameters after -size");

         x0=FIRST_XMIPP_INDEX(sizeX);
         y0=FIRST_XMIPP_INDEX(sizeY);
         z0=FIRST_XMIPP_INDEX(sizeZ);
         xF=LAST_XMIPP_INDEX(sizeX);
         yF=LAST_XMIPP_INDEX(sizeY);
         zF=LAST_XMIPP_INDEX(sizeZ);
         physical_coords=FALSE;
      } else if (check_param(argc,argv,"-r0")) {
         size_mode=FALSE;
         // Get r0
         int i=position_param(argc,argv,"-r0");
         if (i+2>=argc) REPORT_ERROR(1,"Not enough parameters after -r0");
         else {x0=AtoI(argv[i+1]); y0=AtoI(argv[i+2]);}
         if (i+3<argc)
            try {
               z0=AtoI(argv[i+3]);
            } catch (Xmipp_error XE) {z0=1;}

         // Get rF
         i=position_param(argc,argv,"-rF");
         if (i==-1) REPORT_ERROR(1,"-rF not present");
         if (i+2>=argc) REPORT_ERROR(1,"Not enough parameters after -rF");
         else {xF=AtoI(argv[i+1]); yF=AtoI(argv[i+2]);}
         if (i+3<argc)
            try {
               zF=AtoI(argv[i+3]);
            } catch (Xmipp_error XE) {zF=1;}

         physical_coords=check_param(argc,argv,"-physical");
      } else
         REPORT_ERROR(1,"Unknown windowing type");
   }
   
   void show() {
      Prog_parameters::show();
      if (size_mode)
         cout << "New size: (XxYxZ)=" << sizeX << "x" << sizeY << "x"
              << sizeZ << endl;
      else
         cout << "New window: from (z0,y0,x0)=("<< z0 << ","
              << y0 << "," << x0 << ") to (zF,yF,xF)=(" << zF << "," << yF
              << "," << xF << ")\n"
              << "Physical: " << physical_coords << endl;
   }

   void usage() {
      Prog_parameters::usage();
      cerr << " -physical                  : Use physical instead of logical\n"
           << "                              coordinates\n"
	   << " -pad_value                 : value used for padding\n"
	   << " -corner_pad_value          : use the value of the upper\n"
	   << "                              left corner for averaging"
	   << "                              for padding\n"
	   << " -average_pad_value         : use the image average for padding\n"
           << "  [-r0 <x0> <y0> [<z0>]     : Window using window corners\n"
           << "   -rF <xF> <yF> [<zF>]]    : by default indexes are logical\n"
           << "  [-size <sizeX> [<sizeY>] [<sizeZ>]: Window to a new size\n"
           << "                            : if only one is given, the other two\n"
           << "                              are supposed to be the same\n"
	   ;
   }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm) {
   Window_parameters *eprm=(Window_parameters *) prm;
   if (eprm->average_pad==TRUE)
       eprm->init_value=(img()).compute_avg();
   else if (eprm->corner_pad==TRUE)
       eprm->init_value=DIRECT_MAT_ELEM(img(),0,0);
   if (!eprm->physical_coords)
      img().window(eprm->y0,eprm->x0,eprm->yF,eprm->xF,eprm->init_value);
   else img().window(STARTINGY(img())+eprm->y0,STARTINGX(img())+eprm->x0,
           STARTINGY(img())+eprm->yF,STARTINGX(img())+eprm->xF,eprm->init_value);
   return TRUE;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm) {
   Window_parameters *eprm=(Window_parameters *) prm;
   if (eprm->average_pad==TRUE)
       eprm->init_value=(vol()).compute_avg();
   else if (eprm->corner_pad==TRUE)
       eprm->init_value=DIRECT_VOL_ELEM(vol(),0,0,0);
   if (!eprm->physical_coords)
      vol().window(eprm->z0,eprm->y0,eprm->x0,eprm->zF,eprm->yF,eprm->xF,
      eprm->init_value);
   else vol().window(STARTINGZ(vol())+eprm->z0,STARTINGY(vol())+eprm->y0,
      STARTINGX(vol())+eprm->x0,STARTINGZ(vol())+eprm->zF,
      STARTINGY(vol())+eprm->yF,STARTINGX(vol())+eprm->xF,eprm->init_value);
   return TRUE;
}

int main (int argc, char **argv) {
   Window_parameters prm;
   SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Window {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Window/Help/window.html";
      help="Window volumes and images, extract a part of them";
      OPEN MENU menu_window;
      COMMAND LINES {
	+ usual: xmipp_window
               #include "prog_line.mnu"
               $WINDOWING_METHOD
               [-r0 $X0 $Y0 [$Z0] -rf $XF $YF [$ZF]]
               [-size $XDIM [$YDIM] [$ZDIM]]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        $WINDOWING_METHOD {
           label="Window action";
           type=list {
              "Logical corners" {OPT(-r0)=1; OPT(-size)=0;}
              "Centered cube"   {OPT(-r0)=0; OPT(-size)=1;}
           };
        }
        OPT(-r0) {label="Window using logical corners";}
           $X0 {type=float; label="Initial X";}
           $Y0 {type=float; label="Initial Y";}
           $Z0 {type=float; label="Initial Z";}
           $XF {type=float; label="Final X";}
           $YF {type=float; label="Final Y";}
           $ZF {type=float; label="Final Z";}
        OPT(-size) {label="Window within a centered cube";}
           $ZDIM  {type=float; label="Cube Zdim"; by default=$XDIM;}
           $YDIM  {type=float; label="Cube Ydim"; by default=$XDIM;}
           $XDIM  {type=float; label="Cube Xdim";}
      }
   }

   MENU menu_window {
      #include "prog_menu.mnu"
      "Window parameters"
      $WINDOWING_METHOD
      OPT(-r0)
      OPT(-size)
   }
*/
