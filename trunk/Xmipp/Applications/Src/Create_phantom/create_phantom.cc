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

/* ------------------------------------------------------------------------- */
/* Includes                                                                  */
/* ------------------------------------------------------------------------- */
#include <Reconstruction/phantom.hh>

/* ------------------------------------------------------------------------- */
/* Prototypes                                                                */
/* ------------------------------------------------------------------------- */
void Usage(char *argv[]);

/* ------------------------------------------------------------------------- */
/* Program                                                                   */
/* ------------------------------------------------------------------------- */
int main (int argc, char *argv[]) {
   string            fn_phantom;
   string            fn_vol;
   Phantom           phantom;
   VolumeXmipp       vol;

   // Read Parameters ......................................................
   try {
      fn_phantom = get_param(argc, argv, "-i");
      fn_vol     = get_param(argc, argv, "-o");
   }
   catch (Xmipp_error XE) {cout << XE; Usage(argv);}


   try {
      // Read description file .............................................
      phantom.read(fn_phantom);

      // Generate volume and write .........................................
      phantom.draw_in(&vol);
      vol.write(fn_vol);
   }
   catch (Xmipp_error XE) {cout << XE; exit(1);}
   exit(0);
}

/* Usage ------------------------------------------------------------------- */
void Usage (char *argv[]) {
   printf(
      "\n\nPurpose: This program allows you to create phantom XMIPP volumes\n"
      "     from a phantom feature description file.\n"
      );
   printf(
      "\n\nUsage: %s [Parameters]"
      "\nOptions:"
      "\nParameter Values: (note space before value)"
      "\n    -i <description file>"
      "\n    -o <output file>"
      ,argv[0]);
   printf(
      "\nStructure of Parameter file:"
      "\n"
      "# Phantom description file, (generated with phantom help)\n"
      "# General Volume Parameters:\n"
      "#      Xdim      Ydim      Zdim   Background Density [Scale Factor]\n"
      "        65        65        65        0                   [1]\n"
      "# Feature Parameters:\n"
      "#Type +/= Density X_Center Y_Center Z_Center\n"
      " sph   +     1      <x0>     <y0>     <z0>    <radius>\n"
      " blo   +     1      <x0>     <y0>     <z0>    <radius>  <alpha>  <order> \n"
      " cyl   +     1      <x0>     <y0>     <z0>    <radius>            <height>               <rot> <tilt> <psi>\n"
      " dcy   +     1      <x0>     <y0>     <z0>    <radius>            <height>  <separation> <rot> <tilt> <psi>\n"
      " cub   =     1      <x0>     <y0>     <z0>    <xdim>     <ydim>    <zdim>                <rot> <tilt> <psi>\n"
      " ell   =     1      <x0>     <y0>     <z0>    <xradius> <yradius> <zradius>              <rot> <tilt> <psi>\n"
      " con   +     1      <x0>     <y0>     <z0>    <radius>            <height>               <rot> <tilt> <psi>\n"
   );
   printf(
      "\n"
      "Where:\n"
      "sph --> sphere, with center at (x0,y0,z0) and radius in pixels <radius>\n"
      "blo --> blob, with center at (x0,y0,z0), max_radius in pixels <radius>\n"
      "        alpha is the tampering parameter and the order should be self evident\n" 
      "cyl --> cylinder, initially with base at plane XY and height in Z (-h/2 to h/2)\n"
      "        and then it is moved to (x0,y0,z0) (its center) and rotated after\n"
      "        (tilt,rot,psi)\n"
      "dcy --> double cylinder, initially with base at plane XY and height in Z\n"
      "        (-h/2 to h/2), then it is moved to (x0,y0,z0) and rotated after\n"
      "        (tilt,rot,psi)\n"
      "cub --> cube of size (xdim,ydim,zdim) whose center is moved to (x0,y0,z0)\n"
      "        and rotated\n"
      "ell --> ellipsoid with 3 radius and then moved and rotated\n"
      "con --> cone, initially with base at plane XY and height in Z (-h/2 to h/2),\n"
      "        and then the center of the base is moved to (x0,y0,z0) and rotated\n"
      " +  --> the features can overlap and add their densities at overlapping positions\n"
      " =  --> the overlapping positions will keep the last feature density\n"
   );
   exit(1);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Create_phantom {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Create_phantom/Help/create_phantom.html";
      help="Produces a 3D volume from a phantom mathematical description";
      OPEN MENU menu_create_phantom;
      COMMAND LINES {
	+ usual: create_phantom -i $FILE_IN -o $FILE_OUT
      }
      PARAMETER DEFINITIONS {
        $FILE_IN {
	   label="Phantom Description file";
	   help="This file has got a complex structure, better see 
                 the Web help";
	   type=file existing;
	}
        $FILE_OUT {
	   label="Output Volume file";
	   help="Xmipp format";
	   type=file;
	}
      }
   }

   MENU menu_create_phantom {
      "Compulsory variables"
      $FILE_IN
      $FILE_OUT
   }
*/
