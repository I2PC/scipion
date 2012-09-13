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

#include <data/normalize.h>

void process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Normalize_parameters * eprm = (Normalize_parameters *) prm;
    if (eprm->apply_geo) eprm->apply_geo_mask(img);
    eprm->apply(&img);
}

void process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    vol().statistics_adjust(0, 1);
}

int main(int argc, char **argv)
{
    Normalize_parameters prm;
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Colimate menu =========================================================== */
/*Colimate:
   PROGRAM Normalize {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Normalize/Help/normalize.html";
      help="Normalize particle projections";
      OPEN MENU Normalize;
      COMMAND LINES {
         + usual: xmipp_normalize
                   #include "prog_line.mnu"
                   -method $METHOD $METHOD_LIST
                  [-background $SHAPE $SHAPE_LIST $R]
      }
      PARAMETER DEFINITIONS {
         #include "prog_vars.mnu"
         $METHOD {shown=no; type=text;}
         $METHOD_LIST {
            label="Normalizing method";
            type=list {
               "OldXmipp"      {$METHOD="OldXmipp";      OPT($R)=0;}
               "Near_OldXmipp" {$METHOD="Near_OldXmipp"; OPT($R)=1;}
               "NewXmipp"      {$METHOD="NewXmipp";      OPT($R)=1;}
               "Michael"       {$METHOD="Michael";       OPT($R)=1;}
            };
         }
         OPT($R) {label="Background definition";}
         $SHAPE {shown=no; type=text;}
         $SHAPE_LIST {
            label="Background shape";
            type=list {
               "circle" {$SHAPE="circle";}
               "frame"  {$SHAPE="frame";}
            };
         }
         $R     {type=natural; label="Background radius or width";}
      }
   }
   MENU Normalize {
      #include "prog_menu.mnu"
      "Normalizing parameters"
      $METHOD_LIST
      OPT($R)
   }
*/
