/***************************************************************************
 *
 * Authors:     Irene Martinez
 *              Roberto Marabini
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


  /* **********************************************************************

        This file contains some routines for reading and writing
        geometric information in the image's header, using the
        struct "GEO_INFO" defined in "spider.h" as:

		typedef struct {
			float fImami;
			float fFmax;
			float fFmin;
			float fAv;
			float fSig;
			double fGeo_matrix[3][3];
			float fAngle1;
			} GEO_INFO;

   ************************************************************************/

#include <cstdio>
#include <cmath>

#include <sys/types.h>
#include <sys/stat.h>

#include "spider.h"
#include "groe.h"

/***************************************************************
  This routine read or write the struct "infogeo" using a
  3x3 matrix and a pointer to an angle value, depending on
  the flag "rdwr".
 ***************************************************************/
void rdwr_geo(double matriz[3][3],float *angle,GEO_INFO *infogeo,int rdwr)
{
 int i;

 switch(rdwr)
  {
    case WRITING :
      for (i = 0; i < 3; i++)
       {
          infogeo->fGeo_matrix[i][0] = matriz[i][0];
          infogeo->fGeo_matrix[i][1] = matriz[i][1];
          infogeo->fGeo_matrix[i][2] = matriz[i][2];
       }
      infogeo->fAngle1=*angle;
      break;

    case READING :
      for (i = 0; i < 3; i++)
       {
          matriz[i][0] = infogeo->fGeo_matrix[i][0];
          matriz[i][1] = infogeo->fGeo_matrix[i][1];
          matriz[i][2] = infogeo->fGeo_matrix[i][2];
       }
      *angle=infogeo->fAngle1;
      break;
  }

}

/***************************************************************
  This routine read or write the header "cabecero" using the
  struct "infogeo", depending on the flag "rdwr".
 ***************************************************************/
void infog_cabecero(GEO_INFO *infogeo,int rdwr)
{
 int i;

  switch(rdwr)
  {
    case WRITING :
      for (i = 0; i < 3; i++)
       {
          cabecero.fGeo_matrix[i][0] = infogeo->fGeo_matrix[i][0];
          cabecero.fGeo_matrix[i][1] = infogeo->fGeo_matrix[i][1];
          cabecero.fGeo_matrix[i][2] = infogeo->fGeo_matrix[i][2];
       }
      cabecero.fAngle1= infogeo->fAngle1;

      cabecero.fAv   = infogeo->fAv;
      cabecero.fSig  = infogeo->fSig;
      cabecero.fImami= infogeo->fImami;
      cabecero.fFmax = infogeo->fFmax;
      cabecero.fFmin = infogeo->fFmin;
      break;

    case READING :
      for (i = 0; i < 3; i++)
       {
          infogeo->fGeo_matrix[i][0] = cabecero.fGeo_matrix[i][0];
          infogeo->fGeo_matrix[i][1] = cabecero.fGeo_matrix[i][1];
          infogeo->fGeo_matrix[i][2] = cabecero.fGeo_matrix[i][2];
       }
      infogeo->fAngle1=cabecero.fAngle1;

      infogeo->fAv   = cabecero.fAv;
      infogeo->fSig  = cabecero.fSig;
      infogeo->fImami= cabecero.fImami;
      infogeo->fFmax = cabecero.fFmax;
      infogeo->fFmin = cabecero.fFmin;
      break;
  }
}

/**********************************************************/
/**********************************************************/
