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
 ***************************************************************************
 * Modified by Alberto Pascual to add some parameters functionality
 * October 2001
 ****************************************************************************/


  /* **********************************************************************

        This file contains kind routines for the usage checking, and
        for file names handling.

   ************************************************************************/

#include <cstdio>
#include <cstring>
#include <ctime>
#include <cctype>

int check_comm( char **arv)
{

  if (strlen (arv[1]) != 2)
  {
    printf ("Error in the arguments ( %s )",arv[0]);
    puts ("series-name must have two characters");
    return 0;
  }

  if (strlen (arv[2]) != 1)
  {
    printf ("Error in the arguments ( %s )",arv[0]);
    puts ("t|u must have only one character");
    return 0;
  }

  if (tolower(*arv[2]) != 't' && tolower(*arv[2]) != 'u')
  {
    printf ("Error in the arguments ( %s )",arv[0]);
    puts (" t|u must be 't' or 'u'");
    return 0;
  }
 return 1;

}

/* ********************************************************************* */

void getTempoName(char *filNomb)
{
 time_t ltime;

   time(&ltime);

   sprintf(filNomb,"T%d",ltime);
   filNomb[8]='.';
   filNomb[12]='\0';
}
/* ********************************************************************* */
void getDate(char *day, char *hour)
{
 time_t lTiempo;
 struct tm *tmTiempoGreng;

   time(&lTiempo);
   tmTiempoGreng=localtime(&lTiempo);
   tmTiempoGreng->tm_mon++;

   sprintf(hour,"%d%s%d",tmTiempoGreng->tm_hour,":",tmTiempoGreng->tm_min);
   sprintf( day,"%d%s%d%s%d",tmTiempoGreng->tm_mday,"-",
            tmTiempoGreng->tm_mon,"-",tmTiempoGreng->tm_year);
}
/* ********************************************************************* */

/* ********************************************************************* */
/* Gets the complete file name (absolute, not relative to any directory) */

void getNameAbs(char *nameD, char *nameProg)
{
   char *puntero;
   int c= '/';

   puntero = strrchr( nameProg,c ) ;

   if( puntero ==NULL)
     strcpy ( nameD,nameProg);
   else
     strcpy ( nameD,puntero+1);
}
/* ********************************************************************* */
