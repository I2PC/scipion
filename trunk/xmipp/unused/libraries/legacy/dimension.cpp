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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/


/* **********************************************************************

   This file contains the routines that get the image and volume
   dimensions from the header, and other functions for read the
   file names from the selecction files.

 ************************************************************************/

#include <cstring>
#include <cstdio>

#include "spider.h"
#include "groe.h"


/***************************************************************
 This function counts the number of files to process in the
 selecction file pointed by "fil_entrada". A file will be
 processed if its field value in the selecction file is > "flagVal".
 The first file name is copied in "name".
 Finally, it returns the number of counted files.
 ***************************************************************/
int Count_element(FILE *fil_entrada, char *name, int flagVal, int size)
{
    int i;
    int numI = 0;
    int valido = 0;
    char nom_img[16];
    char linea[80];

    while (fgets(linea, 80, fil_entrada) != NULL)
    {
        if ((i = sscanf(linea, "%s %d", nom_img, &valido)) > size || i == 0)
            continue;
        if (valido <= flagVal)
            continue;

        numI++;
        if (numI == 1) strcpy(name, nom_img);

    }

    fseek(fil_entrada, 0L, SEEK_SET);
    return(numI);
}


/***************************************************************
 This function looks in the selecction file for the first file
 name which has a field value > "flagVal".
 The file name is copied to "name" and returns '1'. If no file
 name is found, returns '0'.
 ***************************************************************/
int  First_name(FILE *fil_entrada, char *name, int flagVal, int size)
{
    int i;
    int numI = 0;
    int valido = 0;
    char nom_img[80];
    char linea[80];

    while ((fgets(linea, 80, fil_entrada) != NULL) && numI != 1)
    {
        if ((i = sscanf(linea, "%s %d", nom_img, &valido)) > size || i == 0)
            continue;
        if (valido <= flagVal)
            continue;

        numI++;
    }
    if (numI == 1) strcpy(name, nom_img);

    fseek(fil_entrada, 0L, SEEK_SET);
    return(numI);
}

/**************************************************************************/
/* ***********************   IMAGE DIMENSIONS   ***************************/

/***************************************************************
 This function check the image format and its dimensions calling
 the routine "SPIDheader()" (it's in the file "cabecero.c"), and
 prints an error message if some failure ocur in "SPIDheader()".
 ***************************************************************/
int Get_dimension(char *nom_entrada, int *fil, int *col)
{
    int i;

    ceroCabecero();
    if ((i = SPIDheader(nom_entrada, fil, col)) < 0)
    {
        fprintf(stdout, "\nCouldn't open file %s\n", nom_entrada);
        return(0);
    }
    if (i == 0)
    {
        fprintf(stdout, "\n Error: %s not a SPIDER file\n", nom_entrada);
        return(0);
    }

    return(1);
}

/**************************************************************************/
/* **********************   VOLUME DIMENSIONS   ***************************/

/***************************************************************
 This function check the volume format and its dimensions calling
 the routine "SPIDvolum()" (it's in the file "cabecero.c"), and
 prints an error message if some failure ocur in "SPIDvolum()".
 ***************************************************************/
int Get_Vol_dimension(char *nom_entrada, int *slice, int *row, int *col)
{
    int i;

    ceroCabecero();
    if ((i = SPIDvolum(nom_entrada, slice, row, col)) < 0)
    {
        fprintf(stdout, "\nCouldn't open file %s\n", nom_entrada);
        return(0);
    }
    if (i == 0)
    {
        fprintf(stdout, "\n Error: %s not a SPIDER volume \n", nom_entrada);
        return(0);
    }

    return(1);
}

/****************************************************************************/

