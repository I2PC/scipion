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

        This file contains all the routines to read, the files with
        extensions ".sim" and ".clu".

   ************************************************************************/


#include <cstdio>

#ifdef __STDC__
  int num_componen_sim (FILE *);
  int num_componen_clu (FILE *);
  int lee_lista (FILE *, int , int *);
  int lee_sim(float *, int , char *);
  int lee_clas (float *, int, char *);
  int lee_cen (float *, int , char *);
  int posic_sim(FILE *,int);
  int find_tip(FILE *,char *);
  int posic_clu(FILE *,char *);
  static int old_sim_format(FILE *);
  static int old_clu_format(FILE *);
#else
  int num_componen_sim ();
  int num_componen_clu ();
  int lee_lista ();
  int lee_sim();
  int lee_clas ();
  int lee_cen ();
  int posic_sim();
  int find_tip();
  int posic_clu();
  static int old_sim_format();
  static int old_clu_format();
#endif

/* *************************************************************************/

int num_componen_sim (FILE *fich)
{
 int nsim = 0;
 char linea [2000];
 char aux_name[32];
 int found=0;

 aux_name[0]='\0';
 while ( !found && (fgets (linea, 200, fich) != NULL) )
 {
   if ( sscanf (linea,"%s %d",aux_name,&nsim) )
     if (!strncmp (aux_name,"SYM",3))
        found ++;

 }

 if(!found)  /** try to read the symmetries with the old file format **/
 {
   nsim = old_sim_format(fich);
   if(!nsim )
    {
       fprintf(stderr,"\n Sorry, couldn't read the number of symmetries.");
       fflush(stderr);
    }
 }
 fseek(fich,0L,SEEK_SET);

 return (nsim);
}

/* *************************************************************************/

int num_componen_clu (FILE *fich)
{
 int num_proto = 0;
 int i;
 char linea [200];
 char aux_name[64];
 int found=0;
 size_t len;

 aux_name[0]='\0';

 while ( !found && (fgets (linea, 200, fich) != NULL) )
  {
   if ( sscanf (linea,"%s",aux_name) )
     if (!strncmp (aux_name,"Proto",5))
        found ++;

  }
 if(found)
  {
    len= strlen (linea);
    for(i=0;i<len;i++)
      if(linea[i]==',') num_proto++;
  }
 else  /** try to read the prototypes with the old file format **/
  {
    num_proto = old_clu_format(fich);
    if(!num_proto )
    {
       fprintf(stderr,"\n Sorry, couldn't read the number of prototypes.");
       fflush(stderr);
    }
  }
 fseek(fich,0L,SEEK_SET);
 return (num_proto);
}
/* *************************************************************************/

int  lee_lista (FILE *fich, int num_sim, int *lista_sim)
{
 int i=0, j;
 char linea [2000];
 char formato[2000];
 char aux_name[2000];
 int found=0;
 aux_name[0]=formato[0]='\0';
 while ( !found && (fgets (linea, 2000, fich) != NULL) )
 {
   if ( sscanf (linea,"%s",aux_name) )
     if (!strncmp (aux_name,"Imag",4))
        found ++;
 }
 if(found)
 {
   for (i=0; i < num_sim ; i++)
   {   strcpy (formato , "%*s");
       for (j=0; j<i; j++)
         strcat (formato, " %*c %*d");
       strcat (formato, " %*c %d");
       if (sscanf (linea, formato, &lista_sim[i]) != 1)
       {   fprintf(stdout,"\n Error reading symmetry file ");
           return(0);
       }
   }
 }

 fseek(fich,0L,SEEK_SET);
 return(i);
}

/* *************************************************************************/

int  lee_sim(float *lista, int num_sim, char *linea)
{
 int i, j;
 char formato[2000];

 formato[0]='\0';
 for (i=0; i < num_sim; i++)
  {   strcpy (formato , "%*s");
      for (j=0; j<i; j++)
        strcat (formato, " %*c %*f");
      strcat (formato, " %*c %f");
      if (sscanf (linea, formato, &lista[i]) != 1)
      {   fprintf(stdout,"\n Error reading symmetry file ");
          return(0);
      }
  }
 return(1);
}


/* *************************************************************************/

int  lee_clas (float *lista, int num_clas, char *linea)
{
 int i, j;
 char formato[2000];

 formato[0]='\0';
 for (i=0; i < num_clas; i++)
  {   strcpy (formato , "%*s");
      for (j=0; j<i; j++)
        strcat (formato, " %*f");
      strcat (formato, " %f");
      if (sscanf (linea, formato, &lista[i]) != 1)
      {   fprintf(stdout,"\n : Error reading the classes in cluster file ");
          return(0);
      }
  }
 return(1);
}

/* *************************************************************************/

int  lee_cen (float *lista, int num_cen, char *linea)
{
 int i, j;
 char formato[2000];

 formato[0]='\0';
 for (i=0; i < num_cen; i++)
 {  strcpy (formato , "%*s %*s");
    for (j=0; j<i; j++)
        strcat (formato, " %*f");
    strcat (formato, " %f");
    if (sscanf (linea, formato, &lista[i]) != 1)
    {   fprintf(stdout,"\n Error reading the centres in cluster file ");
        return(0);
    }
    lista[i] *= 100;           /*** cambio de coor. para igualdad con .SIM ***/
  }
 return(1);
}


/* ************************************************************************* */

int posic_sim(FILE *fich,int Num)
{
 char linea [2000];
 char aux_name[2000];
 int found=0,i;

 aux_name[0]='\0';
 while ( !found && (fgets (linea, 2000, fich) != NULL) )
 {
   if ( ( sscanf (linea,"%s",aux_name) ) &&
        (!strncmp (aux_name,"Imag",4) ) )
        found ++;
 }
 if(found)
 {
   for (i=0; i<Num; i++)
     fgets (linea, 2000, fich);
 }
else
 {
    fprintf(stderr,"\n Error reading symmetry file \n");
    fflush(stderr);
 }
 return (found);
}

/* ************************************************************************* */

int find_tip(FILE *fich,char *linea)
{
 int found=0;
 char aux_name[64];

 aux_name[0]='\0';

  while ( !found && (fgets (linea, 2000, fich) != NULL) )
  {
   if ( sscanf (linea,"%s",aux_name) )
     if ((!strncmp (aux_name,"Type",4)) || (!strncmp (aux_name,"Tipo",4)) )
        found ++;
  }
 if(!found)
  {
    fprintf(stdout,"\n Couldn't find 'Type' in cluster file\n");
    fflush(stdout);
  }
 return (found);
}

/* ************************************************************************* */

int posic_clu(FILE *fich,char *linea)
{
 int found=0;
 char aux_name[64];

 aux_name[0]='\0';

  while ( !found && (fgets (linea, 2000, fich) != NULL) )
  {
   if ( sscanf (linea,"%s",aux_name) )
     if (!strncmp (aux_name,"Vect",4) )
        found ++;
  }
 if(!found)
  {
    fprintf(stdout,"\n Couldn't find names in cluster file\n");
    fflush(stdout);
  }
 return (found);
}

/* ************************************************************************* */

static int old_sim_format(FILE *fich)
{

 int i = 0, j;
 char linea [2000];
 char formato[2000];
 char aux_name[64];
 int tempo;
 int found=0;

 fseek(fich,0L,SEEK_SET);
 aux_name[0]='\0';
 formato[0]='\0';
 while ( !found && (fgets (linea, 2000, fich) != NULL) )
 {
   if ( sscanf (linea,"%s",aux_name) )
     if (!strncmp (aux_name,"Imag",4))
        found ++;
 }

 if(found)
 {
  i = 0;
  for (;;)
  {   strcpy (formato , "%*s");
      for (j=0; j<i; j++)
         strcat (formato, " %*c %*d");
      strcat (formato, " %*c %d");
      if (sscanf (linea, formato, &tempo) != 1)
        break;
      i++;
  }
 }

 return i;
}

/* *************************************************************************/

static int old_clu_format(FILE *fich)
{
 int i = 0, j;
 char linea [2000];
 char formato[2000];
 float tempo;
 int found=0;

 fseek(fich,0L,SEEK_SET);
 formato[0]='\0';

 found=posic_clu(fich,linea);

 if(found)
 {
    fgets (linea, 2000, fich);

    i=0;
    for (;;)
     {   strcpy (formato , "%*s");
         for (j=0; j<i; j++)
           strcat (formato, " %*f");
         strcat (formato, " %f");
         if (sscanf (linea, formato, &tempo) != 1)
           break;
         i++;
     }
 }
 return i;
}
/* *************************************************************************/

