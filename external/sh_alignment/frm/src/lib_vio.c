/*********************************************************************
*                            L I B _ V I O                           *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 1998-2003 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    * 
* Auxiliary program to read and write EM maps in Situs format        *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"            

/* function declarations */

void readvol (char *, double *, double *, double *, double *, int *, int *, int *, double **);
void writevol (char *, double, double, double, double, int, int, int, double *);

/*====================================================================*/
void readvol(char *vol_file, double *width, double *gridx, double *gridy, double *gridz, int *extx, int *exty, int *extz, double **phi)

{ unsigned long nvox, count;
  double dgridx, dgridy, dgridz, dwidth;
  double phitest, dtemp;
  FILE *fin;
   
  fin = fopen(vol_file, "r");
  if( fin == NULL ) {
    fprintf(stderr, "lib_vio> Error: Can't open file! %s  [e.c. 13010]\n", vol_file); 
    exit(13010);   
  }

  /* read header and print information */
  int nread=fscanf(fin, "%le %le %le %le %d %d %d", &dwidth, &dgridx, &dgridy, &dgridz, extx, exty, extz);
  *width = dwidth; *gridx = dgridx; *gridy = dgridy; *gridz = dgridz;
  printf ("lib_vio> File %s - Header information: \n", vol_file);
  printf ("lib_vio> Columns, rows, and sections: x=1-%d, y=1-%d, z=1-%d\n",*extx,*exty,*extz);
  printf ("lib_vio> 3D coordinates of first voxel (1,1,1): (%f,%f,%f)\n",*gridx,*gridy,*gridz);  
  printf ("lib_vio> Voxel size in Angstrom: %f \n", *width);
  nvox = *extx * *exty * *extz;

  /* allocate memory and read data */
  printf ("lib_vio> Reading density data... \n");
  *phi = (double *) malloc(nvox * sizeof(double));
  if (*phi == NULL) {
    fprintf(stderr, "lib_vio> Error: Unable to satisfy memory allocation request [e.c. 13020]\n"); 
    exit(13020);
  }

  for (count=0;count<nvox;count++) {
    if (fscanf(fin,"%le", &dtemp) != 1) {
        fprintf(stderr, "lib_vio> Error: file %s is too short (not compatible with header information) or data unreadable [e.c. 13030]\n", vol_file); 
        exit(13030);
    } else *(*phi+count) = dtemp;
  }
  if (fscanf(fin,"%le", &phitest) != EOF) {
    fprintf(stderr, "lib_vio> Error: file %s is too long, not compatible with header information [e.c. 13040]\n", vol_file); 
    exit(13040);
  }
  fclose(fin);  
  printf ("lib_vio> Volumetric data read from file %s\n", vol_file); 
  return;
}



/*====================================================================*/
void writevol(char *vol_file, double width, double gridx, double gridy, double gridz, int extx, int exty, int extz, double *phi)

{ unsigned long nvox, count;
  FILE *fout;
   
  nvox = extx * exty * extz;
  fout = fopen(vol_file, "w");
  if( fout == NULL ) {
    fprintf(stderr, "lib_vio> Error: Can't open file! %s  [e.c. 13210]\n", vol_file); 
    exit(13210);   
  }

  printf ("lib_vio> Writing density data... \n");
  fprintf(fout, "  %f %f %f %f %d %d %d\n", width, gridx, gridy, gridz, extx, exty, extz);
  fprintf(fout, "\n");

  for(count=0;count<nvox;count++) {
    if ((count+1)%10 == 0) fprintf (fout," %10.6f \n",*(phi+count));
    else fprintf (fout," %10.6f ",*(phi+count)); 
  }    
  fclose(fout);  
  printf ("lib_vio> Volumetric data written to file %s \n", vol_file);

  /* header information */
  printf ("lib_vio> File %s - Header information: \n", vol_file);
  printf ("lib_vio> Columns, rows, and sections: x=1-%d, y=1-%d, z=1-%d\n",extx,exty,extz);
  printf ("lib_vio> 3D coordinates of first voxel (1,1,1): (%f,%f,%f)\n",gridx,gridy,gridz);
  printf ("lib_vio> Voxel size in Angstrom: %f \n", width);
  return;
}

