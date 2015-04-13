/*********************************************************************
*                          L I B _ P I O                             *
**********************************************************************
* Library is part of the Situs package (c) Willy Wriggers, 1998-2003 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    * 
* Auxiliary program to read and write files in PDB format.           *
*                                                                    *
**********************************************************************
* Adapted with permission from components of the DOWSER program -    *
* a product of UNC Computational Structural Biology group	     *	
* Authors: Li Zhang, Xinfu Xia, Jan Hermans, Univ. of North Carolina *	
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"

void    readpdb(char *, int *, PDB **);
void	fld2s(char *, char * );
void	fld2i(char *, int, int * );
void	readCol(char *, int, int, char * );
void    writepdb(char *, int, PDB *);
void    appendpdb(char *, int, PDB *);
void    space(int);

int shortLine;
FILE *fp;

void readpdb(char *pdb_file, int *numAtoms, PDB **Atoms)
{
  FILE          *file_ref;
  char	line[101];
  char	field00[7], field01[6], field02[2], field03[3], field04[3],
	field05[2], field06[5], field08[2], field09[5],
	field10[2], field11[4], field12[9], field13[9], field14[9],
	field15[7], field16[7], field17[2], field18[4], field19[3],
	field20[5], field21[3], field22[3];

  char	recdName[7];	/*	 1 -  6	*/

  int	i = 0;
  int	j, lineWidth, nAtom;
  PDB	*atoms, *pureAtoms;

  atoms = (PDB *) malloc( (MAXPDB)*sizeof(PDB) );
  if (atoms == NULL) {
    fprintf(stderr, "lib_pio> Error: Unable to satisfy memory allocation request [e.c. 12010]\n"); 
    exit(12010);
  }

  if( pdb_file == NULL ) {
    fprintf(stderr, "lib_pio> Error: no file name in argument [e.c. 12110]\n");
    exit(12110);   
  } else {
    file_ref = fopen(pdb_file, "r");
    if(file_ref == NULL ) {
      fprintf(stderr, "lib_pio> Error: can't open file! %s  [e.c. 12120]\n", pdb_file); 
      exit(12120);   
    }
  }
  
  for(;;) {
    char* nread=fgets(line, 100, file_ref);
    if (feof(file_ref) || i+1 > MAXPDB) break;
    if (sscanf(line, "%6s", recdName) != 1) continue;
    if (strcmp(recdName,"ATOM")==0 || strcmp(recdName,"HETATM")==0 || strcmp(recdName,"LABEL")==0 ) {
      lineWidth = 80;
      for( j=0; j<80; j++) {
        if( line[j]=='\n' ) { lineWidth = j; break; }
      }
      for( j=lineWidth; j<80; j++) {
        line[j]=' ';
      }

      readCol(line,  1,  6, field00);  fld2s(field00, atoms[i].recdName);
      readCol(line,  7, 11, field01);  fld2i(field01, 5,&(atoms[i].serial));
      readCol(line, 12, 12, field02);
      readCol(line, 13, 14, field03);  fld2s(field03, atoms[i].atomType);
      readCol(line, 15, 16, field04);  fld2s(field04, atoms[i].atomLoc);
      readCol(line, 17, 17, field05);  fld2s(field05, atoms[i].altLoc);
      readCol(line, 18, 21, field06);  fld2s(field06, atoms[i].resName);
      readCol(line, 22, 22, field08);  fld2s(field08, atoms[i].chainID);
      readCol(line, 23, 26, field09);  fld2i(field09, 4,&(atoms[i].resSeq));
      readCol(line, 27, 27, field10);  fld2s(field10, atoms[i].iCode);
      readCol(line, 28, 30, field11);
      readCol(line, 31, 38, field12);  atoms[i].x = (float) atof(field12);
      readCol(line, 39, 46, field13);  atoms[i].y = (float) atof(field13);
      readCol(line, 47, 54, field14);  atoms[i].z = (float) atof(field14);
      readCol(line, 55, 60, field15);  atoms[i].occupancy = (float) atof(field15);
      readCol(line, 61, 66, field16);  atoms[i].tempFactor = (float) atof(field16);
      readCol(line, 67, 67, field17);
      readCol(line, 68, 70, field18);  fld2i(field18, 3,&(atoms[i].ftNote));
      readCol(line, 71, 72, field19);
      readCol(line, 73, 76, field20);  fld2s(field20, atoms[i].segID);
      readCol(line, 77, 78, field21);  fld2s(field21, atoms[i].element);
      readCol(line, 79, 80, field22);  fld2s(field22, atoms[i].charge);
      i++;
    }
  }
  fclose(file_ref);

  if (i==MAXPDB) {
    fprintf(stderr, "lib_pio> Error: Exceeded MAXPDB parameter in file situs.h [e.c. 12210]\n"); 
    exit(12210);
  }

  nAtom = i;
  printf("lib_pio> %d atoms read.\n", nAtom);

  pureAtoms = (PDB *) malloc( nAtom*sizeof(PDB) );
  if (pureAtoms == NULL) {
    fprintf(stderr, "lib_pio> Error: Unable to satisfy memory allocation request [e.c. 12220]\n"); 
    exit(12220);
  }  
  for( i=0; i<nAtom; i++) pureAtoms[i] = atoms[i];
  *numAtoms = nAtom;
  *Atoms = pureAtoms;
  free( atoms );
  return;
}

void readCol(char * line, int i1, int i2, char * field00 )
{
  int i;
  for( i=i1; i<=i2; i++) {
    field00[i-i1] = line[i-1];
  }
  field00[i-i1] = '\0';
}

void fld2s(char * field, char * str)
{
  int i;
  i = sscanf(field, "%s", str);
  if( i < 1 ) { str[0] = '\0'; }
}

void fld2i(char * field, int n, int * num)
{
  int i;
  for( i=0; i<n; i++) {
    if( field[i] != ' ' ) { sscanf(field, "%d", num); return; }
  }
  if( i == n ) { *num = 0;  return; }
}

void writepdb(char *filename, int nAtoms, PDB * atoms)
{
  int i, i1, j, nSpace;

  if( filename == NULL ) {
    fprintf(stderr, "lib_pio> Error: no file name in argument [e.c. 12310]\n");
    exit(12310);   
  }
  
  fp = fopen(filename, "w");
  if(fp == NULL ) {
    fprintf(stderr, "lib_pio> Error: can't open file! %s  [e.c. 12320]\n", filename); 
    exit(12320);   
  }

  for( i=0; i<nAtoms; i++) {
    nSpace = 6 - strlen(atoms[i].recdName);
    fprintf(fp, "%s", atoms[i].recdName);  space(nSpace);

    /* fprintf(fp, "%5d", atoms[i].serial); */
    /* renumber atoms using 5 digits, ignore atoms[i].serial */
    if ((i+1)/100000 == 0) {
      i1 = i+1;
      fprintf(fp, "%5d", i1);
    } else {
      i1=i+1-((i+1)/100000)*100000;
      fprintf(fp, "%05d", i1); 
    }
    space(1);

    nSpace = 2 - strlen(atoms[i].atomType); space(nSpace);
    fprintf(fp, "%s", atoms[i].atomType);

    nSpace = 2 - strlen(atoms[i].atomLoc);
    fprintf(fp, "%s", atoms[i].atomLoc);  space(nSpace);

    space(1);

    j = strlen(atoms[i].resName);
    nSpace = 3 - j;
    if (j<4) {
      space(nSpace);  fprintf(fp, "%s", atoms[i].resName); space(1);
    } else fprintf(fp, "%s", atoms[i].resName);

    nSpace = 1 - strlen(atoms[i].chainID); space(nSpace);
    fprintf(fp, "%s", atoms[i].chainID);

    fprintf(fp, "%4d", atoms[i].resSeq);

    nSpace = 1 - strlen(atoms[i].iCode); space(nSpace);
    fprintf(fp, "%s", atoms[i].iCode);

    space(3);

    fprintf(fp, "%8.3f", atoms[i].x);
    fprintf(fp, "%8.3f", atoms[i].y);
    fprintf(fp, "%8.3f", atoms[i].z);
    fprintf(fp, "%6.2f", atoms[i].occupancy);
    fprintf(fp, "%6.2f", atoms[i].tempFactor);

    space(1);

    fprintf(fp, "%3d",atoms[i].ftNote);

    space(2);

    nSpace = 4 - strlen(atoms[i].segID);
    fprintf(fp, "%s", atoms[i].segID);  space(nSpace);

    nSpace = 2 - strlen(atoms[i].element);
    space(nSpace);  fprintf(fp, "%s", atoms[i].element);

    nSpace = 2 - strlen(atoms[i].charge);
    fprintf(fp, "%s", atoms[i].charge);  space(nSpace);

    fprintf(fp, "\n");
  }
  fclose(fp);
  return;
}

void appendpdb(char *filename, int nAtoms, PDB * atoms)
{
  int i, i1, j, nSpace;

  if( filename == NULL ) {
    fprintf(stderr, "lib_pio> Error: no file name in argument [e.c. 12410]\n");
    exit(12410);   
  }
  
  fp = fopen(filename, "a");
  if(fp == NULL ) {
    fprintf(stderr, "lib_pio> Error: can't open file! %s  [e.c. 12420]\n", filename); 
    exit(12420);   
  }

  for( i=0; i<nAtoms; i++) {
    nSpace = 6 - strlen(atoms[i].recdName);
    fprintf(fp, "%s", atoms[i].recdName);  space(nSpace);

    /* renumber atoms using 5 digits */
    if (atoms[i].serial/100000 == 0) {
      fprintf(fp, "%5d", atoms[i].serial);
    } else {
      i1=atoms[i].serial-(atoms[i].serial/100000)*100000;
      fprintf(fp, "%05d", i1);
    }
    space(1);

    nSpace = 2 - strlen(atoms[i].atomType); space(nSpace);
    fprintf(fp, "%s", atoms[i].atomType);

    nSpace = 2 - strlen(atoms[i].atomLoc);
    fprintf(fp, "%s", atoms[i].atomLoc);  space(nSpace);

    space(1);

    j = strlen(atoms[i].resName);
    nSpace = 3 - j;
    if (j<4) {
      space(nSpace);  fprintf(fp, "%s", atoms[i].resName); space(1);
    } else fprintf(fp, "%s", atoms[i].resName);

    nSpace = 1 - strlen(atoms[i].chainID); space(nSpace);
    fprintf(fp, "%s", atoms[i].chainID);

    fprintf(fp, "%4d", atoms[i].resSeq);

    nSpace = 1 - strlen(atoms[i].iCode); space(nSpace);
    fprintf(fp, "%s", atoms[i].iCode);

    space(3);

    fprintf(fp, "%8.3f", atoms[i].x);
    fprintf(fp, "%8.3f", atoms[i].y);
    fprintf(fp, "%8.3f", atoms[i].z);
    fprintf(fp, "%6.2f", atoms[i].occupancy);
    fprintf(fp, "%6.2f", atoms[i].tempFactor);

    space(1);

    fprintf(fp, "%3d",atoms[i].ftNote);

    space(2);

    nSpace = 4 - strlen(atoms[i].segID);
    fprintf(fp, "%s", atoms[i].segID);  space(nSpace);

    nSpace = 2 - strlen(atoms[i].element);
    space(nSpace);  fprintf(fp, "%s", atoms[i].element);

    nSpace = 2 - strlen(atoms[i].charge);
    fprintf(fp, "%s", atoms[i].charge);  space(nSpace);

    fprintf(fp, "\n");
  }
  fclose(fp);
  return;
}

void space(int nSpace)
{
  int j;
  for( j=0; j<nSpace; j++)
  fprintf(fp, " ");
}


