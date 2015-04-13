/*********************************************************************
*                           L I B _ P W K                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Pablo Chacon and Willy Wriggers, 2001-2003                     *
**********************************************************************
*                                                                    *
* PDB structure manipulation tools.                                  * 
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/


#include "situs.h"

/* external functions */

extern unsigned long gidz_general (int, int, int, unsigned, unsigned);
extern void get_rot_matrix (double [3][3], double, double, double);
extern void zero_vect(double *, unsigned long);
extern void readpdb(char *, int *, PDB **);

/* functions list */

static void copy_atoms (PDB *, PDB *, int, int, int);
void rot_axis (PDB *, PDB *, unsigned, char, double);
void rot_euler (PDB *, PDB *, unsigned, double, double, double);
void translate (PDB *, PDB *, unsigned, double, double, double);
void calc_center (PDB *, unsigned, double *, double *, double *);
double calc_mass (PDB *, unsigned);
double calc_center_mass (PDB *, unsigned, double *, double *, double *);
double calc_sphere (PDB *, unsigned, double, double, double);
void calc_box (PDB *, unsigned, double *, double *, double *, double *, double *, double *);
void project_mass (double **, unsigned long, double, double, double, 
		   unsigned, unsigned, unsigned, PDB *, unsigned, double *);
void read_and_weight_mass (char *, unsigned *, PDB **); 



/*====================================================================*/
static void copy_atoms (PDB *pdb_original, PDB *pdb_duplicate, 
		 int io, int id, int ktot) {
/* copies ktot atoms from *pdb_original to *pdb_duplicate */

  int j, k;
  for (k=0;k<ktot;++k) {
    pdb_duplicate[id+k].weight = pdb_original[io].weight; 
    pdb_duplicate[id+k].x = pdb_original[io].x;   
    pdb_duplicate[id+k].y = pdb_original[io].y;  
    pdb_duplicate[id+k].z = pdb_original[io].z; 
    for (j=0;j<4;++j) pdb_duplicate[id+k].segID[j] = pdb_original[io].segID[j];
    pdb_duplicate[id+k].serial = id+k+1;	
    for (j=0;j<7;++j) pdb_duplicate[id+k].recdName[j] = pdb_original[io].recdName[j];
    for (j=0;j<3;++j) pdb_duplicate[id+k].atomType[j] = pdb_original[io].atomType[j];
    for (j=0;j<3;++j) pdb_duplicate[id+k].atomLoc[j] = pdb_original[io].atomLoc[j];
    for (j=0;j<2;++j) pdb_duplicate[id+k].altLoc[j] = pdb_original[io].altLoc[j];
    for (j=0;j<5;++j) pdb_duplicate[id+k].resName[j] = pdb_original[io].resName[j]; 
    for (j=0;j<2;++j) pdb_duplicate[id+k].chainID[j] = pdb_original[io].chainID[j];
    pdb_duplicate[id+k].resSeq = pdb_original[io].resSeq;
    for (j=0;j<2;++j) pdb_duplicate[id+k].iCode[j] = pdb_original[io].iCode[j]; 
    pdb_duplicate[id+k].occupancy = pdb_original[io].occupancy;   
    pdb_duplicate[id+k].tempFactor = pdb_original[io].tempFactor;
    pdb_duplicate[id+k].ftNote = pdb_original[io].ftNote;  
    for (j=0;j<3;++j) pdb_duplicate[id+k].element[j] = pdb_original[io].element[j];
    for (j=0;j<3;++j) pdb_duplicate[id+k].charge[j] = pdb_original[io].charge[j];
  }
  return;
}


/*====================================================================*/
void rot_axis (PDB *pdb_original, PDB *pdb_rotate, unsigned num_atoms,
		 char axis, double angle ) {
/* structure must be centered, rotates about X, Y, or Z axis */ 

  unsigned id;
  double x,y,z;
  double sint   = sin( angle );
  double cost   = cos( angle );      
  
  switch(axis) { 
  case('X'):
    for (id=0;id<num_atoms;++id) {
      y=pdb_original[id].y;
      z=pdb_original[id].z;
      pdb_rotate[id].x= pdb_original[id].x;
      pdb_rotate[id].y= (cost*y+sint*z);
      pdb_rotate[id].z= (cost*z-sint*y);
    } break;
  case('Y'):
    for (id=0;id<num_atoms;++id) {
      x=pdb_original[id].x;
      z=pdb_original[id].z;
      pdb_rotate[id].y= pdb_original[id].y;
      pdb_rotate[id].x= (cost*x+sint*z);
      pdb_rotate[id].z= (cost*z-sint*x);
    } break;
  case('Z'):
    for (id=0;id<num_atoms;++id) {
      x=pdb_original[id].x;
      y=pdb_original[id].y;
      pdb_rotate[id].z= pdb_original[id].z;
      pdb_rotate[id].x= (cost*x-sint*y);
      pdb_rotate[id].y= (cost*y+sint*x);
    } break;
  }
}


/*====================================================================*/
void rot_euler (PDB *pdb_original, PDB *pdb_rotate, 
		unsigned num_atoms, double psi, double theta, double phi) {
/* structure must be centered, rotates by Euler angles */ 

  unsigned id;
  double rot_matrix[3][3],currx,curry,currz;

  get_rot_matrix(rot_matrix,psi,theta,phi);

  for (id=0;id<num_atoms;++id) {
    currx=pdb_original[id].x; 
    curry=pdb_original[id].y;
    currz=pdb_original[id].z;
    pdb_rotate[id].x = currx * rot_matrix[0][0] +
      curry * rot_matrix[0][1] +
      currz * rot_matrix[0][2];
    pdb_rotate[id].y = currx * rot_matrix[1][0] +
      curry * rot_matrix[1][1] +
      currz * rot_matrix[1][2];
    pdb_rotate[id].z = currx * rot_matrix[2][0] +
      curry * rot_matrix[2][1] +
      currz * rot_matrix[2][2];
  }
}


/*====================================================================*/
void translate (PDB *pdb_original, PDB *pdb_move, 
	       unsigned num_atoms, double x0, double y0, double z0 ) {
/* translates *pdb_original and stores in *pdb_move */
  
  unsigned id;

  for (id=0;id<num_atoms;++id) {
    pdb_move[id].x=pdb_original[id].x+x0; 
    pdb_move[id].y=pdb_original[id].y+y0;
    pdb_move[id].z=pdb_original[id].z+z0;
  }
}


/*====================================================================*/
void calc_center (PDB *pdb0, unsigned num_atoms,
		  double *cx, double *cy, double *cz) {
/* computes geometric center of structure */

  unsigned id;
  
  *cx=0; *cy=0; *cz=0;
  for (id=0;id<num_atoms;++id) {
    *cx+= pdb0[id].x;
    *cy+= pdb0[id].y;
    *cz+= pdb0[id].z;
  }
  if (num_atoms>0) { *cx/=(num_atoms*1.0); *cy/=(num_atoms*1.0);  *cz/=(num_atoms*1.0); }
  else {
    fprintf(stderr, "lib_pwk> Error: Division by zero (atom number) [e.c. 16010]\n"); 
    exit(16010);
  }
}


/*====================================================================*/
double calc_mass (PDB *pdb0, unsigned num_atoms) {
/* returns total mass of structure */
  
  unsigned id;
  double mtot;
  
  mtot=0;
  for (id=0;id<num_atoms;++id) 
    mtot+= pdb0[id].weight; 
  return mtot;
}


/*====================================================================*/
double calc_center_mass (PDB *pdb0, unsigned num_atoms,
		       double *cx, double *cy, double *cz) {
/* computes COM and returns total mass of structure */
  
  unsigned id;
  double mtot;
  
  *cx=0; *cy=0; *cz=0; mtot=0;
  for (id=0;id<num_atoms;++id) {
    mtot += pdb0[id].weight;
    *cx+= pdb0[id].x * pdb0[id].weight;
    *cy+= pdb0[id].y * pdb0[id].weight;
    *cz+= pdb0[id].z * pdb0[id].weight;
  }
  if (mtot>0) { *cx/=mtot; *cy/=mtot; *cz/=mtot; }
  else {
    fprintf(stderr, "lib_pwk> Error: Division by zero (mass) [e.c. 16020]\n"); 
    exit(16020);
  }
  return mtot;
}


/*====================================================================*/
double calc_sphere (PDB *pdb0, unsigned num_atoms,
		   double cx, double cy, double cz) {
/* computes bounding radius of structure relative to input center */
  
  unsigned id;
  double maxradius, currradius;
  
  maxradius=-1e20;
  for (id=0;id<num_atoms;++id) {
    currradius=((pdb0[id].x-cx)*(pdb0[id].x-cx)+
		(pdb0[id].y-cy)*(pdb0[id].y-cy)+
		(pdb0[id].z-cz)*(pdb0[id].z-cz));
    if (currradius > maxradius) maxradius = currradius;
  } 

  if (maxradius >= 0) {
    maxradius = sqrt(maxradius);
    return maxradius;
  } else {
    fprintf(stderr, "lib_pwk> Error: no bounding sphere found [e.c. 16030]\n"); 
    exit(16030);
    return 0;
  }
}


/*====================================================================*/
void calc_box (PDB *pdb0, unsigned num_atoms,
		double *minx, double *miny, double *minz,
		double *maxx, double *maxy, double *maxz) {
/* computes bounding box of structure */
  
  unsigned id;

  *minx = 1e20; *miny = 1e20; *minz = 1e20; *maxx = -1e20; *maxy = -1e20; *maxz = -1e20;
  for (id=0;id<num_atoms;++id) {
    if (*minx > pdb0[id].x) *minx = pdb0[id].x;
    if (*maxx < pdb0[id].x) *maxx = pdb0[id].x;
    if (*miny > pdb0[id].y) *miny = pdb0[id].y;
    if (*maxy < pdb0[id].y) *maxy = pdb0[id].y;
    if (*minz > pdb0[id].z) *minz = pdb0[id].z;
    if (*maxz < pdb0[id].z) *maxz = pdb0[id].z; 
  }
  if (*minx < 1e20 && *miny < 1e20 && *minz < 1e20 && *maxx > -1e20 && *maxy > -1e20 && *maxz > -1e20) return;
  else {
    fprintf(stderr, "lib_pwk> Error: no bounding box found [e.c. 16040]\n"); 
    exit(16040);
  }
}



/*====================================================================*/
void project_mass (double **outmap, unsigned long nvox, 
		   double widthx, double widthy,  double widthz,
		   unsigned extx, unsigned exty, unsigned extz, 
		   PDB *inpdb, unsigned num_atoms, double shift[3]) {
/* projects to lattice using mass-weighting                          */
/* assumes both structure and map are centered (before any shifting) */

  int x0, y0, z0, x1, y1, z1, i;
  double gx, gy, gz;
  double a, b, c;
  double ab, ab1, a1b, a1b1;

  zero_vect(*outmap, nvox);  
  for (i=0;i<num_atoms;++i) {
    /* compute position within grid  */
    gx = extx/2.0+ (inpdb[i].x+shift[0]) / widthx;
    gy = exty/2.0+ (inpdb[i].y+shift[1]) / widthy;
    gz = extz/2.0+ (inpdb[i].z+shift[2]) / widthz;
    x0 = floor (gx); 
    y0 = floor (gy); 
    z0 = floor (gz);
    x1 = x0+1;
    y1 = y0+1; 
    z1 = z0+1; 
    /* interpolate */
    a = x1-gx;
    b = y1-gy;
    c = z1-gz;
    ab= a*b;
    ab1=a * (1-b);
    a1b=(1-a) * b;
    a1b1=(1-a) * (1-b);
    a=(1-c);
    *(*outmap+gidz_general(z0,y0,x0,exty,extx)) += ab   * c * inpdb[i].weight; 
    *(*outmap+gidz_general(z1,y0,x0,exty,extx)) += ab   * a * inpdb[i].weight; 
    *(*outmap+gidz_general(z0,y1,x0,exty,extx)) += ab1  * c * inpdb[i].weight; 
    *(*outmap+gidz_general(z1,y1,x0,exty,extx)) += ab1  * a * inpdb[i].weight; 
    *(*outmap+gidz_general(z0,y0,x1,exty,extx)) += a1b  * c * inpdb[i].weight;
    *(*outmap+gidz_general(z1,y0,x1,exty,extx)) += a1b  * a * inpdb[i].weight; 
    *(*outmap+gidz_general(z0,y1,x1,exty,extx)) += a1b1 * c * inpdb[i].weight;
    *(*outmap+gidz_general(z1,y1,x1,exty,extx)) += a1b1 * a * inpdb[i].weight;
  }
} 







/* note: the following function should be consolidated with readpdb */
/* this would simplify recognition and assignment of mass, and VQ could estimate number of */
/* representative vectors from .weight; num_atoms should be unsigned */
/* in situs.h see colores.h, definition of PDB structure: .weigth */
/* need to support both Situs and colores massweighting */
/* adopt Pablo's element routines in mass weighting (better) */


/*====================================================================*/
void read_and_weight_mass (char *filename, unsigned *num_atoms_rec, PDB **pdb_return) {
/* reads PDB file and massweights all atoms */

  static char *atom_name[] = {
    "H", "HE", "LI", "BE",  "B",  "C",  "N",  "O",  "F", "NE", 
    "NA", "MG", "AL", "SI",  "P",  "S", "CL", "AR",  "K", "CA", 
    "SC", "TI",  "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN"
  };
  static float atom_mass[30] = {
    1.008 ,  4.003,  6.940,  9.012, 10.810, 12.011, 14.007, 16.000,  18.998, 20.170, 
    20.170, 24.305, 26.982, 28.086, 30.974, 32.060, 35.453, 39.948,  39.102, 40.080, 
    44.956, 47.880, 50.040, 51.996, 54.938, 55.847, 58.933, 58.710,  63.546, 65.380
  };
  static int bio_element[6] = {0,5,6,7,14,15};

  PDB *pdb_original, *pdb_recognized;
  int i, j, k;
  int num_atoms;

  readpdb(filename, &num_atoms, &pdb_original);

  /* create the array of recognized atoms */
  
  pdb_recognized = (PDB *) malloc(num_atoms*sizeof(PDB));
  if (pdb_recognized== NULL) {
    fprintf(stderr, "lib_pwk> Error: Unable to satisfy memory allocation request [e.c. 16050]\n"); 
    exit(16050);
  }

  *num_atoms_rec = 0;
  for (i=0;i<num_atoms;++i) {
    for (j=0;j<30;++j) 
      if (strcmp(pdb_original[i].atomType, atom_name[j])==0) {
        pdb_original[i].weight=atom_mass[j]; 
        copy_atoms (pdb_original, pdb_recognized, i, *num_atoms_rec, 1);        
        ++*num_atoms_rec;
        break;
      }   
    
    /* here we look again for the most frequent atom types */  
    if (j==30) {
      for (k=0;k<6;++k) {      
	if (strncmp(pdb_original[i].atomType, atom_name[bio_element[k]],1)==0) {
	  pdb_original[i].weight=atom_mass[bio_element[k]]; 
	  copy_atoms (pdb_original, pdb_recognized, i, *num_atoms_rec, 1);        
	  ++*num_atoms_rec;
	  break;
	}
      }
      if (k==6) printf("lib_pwk> Warning: Unable to identify atom type %s!\n",pdb_original[i].atomType);
    }
  } 
  
  free(pdb_original);
  *pdb_return=pdb_recognized;
}
