/*********************************************************************
*                           L I B _ E U L                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Pablo Chacon and Willy Wriggers, 2001-2003                     *
**********************************************************************
*                                                                    *
* Euler angle - related routines.                                    * 
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

/*

Euler angle order: (psi,theta,phi) <-> (0,1,2)
Euler angle convention: Goldstein or "X-convention":
 
       A=BCD
 
rotation components:

   |  cos_psi  sin_psi 0 |
B= | -sin_psi  cos_psi 0 |
   |    0        0     1 |

   |  1      0          0      |
C= |  0   cos_theta  sin_theta |
   |  0  -sin_theta  cos_theta |

   |  cos_phi  sin_phi 0 |
D= | -sin_phi  cos_phi 0 |
   |    0        0     1 |
 
rotation matrix:
 
  |a_11 a_12 a_13|
A=|a_21 a_22 a_23|
  |a_31 a_22 a_33|

a_11 = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi;
a_12 = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi;
a_13 = sin_psi * sin_theta;
a_21 = -sin_psi * cos_phi  - cos_theta * sin_phi * cos_psi;
a_22 = -sin_psi * sin_phi  + cos_theta * cos_phi * cos_psi;
a_23 =  cos_psi * sin_theta;
a_31 =  sin_theta * sin_phi;
a_32 = -sin_theta * cos_phi;
a_33 =  cos_theta;
 
See http://mathworld.wolfram.com/
*/



#include "situs.h"

/* functions list */

void write_eulers (char *, unsigned long, float *, double [3][2], double);
void read_eulers (char *, unsigned long *, float **);
void remap_eulers (double *,double *,double *,double,double,double,double,double,double); 
void get_rot_matrix (double [3][3], double, double, double);
void eu_spiral (double [3][2], double, unsigned long *, float **); 
void eu_lattman (double [3][2], double, unsigned long *, float **); 
static double deviation_cosine (double [3][3], double [3][3]); 
char similar_eulers(double, double, double, double, double, double);



/*=====================================================================*/
void write_eulers (char *eu_file, unsigned long eu_count, float *eu_store, double eu_range[3][2], double delta) {
/* writes range, step_size, and all angles to specified filename       */

  FILE *out;
  unsigned i;
  
  if (( out = fopen(eu_file, "w")) == NULL) {
    fprintf(stderr, "lib_eul> Error: Unable to open file %s  [e.c. 15010]\n", eu_file); 
    exit(15010);
  }
  
  for (i=0;i<eu_count;i++) 
    fprintf(out, "%10.4f%10.4f%10.4f\n",
	    *(eu_store+3*i+0)/ROT_CONV,*(eu_store+3*i+1)/ROT_CONV,*(eu_store+3*i+2)/ROT_CONV);
  fclose(out);
}



/*=====================================================================*/
void read_eulers (char *in_file, unsigned long *eu_count, float **eu_store) {
/* reads Euler angles from in_file */

  FILE *in;
  unsigned long i;
  double psi,theta,phi;

  /* count number of Euler angle sets */

  if ((in = fopen(in_file, "r")) == NULL) {
    fprintf(stderr, "lib_eul> Error: Unable to open file %s [e.c. 15020]\n", in_file); 
    exit(15020);
  }
  i=0;
  while (!feof(in)) {
    if (fscanf(in,"%lf %lf %lf\n",&psi,&theta,&phi) != 3) break;
    else i++;
  }
  fclose(in);

  /* allocate eu_store */
  *eu_store = (float *) malloc(i * 3 * sizeof(float));
  if (*eu_store == NULL) {
    fprintf(stderr, "lib_eul> Error: Unable to satisfy memory allocation request [e.c. 15040]\n"); 
    exit(15040);
  }
  if ((in = fopen(in_file, "r")) == NULL) {
    fprintf(stderr, "lib_eul> Error: Unable to open file %s  [e.c. 15050]\n", in_file); 
    exit(15050);
  }
  i=0;
  while (!feof(in)) {
    if (fscanf(in,"%lf %lf %lf\n",&psi,&theta,&phi) != 3) break;
    else {
      *(*eu_store+i+0)=psi   * ROT_CONV;
      *(*eu_store+i+1)=theta * ROT_CONV;
      *(*eu_store+i+2)=phi   * ROT_CONV;
      i++;
    }
  }
  fclose(in);
  *eu_count=i; 
  printf("lib_eul> %d Euler angles read from file %s\n", (int)i, in_file);
}


/*=====================================================================*/
void remap_eulers (double *psi_out, double *theta_out, double *phi_out,
		   double psi_in, double theta_in, double phi_in,
		   double psi_ref, double theta_ref, double phi_ref) {
/* remaps Euler angles using rot-matrix invariant transformations  */
/* and attempts to make angles fall in the intervalls:             */
/* [psi_ref,   psi_ref + 2*M_PI]                                   */
/* [theta_ref, theta_ref + M_PI]                                   */
/* [phi_ref,   phi_ref + 2*M_PI]                                   */
/* Note that this is strictly guaranteed only for                  */
/* (psi_ref, theta_ref, phi_ref) == (0,0,0)                        */
/* For non-zero reference we may have for small number of angles   */
/* theta_ref + M_PI <= theta <= theta_ref + 2 * M_PI               */
 
  double curr_psi, curr_theta, curr_phi;
  double new_psi, new_theta, new_phi;

  int i,j;
  double testa,testb;
  double dump[3][3];

  /* bring psi, theta, phi, within 2 PI of reference */

  curr_psi = psi_in - psi_ref; 
  if (curr_psi >= 0) new_psi = fmod(curr_psi,2*M_PI) + psi_ref;
  else new_psi = 2*M_PI - fmod(-curr_psi,2*M_PI) + psi_ref;
 
  curr_theta = theta_in - theta_ref;
  if (curr_theta >= 0) new_theta = fmod(curr_theta,2*M_PI) + theta_ref;
  else new_theta = 2*M_PI - fmod(-curr_theta,2*M_PI) + theta_ref;
  
  curr_phi = phi_in - phi_ref;
  if (curr_phi >= 0) new_phi = fmod(curr_phi,2*M_PI) + phi_ref;
  else new_phi = 2*M_PI - fmod(-curr_phi,2*M_PI) + phi_ref;
  
  /* if theta is not within PI, we use invariant transformations */
  /* and attempt to map to above intervals                       */
  /* this works in most cases even if the reference is not zero  */
  
  if (new_theta - theta_ref > M_PI) { /* theta overflow */

    /* theta -> 2 PI - theta */
    if (new_theta >= 0) curr_theta = fmod(new_theta,2*M_PI);
    else curr_theta = 2*M_PI - fmod(-new_theta,2*M_PI);
    new_theta -= 2 * curr_theta;

    /* remap to [0, 2 PI] interval */
    curr_theta = new_theta - theta_ref;
    if (curr_theta >= 0) new_theta = fmod(curr_theta,2*M_PI) + theta_ref;
    else new_theta = 2*M_PI - fmod(-curr_theta,2*M_PI) + theta_ref;

    /* we have flipped theta so we need to flip psi and phi as well */
    /* to keep rot-matrix invariant                                 */

    /* psi -> psi + PI */
    if (new_psi - psi_ref > M_PI) new_psi -= M_PI;
    else new_psi += M_PI;
    
    /* phi -> phi + PI */
    if (new_phi - phi_ref > M_PI) new_phi -= M_PI;
    else new_phi += M_PI;
  }

  *psi_out = new_psi;
  *theta_out = new_theta;
  *phi_out = new_phi;
}


/*=====================================================================*/
void get_rot_matrix (double dump[3][3], double psi, double theta, double phi) {
/* computes rotation matrix based on Euler angles psi, theta, phi in radians */  

  double sin_psi   = sin( psi );
  double cos_psi   = cos( psi );
  double sin_theta = sin( theta );
  double cos_theta = cos( theta );
  double sin_phi   = sin( phi);
  double cos_phi   = cos( phi );

  /* use Goldstein convention */
  dump[0][0] = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi;
  dump[0][1] = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi;
  dump[0][2] = sin_psi * sin_theta;
  dump[1][0] = -sin_psi * cos_phi  - cos_theta * sin_phi * cos_psi;
  dump[1][1] = -sin_psi * sin_phi  + cos_theta * cos_phi * cos_psi;
  dump[1][2] =  cos_psi * sin_theta;
  dump[2][0] =  sin_theta * sin_phi;
  dump[2][1] = -sin_theta * cos_phi;
  dump[2][2] =  cos_theta;
}


/*=====================================================================*/
void eu_spiral (double eu_range[3][2], double delta, unsigned long *eu_count, float **eu_store) 
/* Spiral algorithm for Euler angle computation in specified range.    */ 
/* The idea behind the algorithm is that one can cut the globe with n  */
/* horizontal planes spaced 2/(n-1) units apart, forming n circles of  */
/* latitude on the sphere, each latitude containing one spiral point.  */
/* To obtain the kth spiral point, one proceeds upward from the        */
/* (k-1)st point (theta(k-1), psi(k-1)) along a great circle to the    */
/* next latitude and travels counterclockwise along it for a fixed     */
/* distance to arrive at the kth point (theta(k), psi(k)).             */
/*                                                                     */
/* Refs: (1) E.B. Saff and A.B.J. Kuijlaars, Distributing Many         */
/*           Points on a Sphere, The Mathematical Intelligencer,       */
/*           19(1), Winter (1997).                                     */
/*       (2) "Computational Geometry in C." Joseph O'Rourke            */

{ unsigned long i,j;
  int phi_steps,n,k;
  double phi, phi_tmp, psi_old, psi_new, psi_tmp, theta, theta_tmp, h;

  /* rough estimate of number of points on sphere that give a surface  */
  /* density that is the squared linear density of angle increments    */
  n=(int)ceil(360.0*360.0/(delta*delta*M_PI));

  /* total nr. points = nr. of points on the sphere * phi increments */
  phi_steps = (int)ceil((eu_range[2][1]-eu_range[2][0])/delta);
  if (phi_steps<1) {
    fprintf(stderr, "lib_eul> Error: negative number of Euler angle steps [e.c. 15080]\n"); 
    exit(15080);
  }

  /* count number of points on the (theta,psi) sphere */
  j=0;

  /* lower pole on (theta,psi) sphere, k=0 */
  theta = M_PI;
  psi_new = (eu_range[0][1]+eu_range[0][0])*0.5*ROT_CONV;
  remap_eulers (&psi_new, &theta, &phi_tmp, psi_new, theta, 0.0, 
		eu_range[0][0]*ROT_CONV, eu_range[1][0]*ROT_CONV, eu_range[2][0]*ROT_CONV);
  if (eu_range[1][0]*ROT_CONV <= theta && eu_range[1][1]*ROT_CONV >= theta ) j++;
  
  /* intermediate sphere latitudes theta */
  psi_old = 0;                 /* longitude */
  for (k=1;k<n-1;k++) {
    h = -1 + 2 * k / (n-1.0); /* linear distance from pole */
    theta = acos(h);
    psi_new = psi_old + 3.6 / (sqrt ((double)n * (1-h*h))); 
    psi_old = psi_new;
    
    remap_eulers (&psi_new, &theta, &phi_tmp, psi_new, theta, 0.0, 
		  eu_range[0][0]*ROT_CONV, eu_range[1][0]*ROT_CONV, eu_range[2][0]*ROT_CONV);
    
    if (eu_range[0][0]*ROT_CONV <= psi_new && eu_range[0][1]*ROT_CONV >= psi_new && eu_range[1][0]*ROT_CONV <= theta && eu_range[1][1]*ROT_CONV >= theta ) j++;
  }
 
  /* upper pole on (theta,psi) sphere, k=n-1 */
  theta = 0.0;
  psi_new = (eu_range[0][1]+eu_range[0][0])*0.5*ROT_CONV;
  remap_eulers (&psi_new, &theta, &phi_tmp, psi_new, theta, 0.0, 
		eu_range[0][0]*ROT_CONV, eu_range[1][0]*ROT_CONV, eu_range[2][0]*ROT_CONV);
  if (eu_range[1][0]*ROT_CONV <= theta && eu_range[1][1]*ROT_CONV >= theta ) j++;
  
  i=phi_steps*j;
  *eu_count=i;  
  printf("lib_eul> Spiral Euler angle distribution, total number %d (delta = %f deg.)\n",(int)i,delta);
  
  /* allocate memory */
  *eu_store = (float *) malloc(i * 3 * sizeof(float));
  if (*eu_store == NULL) {
    fprintf(stderr, "lib_eul> Error: Unable to satisfy memory allocation request for %d orientations [e.c. 15110]\n",(int) i);
    exit(15110);
  }

  j=0;
  /* lower pole on (theta,psi) sphere, k=0 */
  theta = M_PI;
  psi_new = (eu_range[0][1]+eu_range[0][0])*0.5*ROT_CONV;
  remap_eulers (&psi_new, &theta, &phi_tmp, psi_new, theta, 0.0, 
		eu_range[0][0]*ROT_CONV, eu_range[1][0]*ROT_CONV, eu_range[2][0]*ROT_CONV);
  if (eu_range[1][0]*ROT_CONV <= theta && eu_range[1][1]*ROT_CONV >= theta ) {
    for (phi=eu_range[2][0];phi<=eu_range[2][1];phi+=delta) {
      remap_eulers (&psi_tmp, &theta_tmp, &phi_tmp, psi_new, theta, phi*ROT_CONV, 
		    0.0, 0.0, 0.0);
      *(*eu_store+j+0)=psi_tmp;
      *(*eu_store+j+1)=theta_tmp;
      *(*eu_store+j+2)=phi_tmp; 
      j+=3;
    }
  }

  /* intermediate sphere latitudes theta */
  psi_old = 0;   /* longitude */
  for (k=1;k<n-1;k++) {
    h = -1 + 2 * k / (n-1.0); /* linear distance from pole */
    theta = acos(h);
    psi_new = psi_old + 3.6 / (sqrt ((double)n * (1-h*h))); 
    psi_old = psi_new; 
    remap_eulers (&psi_new, &theta, &phi_tmp, psi_new, theta, 0.0, 
		  eu_range[0][0]*ROT_CONV, eu_range[1][0]*ROT_CONV, eu_range[2][0]*ROT_CONV);
    if (eu_range[0][0]*ROT_CONV <= psi_new && eu_range[0][1]*ROT_CONV >= psi_new && eu_range[1][0]*ROT_CONV <= theta && eu_range[1][1]*ROT_CONV >= theta ) {
      for (phi=eu_range[2][0];phi<=eu_range[2][1];phi+=delta) {
	remap_eulers (&psi_tmp, &theta_tmp, &phi_tmp, psi_new, theta, phi*ROT_CONV, 
		      0.0, 0.0, 0.0);
	*(*eu_store+j+0)=psi_tmp;
	*(*eu_store+j+1)=theta_tmp;
	*(*eu_store+j+2)=phi_tmp; 
	j+=3;
      }
    }
  }
 
  /* upper pole on (theta,psi) sphere, k=n-1 */
  theta = 0.0;
  psi_new = (eu_range[0][1]+eu_range[0][0])*0.5*ROT_CONV;
  remap_eulers (&psi_new, &theta, &phi_tmp, psi_new, theta, 0.0, 
		eu_range[0][0]*ROT_CONV, eu_range[1][0]*ROT_CONV, eu_range[2][0]*ROT_CONV);
  if (eu_range[1][0]*ROT_CONV <= theta && eu_range[1][1]*ROT_CONV >= theta ) {
    for (phi=eu_range[2][0];phi<=eu_range[2][1];phi+=delta) {
      remap_eulers (&psi_tmp, &theta_tmp, &phi_tmp, psi_new, theta, phi*ROT_CONV, 
		    0.0, 0.0, 0.0);
      *(*eu_store+j+0)=psi_tmp;
      *(*eu_store+j+1)=theta_tmp;
      *(*eu_store+j+2)=phi_tmp; 
      j+=3;
    }
  }
}


/*=====================================================================*/
void eu_lattman (double eu_range[3][2], double delta, unsigned long *eu_count, float **eu_store) 
/* This subroutine generates a nondegenerate set of Euler angles       */ 
/* in the specified scan range. The angles are first sampled with      */
/* equal spacing in psi, theta, phi. Finally, the distribution is      */
/* sparsed to reduce the density at the poles.                         */ 
/* Reference:                                                          */
/* E.E. Lattman, Optimal Sampling of the Rotation Function, pp 179-185 */
/* The Molecular Replacement Method (Ed. M.G. Rossmann), Gordon and    */
/* Breach, Science Publishers Inc., New York, 1972.                    */
/* See also Gabb et al., J. Mol. Biol. 272:106-120, 1997.              */

{ double matrix1[3][3], matrix2[3][3];
  double psi, psi_tmp, theta, theta_tmp, phi, phi_tmp;
  char *degenerate;
  int psi_steps, theta_steps, phi_steps;
  unsigned long i, j, k;
  unsigned long full_count, curr_count;
  float *curr_eulers;
  double sparsing_cosine = cos(0.5*delta*ROT_CONV); 
  /* suggested Lattman sparsing tolerance: 0.5 delta */

  double deviation_cosine();
  void get_rot_matrix();

  psi_steps   = (int)ceil((eu_range[0][1]-eu_range[0][0])/delta);
  theta_steps = (int)ceil((eu_range[1][1]-eu_range[1][0])/delta);
  phi_steps   = (int)ceil((eu_range[2][1]-eu_range[2][0])/delta);

  if (psi_steps<1 || theta_steps<1 || phi_steps<1) {
    fprintf(stderr, "lib_eul> Error: negative number of Euler angle steps [e.c. 15115]\n"); 
    exit(15115);
  }

  full_count=psi_steps*theta_steps*phi_steps;
  printf("lib_eul> Lattman et al. Euler angles, initial number %d (delta = %f deg.)\n",(int)full_count,delta);

  /* allocate memory */
  curr_eulers = (float *) malloc(full_count * 3 * sizeof(float));
  if (curr_eulers == NULL) {   
    fprintf(stderr, "lib_eul> Error: Unable to satisfy memory allocation request for %d orientations [e.c. 15120]\n", (int)full_count);
    exit(15120);
  }

  degenerate = (char *) malloc(full_count * sizeof(char));
  if (degenerate == NULL) {
    fprintf(stderr, "lib_eul> Error: Unable to satisfy memory allocation request for %d orientations [e.c. 15130]\n", (int)full_count);
    exit(15130);
  }

  j=0;
  for (psi=0;psi<psi_steps*delta;psi+=delta) 
    for (theta=0;theta<theta_steps*delta;theta+=delta) 
      for (phi=0;phi<phi_steps*delta;phi+=delta)
	{	
	  remap_eulers (&psi_tmp, &theta_tmp, &phi_tmp, 
			(eu_range[0][0]+psi)*ROT_CONV, 
			(eu_range[1][0]+theta)*ROT_CONV, 
			(eu_range[2][0]+phi)*ROT_CONV, 
			0.0, 0.0, 0.0);
	  *(curr_eulers+j+0)=psi_tmp;
          *(curr_eulers+j+1)=theta_tmp;
          *(curr_eulers+j+2)=phi_tmp; 
	  j+=3;
	}
  
  for (j=0;j<full_count;j++) *(degenerate+j)=0;  
  
  k=0;
  printf("lib_eul> Searching for redundant orientations ");
  for (j=0;j<full_count;j++) if (*(degenerate+j)==0) {
    if (fmod(j+1,1000)==0.0)  printf(".");
    get_rot_matrix (matrix1,*(curr_eulers+j*3+0),*(curr_eulers+j*3+1),*(curr_eulers+j*3+2));
    for (i=j+1;i<full_count;i++) if (*(degenerate+i)==0) {
      get_rot_matrix (matrix2,*(curr_eulers+i*3+0),*(curr_eulers+i*3+1),*(curr_eulers+i*3+2));
      if (deviation_cosine (matrix1,matrix2) > sparsing_cosine) { /* degenerate orientations */ 
	*(degenerate+i)=1; 
	k++;
      }
    }
  }

  printf("\n");
  *eu_count=(full_count-k); 
  curr_count=(full_count-k);
  
  /* allocate memory */
  *eu_store = (float *) malloc(curr_count * 3 * sizeof(float));
  if (*eu_store == NULL) {
    fprintf(stderr, "lib_eul> Error: Unable to satisfy memory allocation request for %d orientations [e.c. 15210]\n", (int)curr_count);
    exit(15210);
  }

  printf("lib_eul> Lattman et al. Euler angles, final number %d (delta = %f deg.)\n",(int)curr_count,delta);
  i=0;
  for (j=0;j<full_count;j++)
    if (*(degenerate+j)!=1) {
      *(*eu_store+i+0)=*(curr_eulers+j*3+0);
      *(*eu_store+i+1)=*(curr_eulers+j*3+1);
      *(*eu_store+i+2)=*(curr_eulers+j*3+2);
      i+=3; 
    }
  free(curr_eulers);
  free(degenerate);
}



/*=====================================================================*/
static double deviation_cosine (double matrix1[3][3], double matrix2[3][3]) {
/* returns cosine of deviation angle between two rotations             */
/* as defined by Gabb et al., J. Mol. Biol. 272:106-120, 1997.         */

  double trace;
  int i,j;

  trace = 0;
  for (i=0;i<3;++i) for (j=0;j<3;++j) trace += matrix1[i][j]*matrix2[i][j];
  return (trace-1.0)/2.0;
}



/*=====================================================================*/
char similar_eulers(double a1,double a2,double a3,double b1,double b2,double b3) {
/* checks if two orientations are within 1 degree                          */
/* this can be slow if one of the orientations remains the same in a loop  */

  void get_rot_matrix ();
  double deviation_cosine ();
  
  double matrix1[3][3],matrix2[3][3];
  
  get_rot_matrix(matrix2,a1,a2,a3);
  get_rot_matrix(matrix1,b1,b2,b3);
  
  if (deviation_cosine (matrix2,matrix1) <= 0.9998477) return 0;
  else return 1;
}


