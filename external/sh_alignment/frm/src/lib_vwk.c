/*********************************************************************
*                           L I B _ V W K                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Pablo Chacon and Willy Wriggers, 2001-2003                     *
**********************************************************************
*                                                                    *
* Map and kernel creation and manipulation tools.                    * 
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/


#include "situs.h"

/* external functions */

extern void zero_vect (double *, unsigned long);
extern void do_vect (double **, unsigned long);
extern void cp_vect (double **, double **, unsigned long);

/* functions list  */

unsigned long gidz_cube (int, int, int, unsigned);
unsigned long gidz_general (int, int, int, unsigned, unsigned);
void create_padded_map (double **, unsigned *, unsigned *, unsigned *, double *, double *, double *, 
			unsigned long *, double *, unsigned, unsigned, unsigned, double, double, 
			double, double, double, double, unsigned *);
void interpolate_map (double **, unsigned *, unsigned *, unsigned *, double *, double *, double *, 
		      double, double, double, double *, unsigned, unsigned, unsigned, 
		      double, double, double, double, double, double);
void shrink_margin (double **, unsigned *, unsigned *, unsigned *, double *, double *, double *, 
		    unsigned long *, double *, unsigned, unsigned, unsigned, double, double, double, 
		    double, double, double);
void shrink_to_sigma_factor (double **, unsigned *, double *, unsigned, double, double);
double calc_total (double *, unsigned long);
double calc_average (double *, unsigned long);
double calc_sigma (double *, unsigned long);
double calc_norm (double *, unsigned long);
void print_map_info (double *, unsigned long);
void threshold (double *, unsigned long, double);
void normalize (double *, unsigned long, double);
void create_gaussian (double **, unsigned long *, unsigned *, double, double);
void create_identity (double **, unsigned long *, unsigned *);
void create_laplacian (double **, unsigned long *, unsigned *);

void relax_laplacian (double **, unsigned, unsigned, unsigned, unsigned *, double);

void convolve_kernel_inside (double **, double *, unsigned, unsigned, unsigned, double *, unsigned);  
void convolve_kernel_inside_fast (double **, double *, unsigned, unsigned, 
				  unsigned, double *, unsigned, double, unsigned *);
void convolve_kernel_inside_erode (double **, double *, unsigned, unsigned, unsigned, double *, unsigned);  
void convolve_kernel_outside (double **, unsigned *, unsigned *, unsigned *, double *, double *, 
			      double *, double *, unsigned, unsigned, unsigned, double, double, 
			      double, double, double, double, double *, unsigned);


/*====================================================================*/
unsigned long gidz_cube(int k, int j, int i, unsigned ext) {
/* generic index function z>y>x for cubic maps */ 
  return ext*ext*k+ext*j+i; 
}

/*====================================================================*/
unsigned long gidz_general(int k, int j, int i, unsigned ey, unsigned ex) {
/* generic index function z>y>x */ 
  return ex*ey*k+ex*j+i; 
}

/*====================================================================*/
void create_padded_map (double **outmap, unsigned *out_extx, unsigned *out_exty, unsigned *out_extz, 
			double *out_gridx, double *out_gridy, double *out_gridz, unsigned long *out_nvox,
			double *inmap, unsigned in_extx, unsigned in_exty, unsigned in_extz,
			double in_gridx, double in_gridy, double in_gridz,
			double widthx, double widthy, double widthz, 
			unsigned margin[3]) {
/* adds specified margin of zero density voxels to boundary */
/* outmap is allocated and new output map parameters are returned */
  
  unsigned indz, indx, indy;
  unsigned long index_old, index_new;
 
  *out_nvox=(in_extx+margin[0]*2)*(in_exty+margin[1]*2)*(in_extz+margin[2]*2);
  do_vect(outmap,*out_nvox);
  
  for (indz=0;indz<in_extz;indz++)
    for (indy=0;indy<in_exty;indy++)
      for (indx=0;indx<in_extx;indx++) {
	index_old=in_extx*in_exty*indz+in_extx*indy+indx;
	index_new=((in_extx+2*margin[0])*(in_exty+2*margin[1])*
		   (indz+margin[2])+(in_extx+2*margin[0])*
		   (indy+margin[1])+(indx+margin[0])); 
	*(*outmap+index_new) = *(inmap+index_old);
      }
  
  *out_extx = in_extx + margin[0]*2; /* may destroy in_extx if out_extx points to it */
  *out_exty = in_exty + margin[1]*2; /* may destroy in_exty if out_exty points to it */
  *out_extz = in_extz + margin[2]*2; /* may destroy in_extz if out_extz points to it */
  *out_gridx = in_gridx - margin[0] * widthx; /* may destroy in_gridx if out_gridx points to it */
  *out_gridy = in_gridy - margin[1] * widthy; /* may destroy in_gridy if out_gridy points to it */
  *out_gridz = in_gridz - margin[2] * widthz; /* may destroy in_gridz if out_gridz points to it */
  
  /* input variables are no longer used since they may have been overwritten */

  printf ("lib_vwk> Map size expanded from %d x %d x %d to %d x %d x %d by zero-padding.\n",
	  *out_extx-margin[0]*2, *out_exty-margin[1]*2, *out_extz-margin[2]*2,
	  *out_extx, *out_exty, *out_extz);
  printf ("lib_vwk> New map origin (coord of (1,1,1) voxel): (%.3f,%.3f,%.3f)\n",*out_gridx, *out_gridy, *out_gridz);
}


/*====================================================================*/
void interpolate_map (double **outmap, unsigned *out_extx, unsigned *out_exty, unsigned *out_extz, 
		      double *out_gridx, double *out_gridy, double *out_gridz,
		      double out_widthx, double out_widthy, double out_widthz,
		      double *inmap, unsigned in_extx, unsigned in_exty, unsigned in_extz,
		      double in_gridx, double in_gridy, double in_gridz,
		      double in_widthx, double in_widthy, double in_widthz) {
/* change voxel spacings of map and adjust map origin */
/* outmap is allocated and new output map parameters are returned */
  
  int isx_in, isy_in, isz_in;
  int isx_out, isy_out, isz_out;
  int iex_in, iey_in, iez_in;
  int iex_out, iey_out, iez_out;
  double deltax, deltay, deltaz;
  unsigned long out_nvox;
  unsigned indx, indy, indz;
  unsigned sx, sy, sz;
  double xpos, ypos, zpos, gx, gy, gz, a, b, c;
  int x0, y0, z0, x1, y1, z1;


  isx_in  = (int) (in_gridx/in_widthx+0.5);         /* start index rel. to origin, may be neg. */
  isy_in  = (int) (in_gridy/in_widthy+0.5);
  isz_in  = (int) (in_gridz/in_widthz+0.5);

  iex_in  = isx_in+in_extx-1;                       /* end index rel. to origin, may be neg. */
  iey_in  = isy_in+in_exty-1;  
  iez_in  = isz_in+in_extz-1;

  isx_out = ceil (isx_in*in_widthx/out_widthx);     /* start index rel. to origin, may be neg. */
  isy_out = ceil (isy_in*in_widthy/out_widthy);     /* real pos. not before start of inmap */
  isz_out = ceil (isz_in*in_widthz/out_widthz);
  
  iex_out = floor (iex_in*in_widthx/out_widthx);    /* end index rel. to origin, may be neg. */
  iey_out = floor (iey_in*in_widthy/out_widthy);    /* real pos. not after end of inmap */
  iez_out = floor (iez_in*in_widthz/out_widthz);
  
  sx=in_extx; sy=in_exty; sz=in_extz;               /* save input size parameters in case they get over written */
  
  *out_extx = iex_out-isx_out+1;                    /* may destroy in_extx if out_extx points to it */
  *out_exty = iey_out-isy_out+1;                    /* may destroy in_exty if out_exty points to it */
  *out_extz = iez_out-isz_out+1;                    /* may destroy in_extz if out_extz points to it */
  
  if (*out_extx<2 || *out_exty<2 || *out_extz<2) {
    fprintf(stderr, "lib_vwk> Error: interpolation output map size underflow [e.c. 17010]\n"); 
    exit(17010);
  }
  
  deltax = isx_out*out_widthx-isx_in*in_widthx;     /* shift of first voxel */
  deltay = isy_out*out_widthy-isy_in*in_widthy;
  deltaz = isz_out*out_widthz-isz_in*in_widthz;

  out_nvox = (*out_extx) * (*out_exty) * (*out_extz);

  do_vect(outmap,out_nvox);

  /* loop through outmap */

  for (indz=0;indz<*out_extz;indz++)
  for (indy=0;indy<*out_exty;indy++)
  for (indx=0;indx<*out_extx;indx++) {
    
    /* determine position of outmap voxel relative to start of inmap */
    
    xpos = deltax + indx * out_widthx;    
    ypos = deltay + indy * out_widthy;
    zpos = deltaz + indz * out_widthz; 

    /* compute position in index units */
    
    gx = (xpos/in_widthx);
    gy = (ypos/in_widthy);
    gz = (zpos/in_widthz);
 
    /* compute bounding voxels and linear distances */

    x0 = floor (gx); 
    y0 = floor (gy); 
    z0 = floor (gz); 
    x1 = ceil (gx);
    y1 = ceil (gy);
    z1 = ceil (gz);
    a = gx-x0;
    b = gy-y0;
    c = gz-z0;
     
    /* interpolate */
    
    *(*outmap+gidz_general(indz,indy,indx,*out_exty,*out_extx)) =
      a * b * c * *(inmap+gidz_general(z1,y1,x1,sy,sx)) +
      (1-a) * b * c * *(inmap+gidz_general(z1,y1,x0,sy,sx)) +
      a * (1-b) * c * *(inmap+gidz_general(z1,y0,x1,sy,sx)) +
      a * b * (1-c) * *(inmap+gidz_general(z0,y1,x1,sy,sx)) +
      a * (1-b) * (1-c) * *(inmap+gidz_general(z0,y0,x1,sy,sx)) +
      (1-a) * b * (1-c) * *(inmap+gidz_general(z0,y1,x0,sy,sx)) +      
      (1-a) * (1-b) * c * *(inmap+gidz_general(z1,y0,x0,sy,sx)) +
      (1-a) * (1-b) * (1-c) * *(inmap+gidz_general(z0,y0,x0,sy,sx));
  }

  *out_gridx = in_gridx + deltax; /* may destroy in_gridx if out_gridx points to it */
  *out_gridy = in_gridy + deltay; /* may destroy in_gridy if out_gridy points to it */
  *out_gridz = in_gridz + deltaz; /* may destroy in_gridz if out_gridz points to it */

  /* don't use in_gridx, in_gridy, in_gridz below this line  since they may have been overwritten */

  printf ("lib_vwk> Map interpolated from %d x %d x %d to %d x %d x %d.\n",
	  sx, sy, sz, *out_extx, *out_exty, *out_extz);
  printf ("lib_vwk> Voxel spacings changed from (%.3f,%.3f,%.3f) to (%.3f,%.3f,%.3f).\n",
	  in_widthx, in_widthy, in_widthz, out_widthx, out_widthy, out_widthz);
  printf ("lib_vwk> New map origin (coord of (1,1,1) voxel): (%.3f,%.3f,%.3f)\n",*out_gridx, *out_gridy, *out_gridz);
}


/*====================================================================*/
void shrink_margin (double **outmap, unsigned *out_extx, unsigned *out_exty, unsigned *out_extz, 
		    double *out_gridx, double *out_gridy, double *out_gridz, unsigned long *out_nvox,
		    double *inmap, unsigned in_extx, unsigned in_exty, unsigned in_extz,
		    double in_gridx, double in_gridy, double in_gridz,
		    double widthx, double widthy, double widthz) {
/* shrink inmap and replace with smaller *outmap that has odd intervals */
/* outmap is allocated and new output map parameters are returned */
  
  unsigned m,p,q,sx,sy,sz;
  unsigned minx,miny,minz,maxx,maxy,maxz;
  int margin[6];
  unsigned indz,indx,indy;
  unsigned long index_old, index_new;

  maxx=0; maxy=0; maxz=0;
  minx=in_extx-1; miny=in_exty-1; minz=in_extz-1;

  for (q=0;q<in_extz;q++) 
    for (m=0;m<in_exty;m++) 
      for (p=0;p<in_extx;p++) 
	if ( *(inmap+gidz_general(q,m,p,in_exty,in_extx)) > 0) {
	  if (p<=minx) minx=p;
	  if (p>=maxx) maxx=p;
	  if (m<=miny) miny=m;
	  if (m>=maxy) maxy=m;
	  if (q<=minz) minz=q;
	  if (q>=maxz) maxz=q;
	}
  
  if (maxx<minx) {
    fprintf(stderr, "lib_vwk> Error: No positive density found [e.c. 17020]\n"); 
    exit(17020);
  }
  
  margin[1]=in_extx-maxx-1;
  margin[3]=in_exty-maxy-1;
  margin[5]=in_extz-maxz-1;
  margin[0]=minx;
  margin[2]=miny;
  margin[4]=minz;

  /* compute new grid size */
  p=in_extx-(margin[0]+margin[1]);
  m=in_exty-(margin[2]+margin[3]);
  q=in_extz-(margin[4]+margin[5]);

  /* number of map intervals to be odd, if necessary they are increase */
  if (2*(p/2)==p) { p++; if (margin[0]>0)  margin[0]--;  }
  if (2*(m/2)==m) { m++; if (margin[2]>0)  margin[2]--;  }
  if (2*(q/2)==q) { q++; if (margin[4]>0)  margin[4]--;  }

  
  *out_nvox=p*m*q;
  do_vect(outmap,*out_nvox);
  
  for (indz=margin[4];indz<in_extz-margin[5];indz++) 
    for (indy=margin[2];indy<in_exty-margin[3];indy++)
      for (indx=margin[0];indx<in_extx-margin[1];indx++) {
	index_new=p*m*(indz-margin[4])+p*(indy-margin[2])+(indx-margin[0]);
	index_old=in_extx*in_exty*indz+in_extx*indy+indx;  
	*(*outmap+index_new)=*(inmap+index_old);
      }
   
  sx = in_extx;
  sy = in_exty;
  sz = in_extz;
  *out_extx = p;                              /* may destroy in_extx if out_extx points to it */
  *out_exty = m;                              /* may destroy in_exty if out_exty points to it */
  *out_extz = q;                              /* may destroy in_extz if out_extz points to it */
  *out_gridx = in_gridx + margin[0] * widthx; /* may destroy in_gridx if out_gridx points to it */
  *out_gridy = in_gridy + margin[2] * widthy; /* may destroy in_gridy if out_gridy points to it */
  *out_gridz = in_gridz + margin[4] * widthz; /* may destroy in_gridz if out_gridz points to it */
  
  /* input variables are no longer used since they may have been overwritten */

  printf ("lib_vwk> Map size reduced from %d x %d x %d to %d x %d x %d.\n",
	  sx, sy, sz, *out_extx, *out_exty, *out_extz);
  printf ("lib_vwk> New map origin (coord of (1,1,1) voxel): (%.3f,%.3f,%.3f)\n",*out_gridx, *out_gridy, *out_gridz);
}




/*====================================================================*/
double calc_total (double *phi, unsigned long nvox) {
/* returns total density in map */
  unsigned m; double total_density;
  total_density=0;
  for (m=0;m<nvox;m++) total_density+= (*(phi+m));  
  return total_density;
}


/*====================================================================*/
double calc_average (double *phi, unsigned long nvox) { 
/* returns average density */
 
  double total_density; unsigned m,p;
  total_density=0; p=0;
  for (m=0;m<nvox;m++) { 
    total_density+= (*(phi+m)); p++; 
  } 
  return total_density/((double)p);
}


/*====================================================================*/
double calc_sigma (double *phi, unsigned long nvox) {
/* returns sigma */

  unsigned m,p; double ave,varsum;

  ave=calc_average(phi,nvox);
  varsum=0; p=0;
  for (m=0;m<nvox;m++) { 
    varsum+=(*(phi+m)-ave)*(*(phi+m)-ave); 
    p++;
  }
  return sqrt((varsum)/((double)p));
}

/*====================================================================*/
double calc_norm (double *phi, unsigned long nvox) {
/* returns norm */

  unsigned m,p; double varsum;

  varsum=0; p=0;
  for (m=0;m<nvox;m++) { 
    varsum+=(*(phi+m))*(*(phi+m)); 
    p++;
  }
  return sqrt((varsum)/((double)p));
}


/*====================================================================*/
void print_map_info (double *phi, unsigned long nvox) {
/* outputs info about map density distribution */

  double maxdens,mindens,sig,ave;
  unsigned p,m;

  maxdens=-1e20;
  mindens=1e20;
  ave=0;
  p=0;
  for (m=0;m<nvox;m++) {
    if (*(phi+m)>maxdens) maxdens=*(phi+m);
    if (*(phi+m)<mindens) mindens=*(phi+m);
    if (*(phi+m)> 0) { 
	p++;
      ave+=*(phi+m);
    }	
  }
  ave/=(double)p;
  sig=0;
  for (m=0;m<nvox;m++) if (*(phi+m)>0) sig+=(*(phi+m)-ave)*(*(phi+m)-ave);
  sig/=(double)p; 
  printf("lib_vwk> Map density info: max %f, min %f, ave %f, sig %f.\n",maxdens,mindens,ave,sqrt(sig));
}


/*====================================================================*/
void threshold (double *phi, unsigned long nvox, double limit) {
/* set densities below limit to zero */

  unsigned long m;

  if (limit < 0) {
    fprintf(stderr, "lib_vwk> Error: Threshold value negative [e.c. 17120]\n"); 
    exit(17120);
  }
  for (m=0;m<nvox;m++) if (*(phi+m)<limit) *(phi+m)=0.0;
  printf("lib_vwk> Setting density values below %f to zero.\n", limit);
}


/*====================================================================*/
void normalize (double *phi, unsigned long nvox, double factor) {
/* divide density values by factor */

  unsigned i;
  
  if (factor==0) {
    fprintf(stderr, "lib_vwk> Error: Normalization by zero [e.c. 17130]\n"); 
    exit(17130);
  }
  for (i=0;i<nvox;i++) *(phi+i)/=factor;
}


/*====================================================================*/
void create_gaussian (double **phi, unsigned long *nvox, unsigned *ext, 
		      double sigmap, double sigma_factor)  {
/* allocates and generates truncated Gaussian 3D kernel with */
/* sigma1D = sigmap, within sigma_factor*sigmap  */

  int  exth, indx, indy, indz;
  double dsqu;
  double mscale;
  unsigned long count;  
  double bvalue, cvalue;
       
  /* truncate at sigma_factor * sigma1D */
  exth = (int) ceil(sigma_factor*sigmap);
  *ext = 2 * exth - 1;  
  *nvox = *ext * *ext * *ext;
    
  printf ("lib_vwk> Generating Gaussian kernel with %d^3 = %d voxels.\n", (int)*ext, (int)*nvox);
  do_vect(phi, *nvox);   

  /* write Gaussian within sigma_factor * sigma-1D to map */
  bvalue = -1 / (2.0*sigmap*sigmap);
  cvalue = sigma_factor*sigma_factor*sigmap*sigmap; 
    
  mscale=0;
  for (indz=0;indz<*ext;indz++)
    for (indy=0;indy<*ext;indy++)
      for (indx=0;indx<*ext;indx++) {
	dsqu = (indx-exth+1)*(indx-exth+1)+
	  (indy-exth+1)*(indy-exth+1)+
	  (indz-exth+1)*(indz-exth+1);    
	if (dsqu <= cvalue)
	  *(*phi+gidz_cube(indz,indy,indx,*ext)) = exp (dsqu * bvalue);
	mscale+=*(*phi+gidz_cube(indz,indy,indx,*ext));
      }
  for(count=0;count<*nvox;count++) *(*phi+count) /=mscale; 
}


/*====================================================================*/
void shrink_to_sigma_factor (double **outmap, unsigned *out_ext, double *inmap, 
			     unsigned in_ext, double sigmap, double sigma_factor)  {
/* shrink input kernel inmap and replace with smaller *outmap within sigma_factor*sigmap  */
/* that has odd intervals; outmap is allocated and new output map parameters are returned */

  int exth, indx, indy, indz;
  unsigned long nvox, index_old, index_new;
  unsigned margin;
  double cvalue,dsqu;

  /* truncate at sigma_factor * sigma1D */
  exth = (int) ceil(sigma_factor*sigmap);
  *out_ext = 2 * exth - 1;  
  if (*out_ext > in_ext || 2*(in_ext/2)==in_ext) {
    fprintf(stderr, "lib_vwk> Error: Input and output kernels not compatible [e.c. 17030]\n"); 
    exit(17030);
  }
  nvox = *out_ext * *out_ext * *out_ext;
  cvalue = sigma_factor*sigma_factor*sigmap*sigmap; 

  printf ("lib_vwk> Generating kernel with %d^3 = %d voxels.\n",(int)*out_ext,(int)nvox);
  do_vect(outmap,nvox);   
  
  margin = (in_ext- *out_ext)/2;
  for (indz=margin;indz<in_ext-margin;indz++) 
    for (indy=margin;indy<in_ext-margin;indy++)
      for (indx=margin;indx<in_ext-margin;indx++) {


	index_new=*out_ext * *out_ext * (indz-margin) + *out_ext * (indy-margin)+(indx-margin);
	index_old=in_ext*in_ext*indz+in_ext*indy+indx;  
	*(*outmap+index_new)=*(inmap+index_old);
      }

  /* make the kernel spherical */
  for (indz=0;indz<*out_ext;indz++) 
    for (indy=0;indy<*out_ext;indy++) 
      for (indx=0;indx<*out_ext;indx++) {
        index_new=(*out_ext)*(*out_ext)*indz+(*out_ext)*indy+indx;  
	dsqu = (indx-exth+1)*(indx-exth+1)+
	  (indy-exth+1)*(indy-exth+1)+
	  (indz-exth+1)*(indz-exth+1);    
	if (dsqu > cvalue) *(*outmap+index_new)=0.0;
      } 
     

}
/*====================================================================*/
void create_identity (double **phi, unsigned long *nvox, unsigned *ext)  {
/* creates identity 1x1x1 kernel */

  *ext = 1; 
  *nvox = 1; 
  do_vect(phi, *nvox);
  *(*phi+0)=1;
}


/*====================================================================*/
void create_laplacian (double **phi, unsigned long *nvox, unsigned *ext)  {
/* creates Laplacian 3x3x3 kernel */

  unsigned indx, indy, indz;
  double lap_discrete[3][3][3] = 
  { { { 0,  0,      0},  {      0,  1/12.,      0},  { 0,      0, 0} },
    { { 0,  1/12.,  0},  {  1/12., -6/12., 1/12.0},  { 0,  1/12., 0} },
    { { 0,  0,      0},  {      0,  1/12.,      0},  { 0,      0, 0} } };
  
  *ext = 3; 
  *nvox = *ext * *ext * *ext;

  printf ("lib_vwk> Generating Laplacian kernel with %d^3 = %d voxels.\n", (int)*ext, (int)*nvox);
  do_vect(phi, *nvox);

  for (indz=0;indz<*ext;indz++)
    for (indy=0;indy<*ext;indy++)
      for (indx=0;indx<*ext;indx++) 
	*(*phi+gidz_cube(indz,indy,indx,*ext)) = lap_discrete[indx][indy][indz];
}


/*====================================================================*/
void relax_laplacian (double **phi, unsigned extx, unsigned exty, unsigned extz, unsigned ignored[3], double radius)  {
/* relaxes shell of width radius, just outside thresholded area, by the Poisson equation */
/* assumes radius is smaller than previously applied zero padding */
/* assumes that phi has been thresholded and contains some non-zero densities */

  double average[27]={0.,0.,0.,0.,1/6.,0.,0.,0.,0.,0.,1/6.,0.,1/6.,0.,1/6.,0.,1/6.,0.,0.,0.,0.,0.,1/6.,0.,0.,0.,0.};
  double  *nextphi, diff, norm, epsilon;
  unsigned long nvox, indv, indw, threscount, maskcount;
  char *mask; 
  unsigned indx, indy, indz;
  int indx2, indy2, indz2;
  int margx, margy, margz, margin;

  margx = (int)(ignored[0]+radius); 
  margy = (int)(ignored[1]+radius); 
  margz = (int)(ignored[2]+radius);
  margin = (int) ceil(radius);
  printf("lib_vwk> Relaxing %d voxel thick shell about thresholded density... \n",margin);

  /* allocate phi mask */
  nvox = extx * exty * extz;
  mask = (char *) malloc(nvox * sizeof(char)); 
  if (mask == NULL) {
    fprintf(stderr, "lib_vwk> Error: Unable to satisfy memory allocation request [e.c. 17130]\n"); 
    exit(17130);
  }  
  for (indv=0;indv<nvox;indv++) *(mask+indv)=1;

  /* assign phi mask value based on distance to thresholded map */
  for (indz=margz;indz<extz-margz;indz++) 
    for (indy=margy;indy<exty-margy;indy++) 
      for (indx=margx;indx<extx-margx;indx++) {
	indv = gidz_general(indz,indy,indx,exty,extx);
	if (*(*phi+indv)!=0) 
	  for (indz2=-margin;indz2<=margin;indz2++)
	    for (indy2=-margin;indy2<=margin;indy2++) 
	      for (indx2=-margin;indx2<=margin;indx2++) {
		indw = gidz_general(indz+indz2,indy+indy2,indx+indx2,exty,extx);
		if (*(*phi+indw)==0 && indz2*indz2+indy2*indy2+indx2*indx2<radius*radius) *(mask+indw)=0;
	    }
      }
  
  /* compute norm */
  maskcount = 0; threscount = 0; norm = 0;
  for (indv=0;indv<nvox;indv++) {
    if (*(*phi+indv)!=0) {
      ++threscount;
      norm += *(*phi+indv);
    } else if (*(mask+indv)==0) ++maskcount;
  }
  norm /= (double)threscount; /* average density for thresholded volume, assuming threscount>0 */
  norm *= maskcount;          /* density total one would get if mask=0 was filled with average */ 
    
  /* iterate on original lattice, no focusing */
  do_vect(&nextphi,nvox);      
  do {  
    convolve_kernel_inside_fast (&nextphi,*phi,extx,exty,extz,average,3,1.0,ignored);
    diff=0;
    for (indz=ignored[2];indz<extz-ignored[2];indz++) 
      for (indy=ignored[1];indy<exty-ignored[1];indy++) 
	for (indx=ignored[0];indx<extx-ignored[0];indx++) {
	  indv = gidz_general(indz,indy,indx,exty,extx);
	  if (*(mask+indv)==0) {
	    diff += fabs(*(nextphi+indv) - *(*phi+indv));
	    *(*phi+indv) = *(nextphi+indv);
	  }
	}
  } while (diff > 1E-8 * norm); 
  free(nextphi);
  free(mask);
}


/*====================================================================*/
void convolve_kernel_inside (double **outmap, double *inmap, 
			     unsigned in_extx, unsigned in_exty, unsigned in_extz, 
			     double *kernel, unsigned kernel_size)  {
/* convolves inmap with kernel inside, and writes to same-size outmap        */
/* sufficient memory for *outmap must already be allocated                   */
/* allows *outmap and inmap to point to same array - not speed optimized     */

  double dval;
  unsigned long out_nvox;
  unsigned indx,indy,indz;
  int indx2,indy2,indz2,margin;
  double *tmpmap;
 
  if (kernel_size < 1 || 2*((kernel_size+1)/2)-kernel_size-1 != 0) {
    fprintf(stderr, "lib_vwk> Error: Kernel size must be a positive odd number [e.c. 17140]\n"); 
    exit(17140);
  }
  
  margin=(kernel_size-1)/2;
  out_nvox = in_extx * in_exty * in_extz;

  do_vect(&tmpmap,out_nvox);
  cp_vect(&tmpmap,&inmap,out_nvox); /* save inmap in case it gets overwritten */
  zero_vect(*outmap,out_nvox);

  for (indz=margin;indz<in_extz-margin;indz++)
    for (indy=margin;indy<in_exty-margin;indy++)
      for (indx=margin;indx<in_extx-margin;indx++) {
	dval = *(tmpmap+gidz_general(indz,indy,indx,in_exty,in_extx)); 
	if (dval!=0) for (indz2=-margin;indz2<=margin;indz2++)
	  for (indy2=-margin;indy2<=margin;indy2++)
	    for (indx2=-margin;indx2<=margin;indx2++) { 	      
	      *(*outmap+gidz_general(indz+indz2,indy+indy2,indx+indx2,in_exty,in_extx))
		+= *(kernel+gidz_cube(indz2+margin,indy2+margin,indx2+margin,kernel_size)) * dval;
	    }                 
      }       
  free(tmpmap);
}


/*====================================================================*/
void convolve_kernel_inside_fast (double **outmap, double *inmap, 
				  unsigned in_extx, unsigned in_exty, unsigned in_extz, 
				  double *kernel, unsigned kernel_size, double normfac, unsigned ignored[3])
/* convolves inmap with kernel inside "ignored" boundary, and writes to same-size outmap     */
/* memory for *outmap must be allocated, and *outmap and inmap must point to different array */
/* speed-optimized version for colores, no memory allocation, no error checks, norm variable */
{
  double dval;
  unsigned long out_nvox;
  unsigned indx,indy,indz;
  int indx2,indy2,indz2,margin,marginx,marginy,marginz;
  
  margin = (kernel_size-1)/2;
  marginx= margin+ignored[0];
  marginy= margin+ignored[1];
  marginz= margin+ignored[2];
  out_nvox = in_extx * in_exty * in_extz;
 
  zero_vect(*outmap,out_nvox);
  
  for (indz=marginz;indz<in_extz-marginz;indz++)
    for (indy=marginy;indy<in_exty-marginy;indy++)
      for (indx=marginx;indx<in_extx-marginx;indx++) {
	dval = (*(inmap+gidz_general(indz,indy,indx,in_exty,in_extx))) / normfac;
	if (dval!=0) for (indz2=-margin;indz2<=margin;indz2++)
	  for (indy2=-margin;indy2<=margin;indy2++)
	    for (indx2=-margin;indx2<=margin;indx2++) {             
	      *(*outmap+gidz_general(indz+indz2,indy+indy2,indx+indx2,in_exty,in_extx))
		+= *(kernel+gidz_cube(indz2+margin,indy2+margin,indx2+margin,kernel_size)) * dval;
	    }                 
      }
}

/*====================================================================*/
void convolve_kernel_inside_erode (double **outmap, double *inmap, 
				   unsigned in_extx, unsigned in_exty, unsigned in_extz, 
				   double *kernel, unsigned kernel_size)  {
/* convolves inmap with kernel inside, and writes to same-size outmap        */
/* filters out convolutions where kernel hits zero density in *inmap         */ 
/* to avoid cutoff edge effects when using finite difference kernel          */
/* sufficient memory for *outmap must already be allocated                   */
/* allows *outmap and inmap to point to same array - not speed optimized     */
  
  unsigned long out_nvox;
  unsigned indx,indy,indz;
  int indx2,indy2,indz2,margin;
  double *tmpmap;
  double dval, dval2;
  unsigned skip;
  
  if (kernel_size < 1 || 2*((kernel_size+1)/2)-kernel_size-1 != 0) {
    fprintf(stderr, "lib_vwk> Error: Kernel size must be a positive odd number [e.c. 17150]\n"); 
    exit(17150);
  }
  
  margin=(kernel_size-1)/2;
  out_nvox = in_extx * in_exty * in_extz;

  do_vect(&tmpmap,out_nvox);
  cp_vect(&tmpmap,&inmap,out_nvox); /* save inmap in case it gets overwritten */
  zero_vect(*outmap,out_nvox);

  for (indz=margin;indz<in_extz-margin;indz++)
    for (indy=margin;indy<in_exty-margin;indy++)
      for (indx=margin;indx<in_extx-margin;indx++) {
        skip = 0;
        for (indz2=-margin;skip==0&&indz2<=margin;indz2++)
          for (indy2=-margin;skip==0&&indy2<=margin;indy2++)
            for (indx2=-margin;skip==0&&indx2<=margin;indx2++) { /* check if kernel hits zero density */
              dval  = *(tmpmap+gidz_general(indz+indz2,indy+indy2,indx+indx2,in_exty,in_extx));
              dval2 = *(kernel+gidz_cube(margin-indz2,margin-indy2,margin-indx2,kernel_size));
              if (dval == 0 && dval2 != 0) skip = 1;
            }
        if (skip == 0) {
          for (indz2=-margin;indz2<=margin;indz2++)
            for (indy2=-margin;indy2<=margin;indy2++)
              for (indx2=-margin;indx2<=margin;indx2++) {
                dval  = *(tmpmap+gidz_general(indz+indz2,indy+indy2,indx+indx2,in_exty,in_extx));
                dval2 = *(kernel+gidz_cube(margin-indz2,margin-indy2,margin-indx2,kernel_size));
                *(*outmap+gidz_general(indz,indy,indx,in_exty,in_extx)) += dval * dval2;
              }
        }           
      }       
  free(tmpmap);
}


/*====================================================================*/
void convolve_kernel_outside (double **outmap, unsigned *out_extx, unsigned *out_exty, unsigned *out_extz, 
			      double *out_gridx, double *out_gridy, double *out_gridz,
			      double *inmap, unsigned in_extx, unsigned in_exty, unsigned in_extz, 
			      double in_gridx, double in_gridy, double in_gridz,
			      double widthx, double widthy, double widthz, 
			      double *kernel, unsigned kernel_size)  {
/* convolves inmap with kernel overlapping the border, and writes to larger-size outmap */
/* sufficient memory for *outmap must already be allocated                              */
/* allows *outmap and inmap to point to same array - not speed optimized                */

  double dval;
  unsigned long out_nvox;
  unsigned long in_nvox;
  unsigned indx,indy,indz,tmp_extx,tmp_exty;
  int indx2,indy2,indz2,margin;
  double *tmpmap;
  
  if (kernel_size < 1 || 2*((kernel_size+1)/2)-kernel_size-1 != 0) {
    fprintf(stderr, "lib_vwk> Error: Kernel size %d must be a positive odd number [e.c. 17160]\n", kernel_size); 
    exit(17160);
  }
  
  margin=(kernel_size-1)/2;
  
  tmp_extx = in_extx;
  tmp_exty = in_exty;
  *out_extx = kernel_size - 1 + in_extx; /* may overwrite in_extx */
  *out_exty = kernel_size - 1 + in_exty; /* may overwrite in_exty */
  *out_extz = kernel_size - 1 + in_extz; /* may overwrite in_extz */

  in_nvox = in_extx * in_exty * in_extz;
  out_nvox = (*out_extx) * (*out_exty) * (*out_extz);
  
  do_vect(&tmpmap,in_nvox);
  cp_vect(&tmpmap,&inmap,in_nvox); /* save inmap in case it gets overwritten */
  zero_vect(*outmap,out_nvox);
  
  for (indz=margin;indz<(*out_extz)-margin;indz++)
    for (indy=margin;indy<(*out_exty)-margin;indy++)
      for (indx=margin;indx<(*out_extx)-margin;indx++) {
	dval = *(tmpmap+gidz_general(indz-margin,indy-margin,indx-margin,tmp_exty,tmp_extx)); 
	if (dval!=0) for (indz2=-margin;indz2<=margin;indz2++)
	  for (indy2=-margin;indy2<=margin;indy2++)
	    for (indx2=-margin;indx2<=margin;indx2++) {             
	      *(*outmap+gidz_general(indz+indz2,indy+indy2,indx+indx2,(*out_exty),(*out_extx)))
		+= *(kernel+gidz_cube(indz2+margin,indy2+margin,indx2+margin,kernel_size)) * dval;	
	    }                 
      }     
  free(tmpmap);
  *out_gridx = in_gridx - widthx * margin;
  *out_gridy = in_gridy - widthy * margin;
  *out_gridz = in_gridz - widthz * margin;
}


















































