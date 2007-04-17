/* ---------------------------------------------------------------------------

Filename		: WaveletTools.c 

Author			: Daniel Sage

Organization	: EPFL, Biomedical Imaging Group

Date	       	: 21/12/1998  

Purpose			: 

	3D signal splitting using a symmetrical filter 
	3D signal merging using a symmetrical filter 

	The wavelet filters should be symmetrical are loaded before running. 
	The boundary conditions can be choosen: Mirror or Peridic.
	
	There are two external functions:
		WaveletSplit_3D()   x[nx] -> y1[nx/2] & y2[nx/2]   
		WaveletMerge_3D()   y1[ny] & y2[ny] -> x[2*ny]

Convention		:

 	v[nv] : lowpass filter, array of coefficients 
 	w[nw] : highpass filter, array of coefficients 

---------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------- */
/* Declaration of extern procedures                                          */
/* ------------------------------------------------------------------------- */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "../configs.h"
#include "../headers/messagedisplay.h"
#include "../headers/minmax.h"
#include "../headers/getputd.h"
#include "../headers/getput.h"

#include "../headers/wavelettools.h"
#include "../headers/waveletfilters.h"
#include "../headers/waveletfiltersfract.h"

/* ------------------------------------------------------------------------- */
/* Declaration of static procedures                                          */
/* ------------------------------------------------------------------------- */
static int WaveletSplit_1D(	
					double x[], long nx,
					double y1[], double y2[], long *ny,
					double v[], long nv, 
					double w[], long nw,
					short BoundaryConditions);

static int WaveletMerge_1D(
					double y1[],double y2[],long ny,
					double x[],long *nx,
					double v[],long nv,
					double w[],long nw,
					short BoundaryConditions);

static void WaveletSplit_1DMirror(	
					double x[], long nx, 
					double y1[],double y2[], long *ny,
					double v[], long nv, 
					double w[], long nw);
					
static void WaveletSplit_1DPeriodic(
					double x[], long nx,
					double y1[],double y2[], long *ny,
					double v[], long nv, 
					double w[], long nw);

static void WaveletMerge_1DMirror(	
					double y1[],double y2[],long ny,
					double x[],long *nx,
					double v[],long nv,
					double w[],long nw);
				
static void WaveletMerge_1DPeriodic(	
					double y1[],double y2[],long ny,
					double x[],long *nx,
					double v[],long nv,
					double w[],long nw);

/* ----------------------------------------------------------------------------

	Function:	WaveletSplit_3D
	
	Purpose:	Two channel filtering and down-sampling for 3D real signals. 
				Input[nz][ny][nx] -> output[nz][ny][nx]
 				The output is build with four sectors in the classical way to 
 				show the wavelet transform

	Author:		Daniel Sage, EPFL, Biomedical Imaging Group
					
---------------------------------------------------------------------------- */
extern int	WaveletSplit_3D( 
					double	*Input, 
					double	*Output,
					long	nx, long ny, long nz,
					short	Filter,
					short	Order,
					double	Alpha,
					short	BoundaryConditions,
					int		*Status)
{
long	kx, ky, kz, nmax, nout, k, i, j, index, jindex, nx2, ny2, nz2, nxny,
        idxbase;
double	*FilterLowpass,	*FilterHighpass; 
long	NbLowpass[1], NbHighpass[1]; 
double	*Low, *High, *LowLow, *LowHigh, *HighLow, *HighHigh;
double	*x, *y1, *y2;

	nmax = MaxLong(MaxLong(nx, ny), nz);
	nx2  = MaxLong(nx / 2L, 1L);
	ny2  = MaxLong(ny / 2L, 1L);
	nz2  = MaxLong(nz / 2L, 1L);
	nxny = nx*ny;
		
	/* Get the size of the filter */
  	if ( Filter != 5) {
		if ( WaveletFiltersGetSize(Filter, Order, NbLowpass, NbHighpass) == ERROR) {
			MessageDisplay("ERROR - Impossible to get the size of the filter"); 
			*Status = ERROR;
			return( ERROR);
		}
	}
	else {
		if ( WaveletFiltersGetSize_Fract(Alpha, NbLowpass, NbHighpass) == ERROR) {
			MessageDisplay("ERROR - Impossible to get the size of the filter"); 
			*Status = ERROR;
			return( ERROR);
		}
	}

	/* Allocate the memory for the highpass filter */
  	AllocateLineDouble( &FilterHighpass, *NbHighpass, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		*Status = ERROR;
		return( ERROR);
  	}
  	
	/* Allocate the memory for the lowpass filter */
  	AllocateLineDouble( &FilterLowpass, *NbLowpass, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation 22");
		FreeLineDouble(&FilterHighpass);
		*Status = ERROR;
		return( ERROR);
  	}
  	
  	/* Get the coefficient of the filters */ 
  	if ( Filter != 5) {
		if ( WaveletFiltersGetCoef( Filter, Order, FilterLowpass, FilterHighpass) == ERROR) {
			MessageDisplay("ERROR - Impossible to get the coefficient of the filter"); 
			*Status = ERROR;
			return( ERROR);
		}
	}
	else {
		if ( WaveletFiltersGetCoef_Fract( Alpha, FilterLowpass, FilterHighpass) == ERROR) {
			MessageDisplay("ERROR - Impossible to get the coefficient of the filter"); 
			*Status = ERROR;
			return( ERROR);
		}
	}
	
	/* Allocate the temporary vector */
  	AllocateVolumeDouble( &Low, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		return( ERROR);
  	}
  	
  	AllocateVolumeDouble( &High, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		return( ERROR);
  	}

  	AllocateVolumeDouble( &LowLow, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		return( ERROR);
  	}

  	AllocateVolumeDouble( &LowHigh, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		return( ERROR);
  	}

  	AllocateVolumeDouble( &HighLow, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		return( ERROR);
  	}

  	AllocateVolumeDouble( &HighHigh, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		FreeVolumeDouble(&HighLow);
		return( ERROR);
  	}

  	AllocateLineDouble( &x, nmax, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		FreeVolumeDouble(&HighLow);
		FreeVolumeDouble(&HighHigh);
		return( ERROR);
  	}

  	AllocateLineDouble( &y1, nmax, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		FreeVolumeDouble(&HighLow);
		FreeVolumeDouble(&HighHigh);
		FreeLineDouble(&x);
		return( ERROR);
  	}
  	
  	AllocateLineDouble( &y2, nmax, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		FreeVolumeDouble(&HighLow);
		FreeVolumeDouble(&HighHigh);
		FreeLineDouble(&x);
		FreeLineDouble(&y1);
		return( ERROR);
  	}

	/******************************************************************************/
	/* 1D                                                                         */
	/******************************************************************************/
	      
	if (ny == 1 && nz==1) {
		for ( i = 0L; i < nx; i++) x[i] = (double)Input[i];
		WaveletSplit_1D(	x, nx, y1, y2, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		for ( i = 0L; i < nx/2; i++)  Output[i] 	= (double)y1[i];
		for ( i = 0L; i < nx2; i++)  Output[i+nx2] 	= (double)y2[i]; 
	}
	
	/******************************************************************************/
	/* 2D                                                                         */
	/******************************************************************************/
	
	if (ny > 1) {
	  	
	/* y-processing */
	for (kz=0L;kz<nz;kz++)
	for (kx=0L;kx<nx;kx++) {
		GetyDoubleToDouble(Input,nx, ny, nz, kx, 0L, kz, x, ny);
		WaveletSplit_1D(	x, ny, y1, y2, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		PutyDoubleToDouble(Low, nx, ny, nz, kx, 0L, kz, y1, ny); 
		PutyDoubleToDouble(High,nx, ny, nz, kx, 0L, kz, y2, ny); 
	}
	
	
	/* x-processing for lowpass part */
	for (kz=0L;kz<nz;kz++)
	for (ky=0L;ky<ny/2L;ky++) {
		GetxDoubleToDouble(Low, nx, ny, nz, 0L, ky, kz, x, nx); 
		WaveletSplit_1D(	x, nx, y1, y2, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		PutxDoubleToDouble(LowLow, nx, ny, nz, 0L, ky, kz, y1, nx); 
		PutxDoubleToDouble(LowHigh,nx, ny, nz, 0L, ky, kz, y2, nx); 
	}

	/* x-processing for highpass part */
	for (kz=0L;kz<nz;kz++)
	for (ky=0L;ky<ny/2L;ky++) {
		GetxDoubleToDouble(High, nx, ny, nz, 0L, ky, kz, x, nx); 
		WaveletSplit_1D(	x, nx, y1, y2, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		PutxDoubleToDouble(HighLow, nx, ny, nz, 0L, ky, kz, y1, nx); 
		PutxDoubleToDouble(HighHigh,nx, ny, nz, 0L, ky, kz, y2, nx); 
	}
	
	/* Build the output */
	for (kz=0L;kz<nz;kz++)
	for ( i = 0L; i < nx2; i++)
	for ( j = 0L; j < ny2; j++)  {
		index  = kz*nxny+j*nx+i;
		jindex = kz*nxny+(j+ny2)*nx+i;
		Output[index] 		= LowLow[index];
		Output[index+nx2] 	= LowHigh[index];
		Output[jindex] 		= HighLow[index];
		Output[jindex+nx2] 	= HighHigh[index] ;
	}
	
	}	
	
	/******************************************************************************/
	/* 3D                                                                         */
	/******************************************************************************/
	if (nz > 1)
   	   for (i=0L; i<ny; i++)
	   for (j=0L; j<nx; j++) {
	        idxbase=i*nx+j;
		for ( k = 0L; k < nz; k++) {
		   x[k] = (double)Output[k*nxny+idxbase];
		   //if (i==2 && j==2) printf("Analysis x[%d]=%f\n",k,x[k]);
	        }
		WaveletSplit_1D(	x, nz, y1, y2, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		for ( k = 0L; k < nz/2; k++)  Output[k*nxny+idxbase] 	   = (double)y1[k];
		for ( k = 0L; k < nz2; k++)   Output[(k+nz2)*nxny+idxbase] = (double)y2[k];
		//if (i==2 && j==2) for ( k = 0L; k < nz; k++)   printf("Analysis Output[%d]=%f\n",k,Output[k*nxny+idxbase]);
	   }

	FreeLineDouble(&FilterHighpass);
	FreeLineDouble(&FilterLowpass);  
	FreeVolumeDouble(&Low);
	FreeVolumeDouble(&High);
	FreeVolumeDouble(&LowLow);
	FreeVolumeDouble(&LowHigh);
	FreeVolumeDouble(&HighLow);
	FreeVolumeDouble(&HighHigh);
	FreeLineDouble(&x);
	FreeLineDouble(&y1);
	FreeLineDouble(&y2);
	
	*Status = !ERROR;
	return(!ERROR);

}

/* ----------------------------------------------------------------------------

	Function:	WaveletMerge_3D
	
	Purpose:	Two channel filtering and down-sampling for 2D real signals. 
				Input[nx][ny] -> output[nx][ny]
 				The input should be provided with four sectors in the classical 
 				way to show the wavelets transform

	Author:		Daniel Sage, EPFL, Biomedical Imaging Group
					
---------------------------------------------------------------------------- */
extern int	WaveletMerge_3D( 
					double	*Input, 
					double	*Output,
					long	nx, long ny, long nz,
					short	Filter,
					short	Order,
					double	Alpha,
					short	BoundaryConditions,
					int		*Status)
{
long	kx, ky, kz, nmax, nout, k, i, j, index, nx2, ny2, nz2, nxny,
        idxbase;
double	*FilterLowpass,	*FilterHighpass; 
long	NbLowpass[1], NbHighpass[1]; 
double	*Low, *High, *LowLow, *LowHigh, *HighLow, *HighHigh;
double	*x, *y1, *y2;

	nmax = MaxLong(MaxLong(nx, ny), nz);
	nx2 = MaxLong(nx / 2L, 1L);
	ny2 = MaxLong(ny / 2L, 1L);
	nz2 = MaxLong(nz / 2L, 1L);
	nxny=nx*ny;
	
	/* Get the size of the filter */
	if ( Filter != 5) {
		if ( WaveletFiltersGetSize(Filter, (short)(-Order), NbLowpass, NbHighpass) == ERROR) {
			MessageDisplay("ERROR - Impossible to get the size of the filter"); 
			*Status = ERROR;
			return( ERROR);
		}
	}
	else {
		if ( WaveletFiltersGetSize_Fract(Alpha, NbLowpass, NbHighpass) == ERROR) {
			MessageDisplay("ERROR - Impossible to get the size of the filter"); 
			*Status = ERROR;
			return( ERROR);
		}
	}
		
	/* Allocate the memory for the highpass filter */
  	AllocateLineDouble( &FilterHighpass, *NbHighpass, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (FilterHighpass)");
		return( ERROR);
  	}
  	
	/* Allocate the memory for the lowpass filter */
  	AllocateLineDouble( &FilterLowpass, *NbLowpass, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (FilterLowpass)");
		FreeLineDouble(&FilterHighpass);
		return( ERROR);
  	}
  	
  	/* Get the coefficient of the filters */ 
	if ( Filter != 5) {
		if ( WaveletFiltersGetCoef( Filter, (short)(-Order), FilterLowpass, FilterHighpass) == ERROR) {
			MessageDisplay("ERROR - Impossible to get the coefficient of the filter"); 
			*Status = ERROR;
			return( ERROR);
		}
	}
	else {
		if ( WaveletFiltersGetCoef_Fract( Alpha, FilterLowpass, FilterHighpass) == ERROR) {
			MessageDisplay("ERROR - Impossible to get the coefficient of the filter"); 
			*Status = ERROR;
			return( ERROR);
		}
	}
	
	
	/* Allocate the temporary vector */ 
  	AllocateVolumeDouble( &Low, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (Low)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		return( ERROR);
  	}
  	
  	AllocateVolumeDouble( &High, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (High)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		return( ERROR);
  	}

  	AllocateVolumeDouble( &LowLow, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (LowLow)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		return( ERROR);
  	}

  	AllocateVolumeDouble( &LowHigh, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (LowHigh)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		return( ERROR);
  	}

  	AllocateVolumeDouble( &HighLow, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (HighLow)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		return( ERROR);
  	}

  	AllocateVolumeDouble( &HighHigh, nx, ny, nz, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (HighHigh)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		FreeVolumeDouble(&HighLow);
		return( ERROR);
  	}

  	AllocateLineDouble( &x, nmax, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (x)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		FreeVolumeDouble(&HighLow);
		FreeVolumeDouble(&HighHigh);
		*Status = ERROR;
		return( ERROR);
  	}

  	AllocateLineDouble( &y1, nmax, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (y1)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		FreeVolumeDouble(&HighLow);
		FreeVolumeDouble(&HighHigh);
		FreeLineDouble(&x);
		return( ERROR);
  	}
  	
  	AllocateLineDouble( &y2, nmax, Status);
  	if (*Status == ERROR) {
		MessageDisplay( "ERROR - Unable to perform allocation (y2)");
		FreeLineDouble(&FilterHighpass);
		FreeLineDouble(&FilterLowpass);  
		FreeVolumeDouble(&Low);
		FreeVolumeDouble(&High);
		FreeVolumeDouble(&LowLow);
		FreeVolumeDouble(&LowHigh);
		FreeVolumeDouble(&HighLow);
		FreeVolumeDouble(&HighHigh);
		FreeLineDouble(&x);
		FreeLineDouble(&y1);
		return( ERROR);
  	}
  	
	
	/******************************************************************************/
	/* 3D                                                                         */
	/******************************************************************************/
	
	if (nz > 1)
   	   for (i=0L; i<ny; i++)
	   for (j=0L; j<nx; j++) {
	        idxbase=i*nx+j;
		for ( k = 0L; k < nz2; k++)  {
			y1[k] 	= Input[k*nxny+idxbase];
			y2[k] 	= Input[(k+nz2)*nxny+idxbase];
		}
		WaveletMerge_1D(	y1, y2, nz2, x, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		for ( k = 0L; k < nz; k++)
			Input[k*nxny+idxbase] = (double)x[k];
	   }

	/******************************************************************************/
	/* 1D                                                                         */
	/******************************************************************************/
	      
	if (ny == 1 && nz==1) {
		for ( i = 0L; i < nx2; i++)  {
			y1[i] 	= Input[i];
			y2[i] 	= Input[i+nx2];
		}
		WaveletMerge_1D(	y1, y2, nx2, x, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		for ( i = 0L; i < nx; i++)
			Output[i] = (double)x[i];
	}
	
	/******************************************************************************/
	/* 2D                                                                         */
	/******************************************************************************/
	
	if (ny > 1) {
	/* Extract component from the input */
	for (kz=0L;kz<nz;kz++) {
	   idxbase=kz*nxny;
	for ( i = 0L; i < nx2; i++)
	for ( j = 0L; j < ny2; j++)  {
		index = j*nx+i;
		LowLow[idxbase+index] 	= Input[idxbase+index];
		LowHigh[idxbase+index] 	= Input[idxbase+index+nx2];
		HighLow[idxbase+index] 	= Input[idxbase+(j+ny2)*nx+i];
		HighHigh[idxbase+index] = Input[idxbase+i+nx2+(j+ny2)*nx];
	
	}
	}
	
	/* x-processing for lowpass part*/
	for (kz=0L;kz<nz;kz++)
	for (ky=0L;ky<ny/2L;ky++) {
		GetxDoubleToDouble(LowLow,  nx, ny, nz, 0L, ky, kz, y1, nx); 
		GetxDoubleToDouble(LowHigh, nx, ny, nz, 0L, ky, kz, y2, nx); 
		WaveletMerge_1D(	y1, y2, nx2, x, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		PutxDoubleToDouble(Low, nx, ny, nz, 0L, ky, kz, x, nx); 
	}

	/* x-processing for highpass part*/
	for (kz=0L;kz<nz;kz++)
	for (ky=0L;ky<ny/2L;ky++) {
		GetxDoubleToDouble(HighLow,  nx, ny, nz, 0L, ky, kz, y1, nx); 
		GetxDoubleToDouble(HighHigh, nx, ny, nz, 0L, ky, kz, y2, nx); 
		WaveletMerge_1D(	y1, y2, nx2, x, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		PutxDoubleToDouble(High, nx, ny, nz, 0L, ky, kz, x, nx); 
	}
		
	/* y-processing */
	for (kz=0L;kz<nz;kz++)
	for (kx=0L;kx<nx;kx++) {
		GetyDoubleToDouble(Low,  nx, ny, nz, kx, 0L, kz, y1, ny);
		GetyDoubleToDouble(High, nx, ny, nz, kx, 0L, kz, y2, ny);
		WaveletMerge_1D(	y1, y2, ny2, x, &nout,
							FilterLowpass, *NbLowpass,
							FilterHighpass, *NbHighpass,
							BoundaryConditions); 
		PutyDoubleToDouble(Output,nx, ny, nz, kx, 0L, kz, x, ny); 
	}
	}
	
	/* --- Exit --- */
	FreeLineDouble(&FilterHighpass);
	FreeLineDouble(&FilterLowpass);  
	FreeVolumeDouble(&Low);
	FreeVolumeDouble(&High);
	FreeVolumeDouble(&LowLow);
	FreeVolumeDouble(&LowHigh);
	FreeVolumeDouble(&HighLow);
	FreeVolumeDouble(&HighHigh);
	FreeLineDouble(&x);
	FreeLineDouble(&y1);
	FreeLineDouble(&y2);
	
	*Status = !ERROR;
	return(!ERROR);
}


/* ----------------------------------------------------------------------------

	Function:	WaveletSplit_1D
	
	Purpose:	Two channel filtering and down-sampling for 1D real signals. 
				Interface to call the basic routines:
					WaveletSplit_1DMirror() if BoundaryConditions is 1
					WaveletSplit_1DPeriodic() if BoundaryConditions is 2

				lowpass filter:  v[nv]
				highpass filter: w[nw]

				x[nx] -> y1[nx/2] & y2[nx/2]
 

	Author:		Daniel Sage, EPFL, Biomedical Imaging Group
					
---------------------------------------------------------------------------- */
static int WaveletSplit_1D(	
					double x[], long nx,
					double y1[], double y2[], long *ny,
					double v[], long nv, 
					double w[], long nw,
					short BoundaryConditions)
{

	/* --- Check the parameters --- */
	
	if (x == (double *)NULL) {
		MessageDisplay("ERROR - Pointer of input is null");
		return(ERROR);
	}

	if (nx <= 1) {
		MessageDisplay("ERROR - Length of input should be greater than 1");
		return(ERROR);
	}
	
	if ((nx/2L)*2L != nx) {  	/* The input length should be even */
		MessageDisplay("ERROR - Length of input should be even");
		return(ERROR);
	}
	
	if (y1 == (double *)NULL) {
		MessageDisplay("ERROR - Pointer to an output is null");
		return(ERROR);
	}
	
	if (y2 == (double *)NULL) {
		MessageDisplay("ERROR - Pointer to an output is null");
		return(ERROR);
	}
	
	if (ny == (long *)NULL) {
		MessageDisplay("ERROR - Pointer to an output is null");
		return(ERROR);
	}
	
	if (v == (double *)NULL) {
		MessageDisplay("ERROR - Pointer to the lowpass filter is null");
		return(ERROR);
	}
	
	if (nv <= 0 ) {
		MessageDisplay("ERROR - Length of the filter be greater than 0");
		return(ERROR);
	}
	
	if (w == (double *)NULL) {
		MessageDisplay("ERROR - Pointer to the highpass filter is null");
		return(ERROR);
	}
	
	if (nw <= 0 ) {
		MessageDisplay("ERROR - Length of the filter be greater than 0");
		return(ERROR);
	}
			
	/* --- Call the processing function */
	
	switch(BoundaryConditions) {
		case 1:
			WaveletSplit_1DMirror( x, nx,y1, y2, ny, v, nv, w, nw);
			break;
		case 2: 
			WaveletSplit_1DPeriodic( x, nx,y1, y2, ny, v, nv, w, nw);
			break;
		default:	
			MessageDisplay("ERROR - Boundary should be Periodic or Mirror");
			return(ERROR);
	}
		 
	return(!ERROR);
}
					
/* ----------------------------------------------------------------------------

	Function: 	WaveletMerge_1D
	
	Purpose: 	Two channel up-sampling and filtering for 1D real signals. 
				Interface to call the basic routines:
					WaveletSplit_1DMirror() if BoundaryConditions is 1
					WaveletSplit_1DPeriodic() if BoundaryConditions is 2

				lowpass filter:  v[nv]
				highpass filter: w[nw]

				y1[ny] & y2[ny] -> x[2*ny]

	Author:		Daniel Sage, EPFL, Biomedical Imaging Group
	
---------------------------------------------------------------------------- */
static int WaveletMerge_1D(
				double y1[],double y2[],long ny,
				double x[],long *nx,
				double v[],long nv,
				double w[],long nw,
				short BoundaryConditions)
{

	/* --- Check the parameters --- */
	
	if (y1 == (double *)NULL) {
		MessageDisplay("ERROR - Pointer of input is null");
		return(ERROR);
	}
	
	if (y2 == (double *)NULL) {
		MessageDisplay("ERROR - Pointer of input is null");
		return(ERROR);
	}
	
	if (ny <= 1) {
		MessageDisplay("ERROR - Length of input should be greater than 1");
		return(ERROR);
	}
	
	if (nx == (long *)NULL) {
		MessageDisplay("ERROR - Pointer to an output is null");
		return(ERROR);
	}
	
	if (nx == (long *)NULL) {
		MessageDisplay("ERROR - Pointer to an output is null");
		return(ERROR);
	}
	
	if (v == (double *)NULL) {
		MessageDisplay("ERROR - Pointer to the lowpass filter is null");
		return(ERROR);
	}
	
	if (nv <= 0 ) {
		MessageDisplay("ERROR - Length of the filter be greater than 0");
		return(ERROR);
	}
	
	if (w == (double *)NULL) {
		MessageDisplay("ERROR - Pointer to the highpass filter is null");
		return(ERROR);
	}
	
	if (nw <= 0 ) {
		MessageDisplay("ERROR - Length of the filter be greater than 0");
		return(ERROR);
	}
	
	/* --- Call the processing function */
	
	switch(BoundaryConditions) {
		case 1:		
			WaveletMerge_1DMirror( y1, y2, ny, x, nx, v, nv, w, nw);
			break;
			
		case 2:
			WaveletMerge_1DPeriodic( y1, y2, ny, x, nx, v, nv, w, nw);
			break;
			
		default: 
			MessageDisplay("ERROR - Boundary should be Periodic or Mirror");
			return(ERROR);
		 }
		 
	return(!ERROR);

}

/* ----------------------------------------------------------------------------

	Function:	WaveletSplit_1DMirror
	
	Purpose:	Two channel filtering and down-sampling for 1D
				real signals. The filters are symmetrical and the signal
				is extended using mirror extrapolation.

				lowpass filter:  v[nv]
				highpass filter: w[nw]

				x[nx] -> y1[nx/2] & y2[nx/2]
 
	Author:		Michael Unser, EPFL, Biomedical Imaging Group
					
	History:	MU	Oct 1992	tested
				DS	Dec 1998	renaming, double version  

---------------------------------------------------------------------------- */
static void WaveletSplit_1DMirror(	
					double x[], long nx,
					double y1[], double y2[], long *ny,
					double v[], long nv, 
					double w[], long nw)
{	  
double 	pix;
long 	i, j, k, nyy, n, j1, j2, kn;

	nyy = nx / 2L;
	n   = nyy * 2L;					/* signal length */
	kn  = 2L * n - 2L;				/* Global period */

	if ((nv<=1L)||(nw<=1L))			/* Haar transform */
		for (i=0L;i<nyy;i++) {
			j=2*i;
			y1[i]=(+x[j]+x[j+1])/2.;
			y2[i]=(-x[j]+x[j+1])/2.;
	}
	else
		for (i=0L;i<nyy;i++) {
			j=i*2;
			pix=x[j]*v[0];
			for (k=1;k<nv;k++) {
				j1=j-k;
				j2=j+k;
		
				if (j1<0L) {
					j1=(-j1) % kn;
					if (j1>n-1) j1=kn-j1;
				}
				if (j2>n-1) {
					j2=j2 % kn;
					if (j2>n-1) j2=kn-j2;
				}

				pix=pix+v[k]*(x[j1]+x[j2]);
		}
		y1[i]=pix;

		j=j+1L;
		pix=x[j]*w[0];
		for (k=1L;k<nw;k++) {
			j1=j-k;
			j2=j+k;

			if (j1<0L) {
				j1=(-j1) % kn;
				if (j1>n-1) j1=kn-j1;
			}
			if (j2>n-1) {
				j2=j2 % kn;
				if (j2>n-1L) j2=kn-j2;
			}

			pix=pix+w[k]*(x[j1]+x[j2]);
		}
		y2[i]=pix;
	}
	*ny=nyy;
}

/* ----------------------------------------------------------------------------

	Function: 	WaveletSplit_1DPeriodic
	
	Purpose: 	Two channel filtering and down-sampling for 1D
				real signals. The filters are symmetrical and the signal
				is extended using periodic boundary conditions.
				
				lowpass filter:  v[nv]
				highpass filter: w[nw]

				x[nx] -> y1[nx/2] & y2[nx/2]

	Author:		Michael Unser, EPFL, Biomedical Imaging Group
					
	History:	MU	Oct 1992	tested
				DS	Dec 1998	renaming, double version  

---------------------------------------------------------------------------- */
static void WaveletSplit_1DPeriodic(
					double x[], long nx,
					double y1[], double y2[], long *ny,
					double v[], long nv, 
					double w[], long nw)
{	  
double	pix;
long	i, j, k, nyy, n, j1, j2, kn, knn;

	nyy=nx/2;
	n=nyy*2;						/* signal length */
	kn=n;						    /* Global period */
	knn=20*n;

	if ((nv<=1)||(nw<=1))			/* Haar transform */
	for (i=0;i<nyy;i++) {
		j=2*i;
		y1[i]=(+x[j]+x[j+1])/2.;
		y2[i]=(-x[j]+x[j+1])/2.;
	}
	else
	for (i=0;i<nyy;i++) {
		j=i*2;
		pix=x[j]*v[0];
		for (k=1;k<nv;k++) {
			j1=j-k;
			j2=j+k;
			if (j1<0) j1=(knn+j1) % kn;
			if (j2>n-1) j2=j2 % kn;
			pix=pix+v[k]*(x[j1]+x[j2]);
		}
		y1[i]=pix;

		j=j+1;
		pix=x[j]*w[0];
		for (k=1;k<nw;k++) {
			j1=j-k;
			j2=j+k;
			if (j1<0) j1=(knn+j1) % kn;
			if (j2>n-1) j2=j2 % kn;
			pix=pix+w[k]*(x[j1]+x[j2]);
		}
		y2[i]=pix;
	}
	*ny=nyy;
}


/* ----------------------------------------------------------------------------

	Function: 	WaveletMerge_1DMirror
	
	Purpose: 	Two channel up-sampling and filtering for 1D
				real signals. The filters are symmetrical and the
				boundary conditions are mirror.

				lowpass filter:  v[nv]
				highpass filter: w[nw]

				y1[ny] & y2[ny] -> x[2*ny]

 	Author:		Michael Unser, EPFL, Biomedical Imaging Group
					
	History:	MU	Oct 1992	tested
                MU	Apr 1997	corrected boundary conditions)
				DS	Dec 1998	renaming, double version  

---------------------------------------------------------------------------- */
static void WaveletMerge_1DMirror(
				double y1[],double y2[],long ny,
				double x[],long *nx,
				double v[],long nv,
				double w[],long nw)
{	  
double	pix1, pix2;
long	i, j, k, kk, i1, i2, k01, k02, n, kn;

	*nx=ny*2;
	k01=(nv/2)*2-1;
	k02=(nw/2)*2-1;
	n=ny;
	kn=2*n-1;

	if ((nv<=1)||(nw<=1))				/* Haar transform */
		for (i=0;i<ny;i++) {
			j=2*i;
			x[j]=y1[i]-y2[i];
			x[j+1]=y1[i]+y2[i];
		}
	else
	for (i=0;i<ny;i++)	{
	      j=2*i;
	      pix1=v[0]*y1[i];
	      for (k=2;k<nv;k=k+2) {
	        i1=i-(k/2);
	        i2=i+(k/2);
			if (i1<0) {					/* standard boundary */
			  i1=(-i1) % kn;
			  if (i1>n-1) i1=kn-i1;
			}
			if (i2>n-1) {				/* non-standard boundary */
			  i2=(i2) % kn;
			  if (i2>n-1) i2=kn-i2;
			}
			/* NOTE : This is the correct boundary condition assuming that 
			the original image has an even size (2n) and the decimation 
			was performed according to the sequence 1,3,5,É,2n-1.	*/
	        pix1=pix1+v[k]*(y1[i1]+y1[i2]);
	      }

		pix2=0.;
		for (k=-k02;k<nw;k=k+2) {
			kk=abs(k);
			i1=i+(k-1)/2;
			if (i1<0) {					/* non-standard boundary */
				i1=(-i1-1) % kn;
				if (i1>n-1) i1=kn-1-i1;
			}
			if (i1>n-1) {				/* standard boundary */
				i1=i1 % kn;
				if (i1>n-1) i1=kn-1-i1;
			}
			/* NOTE : This is the correct boundary condition assuming that 
			the original image has an even size (2n) and the decimation 
			was performed according to the sequence 2,4,6,É,2n.	*/
			pix2=pix2+w[kk]*y2[i1];
		}

		x[j]=(pix1+pix2);

		j=j+1;
		pix1=0.;
		for (k=-k01;k<nv;k=k+2) {
			kk=abs(k);
			i1=i+(k+1)/2;
			if (i1<0) {					/* standard boundary */
				i1=(-i1) % kn;
				if (i1>n-1) i1=kn-i1;
			}
			if (i1>n-1) {				/* non-standard boundary */
				i1=(i1) % kn;
				if (i1>n-1) i1=kn-i1;
			}
			pix1=pix1+v[kk]*y1[i1];
		}
		pix2=w[0]*y2[i];
		for (k=2;k<nw;k=k+2) {
			i1=i-(k/2);
			i2=i+(k/2);
			if (i1<0) {					/* non-standard boundary */
				i1=(-i1-1) % kn;
				if (i1>n-1) i1=kn-1-i1;
			}
			if (i2>n-1) {				/* standard boundary */
				i2=i2 % kn;
				if (i2>n-1) i2=kn-1-i2;
			}
			pix2=pix2+w[k]*(y2[i1]+y2[i2]);
		}
		x[j]=(pix1+pix2);
	}
}

/* ----------------------------------------------------------------------------

	Function:	WaveletMerge_1DPeriodic
	
	Purpose:	Two channel up-sampling and filtering for 1D
				real signals. The filters are symmetrical and the
				boundary conditions are periodic.
	
				lowpass filter:  v[nv]
				highpass filter: w[nw]

				y1[ny] & y2[ny] -> x[2*ny]
  
	Author:		Michael Unser, EPFL, Biomedical Imaging Group
					
	History:	MU	Oct 1992	tested
				DS	Dec 1998	renaming, double version  


---------------------------------------------------------------------------- */
static void WaveletMerge_1DPeriodic(
				double y1[],double y2[], long ny, 
				double x[], long *nx, 
				double v[],long nv,
				double w[],long nw)
{	  
double 	pix1, pix2;
long	 i, j, k, kk, i1, i2, k01, k02, n, kn, knn;

	*nx=ny*2;
	k01=(nv/2)*2-1;
	k02=(nw/2)*2-1;
	n=ny;
	kn=n;
	knn=20*n;

	if ((nv<=1)||(nw<=1))				/* Haar transform */
		for (i=0;i<ny;i++) {
			j=2*i;
			x[j]=y1[i]-y2[i];
			x[j+1]=y1[i]+y2[i];
		}
	else
	for (i=0;i<ny;i++)	{
		j=2*i;
		pix1=v[0]*y1[i];
		for (k=2;k<nv;k=k+2) {
			i1=i-(k/2);
			i2=i+(k/2);

			if (i1<0) i1=(knn+i1) % kn;
			if (i2>n-1) i2=i2 % kn;

			pix1=pix1+v[k]*(y1[i1]+y1[i2]);
		}

		pix2=0.;
		for (k=-k02;k<nw;k=k+2) {
			kk=abs(k);
			i1=i+(k-1)/2;
			if (i1<0) i1=(knn+i1) % kn;
			if (i1>n-1) i1=i1 % kn;
			pix2=pix2+w[kk]*y2[i1];
		}
		x[j]=(pix1+pix2);

		j=j+1;
		pix1=0.;
		for (k=-k01;k<nv;k=k+2) {
			kk=abs(k);
			i1=i+(k+1)/2;
			if (i1<0) i1=(knn+i1) % kn;
			if (i1>n-1) i1=i1 % kn;
			pix1=pix1+v[kk]*y1[i1];
		}
		pix2=w[0]*y2[i];
		for (k=2;k<nw;k=k+2) {
			i1=i-(k/2);
			i2=i+(k/2);
			if (i1<0) i1=(knn+i1) % kn;
			if (i2>n-1) i2=i2 % kn;
			pix2=pix2+w[k]*(y2[i1]+y2[i2]);
		}
		x[j]=(pix1+pix2);
	}
}



