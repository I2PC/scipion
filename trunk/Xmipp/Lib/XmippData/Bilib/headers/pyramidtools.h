/* ----------------------------------------------------------------------------
	Filename:  	pyramidtools.h
	
	Project:	Biomedical Imaging Library
 
	Author:		Daniel Sage
				Swiss Federal Institute of Technology - Lausanne
				Biomedical Imaging Group
				EPFL/DMT/IOA, BM-Ecublens, CH-1015 Lausanne, Switzerland
				
	Date:		17 March 1999
	
	Purpose:	Header associated to pyramidtools.c
			 
---------------------------------------------------------------------------- */

/**@name Pyramid tools */
//@{
/** Get Pyramid Filter.
	Function: 
		GetPyramidFilter
		
	Purpose:
		Get the coefficients of the filter (reduce and expand filter)
		Return the coefficients in g[ng] and in h[nh]
		
	Convention:
		g[ng] for the reduce filter
		h[nh] for the expansion filter
		
	Parameters:
		Filter is the name of the filter
		
		Order is the order for the filters based on splines
			For the "Spline" filter, Order is 0, 1, 2 or 3
			For the "Spline L2" filter, Order is 0, 1, 3 or 5
			For the "Centered Spline" filter, Order is 0, 1, 2, 3 or 4
			For the "Centered Spline L2" filter, Order is 0, 1, 2, 3 or 4
			
		IsCentered is a return value indicates if the filter is a centered filter
			TRUE if it is a centered filter
			FALSE if it is not a centered filter 
*/
extern int GetPyramidFilter(
				char *Filter, 
				long Order, 
				double g[],long *ng,
				double h[],long *nh, 
				short *FlagCentered		
				);	
					
/** Reduce 2D.
	Function: 
		Reduce_2D
	
	Purpose: 
 		Reduces an image by a factor of two in each dimension.
 		
	Note: 
 		Expects the output array (Out) to be allocated.

	Parameters:
		Input image:  	In[NxIn*NyIn]
		Output image: 	Out[NxIn/2*NyIn/2]
		Filter:			g[ng] coefficients of the filter
*/
extern int Reduce_2D(	
				double *In, long NxIn, long NyIn,
				double *Out,
				double w[], long nw,
				short FlagCentered
				);
				
/** Reduce 3D.
	Function: 
		Reduce_3D
	
	Purpose: 
 		Reduces a volume by a factor of two in each dimension.
 		
	Note: 
 		Expects the output array (Out) to be allocated.

	Parameters:
		Input image:  	In[NxIn*NyIn*NzIn]
		Output image: 	Out[NxIn/2*NyIn/2*NzIn/2]
		Filter:			g[ng] coefficients of the filter
*/
extern int Reduce_3D(	
				double *In, long NxIn, long NyIn, long NzIn,
				double *Out,
				double w[], long nw,
				short FlagCentered
				);
				
/** Expand 2D.
	Function: 
		Expand_2D
	
	Purpose: 
 		Expands an image by a factor of two in each dimension.
 		
	Note: 
 		Expects the output array (Out) to be allocated.

	Parameters:
		Input volume:  	In[NxIn,NyIn]
		Output voulme: 	Out[NxIn*2,NyIn*2]
		Filter coef:	h[nh]
*/
extern int Expand_2D(	
				double *In, long NxIn, long NyIn,
				double *Out,
				double h[], long nh,
				short FlagCentered
				);

/** Expand 3D.
	Function: 
		Expand_3D
	
	Purpose: 
 		Expands a volume by a factor of two in each dimension.
 		
	Note: 
 		Expects the output array (Out) to be allocated.

	Parameters:
		Input volume:  	In[NxIn,NyIn,NzIn]
		Output voulme: 	Out[NxIn*2,NyIn*2,NzIn*2]
		Filter coef:	h[nh]
*/
extern int Expand_3D(	
				double *In, long NxIn, long NyIn, long NzIn,
				double *Out,
				double h[], long nh,
				short FlagCentered
				);

/** Reduce 1D.
	Purpose:	
		Basic function to reduce a 1D signal
	
	Parameters:
		In[NxIn] is the input signal	(NxIn should be greater than 2 and even)
		Out[NxIn/2] is the output signal
		g[ng] is an array that contains the coefficients of the filter
		
	Author:  
		Michael Unser, NIH, BEIP, June 1994
		Daniel Sage, EPFL, Biomedical Imaging Group, April 1999
*/
extern void Reduce_1D(	
				double x[], long nx,
				double y[],
				double w[], long nw,
				short FlagCentered
				);
						
/** Expand 1D.
	Purpose:	
		Basic function to expand a 1D signal
	
	Parameters:
		In[NxIn] is the input signal	(NxIn should be greater than 1)
		Out[NxIn*2] is the output signal
		w[nw] is an array that contains the coefficients of the filter
		
	Author:  
		Michael Unser, NIH, BEIP, June 1994
		Daniel Sage, EPFL, Biomedical Imaging Group, April 1999
*/
extern void Expand_1D(	
				double x[], long nx,
				double y[],
				double w[], long nw,
				short FlagCentered
				);
						
//@}
