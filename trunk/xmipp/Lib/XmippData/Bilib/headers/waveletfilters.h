/* ----------------------------------------------------------------------------
 Function:  WaveletFilters.h
 
 Purpose:	Header file for WaveletFilter.c
---------------------------------------------------------------------------- */

extern int WaveletFiltersGetSize(short Filter, short Order, long *nh, long *ng);
extern int WaveletFiltersGetCoef(short Filter, short Order, double *h, double *g);
extern int WaveletFiltersGetCoef_Fact(double Alpha, double *h, double *g);


