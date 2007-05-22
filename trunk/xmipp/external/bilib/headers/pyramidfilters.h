#ifndef _PYRAMIDFILTERS
#define _PYRAMIDFILTERS
/* ----------------------------------------------------------------------------
 Filename:   pyramidfilters.h

 Project: Biomedical Imaging Library

 Author:  Daniel Sage
    Swiss Federal Institute of Technology - Lausanne
    Biomedical Imaging Group
    EPFL/DMT/IOA, BM-Ecublens, CH-1015 Lausanne, Switzerland

 Date:  17 March 1999

 Purpose: Header associated to pyramidfilters.c

---------------------------------------------------------------------------- */

/**@name Pyramid filters */
//@{
/** PyramidFilterSplinel2.
 Function:  PyramidFilterSplinel2

  Purpose:   Initializes down- and up-sampling filter arrays for
     least squares splines of order 0 to 3.  (little l_2 norm)
      g : reduce filter
      h : expand filter

  Author:  Michael Unser, NIH, BEIP, May 1992
*/
extern void PyramidFilterSplinel2(double g[], long *ng, double *h, long *nh, long Order);

/** PyramidFilterSplineL2.
 Function:  PyramidFilterSplineL2

  Purpose:   Initializes down- and up-sampling filter arrays for
     L2 spline pyramid of order 0 to 5.
      g : reduce filter
      h : expand filter

  Author:  Michael Unser, NIH, BEIP, May 1992
*/
extern void PyramidFilterSplineL2(double g[], long *ng, double *h, long *nh, long Order);

/** PyramidFilterCentered.
 Function: PyramidFilterCentered

 Purpose: Initializes down- and up-sampling filter arrays for
    least squares CENTERED splines of order 0 to 4.  (little l_2 norm)
     g : reduce filter
     h : expand filter

 Note:  filter arrays should be defined as
     double g[20],h[20] filter arrays
     short *ng,*nh; number of taps
     short Order; order of the spline

 Author:  Patrick Brigger, NIH, BEIP May 1996
    Daniel Sage, EPFL, Biomedical Imaging Group, November 1999
*/
extern void PyramidFilterCentered(double g[], long *ng, double h[], long *nh, long Order);

/** PyramidFilterCenteredL2.
 Function: PyramidFilterCenteredL2

 Purpose: Initializes the symmetric down- and up-sampling filter arrays for
     L2 spline pyramid of order 0 to 5 when the downsampled grid is centered.
    These filters have then to be followed by a Haar filter.
     g: reduce filter
     h: expand filter

 Note:  filter arrays should be defined as
     float g[35],h[35] filter arrays
     short *ng,*nh number of taps
     short Order order of the spline

 Author:  Patrick Brigger, NIH, BEIP, April 1996
    Daniel Sage, EPFL, Biomedical Imaging Group, November 1999
*/
extern void PyramidFilterCenteredL2(double g[], long *ng, double h[], long *nh, long Order);

/** PyramidFilterCenteredL2Derivate.
 Function: PyramidFilterCenteredL2Derivate

 Purpose: Initializes the symmetric down- and up-sampling filter arrays for
    L2 DERIVATIVE spline pyramid of order 0 to 5 when the downsampled
    grid is centered.
    These filters have then to be followed by a Derivative Haar filter.
     g : reduce filter
     h : expand filter
  Note:  filter arrays should be defined as
       float g[35],h[35] filter arrays
     short *ng,*nh number of taps
     short Order order of the spline

 Author:  Patrick Brigger, NIH, BEIP, April 1996
    Daniel Sage, EPFL, Biomedical Imaging Group, November 1999
*/
extern void PyramidFilterCenteredL2Derivate(double g[], long *ng, double h[], long *nh, long Order);

//@}
#endif
