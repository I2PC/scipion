/* 
 * Copyright INRIA
 * Author Gregoire Malandain (greg@sophia.inria.fr)
 * Date: June, 9 1998
 */
#include "../recbuffer.h"

static int _VERBOSE_ = 0;

#define EXIT_ON_FAILURE 0
#define EXIT_ON_SUCCESS 1

int RecursiveFilterOnBuffer( void *bufferIn,
			     bufferType typeIn,
			     void *bufferOut,
			     bufferType typeOut,
			     int *bufferDims,
			     int *borderLengths,
			     derivativeOrder *derivatives,
			     float *filterCoefs,
			     recursiveFilterType filterType )
{
  char *proc = "RecursiveFilterOnBuffer";
  register int dimx, dimxXdimy;
  int dimy, dimz;
  register int x, y, z;
  /* 
   *obviously, we need to perform the computation 
   * with float or double values. For this reason,
   * we allocate an auxiliary buffer if the output buffer
   * is not of type float or double.
   */
  void *bufferToBeProcessed = (void*)NULL;
  bufferType typeToBeProcessed = TYPE_UNKNOWN;
  void *bufferResult = (void*)NULL;
  bufferType typeResult = TYPE_UNKNOWN;
  /*
   * lines' lengths
   */
  int lengthX = 0;
  int lengthY = 0;
  int lengthZ = 0;
  int maxLengthline = 0;
  int borderXlength = 0;
  int borderYlength = 0;
  int borderZlength = 0;
  /*
   * 1D arrays for computations.
   */
  double *theLine = (double*)NULL;
  double *resLine = (double*)NULL;
  double *tmpLine = (double*)NULL;
  /*
   * pointers for computations;
   */
  register r32 *r32firstPoint = (r32*)NULL;
  register r64 *r64firstPoint = (r64*)NULL;
  register r32 *r32_pt = (r32*)NULL;
  register r64 *r64_pt = (r64*)NULL;
  register double *dbl_pt1 = (double*)NULL;
  register double *dbl_pt2 = (double*)NULL;
  register double dbl_first = 0.0;
  register double dbl_last = 0.0;
  int offsetLastPoint = 0;
  int offsetNextFirstPoint = 0;
  register r32 *r32firstPointResult = (r32*)NULL;
  register r64 *r64firstPointResult = (r64*)NULL;
  double *theLinePlusBorder = (double*)NULL;
  double *resLinePlusBorder = (double*)NULL;

  
  /* 
   * We check the buffers' dimensions.
   */
  dimx = bufferDims[0];   dimy = bufferDims[1];   dimz = bufferDims[2];
  dimxXdimy = dimx * dimy;
  if ( (dimx <= 0) || (dimy <= 0) || (dimz <= 0) ) {
    if ( _VERBOSE_ > 0 )
      fprintf( stderr, " Fatal error in %s: improper buffer's dimension.\n", proc );
    return( EXIT_ON_FAILURE );
  }
  /*
   * We check the pointers.
   */
  if ( (bufferIn == (void*)NULL) || (bufferOut == (void*)NULL) ) {
    if ( _VERBOSE_ > 0 )
      fprintf( stderr, " Fatal error in %s: NULL pointer on buffer.\n", proc );
    return( EXIT_ON_FAILURE );
  }

  /* 
   * May we use the buffer bufferOut as the bufferResult?
   * If its type is FLOAT or DOUBLE, then yes.
   * If not, we have to allocate an auxiliary buffer.
   */
  if ( (typeOut == FLOAT) || (typeOut == DOUBLE) ) {
    bufferResult = bufferOut;
    typeResult = typeOut;
  } else {
    bufferResult = (void*)malloc( (dimx*dimy*dimz) * sizeof(r32) );
    if ( bufferResult == (void*)NULL ) {
      if ( _VERBOSE_ > 0 )
	fprintf( stderr, " Fatal error in %s: unable to allocate auxiliary buffer.\n", proc );
      return( EXIT_ON_FAILURE );
    }
    typeResult = FLOAT;
  }
  
  /* 
   * May we consider the buffer bufferIn as the bufferToBeProcessed?
   * If its type is FLOAT or DOUBLE, then yes.
   * If not, we convert it into the buffer bufferResult, and this
   * last buffer is the bufferToBeProcessed.
   */
  if ( (typeIn == FLOAT) || (typeIn == DOUBLE) ) {
    bufferToBeProcessed = bufferIn;
    typeToBeProcessed = typeIn;
  } else {
    ConvertBuffer( bufferIn, typeIn, bufferResult, typeResult, (dimx*dimy*dimz) );
    bufferToBeProcessed = bufferResult;
    typeToBeProcessed = typeResult;
  }

  /*
   * Estimation of the lines' length along each direction.
   */
  borderXlength = borderLengths[0];
  borderYlength = borderLengths[1];
  borderZlength = borderLengths[2];
  if ( borderXlength < 0 ) borderXlength = 0;
  if ( borderYlength < 0 ) borderYlength = 0;
  if ( borderZlength < 0 ) borderZlength = 0;
  lengthX = dimx + 2 * borderXlength;
  lengthY = dimx + 2 * borderYlength;
  lengthZ = dimx + 2 * borderZlength;
  maxLengthline = lengthX;
  if ( maxLengthline < lengthY ) maxLengthline = lengthY;
  if ( maxLengthline < lengthZ ) maxLengthline = lengthZ;
  if ( maxLengthline <= 0 ) {
    if ( _VERBOSE_ > 0 )
      fprintf( stderr, " Error in %s: unable to deal with dimensions = 0.\n", proc );
    if ( (typeOut != FLOAT) && (typeOut != DOUBLE) )
      free( bufferResult );
    return( EXIT_ON_FAILURE );
  }
  /*
   * Allocations of work arrays. 
   * We will use them to process each line.
   */
  theLine = (double*)malloc( 3 * maxLengthline * sizeof(double) );
  if ( theLine == (double*)NULL ) {
    if ( _VERBOSE_ > 0 )
      fprintf( stderr, " Fatal error in %s: unable to allocate auxiliary work arrays.\n", proc );
    if ( (typeOut != FLOAT) && (typeOut != DOUBLE) )
      free( bufferResult );
    return( EXIT_ON_FAILURE );
  }
  resLine = theLine + maxLengthline;
  tmpLine = resLine + maxLengthline;

  /*
   * From now,
   * typeToBeProcessed is either FLOAT or DOUBLE
   * so is typeResult.
   */


  /*
   * Processing along X.
   */
  if ( dimx > 4 )
  if (derivatives[0] != NODERIVATIVE)
  if (filterCoefs[0] > 0.0) {
    if ( _VERBOSE_ != 0 )
      fprintf( stderr, " %s: processing along X.\n", proc );
    InitRecursiveCoefficients( (double)filterCoefs[0], filterType, derivatives[0] );
    
    r64firstPoint = (r64*)bufferToBeProcessed;
    r32firstPoint = (r32*)bufferToBeProcessed;

    r64firstPointResult = (r64*)bufferResult;
    r32firstPointResult = (r32*)bufferResult;

    offsetLastPoint = borderXlength + dimx - 1;

    theLinePlusBorder = theLine + borderXlength;
    resLinePlusBorder = resLine + borderXlength;

    /*
     * There are dimz*dimy X lines to be processed.
     */
    for ( z=0; z<dimz; z++ )
    for ( y=0; y<dimy; y++ ) {
      /*
       * Acquiring a X line.
       */ 
      dbl_pt1 = theLinePlusBorder;
      switch ( typeToBeProcessed ) {
      case DOUBLE :
	(void)memcpy( (void*)dbl_pt1, (void*)r64firstPoint, dimx );
	r64firstPoint += dimx;
	break;
      case FLOAT :
      default :
	for ( x=0; x<dimx; x++, dbl_pt1++, r32firstPoint++ ) *dbl_pt1 = *r32firstPoint;
      }
      /*
       * Adding points at both ends of the line.
       */
      if ( borderXlength > 0 ) {
	dbl_pt1 = theLine + borderXlength;   dbl_first = *dbl_pt1;
	dbl_pt2 = theLine + offsetLastPoint; dbl_last  = *dbl_pt2;
	for ( x=0; x<borderXlength; x++ ) {
	  *--dbl_pt1 = dbl_first;
	  *++dbl_pt2 = dbl_last;
	}
      }
      /*
       * Processing the line.
       */
      if ( RecursiveFilter1D( theLine, resLine, tmpLine, resLine, lengthX ) == 0 ) {
	if ( _VERBOSE_ != 0 ) 
	  fprintf(stderr," Error in %s: unable to process X line (y=%d,z=%d).\n", proc, y, z);
	if ( (typeOut != FLOAT) && (typeOut != DOUBLE) )
	  free( bufferResult );
	free( (void*)theLine );
	return( EXIT_ON_FAILURE );
      }
      /*
       * Copy the result into the buffer bufferResult.
       */
      dbl_pt1 = resLinePlusBorder;
      switch ( typeResult ) {
      case DOUBLE :
	(void)memcpy( (void*)r64firstPointResult, (void*)dbl_pt1, dimx );
	r64firstPointResult += dimx;
	break;
      case FLOAT :
      default :
	for ( x=0; x<dimx; x++, dbl_pt1++, r32firstPointResult++ )
	  *r32firstPointResult = (r32)(*dbl_pt1);
      }
    }
    
    /*
     * The next buffer to be processed is the buffer
     * bufferResult.
     */
    bufferToBeProcessed = bufferResult;
    typeToBeProcessed = typeResult;
  
  } /* end of Processing along X. */
  
  /*
   * Processing along Y.
   */
  if ( dimy > 4 )
  if (derivatives[1] != NODERIVATIVE)
  if (filterCoefs[1] > 0.0) {
    if ( _VERBOSE_ != 0 )
      fprintf( stderr, " %s: processing along Y.\n", proc );
    InitRecursiveCoefficients( (double)filterCoefs[1], filterType, derivatives[1] );
    
    r64firstPoint = (r64*)bufferToBeProcessed;
    r32firstPoint = (r32*)bufferToBeProcessed;

    r64firstPointResult = (r64*)bufferResult;
    r32firstPointResult = (r32*)bufferResult;

    offsetLastPoint = borderYlength + dimy - 1;
    offsetNextFirstPoint = dimx * dimy - dimx;

    theLinePlusBorder = theLine + borderYlength;
    resLinePlusBorder = resLine + borderYlength;

    /*
     * There are dimz*dimx Y lines to be processed.
     */
    for ( z=0; z<dimz; z++ ) {
      for ( x=0; x<dimx; x++ ) {
      /*
       * Acquiring a Y line.
       */ 
	dbl_pt1 = theLinePlusBorder;
	switch ( typeToBeProcessed ) {
	case DOUBLE :
	  r64_pt = r64firstPoint;
	  for ( y=0; y<dimy; y++, dbl_pt1++, r64_pt += dimx ) *dbl_pt1 = *r64_pt;
	  /*
	   * Going to the first point of the next Y line
	   */
	  r64firstPoint ++;
	  break;
	case FLOAT :
	default :
	  r32_pt = r32firstPoint;
	  for ( y=0; y<dimy; y++, dbl_pt1++, r32_pt += dimx ) *dbl_pt1 = *r32_pt;
	  r32firstPoint ++;
	}
	/*
	 * Adding points at both ends of the line.
	 */
	if ( borderYlength > 0 ) {
	  dbl_pt1 = theLine + borderYlength;   dbl_first = *dbl_pt1;
	  dbl_pt2 = theLine + offsetLastPoint; dbl_last  = *dbl_pt2;
	  for ( y=0; y<borderYlength; y++ ) {
	    *--dbl_pt1 = dbl_first;
	    *++dbl_pt2 = dbl_last;
	  }
	}
	/*
	 * Processing the line.
	 */
	if ( RecursiveFilter1D( theLine, resLine, tmpLine, resLine, lengthY ) == 0 ) {
	  if ( _VERBOSE_ != 0 ) 
	    fprintf(stderr," Error in %s: unable to process Y line (x=%d,z=%d).\n", proc, x, z);
	  if ( (typeOut != FLOAT) && (typeOut != DOUBLE) )
	    free( bufferResult );
	  free( (void*)theLine );
	  return( EXIT_ON_FAILURE );
	}
	/*
	 * Copy the result into the buffer bufferResult.
	 */
	dbl_pt1 = resLinePlusBorder;
	switch ( typeResult ) {
	case DOUBLE :
	  r64_pt = r64firstPointResult;
	  for ( y=0; y<dimy; y++, dbl_pt1++, r64_pt += dimx ) *r64_pt = *dbl_pt1;
	  r64firstPointResult ++;
	  break;
	case FLOAT :
	default :
	  r32_pt = r32firstPointResult;
	  for ( y=0; y<dimy; y++, dbl_pt1++, r32_pt += dimx ) *r32_pt = *dbl_pt1;
	  r32firstPointResult ++;
	}
      }
      /*
       * Going to the first point of the next Y line
       * which is the first Y line of the next slice.
       *
       * The pointer r[32,64]firstPoint[Result] has
       * already been increased by dimx. To reach
       * the first point of the next slice, we
       * have to increase it by (dimx*dimy)-dimx.
       */
      switch ( typeToBeProcessed ) {
      case DOUBLE :
	r64firstPoint += offsetNextFirstPoint;
	break;
      case FLOAT :
      default :
	r32firstPoint += offsetNextFirstPoint;
      }
      switch ( typeResult ) {
      case DOUBLE :
	r64firstPointResult += offsetNextFirstPoint;
	break;
      case FLOAT :
      default :
	r32firstPointResult += offsetNextFirstPoint;
      }
    }
    
    /*
     * The next buffer to be processed is the buffer
     * bufferResult.
     */
    bufferToBeProcessed = bufferResult;
    typeToBeProcessed = typeResult;
  
  } /* end of Processing along Y. */
  

  /*
   * Processing along Z.
   */
  if ( dimz > 4 )
  if (derivatives[2] != NODERIVATIVE)
  if (filterCoefs[2] > 0.0) {
    if ( _VERBOSE_ != 0 )
      fprintf( stderr, " %s: processing along Z.\n", proc );
    InitRecursiveCoefficients( (double)filterCoefs[2], filterType, derivatives[2] );
    
    r64firstPoint = (r64*)bufferToBeProcessed;
    r32firstPoint = (r32*)bufferToBeProcessed;

    offsetLastPoint = borderZlength + dimz - 1;

    r64firstPointResult = (r64*)bufferResult;
    r32firstPointResult = (r32*)bufferResult;

    offsetLastPoint = borderZlength + dimz - 1;

    theLinePlusBorder = theLine + borderYlength;
    resLinePlusBorder = resLine + borderYlength;

    /*
     * There are dimy*dimx Z lines to be processed.
     */
    for ( y=0; y<dimy; y++ )
    for ( x=0; x<dimx; x++ ) {
      /*
       * Acquiring a Z line.
       */ 
      dbl_pt1 = theLinePlusBorder;
      switch ( typeToBeProcessed ) {
      case DOUBLE :
	r64_pt = r64firstPoint;
	for ( z=0; z<dimz; z++, dbl_pt1++, r64_pt += dimxXdimy ) *dbl_pt1 = *r64_pt;
	/*
	 * Going to the first point of the next Z line
	 */
	r64firstPoint ++;
	break;
      case FLOAT :
      default :
	r32_pt = r32firstPoint;
	for ( z=0; z<dimz; z++, dbl_pt1++, r32_pt += dimxXdimy ) *dbl_pt1 = *r32_pt;
	r32firstPoint ++;
      }
      /*
       * Adding points at both ends of the line.
       */
      if ( borderZlength > 0 ) {
	dbl_pt1 = theLine + borderZlength;   dbl_first = *dbl_pt1;
	dbl_pt2 = theLine + offsetLastPoint; dbl_last  = *dbl_pt2;
	for ( z=0; z<borderZlength; z++ ) {
	  *--dbl_pt1 = dbl_first;
	  *++dbl_pt2 = dbl_last;
	}
      }
      /*
       * Processing the line.
       */
      if ( RecursiveFilter1D( theLine, resLine, tmpLine, resLine, lengthZ ) == 0 ) {
	if ( _VERBOSE_ != 0 ) 
	  fprintf(stderr," Error in %s: unable to process Z line (x=%d,y=%d).\n", proc, x, y);
	if ( (typeOut != FLOAT) && (typeOut != DOUBLE) )
	  free( bufferResult );
	free( (void*)theLine );
	return( EXIT_ON_FAILURE );
      }
      
      /*
       * Copy the result into the buffer bufferResult.
       */
      dbl_pt1 = resLinePlusBorder;
      switch ( typeResult ) {
      case DOUBLE :
	r64_pt = r64firstPointResult;
	for ( z=0; z<dimz; z++, dbl_pt1++, r64_pt += dimxXdimy ) *r64_pt = *dbl_pt1;
	r64firstPointResult ++;
	break;
      case FLOAT :
      default :
	r32_pt = r32firstPointResult;
	for ( z=0; z<dimz; z++, dbl_pt1++, r32_pt += dimxXdimy ) *r32_pt = *dbl_pt1;
	r32firstPointResult ++;
      }
    }
  } /* end of Processing along Z. */
  



  /*
   * From bufferResult to bufferOut
   */
  ConvertBuffer( bufferResult, typeResult, bufferOut, typeOut, (dimx*dimy*dimz) );

  /*
   * Releasing the buffers.
   */
  if ( (typeOut != FLOAT) && (typeOut != DOUBLE) )
    free( bufferResult );
  free( (void*)theLine );
  
  return( EXIT_ON_SUCCESS );
}






void Recbuffer_verbose ( )
{
  _VERBOSE_ = 1;
  Recline_verbose ( );
}





void Recbuffer_noverbose ( )
{
  _VERBOSE_ = 0;
  Recline_noverbose ( );
}
