#include "../extrema.h"

static int _VERBOSE_ = 0;

/* 
 * epsilon value to select gradient extrema candidates
 */
static double _EPSILON_NORM_ = 0.5;

/*
 * epsilon value to decide of the interpolation type.
 * If one derivative's absolute value is larger than this
 * epsilon (close to one), then we use the nearest value
 * else we perform a [bi,tri]linear interpolation.
 */
static double _EPSILON_DERIVATIVE_ = 0.95;


#define EXIT_ON_FAILURE 0
#define EXIT_ON_SUCCESS 1
int partial_derivative_3D( void *bufferIn,
				bufferType typeIn,
				void *bufferOut,
				bufferType typeOut,
				int *bufferDims,
				int *borderLengths,
				float *filterCoefs,
				recursiveFilterType filterType,
				derivativeOrder Zgradient[3])
{
  char *proc="Extract_Gradient_Maxima_3D";
  /*
   * auxiliary buffer
   */ 
  float *tmpBuffer = (float*)NULL;
  float *bufferZsmoothed = (float*)NULL;
  float *bufferZderivated = (float*)NULL;
  /*
   * Pointers
   */
  /* 
   * gx[0] points toward the X gradient of the current slice
   * gx[0] points toward the X gradient of the next slice
   */
  float *gx[2] = { (float*)NULL, (float*)NULL };
  /*
   * gy: idem gx but for the Y gradient
   */
  float *gy[2] = { (float*)NULL, (float*)NULL };
  float *gz = (float*)NULL;
  /*
   * norme[0] points toward the gradient modulus of the previous slice
   * norme[1] points toward the gradient modulus of the current slice
   * norme[2] points toward the gradient modulus of the next slice
   */
  float *norme[3] = { (float*)NULL, (float*)NULL, (float*)NULL }; 
  float *sliceZsmoothed = (float*)NULL;
  float *pt = (float*)NULL;
  /*
   * additional parameters for recursive filtering
   */
  derivativeOrder Xgradient[3] = { DERIVATIVE_1_EDGES, SMOOTHING, NODERIVATIVE };
  derivativeOrder Ygradient[3] = { SMOOTHING, DERIVATIVE_1_EDGES, NODERIVATIVE };
  /*derivativeOrder Zgradient[3] = { NODERIVATIVE, DERIVATIVE_1,
  NODERIVATIVE};*/
  derivativeOrder Zsmoothing[3] = { NODERIVATIVE, NODERIVATIVE, SMOOTHING };
  int sliceDims[3];
  /*
   *
   */
  int z, dimxXdimy;

  /* 
   * We check the buffers' dimensions.
   */
  if ( (bufferDims[0] <= 0) || (bufferDims[1] <= 0) || (bufferDims[2] <= 0) ) {
    if ( _VERBOSE_ > 0 )
      fprintf( stderr, " Fatal error in %s: improper buffer's dimension.\n", proc );
    return( EXIT_ON_FAILURE );
  }

  /*
   * May we perform a 3D edge detection?
   */
  if ( bufferDims[2] <= 4 ) {
    if ( _VERBOSE_ > 0 )
      fprintf( stderr, " Warning in %s: switch to 2D edge extraction.\n", proc );
    return( Extract_Gradient_Maxima_2D( bufferIn, typeIn,
					bufferOut, typeOut,
					bufferDims, borderLengths,
					filterCoefs, filterType ) );
  }

  /*
   *
   */
  dimxXdimy = bufferDims[0] * bufferDims[1];
  sliceDims[0] = bufferDims[0];
  sliceDims[1] = bufferDims[1];
  sliceDims[2] = 1;
  
  /*
   * test of the coefficients
   */
  if ( (filterCoefs[0] < 0.0) || (filterCoefs[1] < 0.0) ||
       (filterCoefs[2] < 0.0) ) {
    if ( _VERBOSE_ > 0 )
      fprintf( stderr, " Error in %s: negative coefficient's value.\n", proc );
    return( EXIT_ON_FAILURE );
  }
  
  /* 
   * Allocation of auxiliary buffers.
   *
   * We need a 3D buffer for the Z component of the
   * gradient, plus a 3D buffer for the 3D buffer 
   * smoothed along Z, plus 7 2D buffers for the
   * X component of the gradient in both the current
   * and the next slices, idem for the Y component,
   * idem for the modulus plus one 2D buffer for
   * the modulua in the previous slice.
   *
   * If the buffer bufferOut is of type FLOAT,
   * we use it as the Z component of the gradient.
   *
   * This Z component will be used to stored the
   * extrema of the gradient.
   */
  tmpBuffer = (float*)malloc( 7 * dimxXdimy * sizeof( float ) );
  if ( tmpBuffer == (float*)NULL ) {
    if ( _VERBOSE_ > 0 ) {
      fprintf( stderr, " Fatal error in %s:", proc );
      fprintf( stderr, " unable to allocate auxiliary buffer.\n" );      
    }
    return( EXIT_ON_FAILURE );
  }
  gx[0] = tmpBuffer;
  gx[1] = gx[0] + dimxXdimy;
  gy[0] = gx[1] + dimxXdimy;
  gy[1] = gy[0] + dimxXdimy;
  norme[0] = gy[1] + dimxXdimy;
  norme[1] = norme[0] + dimxXdimy;
  norme[2] = norme[1] + dimxXdimy;
  
  bufferZsmoothed = (float*)malloc( bufferDims[2] * dimxXdimy * sizeof( float ) );
  if ( bufferZsmoothed == (float*)NULL ) {
    if ( _VERBOSE_ > 0 ) {
      fprintf( stderr, " Fatal error in %s:", proc );
      fprintf( stderr, " unable to allocate auxiliary first 3D buffer.\n" );
    }
    free( tmpBuffer );
    return( EXIT_ON_FAILURE );
  }
  
  if ( typeOut == FLOAT ) {
  /* Daba warning porque c++ prohibe la asignacion de void * */
  /*  bufferZderivated = bufferOut;*/
  bufferZderivated = (float *)bufferOut;
  } else {
    bufferZderivated = (float*)malloc( bufferDims[2] * dimxXdimy * sizeof( float ) );
    if ( bufferZderivated == (float*)NULL ) {
      if ( _VERBOSE_ > 0 ) {
	fprintf( stderr, " Fatal error in %s:", proc );
	fprintf( stderr, " unable to allocate auxiliary first 3D buffer.\n" );
      }
      free( tmpBuffer );
      free( bufferZsmoothed );
      return( EXIT_ON_FAILURE );
    }
  }
  
  /* 
   * Computation of the Z component of the gradient.
   * Computation of the input buffer smoothed along Z.
   */
  if ( RecursiveFilterOnBuffer( bufferIn, typeIn,
				bufferZderivated, FLOAT,
				bufferDims, borderLengths,
				Zgradient, filterCoefs,
				filterType ) == 0 ) {
    if ( _VERBOSE_ > 0 ) {
      fprintf( stderr, " Fatal error in %s:", proc );
      fprintf( stderr, " unable to compute Z gradient.\n" );
    }
    free( tmpBuffer );
    free( bufferZsmoothed );
    if ( typeOut != FLOAT ) free( bufferZderivated );
    return( EXIT_ON_SUCCESS );
  }
				     
  /*
   * conversion of the buffer bufferZderivated of type FLOAT
   * into the buffer bufferOut.
   */
  
  if (typeOut != FLOAT ) {
    ConvertBuffer( bufferZderivated, FLOAT, 
		   bufferOut, typeOut, bufferDims[2]*dimxXdimy);
  }

  free( tmpBuffer );
  free( bufferZsmoothed );
  if ( typeOut != FLOAT ) free( bufferZderivated );
  return( EXIT_ON_SUCCESS );
}
