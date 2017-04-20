#include <time.h>
#include <sys/time.h>


/** Native implementation when interpolation is either nearest neighbor or linear*/
void cuda_rotate_image_linear_nearest(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim,
		double *ang, int interp, int wrap, int first_call, struct ioTime* time);

