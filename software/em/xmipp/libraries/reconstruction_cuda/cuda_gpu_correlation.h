
#include <time.h>
#include <sys/time.h>



void cuda_fft(float *image, float *recovered_image, float *module, float *angle, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch);

void cuda_cart2polar(float *image, float *polar_image, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch);

