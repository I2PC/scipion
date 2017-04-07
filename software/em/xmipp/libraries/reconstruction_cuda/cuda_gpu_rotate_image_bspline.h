
/** Bspline interpolation as implemented in http://www.dannyruijters.nl/cubicinterpolation/ */
void cuda_rotate_image_bspline(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim, double *ang, int wrap, int first_call);

