
void freeRotatedImage (float **rotated_image);

/** Native implementation when interpolation is either nearest neighbor or linear*/
float** cuda_rotate_image_linear_nearest_parallel(float **image, size_t *Xdim, size_t *Ydim, size_t *Zdim, double **ang, int interp, int wrap, int first_call);

