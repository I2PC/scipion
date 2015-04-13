%module swig_frm 
%{
#define SWIG_FILE_WITH_INIT
int frm(double *f, int dim1, double *g, int dim2, double *res, int dim3);
int frm_corr(double *f, int dim1, double *g, int dim2, double *coef_corr, int dim3);
int frm_fourier_corr(double *fr, int dim1, double *fi, int dim2, double *gr, int dim3, double *gi, int dim4, double *coef_corr, int dim5);
int find_topn_angles(double *coef_corr, int dim1, int bw, double *res, int dim2, double dist_cut);
int enlarge2(double *in, int dim1, int nx, int ny, int nz, double *out, int dim2);
int get_constraint_vol(double *cv, int dim1, int bw, double phi, double psi, double the, double nearby);
%}
%include "numpy.i"
%init %{
    import_array();
%}

int frm(double *IN_ARRAY1, int DIM1, double *IN_ARRAY1, int DIM1, double *INPLACE_ARRAY1, int DIM1);
int frm_corr(double *IN_ARRAY1, int DIM1, double *IN_ARRAY1, int DIM1, double *INPLACE_ARRAY1, int DIM1);
int frm_fourier_corr(double *IN_ARRAY1, int DIM1, double *IN_ARRAY1, int DIM1, double *IN_ARRAY1, int DIM1, double *IN_ARRAY1, int DIM1, double *INPLACE_ARRAY1, int DIM1);
int find_topn_angles(double *IN_ARRAY1, int DIM1, int bw, double *INPLACE_ARRAY1, int DIM1, double dist_cut);
int enlarge2(double *IN_ARRAY1, int DIM1, int nx, int ny, int nz, double *INPLACE_ARRAY1, int DIM1);
int get_constraint_vol(double *INPLACE_ARRAY1, int DIM1, int bw, double phi, double psi, double the, double nearby);
