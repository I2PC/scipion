/**@name Swap axes in a volume */
//@{
/*--------------------------------------------------------------------------*/
/** Swap the x-axis and the y-axis of a volume.
    Input VolumeSource is a (float)volume of size (Nxy x Nyx x Nz).
    Output VolumeDestination is a (float)volume of size (Nyx x Nxy x Nz).

    success: return(!ERROR); failure: return(ERROR); */
extern int  SwapXyVolumeFloat
    (
        float *VolumeSource,  /* float input data */
        float *VolumeDestination, /* float output data */
        long Nxy,    /* input width, output height */
        long Nyx,    /* input height, output width */
        long Nz     /* depth */
    );

/*--------------------------------------------------------------------------*/
/** Swap the y-axis and the z-axis of a volume.
    Input VolumeSource is a (float)volume of size (Nx x Nyz x Nzy).
    Output VolumeDestination is a (float)volume of size (Nx x Nzy x Nyz).

    success: return(!ERROR); failure: return(ERROR); */
extern int  SwapYzVolumeFloat
    (
        float *VolumeSource,  /* float input data */
        float *VolumeDestination, /* float output data */
        long Nx,     /* width */
        long Nyz,    /* input height, output depth */
        long Nzy     /* input depth, output height */
    );

/*--------------------------------------------------------------------------*/
/** Swap the z-axis and the x-axis of a volume.
    Input VolumeSource is a (float)volume of size (Nxz x Ny x Nzx).
    Output VolumeDestination is a (float)volume of size (Nzx x Ny x Nxz).

    success: return(!ERROR); failure: return(ERROR); */
extern int  SwapZxVolumeFloat
    (
        float *VolumeSource,  /* float input data */
        float *VolumeDestination, /* float output data */
        long Nxz,    /* input width, output depth */
        long Ny,     /* height */
        long Nzx     /* input depth, output width */
    );
//@}
