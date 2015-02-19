/*--------------------------------------------------------------------------*/
extern int  AllocateLineDouble
    (
        double *(Line[]),   /* double output pointer */
        long LineLength,   /* length of the line */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
extern int  AllocateVolumeFloat
    (
        float **Volume,   /* float output pointer */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
extern int  AllocateVolumeShort
    (
        short **Volume,   /* short output pointer */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
extern int  CopyFloatToFloat
    (
        float *VolumeSource,  /* float input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        float *VolumeDestination, /* float output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        long NxCopy,    /* width of the block to copy */
        long NyCopy,    /* height of the block to copy */
        long NzCopy    /* depth of the block to copy */
    );

/*--------------------------------------------------------------------------*/
extern int  CopyFloatToShort
    (
        float *VolumeSource,  /* float input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        short *VolumeDestination, /* short output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        long NxCopy,    /* width of block to copy */
        long NyCopy,    /* height of block to copy */
        long NzCopy    /* depth of the block to copy */
    );

/*--------------------------------------------------------------------------*/
extern int  CopyShortToFloat
    (
        short *VolumeSource,  /* short input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        float *VolumeDestination, /* float output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        long NxCopy,    /* width of block to copy */
        long NyCopy,    /* height of block to copy */
        long NzCopy    /* depth of the block to copy */
    );

/*--------------------------------------------------------------------------*/
extern int  CopyShortToShort
    (
        short *VolumeSource,  /* short input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        short *VolumeDestination, /* short output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        long NxCopy,    /* width of block to copy */
        long NyCopy,    /* height of block to copy */
        long NzCopy    /* depth of the block to copy */
    );

/*--------------------------------------------------------------------------*/
extern int  FreeLineDouble
    (
        double *(Line[])   /* 1D double array */
    );

/*--------------------------------------------------------------------------*/
extern int  FreeVolumeFloat
    (
        float **Volume   /* 3D float array */
    );

/*--------------------------------------------------------------------------*/
extern int  FreeVolumeShort
    (
        short **Volume   /* 3D short array */
    );

/*--------------------------------------------------------------------------*/
extern int  GetxFloatToDouble
    (
        float *VolumeSource,  /* float input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        double RowDestination[], /* double output data */
        long NxCopy    /* length of the output */
    );

/*--------------------------------------------------------------------------*/
extern int  GetxShortToDouble
    (
        short *VolumeSource,  /* short input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        double RowDestination[], /* double output data */
        long NxCopy    /* length of the output */
    );

/*--------------------------------------------------------------------------*/
extern int  GetyFloatToDouble
    (
        float *VolumeSource,  /* float input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        double ColumnDestination[],/* double output data */
        long NyCopy    /* length of the output */
    );

/*--------------------------------------------------------------------------*/
extern int  GetyShortToDouble
    (
        short *VolumeSource,  /* short input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        double ColumnDestination[],/* double output data */
        long NyCopy    /* length of the output */
    );

/*--------------------------------------------------------------------------*/
extern int  GetzFloatToDouble
    (
        float *VolumeSource,  /* float input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        double PillarDestination[],/* double output data */
        long NzCopy    /* length of the output */
    );

/*--------------------------------------------------------------------------*/
extern int  GetzShortToDouble
    (
        short *VolumeSource,  /* short input data */
        long NxSource,   /* width of the input */
        long NySource,   /* height of the input */
        long NzSource,   /* depth of the input */
        long XSource,   /* x coordinate to get from */
        long YSource,   /* y coordinate to get from */
        long ZSource,   /* z coordinate to get from */
        double PillarDestination[],/* double output data */
        long NzCopy    /* length of the output */
    );

/*--------------------------------------------------------------------------*/
extern int  PutxDoubleToFloat
    (
        float *VolumeDestination, /* float output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        double RowSource[],  /* double input data */
        long NxCopy    /* length of the input */
    );

/*--------------------------------------------------------------------------*/
extern int  PutxDoubleToShort
    (
        short *VolumeDestination, /* short output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        double RowSource[],  /* double input data */
        long NxCopy    /* length of the input */
    );

/*--------------------------------------------------------------------------*/
extern int  PutyDoubleToFloat
    (
        float *VolumeDestination, /* float output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        double ColumnSource[],  /* double input data */
        long NyCopy    /* length of the input */
    );

/*--------------------------------------------------------------------------*/
extern int  PutyDoubleToShort
    (
        short *VolumeDestination, /* short output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        double ColumnSource[],  /* double input data */
        long NyCopy    /* length of the input */
    );

/*--------------------------------------------------------------------------*/
extern int  PutzDoubleToFloat
    (
        float *VolumeDestination, /* float output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        double PillarSource[],  /* double input data */
        long NzCopy    /* length of the input */
    );

/*--------------------------------------------------------------------------*/
extern int  PutzDoubleToShort
    (
        short *VolumeDestination, /* short output data */
        long NxDestination,  /* width of the output */
        long NyDestination,  /* height of the output */
        long NzDestination,  /* depth of the output */
        long XDestination,  /* x coordinate to put into */
        long YDestination,  /* y coordinate to put into */
        long ZDestination,  /* z coordinate to put into */
        double PillarSource[],  /* double input data */
        long NzCopy    /* length of the input */
    );

