/*--------------------------------------------------------------------------*/
extern int  AllocateVolumeDouble
(
    double **Volume,   /* double output pointer */
    long Nx,     /* width of the volume */
    long Ny,     /* height of the volume */
    long Nz,     /* depth of the volume */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
extern int  CopyDoubleToFloat
(
    double *VolumeSource,  /* double input data */
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
extern int CopyDoubleToDouble
(
    double  *VolumeSource,      /* double input data */
    long    NxSource,           /* width of the input */
    long    NySource,           /* height of the input */
    long    NzSource,           /* depth of the input */
    long    XSource,            /* x coordinate to get from */
    long    YSource,            /* y coordinate to get from */
    long    ZSource,            /* z coordinate to get from */
    double  *VolumeDestination, /* double output data */
    long    NxDestination,      /* width of the output */
    long    NyDestination,      /* height of the output */
    long    NzDestination,      /* depth of the output */
    long    XDestination,       /* x coordinate to put into */
    long    YDestination,       /* y coordinate to put into */
    long    ZDestination,       /* z coordinate to put into */
    long    NxCopy,             /* width of the block to copy */
    long    NyCopy,             /* height of the block to copy */
    long    NzCopy              /* depth of the block to copy */
);

/*--------------------------------------------------------------------------*/
extern int  CopyFloatToDouble
(
    float *VolumeSource,  /* float input data */
    long NxSource,   /* width of the input */
    long NySource,   /* height of the input */
    long NzSource,   /* depth of the input */
    long XSource,   /* x coordinate to get from */
    long YSource,   /* y coordinate to get from */
    long ZSource,   /* z coordinate to get from */
    double *VolumeDestination, /* double output data */
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
extern int  FreeVolumeDouble
(
    double **Volume   /* 3D double array */
);

/*--------------------------------------------------------------------------*/
extern int  GetxDoubleToDouble
(
    double *VolumeSource,  /* double input data */
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
extern int  GetyDoubleToDouble
(
    double *VolumeSource,  /* double input data */
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
extern int  GetzDoubleToDouble
(
    double *VolumeSource,  /* double input data */
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
extern int  PutxDoubleToDouble
(
    double *VolumeDestination, /* double output data */
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
extern int  PutyDoubleToDouble
(
    double *VolumeDestination, /* double output data */
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
extern int  PutzDoubleToDouble
(
    double *VolumeDestination, /* double output data */
    long NxDestination,  /* width of the output */
    long NyDestination,  /* height of the output */
    long NzDestination,  /* depth of the output */
    long XDestination,  /* x coordinate to put into */
    long YDestination,  /* y coordinate to put into */
    long ZDestination,  /* z coordinate to put into */
    double PillarSource[],  /* double input data */
    long NzCopy    /* length of the input */
);

