

extern int WaveletSplit_3D(
    double *Input,
    double *Output,
    long nx, long ny, long nz,
    short Filter,
    short Order,
    double Alpha,
    short BoundaryConditions,
    int  *Status);

extern int WaveletMerge_3D(
    double *Input,
    double *Output,
    long nx, long ny, long nz,
    short Filter,
    short Order,
    double Alpha,
    short BoundaryConditions,
    int  *Status);

