/*--------------------------------------------------------------------------*/
struct TTraceLine
{
    long k[3];    /* index */
    int  Quadrant;   /* direction of propagation [0..26] \ {13} */
    int  ModificationCode; /* which component will change */
    double Entry;    /* intersection of the line entering the current voxel */
    double Exit;    /* intersection of the line exiting the current voxel */
    double P[3];    /* point on line */
    double q[3];    /* direction of line */
};

