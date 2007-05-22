/*--------------------------------------------------------------------------*/
extern int  GetFoldedIndex
    (
        long InputIndex,   /* index to fold back */
        long *OutputIndex,  /* folded index */
        long SignalLength,  /* length of the signal */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
extern int  GetFoldedValueDouble
    (
        double Signal[],   /* double input data */
        long InputIndex,   /* index to fold back */
        double *OutputValue,  /* double output value */
        long SignalLength,  /* length of the input data */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
extern int  GetFoldedValueShort
    (
        short Signal[],   /* short input data */
        long InputIndex,   /* index to fold back */
        short *OutputValue,  /* short output value */
        long SignalLength,  /* length of the input data */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

