/**@name Signal processing windows */
//@{
/*--------------------------------------------------------------------------*/
/** Returns the value of a Bartlet (triangular) window evaluated at Argument.
    The width of the symmetric window is (2 * H) */
extern double Bartlet
    (
        double Argument,   /* input */
        long WindowHalfLength /* half size H of the window */
    );

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blackman window evaluated at Argument.
    The width of the symmetric window is (2 * H) */
extern double Blackman
    (
        double Argument,   /* input */
        long WindowHalfLength /* half size H of the window */
    );

/*--------------------------------------------------------------------------*/
/** Returns the value of a Dirichlet (rectangular) window evaluated at
    Argument. The width of the symmetric window is (2 * H) */
extern double Dirichlet
    (
        double Argument,   /* input */
        long WindowHalfLength /* half size H of the window */
    );

/*--------------------------------------------------------------------------*/
/** Returns the value of a Hamming window evaluated at Argument.
    Classic Hamming weights are used.
    The width of the symmetric window is (2 * H) */
extern double HammingClassic
    (
        double Argument,   /* input */
        long WindowHalfLength /* half size H of the window */
    );

/*--------------------------------------------------------------------------*/
/** Returns the value of a Hamming window evaluated at Argument.
    Optimal (non-classic) weights are used.
    The width of the symmetric window is (2 * H) */
extern double HammingExact
    (
        double Argument,   /* input */
        long WindowHalfLength /* half size H of the window */
    );

/*--------------------------------------------------------------------------*/
/** Returns the value of a Hanning window evaluated at Argument.
    The width of the symmetric window is (2 * H) */
extern double Hanning
    (
        double Argument,   /* input */
        long WindowHalfLength /* half size H of the window */
    );
//@}
