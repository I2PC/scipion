/*--------------------------------------------------------------------------*/
extern int  GetTimeStamp
    (
        struct TTimeStamp
        *TimeStamp   /* output time stamp */
    );

/*--------------------------------------------------------------------------*/
extern int  GetElapsedTime
    (
        struct TTimeStamp
        OldTimeStamp,  /* oldest time stamp */
        struct TTimeStamp
        NewTimeStamp,  /* newest time stamp */
        double *SystemTimeDifference,
        /* elapsed system time in seconds */
        double *RealTimeDifference /* true elapsed time in seconds */
    );

