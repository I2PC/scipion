/**@name Exception configuration
   Xmipp admits 2 kinds of error control systems: the first one is very simple,
   every time an illegal operation is performed, the running program is
   aborted and an exit error code and the corresponding message are shown;
   the second method makes use of C++ exception handling capability, then
   at every illegal operation an exception of type Xmipp_error is produced,
   showing at that time the same error code and message as in the first
   method, but instead of aborting the program, there are routines capable
   to handle this exceptions and that may take actions to solve the problem
   without aborting (for instance, when reading a document file imagine that
   there is a line which doesn't conform the file format, the routine reading
   lines would throw an exception that would be caught by the routine reading
   the whole file and this latter routine decides that it is not such a big
   problem to abort the program, it informs about the error in the line, 
   and skip it). As you see the exception handling is a very powerful
   method to save your data before the catastrophe, to ignore not so
   important errors, or whatever. Many functions have got exceptions
   associated, you may catch the errors in your programs and take the
   convenient actions to avoid them.
   
   The problem with exceptions is that not all compilers manage exceptions
   (g++ do), and then some compilation problems would be produced and
   you couldn't compile Xmipp at all. So this option allows you to
   switch between method one and method two. If _NO_EXCEPTION is defined
   then method 1 is selected, otherwise you will work with the exception
   error handling.
   
   @see Xmipp_error
*/
//#define _NO_EXCEPTION

/**@name Time measurements
   Xmipp is also prepared to make some time measurements of the elapsed time
   of some parts of the algorithms, some estimation about the remaining time
   to go, ... This part are not essential in Xmipp algorithms and depends
   a lot on the specific machine Xmipp is running in. The included header
   files for that purposes might be different in different platforms, even
   the data structures for time may differ. Xmipp has been proved to work
   well in SGI IRIX, and Linux. If you experience some compilation troubles 
   due to the time functions you can disable them (all functions are empty,
   they are called but they do nothing) with this option. Define _NO_TIME
   if you don't want the time functions to be active.
*/
//#define _NO_TIME
