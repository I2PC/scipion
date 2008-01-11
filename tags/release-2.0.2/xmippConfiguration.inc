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
