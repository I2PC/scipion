/*
        Base on rwMRC.h
        Header file for reading and writing MRC files
        Format: 3D crystallographic image file format for the MRC package
        Author: Bernard Heymann
        Created: 19990321       Modified: 20030723
*/

#ifndef RWMRC_H
#define RWMRC_H

///@defgroup MRC MRC File format
///@ingroup ImageFormats

// I/O prototypes
/** MRC Reader
  * @ingroup MRC
*/
int readMRC(int img_select, bool isStack=false);

/** MRC Writer
  * @ingroup MRC
*/
int writeMRC(int img_select, bool isStack=false, int mode=WRITE_OVERWRITE, std::string bitDepth="", bool adjust=false);

#endif
