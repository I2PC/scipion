/*
 * Copyright INRIA
 * Author Gregoire Malandain (greg@sophia.inria.fr)
 * Date: June, 9 1998
 */

#ifndef _typedefs_h_
#define _typedefs_h_

#ifdef __cplusplus
extern "C"
{
#endif

    /* Differents type coding for images and buffers.
     */
    typedef enum {
        TYPE_UNKNOWN /* unknown type */,
        UCHAR  /* unsigned char */,
        SCHAR  /* signed char */,
        USHORT /* unsigned short int */,
        SSHORT /* signed short int */,
        INT    /* signed int */,
        ULINT  /* unsigned long int */,
        FLOAT  /* float */,
        DOUBLE  /* double */
    } ImageType, bufferType;

    typedef short int          s8;
    /*Quite esto porque daba un warning en convert.c line:319 (asig -128 a char)*/
    /*typedef char               s8;*/
    typedef unsigned char      u8;
    typedef short int          s16;
    typedef unsigned short int u16;
    typedef int                i32;
    typedef int                s32;
    typedef unsigned long int  u64;
    typedef float              r32;
    typedef double             r64;



    /* Typedef Booleen
     */
    typedef enum {
        False = 0,
        True = 1
    } typeBoolean;


#ifdef __cplusplus
}
#endif

#endif
