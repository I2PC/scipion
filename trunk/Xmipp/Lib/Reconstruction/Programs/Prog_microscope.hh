/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2001 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/
#ifndef _PROG_MICROSCOPE_HH
#  define _PROG_MICROSCOPE_HH

#ifdef _HAVE_VTK
#include "Prog_FourierFilter.hh"
#include <XmippData/xmippProgs.hh>

/**@name Microscope program */
//@{
/* Microscope Program Parameters ------------------------------------------- */
/** Parameter class for the project program */
class Prog_Microscope_Parameters: public Prog_parameters {
public:
   /// Filename with the CTF
   FileName fn_ctf;
   /// Total noise power
   double   sigma;
   /// Low pass frequency before CTF
   double   low_pass_before_CTF;
   /// Filename with the root squared spectrum for noise after CTF
   FileName fn_after_ctf;
   /// Defocus change (%)
   double   defocus_change;
public:
   /// CTF
   FourierMask ctf;
   /// Low pass filter, if it is 0 no lowpass filter is applied
   FourierMask lowpass;
   /// After CTF noise root squared spectrum
   FourierMask after_ctf;
   /// Noise power before CTF
   double   sigma_before_CTF;
   /// Noise power after CTF
   double   sigma_after_CTF;
   /// Input image Xdim
   int Xdim;
   /// Input image Ydim
   int Ydim;
public:
   /** Read from a command line.
       An exception might be thrown by any of the internal conversions,
       this would mean that there is an error in the command line and you
       might show a usage message. */
   void read(int argc, char **argv) _THROW;

   /** Usage message.
       This function shows the way of introducing this parameters. */
   void usage();

   /** Show parameters. */
   void show();

   /** Produce side information. */
   void produce_side_info();

   /** Apply to a single image. The image is modified.
       If the CTF is randomly selected then a new CTF is generated
       for each image */
   void apply(matrix2D<double> &I);
};
//@}
#endif
#endif
