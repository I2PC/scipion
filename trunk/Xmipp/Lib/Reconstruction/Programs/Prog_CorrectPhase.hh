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

#ifndef _CORRECTPHASE_HH
   #define _CORRECTPHASE_HH

#ifdef _HAVE_VTK

#include "Prog_FourierFilter.hh"
#include <XmippData/xmippSelFiles.hh>

/**@name Correct Phase */
//@{
/// Correct Phase parameters
class CorrectPhase_Params {
public:
   /// Name of the CTF file
   FileName fn_ctf;

   /// CTF description
   bool CTF_description_file;

   /// Epsilon=minimum CTF value to correct
   double epsilon;

   #define CORRECT_SETTING_SMALL_TO_ZERO 0
   #define CORRECT_LEAVING_SMALL         1
   #define CORRECT_AMPLIFYING_NOT_SMALL  2
   /** Correcting method. Valid methods are:
      \\CORRECT_SETTING_SMALL_TO_ZERO: Values where the CTF<epsilon are
         set to 0
      \\CORRECT_LEAVING_SMALL: Values where the CTF<epsilon are
         left as they are
      \\CORRECT_AMPLIFYING_NOT_SMALL: Values where the ABS(CTF)>epsilon are
         divided by the CTF
   */
   int method;

   /// Side Info: CTF
   FourierMask ctf;
   bool        multiple_CTFs;
   SelFile     SF_CTF; // In case of multiple CTFs
public:
   /** Empty constructor */
   CorrectPhase_Params(): epsilon(0), method(0) {}

   /** Read parameters from command line. */
   void read(int argc, char **argv);

   /** Show. */
   void show();

   /** Usage. */
   void usage();

   /** Produce side information.
       The CTF file is read. */
   void produce_side_info();

   /** Exists CTF for this file.
       Given a filename it returns the CTF file if there exists a CTF file called
       ctf-<root><num> disregarding the path. This implies that there cannot be
       two files in different directories with the same filename.

       If no CTF file is found then "" is returned */
   FileName CTF_filename(const FileName &fn);

   /** Correct phase of an image.
       An exception is thrown if the input image is not of the same size
       as the ctf or if the CTF is not real */
   void correct(vtkImageData *v) _THROW;

   /** Correct phase of a set of images.*/
   void correct(SelFile &SF);
};
//@}
#endif
#endif
