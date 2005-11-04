/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or   
 * (at your option) any later version.                                 
 *                                                                     
 * This program is distributed in the hope that it will be useful,     
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
 * GNU General Public License for more details.                        
 *                                                                     
 * You should have received a copy of the GNU General Public License   
 * along with this program; if not, write to the Free Software         
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA            
 * 02111-1307  USA                                                     
 *                                                                     
 *  All comments concerning this program package may be sent to the    
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

#ifndef _XMIPP_ERROR
   #define _XMIPP_ERROR

// Error handling ----------------------------------------------------------
/**@name Error handling
   The error handling is performed in two different ways depending on the
   configuration selected for Xmipp in the file xmippConfiguration: a
   simple error management and a second method based on C++ exceptions.
   The first method aborts the program with an error code (different for
   each error) while the second throws an exception which might be caught
   by an external routine or program.
   \\
   \\The prototype definitions in both cases are the same as they are
   based on some macros which change with the configuration. Here goes
   a programming example considering both implementations.
   \\
   \begin{verbatim}
   // Class definition
   class ReconstructingVolume:
   {
      ...
      void write(const FileName &fn) const;
      ...
   }
   
   // Class implementation
   void ReconstructingVolume::write(const FileName &fn) const {
      ...
      if (...) REPORT_ERROR(6001,"Volume too small to be stored");
      ...
   }
   
   // Use of this class in an external program
      ...
      #ifndef _NO_EXCEPTION
         try {vol_blobs.write(fn_blobs);}
         catch (Xmipp_error XE) {
            cout << XE;
            cout << "The reconstructed volume is too small to be saved in blobs\n"
                 << "So, there is no blob version of it at this iteration\n"
                 << "I go on processing" << endl;
         }
      #else
         vol_blobs.write(fn_blobs);
      #endif
      ...
   \end{verbatim}   
   You see that the routine implementation is the same in both cases but
   the external program varies from one to the other as in the exception
   case we can catch the exception and go on processing, while in the
   exit mode, the program always is aborted. If you don't put the routine
   in a try-catch structure and an exception is thrown then a core is
   generated and the program is automatically aborted.
   
   \\ See \URL[Configuration]{../../Extra_Docs/Configuration.html}
   \\ See \URL[Error Codes]{../../Extra_Docs/error_codes.html} */
//@{
/** Show message and exit.
    This macro shows the given message and exits with the error code.
    \\ Ex: if (...) EXIT_ERROR(1,"Error 1"); */
#define EXIT_ERROR(nerr,ErrormMsg) _Xmipp_error(nerr,ErrormMsg) 

void _Xmipp_error (const int nerr, const string &what);

#define REPORT_ERROR(nerr,ErrormMsg) throw Xmipp_error(nerr,ErrormMsg)
/** Exception class.
    This is the class type for the errors thrown by the routines when the
    exception handling mode is active (see Xmipp Configuration for details
    about enabling the exception handling).
    @see Configuration
    @see Error Codes*/
   class Xmipp_error {
      public:
        /** Error code.*/
        int       __errno;
        
        /** Message shown.*/
        string    msg;
        
        Xmipp_error(const int nerr, const string &what);
        friend ostream& operator << (ostream& o, Xmipp_error &XE);
   };
//@}
#endif
