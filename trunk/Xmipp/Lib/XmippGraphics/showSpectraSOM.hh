/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Alberto Pascual (pascual@cnb.uam.es)
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

#ifndef SHOWSPECTRASOM_H
#define SHOWSPECTRASOM_H

#include "showSpectra.hh"
#include "showSOM.hh"

/* COSS: This code is exactly the same as the one of ShowSOM but I haven't
         been capable of writing a class derived from ShowSOM and ShowSpectra.
*/

/*@name Show Spectra SOM. */
//@{
/** Class to show Spectra SOM. */
class ShowSpectraSOM: public ShowSpectra {
   Q_OBJECT;
public:
   // Apply geometric transformation when showing average & stddev of original 
   // images and to compute the error image
   bool        apply_geo;
protected:
   // Original data
   xmippCTVectors *Vdat;
   // Array of list of images, a list for each codevector
   vector<string> *SFcv;
   // Array of list of spectra, a list for each codevector
   vector<int>    *SFcvs;
   // Histogram of images assigned
   string  *hisAssigned;
   // Errors within group
   string  *cv_errors;
   // Info file
   QString infStr;
   
    /* Initialization.
       Sets SFcv = NULL; and then calls to ShowSel::init() */
    virtual void init();
    /* Clear everything */
    virtual void clear();
    /* Read SOM files and Dat file. */
    virtual void readFile(const FileName &_fn,
       double _minGray=0, double _maxGray=0) _THROW;
    /* Read only SOM files. */
    virtual void readSOMFiles(const FileName &_fn_root) _THROW;

    /* Extract the represented images to a list.
       The SelFile is not cleared at the beginning. */
    virtual void extractRepresented(SelFile &SF_represented);

    /* Extract the represented spectra to a list.
       The xmippCTVector is not cleared at the beginning. */
    virtual void extractRepresented(xmippCTVectors &v_represented);

    /* Initialize right click menubar */
    virtual void initRightclickMenubar();
    /* Cell label depending on labeltype */
    virtual const char* cellLabel(int i) const;
    /* Some menu items for the options */
    int mi_hisAsLabels, mi_errAsLabels;
       #define Histogram_LABEL 2
       #define Err_LABEL       3
    /* Change the labels */
    virtual void changeLabel(int _clicked_mi);

protected slots:
    /* Save assigned images in a new selfile */
    virtual void saveAssigned();
    /* Different kind of labels */
    virtual void changeLabelsToImg();
    virtual void changeLabelsToSel();
    virtual void changeLabelsToHis();
    virtual void changeLabelsToErr();

    /* Show Average and SD of the images represented */
    virtual void showRepresentedImagesStats();
    /* Show Average and SD of the spectra represented */
    virtual void showRepresentedSpectraStats();
    /* Show represented images */
    virtual void showRepresentedSel();
    /* Show represented spectra */
    virtual void showRepresentedSpectra();
    /* Show the difference between the current code vector and the
       average of the assigned images. */
    virtual void showErrorSpectrum();
    /* Show algorithm information */
    virtual void showAlgoInfo();
public:
   /** Init with SOM files.
       The given parameter is the rootname and the following files are needed:
   */
   void initWithFile(const FileName &_fn_root, const FileName &_fn_dat);

};
//@}
#endif
