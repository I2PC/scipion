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

#ifndef SHOWSEL_H
#define SHOWSEL_H

#include "showTable.hh"
#include <XmippData/xmippSelFiles.hh>

/**@name Show Selection File. */
//@{
/** Class to show selfiles.
    Example of use:
    \begin{verbatim}
       ShowSel *showsel = new ShowSel;
       showsel->initWithFile(numRows, numCols, argv[i]);
       showsel->show();
    \end{verbatim}
*/
class ShowSel: public ShowTable {
    Q_OBJECT;
protected:
    // Filenames within the Selfile
    FileName   *imgnames;
    // Selstatus of those files
    bool       *selstatus;
    // Kind of label: there are two types
    // SFLabel_LABEL  --> use the Selfile status as label
    // Filename_LABEL --> use the Filename       as label
    int         labeltype;
		// Valid types
		#define SFLabel_LABEL  0
		#define Filename_LABEL 1

    // Some menu items to enable/disable them
    int mi_Individual_norm, mi_Global_norm,
        mi_showLabel,       mi_hideLabel,
        mi_selAsLabels,     mi_imgAsLabels;

    // Reimplementation of the insert_in queue to account for
    // the deletion of unused cells
    void insert_content_in_queue(int i);

    /* Initialization.
       Sets labeltype = SFLabel_LABEL; imgnames, selstatus = NULL; and then
       calls to ShowTable::init() */
    virtual void init();
    /* Clear everything */
    virtual void clear();
    /* Configures the QTable widget */
    virtual void initTable();
    /* Compute the minimum and maximum pixel within selfile */
    virtual void compute_global_normalization_params();
    /* Write SelFile, opens a GUI for the file selection and checks if the
       output filename exists. If it does, it asks for confirmation.*/
    virtual void writeSelFile(SelFile &_SF);
    /* Cell label depending on labeltype */
    virtual const char* cellLabel(int i) const;
    /* Initialize right click menubar */
    virtual void initRightclickMenubar();
    /* Produce Pixmap for cell i */
    virtual void producePixmapAt(int i);
    /* Open a new file. Old parameters must be discarded */
            void openNewFile (const FileName &);
    /* Read a Selfile */
    void readFile(const FileName &_fn) _THROW;
protected slots:
    // These slots are related with the right click menubar ---------------- */
    /* Unselect all cells */
    virtual void unSelectAll();
    /* Select all cells */
    virtual void SelectAll();
    /* Save current SelFile with marked projections as discarded */
    virtual void saveSelFileDiscarded();
    /* Save current SelFile with marked projections as active, and the
       rest as discarded */
    virtual void saveSelFileActive();
    /* Save marked projections in a new Selfile as Active */
    virtual void saveSelFileNew();
    /* Show Average and Stddev images */
    virtual void showStats();
    /* Show statistics of the SelFile */
    virtual void showSelStats();
    /* Change the global/individual normalization status */
    virtual void changeNormalize();
    /* Change the show/hide labels status */
    virtual void changeShowLabels();
    /* Change the kind of labels.*/
    virtual void changeLabels();
public:
    /** Read selfile for the first time.
        If you need to update the volume representation, use openNewFile(()*/
    virtual void initWithFile(int _numRows, int _numCols,
       const FileName &_fn);
};
//@}
#endif
