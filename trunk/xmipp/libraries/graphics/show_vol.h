/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Alberto Pascual (pascual@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef SHOWVOL_H
#define SHOWVOL_H

#include "show_table.h"

#include <data/image.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QMouseEvent>
#endif

/**@defgroup ShowVolumes Show Volumes
   @ingroup GraphicsLibrary */
//@{
/** Class to show volumes.
    Example of use:
    @code
       ShowVol *showvol = new ShowVol;
       if (poll) showvol->setPoll();
       showvol->initWithFile(numRows, numCols, argv[i]);
       showvol->show();
    @endcode
*/
class ShowVol: public ShowTable
{
    Q_OBJECT;
protected:
    // Volume to be represented
    Image<double> V;
    // Axis for slistd::cing. Valid values 'X', 'Y' or 'Z'
    char slices;

    /* Init table widget.
       Current volume must have been already read */
    virtual void initTable();
    /* Init right click menu bar */
    virtual void initRightclickMenubar();

    /* Read a volume file */
    virtual void readFile(const FileName &_fn,
                          double _minGray = 0, double _maxGray = 0);
    /* Open a new file.
       The old window and volume parameters must be discarded */
    virtual void openNewFile(const FileName &);
    /* Reopen the file because its time has changed */
    virtual void reOpenFile();

    /* Produce pixmap at cell i */
    virtual void producePixmapAt(int i);
    /* Update the label with voxel value for physical coordinate (k,i,j)*/
    virtual void updateStatus(int k, int i, int j);
private slots:
    virtual void contentsMouseMoveEvent(QMouseEvent *);
public:
    /** Read volume for the first time.
        If you need to update the volume representation, use setPoll()
    and change the saved volume or use openNewFile().
    The gray limits are used for common normalization. Set them to 0
    if you don't want to use this option. */
    void initWithFile(int _numRows, int _numCols, const FileName &_fn,
                      double _minGray = 0, double _maxGray = 0);
};
//@}
#endif
