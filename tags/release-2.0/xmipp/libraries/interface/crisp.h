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
/*****************************************************************************/
/* INTERACTION WITH CRISP                                                    */
/*****************************************************************************/

#ifndef _XMIPP_CRISP_HH
#define _XMIPP_CRISP_HH

#include <data/funcs.h>
#include <data/volume.h>

/**@defgroup Crisp Crisp
   @ingroup InterfaceLibrary */
//@{
/** Crisp volumes */
class CrispVolume
{
public:
    /// Flags
    short  flags;
    /// real dimension of density map along A axis in steps
    short  Adata;
    /// real dimension of density map along B axis in steps
    short  Bdata;
    /// real dimension of density map along C axis in steps
    short  Cdata;
    /// number of steps along A axis in file
    short  Afile;
    /// number of steps along B axis in file
    short  Bfile;
    /// number of steps along C axis in file
    short  Cfile;
    /// cell size along A axis in angstroms
    short  Asize;
    /// cell size along B axis in angstroms
    short  Bsize;
    /// cell size along C axis in angstroms
    short  Csize;
    /// angle between A and B axes
    short  Gamma;
    /// reserved
    short  dummy0[21];
    /// name of density map (for reference)
    char Name[32];
    /// reserved
    short  dummy1[208];

    /// Filename
    FileName name;
    /// Volume
    Volume V;

public:
    /// Read from file
    void read(const FileName &fn);
    /// Show header
    friend ostream & operator << (ostream &out, const CrispVolume &cv);
    /// Write as spider
    void write_as_spider(const FileName &fn);
};

//@}
#endif
