/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (1999)
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

#include <data/volume_segment.h>

void Usage(const Prog_segment_prm &prm);

int main(int argc, char **argv)
{
    Prog_segment_prm prm;

    // Read arguments
    try
    {
        prm.read(argc, argv);
        cout << prm;
    }
    catch (Xmipp_error Xe)
    {
        cout << Xe;
        Usage(prm);
        exit(1);
    }

    // Really segment
    try
    {
        VolumeXmipp mask;
        prm.produce_side_info();
        prm.segment(mask);
    }
    catch (Xmipp_error Xe)
    {
        cout << Xe;
    }
}

/* Usage ------------------------------------------------------------------- */
void Usage(const Prog_segment_prm &prm)
{
    cerr << "Usage: segment [Parameters]\n";
    prm.usage();
}

