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

#include <data/phantom.h>
#include "volume_labelling.h"

/* Label a volume from a description file ---------------------------------- */
int Labels_from_description(Phantom &phantom, Volume *vol_label)
{
    vol_label->adapt_to_size(phantom.zdim, phantom.ydim, phantom.xdim);
    for (int f = 0; f < phantom.FeatNo(); f++)
        phantom.VF[f]->draw_in(vol_label, 1, f + 1);
    return phantom.FeatNo();
}

/* Label a volume from growing voxels -------------------------------------- */
// Esta funci�n hay que acabarla
int Labels_from_growing(Volume *phantom, Volume *vol_label)
{
    vol_label->adapt_to_size((*phantom)().zdim, (*phantom)().ydim,
                             (*phantom)().xdim);
    return 0;
}

#ifdef NUNCADEFINIDO
/* ------------------------------------------------------------------------- */
/* Label by growing voxels                                                   */
/* ------------------------------------------------------------------------- */
void Labels_from_growing(Matrix3D<float> *vol, Matrix3D<short> *vol_label,
                         short max_feat, short *num_feat)
{

    int     i, j, k;                      /* X,Y,Z indices inside volume        */
    short   next_label = 1;               /* Next label to use when labelling   */

    for (k = 0; k < vol->zdim; k++)
        for (j = 0; j < vol->ydim; j++)
            for (i = 0; i < vol->xdim; i++)
            {
                /* If voxel is background or already labelled skip */
                if (vol->m[k][j][i] == 0.0)     continue;
                if (vol_label->m[k][j][i] != 0) continue;

                /* Else, assign a new label to this pixel and grow it */
                if (next_label < max_feat)
                {
                    grow_voxel(vol_label, vol, k, j, i, next_label, 0);
                    next_label++;
                }
                else
                {
                    printf("Overriden maximum number of features=%d\n", max_feat);
                    exit(1);
                }

            } /* for every voxel */
    *num_feat = next_label - 1;
}

/* ------------------------------------------------------------------------- */
/* Grow Voxel                                                                */
/* ------------------------------------------------------------------------- */
/* This function takes a voxel, labels it and studies the labels of the
   voxels in its surroundings. A prerrequisite to this routine is that
   the entering voxel is greater than 0 (a real density) and it has not
   been already labelled */

/* *** LA FUNCION ESTA HECHA COMO RECURSIVA, SE GANARIA MUCHO EN VELOCIDAD
   Y ESPACIO DE CALCULO SI SE HICIERA CON UNA COLA *** */

void grow_voxel(Matrix3D<short> *vol_label, Matrix3D<float> *vol_phantom,
                int k0, int j0, int i0, short label, short level)
{

    /* Variables ------------------------------------------------------------ */
    int i, j, k;           /* indexes for matrices */
    int ileft, iright;    /* limits for x-axis */
    int jleft, jright;    /* limits for y-axis */
    int kupper, klower;   /* limits for z-axis */

    /* Assign label ...................................................... */
    (*vol_label).m[k0][j0][i0] = label;
#ifdef DEBUG_NEIGHBOURS
    {
        int n;
        for (n = 0; n < level; n++) printf(" ");
        printf("%3d Asignado [%d][%d][%d] -> %d\n", level, k0, j0, i0, label);
    }
#endif

    /* Choose limits for growing ......................................... */
    ileft  = max(0, i0 - 1);
    iright = min(i0 + 1, (*vol_phantom).xdim - 1);
    jleft  = max(0, j0 - 1);
    jright = min(j0 + 1, (*vol_phantom).ydim - 1);
    klower = max(0, k0 - 1);
    kupper = min(k0 + 1, (*vol_phantom).zdim - 1);

    /* Grow every voxel in that cube ..................................... */
    for (k = klower; k <= kupper; k++)
        for (j = jleft; j <= jright; j++)
            for (i = ileft; i <= iright; i++)
            {
                /* If voxel is background or already labelled skip */
                if ((*vol_phantom).m[k][j][i] == 0.0) continue;
                if ((*vol_label).m[k][j][i] != 0)     continue;
                grow_voxel(vol_label, vol_phantom, k, j, i, label, level + 1);
            }

}

#endif
