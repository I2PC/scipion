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

#include "adjust_surface.h"
#include "surface.h"

#include <interface/inria.h>
#include <data/args.h>

/* Read from command line ================================================== */
void Prog_Adjust_Surface_Parameters::read(int argc, char **argv)
{
    tell = 0;
    fn_in_surface  = get_param(argc, argv, "-i");
    fn_out_surface = get_param(argc, argv, "-o", "");
    fn_vol         = get_param(argc, argv, "-vol");
    exhaustive     = check_param(argc, argv, "-exhaustive");
    apply          = check_param(argc, argv, "-apply");
    given_ztop     = check_param(argc, argv, "-ztop");
    if (given_ztop)
    {
        int i = position_param(argc, argv, "-ztop");
        if (i + 3 >= argc) REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::read:"
                                            " Not enough parameters behind -ztop");
        ztop0 = AtoI(argv[i+1]);
        ztopF = AtoI(argv[i+2]);
        ztop_step = AtoI(argv[i+3]);
    }
    given_zbottom  = check_param(argc, argv, "-zbottom");
    if (given_zbottom)
    {
        int i = position_param(argc, argv, "-zbottom");
        if (i + 3 >= argc) REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::read:"
                                            " Not enough parameters behind -zbottom");
        zbottom0 = AtoI(argv[i+1]);
        zbottomF = AtoI(argv[i+2]);
        zbottom_step = AtoI(argv[i+3]);
    }
    if (check_param(argc, argv, "-ang"))
    {
        int i = position_param(argc, argv, "-ang");
        if (i + 3 >= argc) REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::read:"
                                            " Not enough parameters behind -ang");
        angle0 = AtoF(argv[i+1]);
        angleF = AtoF(argv[i+2]);
        angle_step = AtoF(argv[i+3]);
    }
    else
    {
        angle0 = angleF = 0;
        angle_step = 1;
    }
    if (check_param(argc, argv, "-shift"))
    {
        int i = position_param(argc, argv, "-shift");
        if (i + 4 >= argc) REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::read:"
                                            " Not enough parameters behind -shift");
        shiftX0 = AtoF(argv[i+1]);
        shiftY0 = AtoF(argv[i+2]);
        shift_dist = AtoF(argv[i+3]);
        shift_step = AtoI(argv[i+4]);
    }
    else
    {
        shiftX0 = shiftY0 = shift_dist = 0;
        shift_step = 1;
    }
    if (check_param(argc, argv, "-scaleX"))
    {
        int i = position_param(argc, argv, "-scaleX");
        if (i + 3 >= argc) REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::read:"
                                            " Not enough parameters behind -scaleX");
        scaleX0 = AtoF(argv[i+1]);
        scaleXF = AtoF(argv[i+2]);
        scaleX_step = AtoF(argv[i+3]);
    }
    else
    {
        scaleX0 = scaleXF = 1;
        scaleX_step = 1;
    }
    if (check_param(argc, argv, "-scaleY"))
    {
        int i = position_param(argc, argv, "-scaleY");
        if (i + 3 >= argc) REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::read:"
                                            " Not enough parameters behind -scaleY");
        scaleY0 = AtoF(argv[i+1]);
        scaleYF = AtoF(argv[i+2]);
        scaleY_step = AtoF(argv[i+3]);
    }
    else
    {
        scaleY0 = scaleYF = 1;
        scaleY_step = 1;
    }
    if (check_param(argc, argv, "-bottom_surface")) direction = BOTTOM2TOP;
    else direction = TOP2BOTTOM;
    if (check_param(argc, argv, "-corr_2D"))      tell |= CORR_2D;
    if (check_param(argc, argv, "-corr_grad"))    tell |= CORR_GRAD;
    if (check_param(argc, argv, "-manual_order")) tell |= MANUAL_ORDER;
    phantom = check_param(argc, argv, "-phantom");
}

/* Usage =================================================================== */
void Prog_Adjust_Surface_Parameters::usage() const
{
    cout << "\nUsage:\n";
    cout << "adjust_surface\n"
    << "   -i <Input surface>             : Xmipp image\n"
    << "  [-o <Output surface>]           : Output scaled surface\n"
    << "   -vol <volume to fit>           : volume to be fitted\n"
    << "  [-exhaustive]                   : Do exhaustive search\n"
    << "  [-apply]                        : Apply best combination to surface\n"
    << "  [-ztop <ztop0> <ztopF> <step>]  : Range to search top value\n"
    << "  [-zbottom <zbottom0> <zbottomF> <step>]: Range to search bottom value\n"
    << "  [-ang <ang0=0> <angF=0> <step=1>]: Range to search angle\n"
    << "  [-scaleX <sc0=0> <scF=0> <step=1>]: Range to search scaleX\n"
    << "  [-scaleY <sc0=0> <scF=0> <step=1>]: Range to search scaleY\n"
    << "  [-shift <x0=0><y0=0><dist=0><step=1>]: maximum allowed shift\n"
    << "  [-bottom_surface]               : by default, a top surface\n"
    << "                                    is adjusted\n"
    << "  [-corr_2D]                      : Apply 2D correlation method\n"
    << "  [-corr_grad]                    : Apply Gradient correlation method\n"
    << "  [-phantom]                      : Surface is coming from a phantom\n"
    << "  [-manual_order]                 : The heights are given manually\n"
    ;
}

/* Produce side information ================================================ */
void Prog_Adjust_Surface_Parameters::produce_Side_Info()
{
    surface.read(fn_in_surface);
    surface().setXmippOrigin();
    V.read(fn_vol);
    V().setXmippOrigin();
    if (!given_ztop)
    {
        ztop0    = STARTINGZ(V());
        ztopF    = FINISHINGZ(V());
    }
    if (!given_zbottom)
    {
        zbottom0 = STARTINGZ(V());
        zbottomF = FINISHINGZ(V());
    }

    if (ztop0 < STARTINGZ(V()) || ztopF > FINISHINGZ(V()) || ztop0 > ztopF)
        REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::produce_Side_Info: "
                     "ztop out of volume range");
    if (zbottom0 < STARTINGZ(V()) || zbottomF > FINISHINGZ(V()) ||
        zbottom0 > zbottomF)
        REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::produce_Side_Info: "
                     "zbottom out of volume range");
    if (tell & CORR_GRAD)
    {
#ifdef _HAVE_INRIA
        compute_gradient(V(), V_grad);
#else
        REPORT_ERROR(1, "Prog_Adjust_Surface_Parameters::produce_Side_Info: "
                     "You don't have INRIA library, so this option is unavailable");
#endif
    }

    // Compute surface range
    min_val = max_val = DIRECT_IMGPIXEL(surface, 0, 0);
#define prmsij IMGPIXEL(surface,i,j)
#define ZDIM   ZSIZE(V())
    FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(surface))
    {
        if (!(prmsij == ZDIM && phantom))
        {
            if (prmsij < min_val) min_val = prmsij;
            if (prmsij > max_val) max_val = prmsij;
        }
    }

    // Compute valid shift range
    shift_mask.resize(2*CEIL(shift_dist) + 1, 2*CEIL(shift_dist) + 1);
    shift_mask.setXmippOrigin();
    STARTINGX(shift_mask) += ROUND(shiftX0);
    STARTINGY(shift_mask) += ROUND(shiftY0);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(shift_mask)
    shift_mask(i, j) =
        ((i - shiftY0) * (i - shiftY0) + (j - shiftX0) * (j - shiftX0) <=
         shift_dist * shift_dist);
}
#undef prmsij
#undef ZDIM

/* Create surface mask ===================================================== */
void create_surface_mask(const Image *surf, const Volume *V, Volume *Vsurf,
                         int direction)
{
    (*Vsurf)().initZeros((*V)());
    switch (direction)
    {
    case TOP2BOTTOM:
        create_surface_mask(surf, NULL, ZSIZE((*V)()), Vsurf);
        break;
    case BOTTOM2TOP:
        create_surface_mask(NULL, surf, ZSIZE((*V)()), Vsurf);
        break;
    default:
        REPORT_ERROR(1, "create_surface_mask: Unknown direction");
    }

    // Invert mask
    FOR_ALL_ELEMENTS_IN_MATRIX3D(VOLMATRIX(*Vsurf))
    VOLVOXEL(*Vsurf, k, i, j) = 1 - VOLVOXEL(*Vsurf, k, i, j);
}

/* Correlation surface - volume ============================================ */
// 3D correlation mode
#define VOL VOLMATRIX(*V)
double correlate_surface_and_volume_3D(const Image *surf, const Volume *V,
                                       Volume *Vsurf, int ktop, int kbottom, int direction, int tell)
{
    double retval = 0;
    int   N = 0;

    // Create surface mask.
    create_surface_mask(surf, V, Vsurf, direction);

    // Compute mean values within the interest area
    if (direction == TOP2BOTTOM) kbottom--;
    else                       ktop++;
    double V_mean = 0, Vsurf_mean = 0;
    for (int k = ktop; k <= kbottom; k++)
        for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
            for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
            {
                N++;
                V_mean     += VOLVOXEL(*V, k, i, j);
                Vsurf_mean += VOLVOXEL(*Vsurf, k, i, j);
            }
    if (N != 0)
    {
        V_mean /= N ;
        Vsurf_mean /= N;
    }
    else      return 0;

    // Compute Power of both functions
    double V_power = 0, Vsurf_power = 0;
    for (int k = ktop; k <= kbottom; k++)
        for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
            for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
            {
                V_power += (VOLVOXEL(*V, k, i, j) - V_mean) *
                           (VOLVOXEL(*V, k, i, j) - V_mean);
                Vsurf_power += (VOLVOXEL(*Vsurf, k, i, j) - Vsurf_mean) *
                               (VOLVOXEL(*Vsurf, k, i, j) - Vsurf_mean);
            }
    V_power = sqrt(V_power);
    Vsurf_power = sqrt(Vsurf_power);

    // Correlate
    for (int k = ktop; k <= kbottom; k++)
        for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
            for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
            {
                double add_val = (VOLVOXEL(*V, k, i, j) - V_mean) *
                                 (VOLVOXEL(*Vsurf, k, i, j) - Vsurf_mean);
                retval += add_val;
            }

    if (tell & MANUAL_ORDER)
    {
        VolumeXmipp save;
        save().initZeros((*Vsurf)());
        for (int k = ktop; k <= kbottom; k++)
            for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
                for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
                    VOLVOXEL(save, k, i, j) = VOLVOXEL(*Vsurf, k, i, j);
        save.write("PPPsurface.vol");
        for (int k = ktop; k <= kbottom; k++)
            for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
                for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
                    VOLVOXEL(save, k, i, j) =
                        SGN(VOLVOXEL(*Vsurf, k, i, j) - Vsurf_mean) !=
                        SGN(VOLVOXEL(*V, k, i, j) - V_mean);
        save.write("PPPsign.vol");
    }

    // Return
    return retval / (V_power*Vsurf_power);
}

/* Correlation surface - volume ============================================ */
// 2D correlation mode (Eva's method)
double correlate_surface_and_volume_2D(const Image *surf, const Volume *V,
                                       Volume *Vsurf, int ktop, int kbottom, int direction, int tell)
{
    double retval = 0;
    int   N = XSIZE(IMGMATRIX(*surf)) * YSIZE(IMGMATRIX(*surf));

    // Project Volume
    (*Vsurf)().initZeros(1, YSIZE(VOL), XSIZE(VOL));
    STARTINGZ((*Vsurf)()) = 0;
    STARTINGY((*Vsurf)()) = STARTINGY(VOL);
    STARTINGX((*Vsurf)()) = STARTINGX(VOL);

    for (int k = ktop; k <= kbottom; k++)
        for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
            for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
                VOLVOXEL(*Vsurf, 0, i, j) += VOLVOXEL(*V, k, i, j);

    // Compute mean values within the interest area
    double V_mean = 0, Vsurf_mean = 0;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(*surf))
    {
        if (direction == BOTTOM2TOP) V_mean += IMGPIXEL(*surf, i, j);
        else                       V_mean -= IMGPIXEL(*surf, i, j);
        Vsurf_mean += VOLVOXEL(*Vsurf, 0, i, j);
    }
    V_mean /= N;
    Vsurf_mean /= N;

    // Compute Power of both functions
    double V_power = 0, Vsurf_power = 0;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(*surf))
    {
        if (direction == BOTTOM2TOP)
            V_power += (IMGPIXEL(*surf, i, j) - V_mean) * (IMGPIXEL(*surf, i, j) - V_mean);
        else
            V_power += (-IMGPIXEL(*surf, i, j) - V_mean) * (-IMGPIXEL(*surf, i, j) - V_mean);
        Vsurf_power += (VOLVOXEL(*Vsurf, 0, i, j) - Vsurf_mean) *
                       (VOLVOXEL(*Vsurf, 0, i, j) - Vsurf_mean);
    }
    V_power = sqrt(V_power);
    Vsurf_power = sqrt(Vsurf_power);

    // Correlate
    for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
        for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
            if (direction == BOTTOM2TOP)
                retval += (IMGPIXEL(*surf, i, j) - V_mean) *
                          (VOLVOXEL(*Vsurf, 0, i, j) - Vsurf_mean);
            else
                retval += (-IMGPIXEL(*surf, i, j) - V_mean) *
                          (VOLVOXEL(*Vsurf, 0, i, j) - Vsurf_mean);

    if (tell & MANUAL_ORDER)
    {
        ImageXmipp save;
        save() = (*Vsurf)().getSlice(0);
        save.write("PPPsurface.img");
        for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
            for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
                if (direction == BOTTOM2TOP)
                    IMGPIXEL(save, i, j) =
                        SGN(VOLVOXEL(*Vsurf, 0, i, j) - Vsurf_mean) !=
                        SGN(IMGPIXEL(*surf, i, j) - V_mean);
                else
                    IMGPIXEL(save, i, j) =
                        SGN(VOLVOXEL(*Vsurf, 0, i, j) - Vsurf_mean) !=
                        SGN(-IMGPIXEL(*surf, i, j) - V_mean);
        save.write("PPPsign.img");
    }
    return retval / (V_power*Vsurf_power);
}

/* Correlation surface - Gradient method =================================== */
//#define DEBUG
double correlate_surface_and_volume_gradients(const Image *surf,
        const Volume *V, Volume *Vsurf,
        const Vectorial_matrix3D &V_grad, Vectorial_matrix3D & Vsurf_grad,
        int ktop, int kbottom, int direction, int tell)
{
    double retval = 0;
#ifdef _HAVE_INRIA
    // Create surface mask
    create_surface_mask(surf, V, Vsurf, direction);
    (*Vsurf)().window(MAX(STARTINGZ((*Vsurf)()), ktop - 5),
                      STARTINGY((*Vsurf)()), STARTINGX((*Vsurf)()),
                      MIN(FINISHINGZ((*Vsurf)()), kbottom + 5),
                      FINISHINGY((*Vsurf)()), FINISHINGX((*Vsurf)()));

    // Create gradient of mask
    compute_gradient((*Vsurf)(), Vsurf_grad);

    if (tell&MANUAL_ORDER)
    {
        VolumeXmipp aux;
        aux = *Vsurf;
        aux.write("PPPVsurf.vol");
        Vsurf_grad.module(aux());
        aux.write("PPPVsurf_grad_modules.vol");
        Vsurf_grad.write("PPPVsurf_grad.vol");
        V_grad.module(aux());
        aux.write("PPPV_grad_modules.vol");
        V_grad.write("PPPV_grad.vol");
    }

    // Correlate
    Matrix1D<double> grad_in_V(3), grad_in_Vsurf(3);
    for (int i = STARTINGY(VOL); i <= FINISHINGY(VOL); i++)
        for (int j = STARTINGX(VOL); j <= FINISHINGX(VOL); j++)
        {
            bool surface_found;
            surface_found = false;
            int k;
            for (k = ktop; k <= kbottom; k++)
                if (VOLVOXEL(*Vsurf, k, i, j) != 0)
                {
                    surface_found = true;
                    break;
                }
            if (surface_found)
            {
                V_grad.vector_at(k, i, j, grad_in_V);
                Vsurf_grad.vector_at(k, i, j, grad_in_Vsurf);
                //grad_in_Vsurf.selfNormalize();
                retval += ABS(dotProduct(grad_in_V, grad_in_Vsurf));
#ifdef DEBUG
                cout << "(" << k << "," << i << "," << j << ") in V="
                << grad_in_V.transpose() << " in surf "
                << grad_in_Vsurf.transpose() << endl;
#endif
            }
        }
#endif
    return retval;
}
#undef DEBUG
#undef VOL

/* Global variables used for communicating with the main program =========== */
Prog_Adjust_Surface_Parameters *gasprm;

/* Evaluate correlation ==================================================== */
// p[1]=shiftx
// p[2]=shifty
// p[3]=ang
// p[4]=ktop
// p[5]=kbottom
#define prmsij gasprm->surface(i,j)
#define    sij gasprm->wsurface(i,j)
#define ZDIM   ZSIZE(gasprm->V())
double eval_surface(double *p)
{
    // Check valid ranges
    Matrix1D<double> shift(2);
    VECTOR_R2(shift, ROUND(p[1]), ROUND(p[2]));
    if (gasprm->shift_mask.outside(shift))               return 1e10;
    if (!gasprm->shift_mask(shift))                      return 1e10;
    if (p[3] < gasprm->angle0   || p[3] > gasprm->angleF) return 1e10;
    if (p[4] < gasprm->ztop0    || p[4] > gasprm->ztopF) return 1e10;
    if (p[5] < gasprm->zbottom0 || p[5] > gasprm->zbottomF) return 1e10;
    if (p[4] + 1 >= p[5])                                    return 1e10;
    if (p[6] < gasprm->scaleX0  || p[6] > gasprm->scaleXF) return 1e10;
    if (p[7] < gasprm->scaleY0  || p[7] > gasprm->scaleYF) return 1e10;

    // Adjust surface range
    gasprm->wsurface().resize(gasprm->surface());
    double a = (p[5] - p[4]) / (gasprm->max_val - gasprm->min_val);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(gasprm->surface))
    if (!(prmsij == ZDIM && gasprm->phantom))
        sij = p[4] + a * (prmsij - gasprm->min_val);
    else sij = ZDIM;
    if (p[3] != 0) gasprm->wsurface().rotate(p[3], WRAP);
    if (XX(shift) != 0 || YY(shift) != 0)
        gasprm->wsurface().translate(shift, WRAP);
    if (p[6] != 0 || p[7] != 0)
    {
        matrix2D<double> scale_m;
        scale_m.init_identity(3);
        matrix2D<double> aux;
        scale_m(0, 0) = p[6];
        scale_m(1, 1) = p[7];
        apply_geom(aux, scale_m, gasprm->wsurface(), IS_NOT_INV, WRAP);
        gasprm->wsurface() = aux;
    }

    // Correlate
    double corr;
    if (gasprm->tell & CORR_2D)
        corr = correlate_surface_and_volume_2D((Image *) & (gasprm->wsurface),
                                               (Volume *) & (gasprm->V), (Volume *) & (gasprm->Vsurf),
                                               ROUND(p[4]), ROUND(p[5]), gasprm->direction, gasprm->tell);
    else if (gasprm->tell & CORR_GRAD)
        corr = correlate_surface_and_volume_gradients((Image *) & (gasprm->wsurface),
                (Volume *) & (gasprm->V), (Volume *) & (gasprm->Vsurf),
                gasprm->V_grad, gasprm->Vsurf_grad,
                ROUND(p[4]), ROUND(p[5]), gasprm->direction, gasprm->tell);
    else
        corr = correlate_surface_and_volume_3D((Image *) & (gasprm->wsurface),
                                               (Volume *) & (gasprm->V), (Volume *) & (gasprm->Vsurf),
                                               ROUND(p[4]), ROUND(p[5]), gasprm->direction, gasprm->tell);

    // Return
    return -corr;
}
#undef prmsij
#undef    sij
#undef ZDIM

/* Main routine ============================================================ */
#define prmsij prm.surface(i,j)
#define    sij prm.wsurface(i,j)
#define ZDIM   ZSIZE(prm.V())
#define VOL VOLMATRIX(prm.V)
void ROUT_adjust_surface(Prog_Adjust_Surface_Parameters &prm)
{
    int best_ktop, best_kbottom;
    double best_corr = 0, corr, best_shiftx, best_shifty, best_ang, best_scaleX,
                       best_scaleY;

    gasprm = &prm;
    prm.p.resize(7);

    // Search for all possible combinations .................................
    cerr << "Correlating surface and volume ...\n";
    if (prm.exhaustive)
    {
        double *p_aux = prm.p.adaptForNumericalRecipes();
        cout << "# shiftx shifty angle ktop kbottom corr\n";
        int act_corr = 0;
        for (prm.p(1) = STARTINGY(prm.shift_mask); prm.p(1) <= FINISHINGY(prm.shift_mask); prm.p(1) += prm.shift_step)
            for (prm.p(0) = STARTINGX(prm.shift_mask); prm.p(0) <= FINISHINGX(prm.shift_mask); prm.p(0) += prm.shift_step)
            {
                if (!prm.shift_mask((int)prm.p(1), (int)prm.p(0))) continue;
                for (prm.p(2) = prm.angle0; prm.p(2) <= prm.angleF; prm.p(2) += prm.angle_step)
                    for (prm.p(3) = prm.ztop0; prm.p(3) <= prm.ztopF; prm.p(3) += prm.ztop_step)
                        for (prm.p(4) = prm.zbottom0; prm.p(4) <= prm.zbottomF; prm.p(4) += prm.zbottom_step)
                            for (prm.p(5) = prm.scaleX0; prm.p(5) <= prm.scaleXF; prm.p(5) += prm.scaleX_step)
                                for (prm.p(6) = prm.scaleY0; prm.p(6) <= prm.scaleYF; prm.p(6) += prm.scaleY_step)
                                {
                                    if (prm.tell & MANUAL_ORDER)
                                    {
                                        cout << "shiftx =";
                                        cin >> prm.p(0);
                                        cout << "shifty =";
                                        cin >> prm.p(1);
                                        cout << "angle =";
                                        cin >> prm.p(2);
                                        cout << "ktop (" << prm.ztop0 << "..." << prm.ztopF << ")=";
                                        cin >> prm.p(3);
                                        cout << "kbottom (" << prm.zbottom0 << "..." << prm.zbottomF << ")=";
                                        cin >> prm.p(4);
                                        cout << "scaleX =";
                                        cin >> prm.p(5);
                                        cout << "scaleY =";
                                        cin >> prm.p(6);
                                    }

                                    corr = eval_surface(p_aux);
                                    // Best correlation?
                                    if (corr < best_corr)
                                    {
                                        best_corr = corr;
                                        best_ktop = ROUND(prm.p(3));
                                        best_kbottom = ROUND(prm.p(4));
                                        best_ang = prm.p(2);
                                        best_shiftx = prm.p(0);
                                        best_shifty = prm.p(1);
                                        best_scaleX = prm.p(5);
                                        best_scaleY = prm.p(6);
                                        cout << "********    ";
                                    }

                                    cout << prm.p(0) << " " << prm.p(1) << " " << prm.p(2) << " "
                                    << prm.p(3) << " " << prm.p(4) << " " << prm.p(5) << " "
                                    << prm.p(6) << " " << corr << endl;
                                    if (prm.tell & MANUAL_ORDER)
                                    {
                                        cout << "Press any key\n";
                                        char c;
                                        cin >> c;
                                    }
                                }
            }
        prm.p.killAdaptationForNumericalRecipes(p_aux);
    }
    else
    {
        // Search with Powell ...................................................
        prm.p.initZeros();
        prm.p(0) = prm.shiftX0;
        prm.p(1) = prm.shiftY0;
        prm.p(2) = prm.angle0;
        prm.p(3) = prm.ztop0;
        prm.p(4) = prm.zbottomF;
        prm.p(5) = (prm.scaleX0 + prm.scaleXF) / 2;
        prm.p(6) = (prm.scaleY0 + prm.scaleYF) / 2;
        Matrix1D<double> steps(7);
        steps.init_constant(1);
        int iter = 0;
        powellOptimizer(prm.p, 1, 7, &eval_surface,
                         0.01, corr, iter, steps, true);
        best_corr = corr;
        best_ktop = ROUND(prm.p(3));
        best_kbottom = ROUND(prm.p(4));
        best_ang = prm.p(2);
        best_shiftx = prm.p(0);
        best_shifty = prm.p(1);
        best_scaleX = prm.p(5);
        best_scaleY = prm.p(6);
    }

    // Apply best correlation
    cout << "Best correlation found for\n"
    << "  shiftx = " << best_shiftx  << endl
    << "  shifty = " << best_shifty  << endl
    << "  ang    = " << best_ang     << endl
    << "  ktop   = " << best_ktop    << endl
    << "  kbottom= " << best_kbottom << endl
    << "  scaleX = " << best_scaleX  << endl
    << "  scaleY = " << best_scaleY  << endl;
    if (prm.apply)
    {
        prm.wsurface().resize(prm.surface());
        double a = (best_kbottom - best_ktop) / (prm.max_val - prm.min_val);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(IMGMATRIX(prm.wsurface))
        if (!(prmsij == ZDIM && prm.phantom)) sij = best_ktop + a * (prmsij - prm.min_val);
        else                                sij = ZDIM;
        if (best_ang != 0) prm.wsurface().rotate(best_ang, WRAP);
        if (best_shiftx != 0 || best_shifty != 0)
            prm.wsurface().translate(vectorR2(best_shiftx, best_shifty), WRAP);
        if (best_scaleX != 0 || best_scaleY != 0)
        {
            matrix2D<double> scale_m;
            scale_m.init_identity(3);
            matrix2D<double> aux;
            scale_m(0, 0) = best_scaleX;
            scale_m(1, 1) = best_scaleY;
            apply_geom(aux, scale_m, prm.wsurface(), IS_NOT_INV, WRAP);
            prm.wsurface() = aux;
        }

        // Save results
        if (prm.fn_out_surface == "") prm.wsurface.write(prm.surface.name());
        else                        prm.wsurface.write(prm.fn_out_surface);
    }
}
