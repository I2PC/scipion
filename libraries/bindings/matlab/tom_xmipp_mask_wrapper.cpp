/*=================================================================
 *
 * tom_xmipp_mask is a wrapper to xmipp_mask
 *
 * The calling syntax is:
 *
 *		mask = tom_xmipp_mask_wrapper(size,type,mode,R1,R2,pix_width,H,sigma,omega,[rectdim_x, rectdim_y, rectdim_z],[x0 y0 z0],[smin smax]);
 *
 * Electron Tomography toolbox of the
 * Max-Planck-Institute for Biochemistry
 * Dept. Molecular Structural Biology
 * 82152 Martinsried, Germany
 * http://www.biochem.mpg.de
 *
 * and
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 *
 * created: 11/10/2007
 * by: Andreas Korinek
 *
 *=================================================================*/

/*xmipp includes */
#include "mask.h"
#include "tom_xmipp_helpers.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
    Mask maskParams;
    maskParams.fn_mask = "";
    
    /** Mask Type
     *
     * The only valid types are BINARY_CIRCULAR_MASK, BINARY_CROWN_MASK,
     * BINARY_CYLINDER_MASK, BINARY_FRAME_MASK, GAUSSIAN_MASK,
     * RAISED_COSINE_MASK, BLACKMAN_MASK, SINC_MASK, SINC_BLACKMAN_MASK,
     * READ_MASK, RAISED_CROWN_MASK, BINARY_CONE_MASK, BINARY_WEDGE_MASK
     */
    maskParams.type = (int) mxGetScalar(prhs[1]);
    
    /** Mode
     * The valid modes are INNER_MASK and OUTSIDE_MASK.
     */
    maskParams.mode = (int) mxGetScalar(prhs[2]);
    
    /** Radius 1
     * Radius for Circular and Cylinder masks and R1 for crowns and raised
     * cosine.
     */
    maskParams.R1 = (double) mxGetScalar(prhs[3]);
    
    /** Radius 2
     * R2 for crowns and raised cosine.
     */
    maskParams.R2 = (double) mxGetScalar(prhs[4]);
    
    /** Pixel width
     * For raised crowns.
     */
    maskParams.pix_width = (double) mxGetScalar(prhs[5]);
    
    /** Height
     * Height for cylinders.
     */
    maskParams.H = (double) mxGetScalar(prhs[6]);
    
    /** Sigma
     * Sigma for gaussians.
     */
    maskParams.sigma = (double) mxGetScalar(prhs[7]);
    
    /** Omega
     * Frequency for sincs
     */
    maskParams.omega = (double) mxGetScalar(prhs[8]);
    
    /** Rectangular dimensions
     */
    const int *p_dim=(const int *) mxGetData(prhs[9]);
    maskParams.Xrect = (int)p_dim[0];
    maskParams.Yrect = (int)p_dim[1];
    maskParams.Zrect = (int)p_dim[2];
    
    /** Z origin */
    const double *p_origin=mxGetPr(prhs[10]);
    maskParams.x0 = (double)p_origin[0];
    maskParams.y0 = (double)p_origin[1];
    maskParams.z0 = (double)p_origin[2];
    
    /** Minimum scale for DWT masks
     */
    const int *p_dwtscale=(const int*) mxGetData(prhs[11]);
    maskParams.smin = (int)p_dwtscale[0];
    maskParams.smax = (int)p_dwtscale[1];

    mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
    if (ndims == 2)  
    {
        MultidimArray<double> image;
        getMatrix2D(prhs[0],image);
        maskParams.generate_mask(image);

        if (maskParams.datatype() == DOUBLE_MASK)
        {
            setMatrix2D(maskParams.get_cont_mask(), plhs[0]);
        }
        else
        {            
            setMatrix2D(maskParams.get_binary_mask(), plhs[0]);
        }
    }
    else
    {
    	MultidimArray<double> volume;
        getMatrix3D(prhs[0],volume);
        maskParams.generate_mask(volume);

        if (maskParams.datatype() == DOUBLE_MASK)
        {
            setMatrix3D(maskParams.get_cont_mask(), plhs[0]);
        }
        else
        {
            setMatrix3D(maskParams.get_binary_mask(), plhs[0]);
        }
    }

}	

