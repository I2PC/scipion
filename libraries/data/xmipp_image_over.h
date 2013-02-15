/***************************************************************************
 *
 * Authors:  Alberto Pascual Montano (pascual@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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

#ifndef IMAGEOVER_H
#define IMAGEOVER_H

#include "xmipp_image.h"
/** @ingroup Images
 *
 * The oversampled images are images which allow a more accurate treatment of
 * information by oversampling all pixels. The idea of this class is to have an
 * image with a logical size smaller than its physical one, for this reason you
 * could use non integer logical positions to access to different samples of the
 * "same" pixel. Let's set an example, blob footprints are of size 5x5,
 * for instance, with the center at physical position 3x3. It's convenient to
 * have this footprint defined between -2 and 2 with the origin at the center of
 * the image. But it's also convenient to have more than 25 values of the
 * footprint, this could be done by sampling finer each pixel, say 51 samples
 * per pixel, this would result in an image of 255x255 samples (with the center
 * at [127][127]). But we don't to access to the pixel at position (-120,-120)
 * but to the one at logical position (-2.35,-2.35) which we know is close to
 * the border. The oversampled image class allows you this kind of image access.
 * You have to say from where to where your image is defined for each axis (in
 * this case -2...2 for both) and which are the sampling rates for both axes
 * (51 for both). This is the initialisation process. From now on you can work
 * with your image in the way formerly described.
 *
 * Pay attention to two points:
 *
 * * The sampling rate must be an odd number, this is because so the logical
 * oversampled pixels (-2,-2), (-2, -1) ... are exactly on one cell of the
 * underlying 2D matrix.
 *
 * * As the oversampled pixels are centered with respect to the "superpixels"
 * defined by the 2D matrix, the extent of the oversampled image is a little
 * larger than from (-2,-2) to (2,2), ie, from (-2.5,-2.5) to (2.5,2.5)
 *
 * Oversampled images are normal images except for pixel access which can be
 * done in two fashions: either using the normal Image procedures (in this case,
 * you are restricted to the macro IMGPIXEL), or using the fractional indexes
 * defined for this class. Have a look at the following example of how to
 * compute the blob footprint with this kind of images. Pay special attention to
 * the pixel access at the end of the loop. Notice that x and y moves at values
 * between -2 and 2, which keep logical meaning of what we are doing while all
 * the burden of computing where to put these values at the image is done by the
 * library.
 *
 * @code
 * void footprint_blob
 * (
 *     ImageOver& blobprint, // blob foot_print table
 *     const struct blobtype& blob, // blob description
 *     int istep, // number of foot-print samples per one sample
 *                // on projection plane in u,v directions
 *     int   normalise) // set to 1 if you want normalise. Usually
 *                      // you may omit it and no normalisation is performed
 * {
 *     // Resize output image and redefine the origin of it
 *     int footmax = CEIL(blob.radius);
 *     blobprint.init(-footmax, footmax, istep, -footmax, footmax, istep);
 *
 *     // Run for indexes in the Image class
 *     for (int i = STARTINGY(blobprint()); i <= FINISHINGY(blobprint()); i++)
 *         for (int j = STARTINGX(blobprint()); j <= FINISHINGX(blobprint());
 *              j++)
 *         {
 *             // Compute oversampled index, and then its value
 *             double vi, ui;
 *             IMG2OVER(blobprint, i, j, vi, ui);
 *
 *             double r = sqrt(vi*vi + ui*ui);
 *             IMGPIXEL(blobprint, i, j) = blob_val(r, blob);
 *         }
 *
 *     // Adjust the footprint structure
 *     if (normalise)
 *         blobprint() = blobprint() / blobprint().sum();
 * }
 * @endcode
 *
 * Note: for this class the X axis has been renamed as U, Y as V, and Z as W.
 */
class ImageOver : public Image<double>
{
public:
    int uistep, vistep, wistep; // number of samples per
    // one sample on normal image (50)
    // in u,v,w directions

    int overumax, overvmax, overwmax; // table borders in normal units (-2,2)
    int overumin, overvmin, overwmin; // table borders in normal units (-2,2)
    // They should be an integer number

public:
    /** Prepare the Oversampled Image for work
     *
     * This function defines the oversampled image, it's very important to call
     * it before starting to work with the image. If the sampling rate in any of
     * the directions is an even number, then it is substituted by the next odd
     * number (+1). This is done in order to force the logical oversampled
     * pixels to be exactly in one cell of the underlying 2D matrix.
     *
     * @code
     * IO.init(-2, 2, 51, -2, 2, 51);
     * @endcode
     */
    void init(int _umin, int _umax, int _uistep,
              int _vmin, int _vmax, int _vistep,
              int _wmin=0, int _wmax=0, int _wistep=1);

    /** Initialize from Image
     *
     * This functions sets the origins maximum values and steps from image arguments
     * and values passed.
     */
    void init(MultidimArray<double> &im, int _uistep=1, int _vistep=0, int _wistep=0);

    /** Window
     *
     * Set a selfWindow in the logical space.
     */
    void window(int _v0, int _u0, int _vF, int _uF);

    /** Empty image
     *
     * A 0x0 image with no information about oversampling is produced. Remember
     * to re-initialise the oversampled image before starting to use it again.
     *
     * @code
     * IO.clear();
     * @endcode
     */
    void clear();

    /** Returns the exact position of a non-integer location
     *
     * The existing indexes (iu,iv) within the overimage are returned according
     * to a certain (u,v) and the Oversampled Image definition. Usually this
     * function is not needed in normal oversampled images operation as you can
     * access to pixels with fractional indexes. In our example of footprints,
     * this function would make the conversion between the oversampled position
     * (-1,1) to the image index (-51,51) (what is returned). Don't be mislead
     * by the physical position of the pixel (-51,51) which is (76,178).
     *
     * @code
     * IO.over2img(y, x, iy, ix);
     * @endcode
     */
    void over2img(double v, double u, int& iv, int& iu) const
    {
        if (v < overvmin || v > overvmax)
            REPORT_ERROR(ERR_VALUE_INCORRECT, "ImgeOver::over2img: v out of range");

        if (u < overumin || u > overumax)
            REPORT_ERROR(ERR_VALUE_INCORRECT, "ImgeOver::over2img: u out of range");

        iu = (int) ROUND((((u) - overumin) * uistep));
        iv = (int) ROUND((((v) - overvmin) * vistep));
    }

    /** Speed up pixel index macro
     *
     * This macro is the same as the function over2img but faster due to no
     * function call is performed.
     *
     * @code
     * OVER2IMG(IO, y, x, iy, ix);
     * @endcode
     */
#define OVER2IMG(IO, v, u, iv, iu) \
    iu = (int) round((((u)-(IO).overumin) * (IO).uistep)); \
    iv = (int) round((((v)-(IO).overvmin) * (IO).vistep));

#define OVER2IMG_Z(IO, w, iw) \
	    iw = (int) round((((w)-(IO).overwmin) * (IO).wistep));

    /** Returns the logical position of a "physical" location
     *
     * This function makes the conversion from image logical location to
     * oversampled logical position. In our example of footprints the conversion
     * between the oversampled position (-51,51) to the image index (-1,1) (what
     * is returned). Notice that the "physical" position of the pixel (-51,51)
     * is not the true physical one, that would be (76,178).
     *
     * @code
     * IO.img2over(iv, iu, v, u);
     * @endcode
     */
    void img2over(size_t iv, size_t iu, double & v, double & u) const
    {
        if (iu < 0 || iu > XSIZE(data))
            REPORT_ERROR(ERR_VALUE_INCORRECT, "ImageOver::img2over: iu out of range");
        if (iv < 0 || iv > YSIZE(data))
            REPORT_ERROR(ERR_VALUE_INCORRECT, "ImageOver::img2over: iv out of range");
        u = (double)(overumin) + iu / (double)(uistep);
        v = (double)(overvmin) + iv / (double)(vistep);
    }

    /** Speed up pixel index macro
     *
     * This macro is the same as the function img2over but faster due to no
     * function call is performed.
     *
     * @code
     * IMG2OVER(IO, iy, ix, y, x);
     * @endcode
     */
#define IMG2OVER(IO, iv, iu, v, u) \
    u = (double) (IO).overumin + (iu) / (double) ((IO).uistep); \
    v = (double) (IO).overvmin + (iv) / (double) ((IO).vistep);

    /** Constant pixel access with fractional indexes
     *
     * This function returns you a constant pixel access referred with a
     * fractional index. If you want to access to the pixels in the classic
     * Image fashion then you should use the macro IMGPIXEL. An exception is
     * thrown if you try to access an image outside the image.
     *
     * @code
     * std::cout << IO(1.34,-0.56) << std::endl;
     * @endcode
     */
    double operator ()(double v, double u) const
    {
        if (v < overvmin || v > overvmax)
            REPORT_ERROR(ERR_VALUE_INCORRECT, "ImgeOver::over2img: v out of range");
        if (u < overumin || u > overumax)
            REPORT_ERROR(ERR_VALUE_INCORRECT, "ImgeOver::over2img: u out of range");
        int iu, iv;
        OVER2IMG(*this, v, u, iv, iu);
        return A2D_ELEM(data, iv, iu);
    }

    /** Pixel access with fractional indexes
     *
     * This function allows you to set pixels referred with a fractional index.
     * If you want to access to the pixels in the classic Image fashion then you
     * should use the macro IMGPIXEL.
     *
     * @code
     * IO(1.34, -0.56) = 1;
     * @endcode
     */
    double & operator ()(double v, double u)
    {
        int iu, iv;
        OVER2IMG(*this, v, u, iv, iu);
        return A2D_ELEM(data, iv, iu);
    }

    /** Speed up pixel access macro
     *
     * This macro is a speeded up version of the pixel access routines for
     * oversampled images.
     *
     * @code
     * OVERPIXEL(IO, 1.34, -0.56) = 1;
     * @endcode
     */
#define OVERPIXEL(IO, y, x) IMGPIXEL((IO), \
                                     (int) ROUND(((u) * (IO).uistep)), \
                                     (int) ROUND(((v) * (IO).vistep)))

    // The following two functions have been redefined (they should be
    // inherited from Image) because the compiler doesn't admit this
    // inheritance and the redefinition of the fractional index pixel
    // access.

    /** Matrix access
     *
     * This function allows you to access the Matrix2D over which the complex
     * data type is based.
     *
     * @code
     * std::cout << IO();
     * @endcode
     */
    MultidimArray< double >& operator()()
    {
        return data;
    }

    const MultidimArray< double >& operator()() const
    {
        return data;
    }

    /** Maximum value in not sampled units on U axis
     *
     * In our example: 2.
     *
     * @code
     * int Umax = IO.umax();
     * @endcode
     */
    int umax() const
    {
        return overumax;
    }

    /** Maximum value in not sampled units on V axis
     *
     * In our example: 2.
     *
     * @code
     * int Vmax = IO.vmax();
     * @endcode
     */
    int vmax() const
    {
        return overvmax;
    }

    /** Maximum value in not sampled units on W axis
     *
     * In our example: 2.
     *
     * @code
     * int Wmax = IO.wmax();
     * @endcode
     */
    int wmax() const
    {
        return overwmax;
    }

    /** Minimum value in not sampled units on U axis
     *
     * In our example: -2.
     *
     * @code
     * int Umin = IO.umin();
     * @endcode
     */
    int umin() const
    {
        return overumin;
    }

    /** Minimum value in not sampled units on V axis
     *
     * In our example: -2.
     *
     * @code
     * int Vmin = IO.vmin();
     * @endcode
     */
    int vmin() const
    {
        return overvmin;
    }

    /** Sampling rate in U axis
     *
     * In our example: 51.
     *
     * @code
     * int Usampling = IO.Ustep();
     * @endcode
     */
    int Ustep() const
    {
        return uistep;
    }

    /** Sampling rate in V axis
     *
     * In our example: 51.
     *
     * @code
     * int Vsampling = IO.Vstep();
     * @endcode
     */
    int Vstep() const
    {
        return vistep;
    }

    /** Sampling rate in W axis
     *
     * In our example: 51.
     *
     * @code
     * int Wsampling = IO.Wstep();
     * @endcode
     */
    int Wstep() const
    {
        return wistep;
    }

    /** Generate the normal image by averaging
     *
     * This function passes from an oversampled image to a normal one, in our
     * example from the 250x250 image to a 5x5. A pointer to the normal image
     * must be given.
     *
     * @code
     * IO.down_sample(&I);
     * @endcode
     */
    void downsample(Image< double >* I) const;
};
#endif
