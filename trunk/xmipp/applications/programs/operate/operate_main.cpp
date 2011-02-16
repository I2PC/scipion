/***************************************************************************
 *
 * Authors: J.M. de la Rosa Trevin   (jmdelarosa@cnb.csic.es) x 0.95
 *          Joaquin Oton             (joton@cnb.csic.es)      x 0.05
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

#include <data/args.h>
#include <data/filters.h>
#include <data/image.h>
#include <data/funcs.h>
#include <algorithm>

#include "data/program.h"
#include "data/metadata_extension.h"

//This define the prototype of binary operations on images
//the result will be left in op1
typedef void ImageBinaryOperator(Image<double> &op1, const Image<double> &op2);

//This define the prototype of unary operations on images
//the result will be left in op
typedef double ImageUnaryOperator(Image<double> &op);

///==== Operations ======
void minus(Image<double> &op1, const Image<double> &op2)
{
    op1() -= op2();
}

void plus(Image<double> &op1, const Image<double> &op2)
{
    op1() += op2();
}

void mult(Image<double> &op1, const Image<double> &op2)
{
    op1() *= op2();
}

void divide(Image<double> &op1, const Image<double> &op2)
{
    op1() /= op2();
}

void min(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = XMIPP_MIN(dAi(op1(), n), dAi(op2(), n));
    }
}

void max(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = XMIPP_MAX(dAi(op1(), n), dAi(op2(), n));
    }
}

void compare(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) == dAi(op2(), n) ? 0 : (dAi(op1(), n) < dAi(op2(), n) ? -1 : 1 );
    }
}

///Be carefull with integer images for relational operations...due to double comparisions
void eq(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) == dAi(op2(), n) ? 1 : 0;
    }
}

void ne(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) != dAi(op2(), n) ? 1 : 0;
    }
}

void lt(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) < dAi(op2(), n) ? 1 : 0;
    }
}

void le(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) <= dAi(op2(), n) ? 1 : 0;
    }
}

void gt(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) > dAi(op2(), n) ? 1 : 0;
    }
}

void ge(Image<double> &op1, const Image<double> &op2)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op1())
    {
        dAi(op1(), n) = dAi(op1(), n) >= dAi(op2(), n) ? 1 : 0;
    }
}


double sqrt(Image<double> &op)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op())
    {
        dAi(op(), n) = sqrt(dAi(op(), n));
    }
}

double log10(Image<double> &op)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op())
    {
        dAi(op(), n) = log10(dAi(op(), n));
    }
}

double powerExp = 2;
double power(Image<double> &op)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op())
    {
        dAi(op(), n) = pow(dAi(op(), n), powerExp);
    }
}

int nSlice;
char axis;
double getSlice(Image<double> &op)
{
    MultidimArray<double> imAux;

    op().getSlice(nSlice, imAux, axis);
    op() = imAux;
}


FileName fnOut;
double radialAvg(Image<double> &op)
{
    op().setXmippOrigin();
    Matrix1D<int> center(3);
    center.initZeros();
    MultidimArray<double> radial_mean;
    MultidimArray<int> radial_count;
    radialAverage(op(), center, radial_mean, radial_count);
    radial_mean.write((fnOut.withoutExtension()).addExtension("txt"));

    int my_rad;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(op())
    {
        my_rad = (int)floor(sqrt((double)(i * i + j * j + k * k)));
        op(k, i, j) = radial_mean(my_rad);
    }
}

/** Some basic operation on images, or selfiles */
class ProgOperate: public XmippMetadataProgram
{
private:
    ///Functions pointers to selected operation
    ImageBinaryOperator * binaryOperator;
    ImageUnaryOperator * unaryOperator;

    FileName fn2;
    MetaData md2;
    MDIterator md2Iterator;
    Image<double> img2;
    bool isValue;
    double value;
    String file_or_value;

protected:
    void defineParams()
    {
        each_image_produces_an_output = true;
        addUsageLine("A simple Xmipp images calculator. Binary and unary operations");
        XmippMetadataProgram::defineParams();
        addParamsLine("== Binary operations: ==");
        addParamsLine("   -plus <file_or_value>    :Sums two images, volumes or adds a numerical value to an image");
        addParamsLine("or -minus <file_or_value>   :Substracts two images, volumes or substracts a numerical value to an image");
        addParamsLine("or -mult <file_or_value>    :Multiplies two images, volumes, or multiplies per a given number");
        addParamsLine("or -divide <file_or_value>  :Divides two images, volumes, or divides per a given number");
        addParamsLine("or -min <file_or_value>     :Minimum of two images, volumes, or number (pixel-wise)");
        addParamsLine("or -max <file_or_value>     :Maximum of two images, volumes, or number (pixel-wise)");
        addParamsLine("or -compare <file_or_value> :Returns -1 if the left value is less, 0 if are equal or 1 if greater.(pixel-wise)");
        addParamsLine("==+ Relational operations: ==");
        addParamsLine("or -eq <file_or_value>      :Returns 1 if the pixels values are equal, 0 otherwise (pixel-wise)");
        addParamsLine("or -le <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");
        addParamsLine("or -lt <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");
        addParamsLine("or -ge <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");
        addParamsLine("or -gt <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");
        addParamsLine("or -ne <file_or_value>      :Returns 1 if the pixels values are equal  less, 0 otherwise (pixel-wise)");

        addParamsLine("== Unary operations: ==");
        addParamsLine("or -log10                   :Computes the logarithm of an image");
        addParamsLine("or -sqrt                    :Computes the square root of an image");
        addParamsLine("or -pow <value=2>           :Computes the power of an image");
        addParamsLine("or -slice  <value>          :Extracts a given slice from a volume");
        addParamsLine("or -column <value>          :Extracts a given column from a image or volume");
        addParamsLine("or -row    <value>          :Extracts a given row from a image or volume");
        addParamsLine("or -radial_avg              :Compute the radial average of an image");
    }

    void readParams()
    {
        XmippMetadataProgram::readParams();
        binaryOperator = NULL;
        unaryOperator = NULL;
        isValue = false;
        // Check operation to do
        //Binary operations
        if (checkParam("-plus"))
        {
            file_or_value = getParam("-plus");
            binaryOperator = plus;
        }
        else if (checkParam("-minus"))
        {
            file_or_value = getParam("-minus");
            binaryOperator = minus;
        }
        else if (checkParam("-mult"))
        {
            file_or_value = getParam("-mult");
            binaryOperator = mult;
        }
        else if (checkParam("-divide"))
        {
            file_or_value = getParam("-divide");
            binaryOperator = divide;
        }
        else if (checkParam("-max"))
        {
            file_or_value = getParam("-max");
            binaryOperator = max;
        }
        else if (checkParam("-min"))
        {
            file_or_value = getParam("-min");
            binaryOperator = min;
        }
        else if (checkParam("-compare"))
        {
            file_or_value = getParam("-compare");
            binaryOperator = compare;
        }
        ///Relational operations
        else if (checkParam("-eq"))
        {
            file_or_value = getParam("-eq");
            binaryOperator = eq;
        }
        else if (checkParam("-ne"))
        {
            file_or_value = getParam("-ne");
            binaryOperator = ne;
        }
        else if (checkParam("-lt"))
        {
            file_or_value = getParam("-lt");
            binaryOperator = lt;
        }
        else if (checkParam("-le"))
        {
            file_or_value = getParam("-le");
            binaryOperator = le;
        }
        else if (checkParam("-gt"))
        {
            file_or_value = getParam("-gt");
            binaryOperator = gt;
        }
        else if (checkParam("-ge"))
        {
            file_or_value = getParam("-ge");
            binaryOperator = ge;
        }
        ///Unary operations
        else if (checkParam("-log10"))
            unaryOperator = log10;
        else if (checkParam("-sqrt"))
            unaryOperator = sqrt;
        else if (checkParam("-pow"))
        {
            powerExp = getDoubleParam("-pow");
            unaryOperator = power;
        }
        else if (checkParam("-slice"))
        {
            axis = 'Z';
            nSlice = getIntParam("-slice");
            unaryOperator = getSlice;
        }
        else if (checkParam("-column"))
        {
            axis = 'X';
            nSlice = getIntParam("-column");
            unaryOperator = getSlice;
        }
        else if (checkParam("-row"))
        {
            axis = 'Y';
            nSlice = getIntParam("-row");
            unaryOperator = getSlice;
        }
        else if (checkParam("-radial_avg"))
        {
            fnOut = fn_out;
            unaryOperator = radialAvg;
        }
        else if (checkParam("-forcePositive"))
            unaryOperator = log10;
        else
            REPORT_ERROR(ERR_VALUE_INCORRECT, "No valid operation specified");
        if (binaryOperator != NULL)
        {
            if (!exists(file_or_value))
            {
                isValue = true;
                value = textToFloat(file_or_value);
                int xdim, ydim, zdim;
                unsigned long ndim;
                if (fn_in.isMetaData())
                    ImgSize(mdIn, xdim, ydim, zdim, ndim);
                else
                    SingleImgSize(fn_in, xdim, ydim, zdim, ndim);
                img2().resizeNoCopy(zdim, ydim, xdim);
                img2().initConstant(value);
            }
            else
            {
                isValue = false;
                fn2 = file_or_value;
                if (fn2.isMetaData())
                {
                    md2.read(fn2);
                    if (mdIn.size() != md2.size())
                        REPORT_ERROR(ERR_MD, "Both metadatas operands should be of same size.");
                    md2Iterator = MDIterator(md2);
                }
                else
                {
                    isValue = true;
                    img2.read(fn2);
                }
            }
        }
    }

    void processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
    {
        Image<double> img;
        img.readApplyGeo(fnImg, mdIn, objId);

        if (unaryOperator != NULL)
            unaryOperator(img);
        else
        {
            if (!isValue)
            {
                img2.readApplyGeo(md2, md2Iterator.objId);
                md2Iterator.moveNext();
            }
            binaryOperator(img, img2);
        }
        img.write(fnImgOut);
    }
}
;//end of class ProgOperate


int main(int argc, char **argv)
{
    ProgOperate program;
    program.read(argc, argv);
    program.tryRun();
}
