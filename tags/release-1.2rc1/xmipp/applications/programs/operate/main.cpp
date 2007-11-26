/* author: Javier Velazquez
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

#include <data/args.h>
#include <data/image.h>
#include <data/volume.h>
#include <data/funcs.h>

void Usage(void);
bool check_for_operation(int argc, char **argv, char *operation,
                         FileName &fn, int &operand_type);
void compute(int operation, int operand_type1, int operand_type2, FileName &fn_1,
             FileName &fn_2, FileName &fn_out);
void operate_plus(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out);
void operate_minus(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out);
void multiplication(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out);
void division(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out);
void log10(int operand_type1, FileName &fn_1, FileName &fn_out);
void sqrt(int operand_type1, FileName &fn_1, FileName &fn_out);
void extract_row(int operand_type1, int operand_type2,
                 FileName &fn_1, FileName &fn_2, FileName &fn_out);
void extract_column(int operand_type1, int operand_type2,
                    FileName &fn_1, FileName &fn_2, FileName &fn_out);
void extract_slice(int operand_type1, int operand_type2,
                   FileName &fn_1, FileName &fn_2, FileName &fn_out);
void radial_avg(int operand_type1, FileName &fn_1, FileName &fn_out);


// Operations suported
#define    OPERATE_PLUS   1
#define    OPERATE_MINUS  2
#define    MULTIPLICATION 3
#define    DIVISION       4
#define    LOG10          5
#define    SQRT           6
#define    SLICE    7
#define    COLUMN    8
#define    ROW   9
#define    RADIAL_AVG    10

// types supported
#define    VOLUME    1
#define    IMAGE     2
#define    NUMBER    3

/**************************************************************************

   NAME:          main

   DESCRIPTION:

   DATE:          22-8-2001

/**************************************************************************/
int main(int argc, char **argv)
{
    FileName  fn_1, fn_2, fn_out;
    int          operation = 0, operand_type1 = 0, operand_type2 = 0, operand_type_input2 = 0;


    // Obtain command line parameters form the program
    try
    {
        // Read the output image filename
        fn_out  = getParameter(argc, argv, "-o");
        // Check for #, which indicates a binary operation
        check_for_operation(argc, argv, "-i", fn_1, operand_type1);

        // Check the operations supported
        if (check_for_operation(argc, argv, "-plus", fn_2, operand_type2)) operation = OPERATE_PLUS;
        else if (check_for_operation(argc, argv, "-minus", fn_2, operand_type2)) operation = OPERATE_MINUS;
        else if (check_for_operation(argc, argv, "-mult", fn_2, operand_type2)) operation = MULTIPLICATION;
        else if (check_for_operation(argc, argv, "-divide", fn_2, operand_type2)) operation = DIVISION;
        else if (checkParameter(argc, argv, "-log10")) operation = LOG10;
        else if (checkParameter(argc, argv, "-sqrt")) operation = SQRT;
        else if (check_for_operation(argc, argv, "-slice", fn_2, operand_type2)) operation = SLICE;
        else if (check_for_operation(argc, argv, "-column", fn_2, operand_type2)) operation = COLUMN;
        else if (check_for_operation(argc, argv, "-row", fn_2, operand_type2)) operation = ROW;
        else if (checkParameter(argc, argv, "-radial_avg")) operation = RADIAL_AVG;
        else
            REPORT_ERROR(1, "No valid operation specified");

        compute(operation, operand_type1, operand_type2, fn_1, fn_2, fn_out);

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage();
        exit(1);
    }

}


bool check_for_operation(int argc, char **argv, char *operation, FileName &fn, int &operand_type)
{

    if (checkParameter(argc, argv, operation))
    {
        fn  = getParameter(argc, argv, operation);
        // If the file exist, tell if it is an image or a volume
        if (exists(fn))
        {
            if (Is_ImageXmipp(fn))
                operand_type = IMAGE;
            else if (Is_VolumeXmipp(fn))
                operand_type = VOLUME;
            else
                REPORT_ERROR(1502, "One of the files is not a Xmipp image or volume\n");
        }
        // In other case check if it's a number
        else
        {
            // check if we have a number using textToFloat
            // If a problem exist, textToFloat will throw an exception, catched by the main function
            double dummy = textToFloat(fn, 2101, "One of the parameters is neither a number nor a file\n\n");
            operand_type = NUMBER;
        }

        return true;
    }

    return false;

}

// Insert new functions here and create its operations as in plus() or log10(), for example
void compute(int operation, int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out)
{
    switch (operation)
    {
    case    OPERATE_PLUS:
        operate_plus(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    OPERATE_MINUS:
        operate_minus(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    MULTIPLICATION:
        multiplication(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    DIVISION:
        division(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    LOG10:
        log10(operand_type1, fn_1, fn_out);
        break;
    case    SQRT:
        sqrt(operand_type1, fn_1, fn_out);
        break;
    case    SLICE:
        extract_slice(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    COLUMN:
        extract_column(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    ROW:
        extract_row(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    RADIAL_AVG:
        radial_avg(operand_type1, fn_1, fn_out);
        break;
    }
}

void operate_plus(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        ImageXmipp out;
        out.read(fn_2, false, false, true);
        double number1 = textToFloat(fn_1);
        out() = out() + number1;
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        ImageXmipp out;
        out.read(fn_1, false, false, true);
        double number2 = textToFloat(fn_2);
        out() = out() + number2;
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        ImageXmipp Op1, out;
        Op1.read(fn_1, false, false, true);
        out.read(fn_2, false, false, true);
        out() = Op1() + out();
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == NUMBER && operand_type2 == VOLUME)
    {
        VolumeXmipp out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        out() = out() + number1;
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME && operand_type2 == NUMBER)
    {
        VolumeXmipp out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        out() = out() + number2;
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME && operand_type2 == VOLUME)
    {
        VolumeXmipp Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        out() = Op1() + out();
        out.write(fn_out);
    }

}

void operate_minus(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        ImageXmipp out;
        out.read(fn_2, false, false, true);
        double number1 = textToFloat(fn_1);
        out() = number1 - out();
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        ImageXmipp out;
        out.read(fn_1, false, false, true);
        double number2 = textToFloat(fn_2);
        out() = out() - number2;
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        ImageXmipp Op1, out;
        Op1.read(fn_1, false, false, true);
        out.read(fn_2, false, false, true);
        out() = Op1() - out();
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == NUMBER && operand_type2 == VOLUME)
    {
        VolumeXmipp out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        out() = out() - number1;
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME && operand_type2 == NUMBER)
    {
        VolumeXmipp out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        out() = out() - number2;
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME && operand_type2 == VOLUME)
    {
        VolumeXmipp Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        out() = Op1() - out();
        out.write(fn_out);
    }

}

void multiplication(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        ImageXmipp out;
        out.read(fn_2, false, false, true);
        double number1 = textToFloat(fn_1);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(out())
        out(i, j) *= number1;
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        ImageXmipp out;
        out.read(fn_1, false, false, true);
        double number2 = textToFloat(fn_2);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(out())
        out(i, j) *= number2;
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        ImageXmipp Op1, out;
        Op1.read(fn_1, false, false, true);
        out.read(fn_2, false, false, true);
        multiplyElements(Op1(), out(), out());
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == NUMBER && operand_type2 == VOLUME)
    {
        VolumeXmipp out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(out())
        out(k, i, j) *= number1;
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME && operand_type2 == NUMBER)
    {
        VolumeXmipp out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(out())
        out(k, i, j) *= number2;
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME && operand_type2 == VOLUME)
    {
        VolumeXmipp Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(out())
        out(k, i, j) = Op1(k, i, j) * out(k, i, j);
        out.write(fn_out);
    }

}

void division(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        ImageXmipp out;
        out.read(fn_2, false, false, true);
        double number1 = textToFloat(fn_1);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(out())
        out(i, j) = number1 / out(i, j);
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        ImageXmipp out;
        out.read(fn_1, false, false, true);
        double number2 = textToFloat(fn_2);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(out())
        out(i, j) = out(i, j) / number2;
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        ImageXmipp Op1, out;
        Op1.read(fn_1, false, false, true);
        out.read(fn_2, false, false, true);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(out())
        out(i, j) = Op1(i, j) / out(i, j);
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == NUMBER && operand_type2 == VOLUME)
    {
        VolumeXmipp out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(out())
        out(k, i, j) = number1 / out(k, i, j);
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME && operand_type2 == NUMBER)
    {
        VolumeXmipp out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(out())
        out(k, i, j) = out(k, i, j) / number2;
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME && operand_type2 == VOLUME)
    {
        VolumeXmipp Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(out())
        out(k, i, j) = Op1(k, i, j) / out(k, i, j);
        out.write(fn_out);
    }

}


void log10(int operand_type1, FileName &fn_1, FileName &fn_out)
{
    if (operand_type1 == IMAGE)
    {
        ImageXmipp out;
        out.read(fn_1, false, false, true);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(out())
        out(i, j) = log10(1 + out(i, j));
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME)
    {
        VolumeXmipp out;
        out.read(fn_1);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(out())
        out(k, i, j) = log10(1 + out(k, i, j));
        out.write(fn_out);
    }
}

void sqrt(int operand_type1, FileName &fn_1, FileName &fn_out)
{
    if (operand_type1 == IMAGE)
    {
        ImageXmipp out;
        out.read(fn_1, false, false, true);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(out())
        out(i, j) = sqrt(out(i, j));
        out.set_originOffsets(0., 0.);
        out.set_eulerAngles(0., 0., 0.);
        out.write(fn_out);
    }
    else if (operand_type1 == VOLUME)
    {
        VolumeXmipp out;
        out.read(fn_1);
        FOR_ALL_ELEMENTS_IN_MATRIX3D(out())
        out(k, i, j) = sqrt(out(k, i, j));
        out.write(fn_out);
    }
}


void extract_slice(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out)
{
    if (operand_type1 == VOLUME && operand_type2 == NUMBER)
    {
        VolumeXmipp Op1;
        ImageXmipp  out;
        Op1.read(fn_1);
        int number2 = textToInteger(fn_2);
        // If the slice requested exists
        if (number2 >= STARTINGZ(Op1()) && number2 <= FINISHINGZ(Op1()))
        {
            // Resize image
            out().resize(Op1().rowNumber(), Op1().colNumber());
            // Copy
            for (int i = STARTINGY(Op1());i <= FINISHINGY(Op1());i++)
                for (int j = STARTINGX(Op1());j <= FINISHINGX(Op1());j++)
                    out(i, j) = Op1(number2, i, j);
            // Save
            out.write(fn_out);
        }
    }
}

void extract_column(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out)
{
    if (operand_type1 == VOLUME && operand_type2 == NUMBER)
    {
        VolumeXmipp Op1;
        ImageXmipp  out;
        Op1.read(fn_1);
        int number2 = textToInteger(fn_2);
        // If the column requested exists
        if (number2 >= STARTINGX(Op1()) && number2 <= FINISHINGX(Op1()))
        {
            // Resize image
            out().resize(Op1().sliceNumber(), Op1().rowNumber());
            // Copy
            for (int k = STARTINGZ(Op1());k <= FINISHINGZ(Op1());k++)
                for (int i = STARTINGY(Op1());i <= FINISHINGY(Op1());i++)
                    out(k, i) = Op1(k, i, number2);
            // Save
            out.write(fn_out);
        }
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        ImageXmipp  Op1;
        ImageXmipp  out;
        Op1.read(fn_1, false, false, true);
        int number2 = textToInteger(fn_2);
        // If the column requested exists
        if (number2 >= STARTINGX(Op1()) && number2 <= FINISHINGX(Op1()))
        {
            // Resize image
            out().resize(Op1().rowNumber(), 1);
            // Copy
            for (int i = STARTINGY(Op1());i <= FINISHINGY(Op1());i++)
                out(i, 0) = Op1(i, number2);
            // Save
            out.write(fn_out);
        }
    }
}


void extract_row(int operand_type1, int operand_type2, FileName &fn_1, FileName &fn_2, FileName &fn_out)
{
    if (operand_type1 == VOLUME && operand_type2 == NUMBER)
    {
        VolumeXmipp Op1;
        ImageXmipp  out;
        Op1.read(fn_1);
        int number2 = textToInteger(fn_2);
        // If the column requested exists
        if (number2 >= STARTINGY(Op1()) && number2 <= FINISHINGY(Op1()))
        {
            // Resize image
            out().resize(Op1().sliceNumber(), Op1().colNumber());
            // Copy
            for (int k = STARTINGZ(Op1());k <= FINISHINGZ(Op1());k++)
                for (int j = STARTINGX(Op1());j <= FINISHINGX(Op1());j++)
                    out(k, j) = Op1(k, number2, j);
            // Save
            out.write(fn_out);
        }
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        ImageXmipp  Op1;
        ImageXmipp  out;
        Op1.read(fn_1, false, false, true);
        int number2 = textToInteger(fn_2);
        // If the column requested exists
        if (number2 >= STARTINGY(Op1()) && number2 <= FINISHINGY(Op1()))
        {
            // Resize image
            out().resize(1, Op1().colNumber());
            // Copy
            for (int j = STARTINGX(Op1());j <= FINISHINGX(Op1());j++)
                out(0, j) = Op1(number2, j);
            // Save
            out.write(fn_out);
        }
    }
}

void radial_avg(int operand_type1, FileName &fn_1, FileName &fn_out)
{
    if (operand_type1 == IMAGE)
    {
        ImageXmipp input;
        input.read(fn_1, false, false, true);
        input().setXmippOrigin();
        Matrix1D<int> center(2);
        center.initZeros();
        Matrix1D<double> radial_mean;
        Matrix1D<int> radial_count;
        radialAverage(input(), center, radial_mean, radial_count);
        radial_mean.write(fn_out);


        int my_rad;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(input())
        {
//              my_rad=(int)ROUND(sqrt((double)(i*i+j*j))); // bug corrected by JA Velazquez-Muriel
            my_rad = (int)floor(sqrt((double)(i * i + j * j)));
            input(i, j) = radial_mean(my_rad);
        }
        input.write(fn_out + ".img");

    }
    else if (operand_type1 == VOLUME)
    {
        VolumeXmipp input;
        input.read(fn_1);
        input().setXmippOrigin();
        Matrix1D<int> center(2);
        center.initZeros();
        Matrix1D<double> radial_mean;
        Matrix1D<int> radial_count;
        radialAverage(input(), center, radial_mean, radial_count);
        radial_mean.write(fn_out);

        int my_rad;
        FOR_ALL_ELEMENTS_IN_MATRIX3D(input())
        {
//              my_rad=(int)ROUND(sqrt((double)(i*i+j*j+k*k))); // bug corrected by JA Velazquez-Muriel
            my_rad = (int)floor(sqrt((double)(i * i + j * j + k * k)));
            input(k, i, j) = radial_mean(my_rad);
        }
        input.write(fn_out + ".img");
    }
}

/**************************************************************************

   NAME:          usage

   DESCRIPTION:   This function displays how to use the program

   DATE:          19-1-2001

/**************************************************************************/
void Usage()
{
    cout  << " A simple Xmipp images calculator. Binary and unary operations\n"
    << " Parameters:\n"
    << " -i xmipp image or volume. This is the input for the program. \n"
    << " -o output of the program. An xmipp image or volue. \n"
    << "\n"
    << " CURRENTLY SUPPORTED OPERATIONS \n"
    << "================================\n"
    << " -plus <file or value>    Sums two images, volumes or adds a numerical value to an image\n"
    << " -minus <file or value>   Substracts two images, volumes or substracts a numerical value to an image\n"
    << " -mult <file or value>    Multiplies two images, volumes, or multiplies per a given number\n"
    << " -divide <file or value>  Divides two images, volumes, or divides per a given number\n"
    << " -log10                   Computes the logarithm of an image\n"
    << " -sqrt                    Computes the square root of an image\n"
    << " -slice  <value>          Extracts a given slice from a volume\n"
    << " -column <value>          Extracts a given column from a image or volume\n"
    << " -row    <value>          Extracts a given row from a image or volume\n"
    << " -radial_avg              Compute the radial average of an image\n"
    << " EXAMPLES \n"
    << "==========\n"
    << "operate -i image1 -plus image2 -o image3\n"
    << "operate -i image1 -mult image2 -o image3\n"
    << "operate -i image1 -divide 2 -o image3\n"
    << "operate -i image1 -sqrt -o image3\n"
    ;
}
