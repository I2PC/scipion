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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/image.h>
#include <data/funcs.h>
#include <algorithm>

void Usage(void);
bool check_for_operation(int argc, char **argv, char *operation,
							FileName &fn, int &operand_type);
void compute(int operation, int operand_type1, int operand_type2,
			   const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void operate_plus(int operand_type1, int operand_type2,
					const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void operate_minus(int operand_type1, int operand_type2,
                     const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void multiplication(int operand_type1, int operand_type2,
                      const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void division(int operand_type1, int operand_type2,
               const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void operate_min(int operand_type1, int operand_type2,
                   const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void operate_max(int operand_type1, int operand_type2,
                   const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void operate_which_min(int operand_type1, int operand_type2,
                          const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void operate_which_max(int operand_type1, int operand_type2,
                          const FileName &fn_1, const FileName &fn_2, const FileName &fn_out);
void log10(int operand_type1, const FileName &fn_1,const FileName &fn_out);
void forcePositive(int operand_type1, const FileName &fn_1,const FileName &fn_out);
void sqrt(int operand_type1, const FileName &fn_1, const FileName &fn_out);
void extract_row(int operand_type1, int operand_type2,
                 const FileName &fn_1, const FileName &fn_2,
                 const FileName &fn_out);
void extract_column(int operand_type1, int operand_type2,
                    const FileName &fn_1, const FileName &fn_2,
                    const FileName &fn_out);
void extract_slice(int operand_type1, int operand_type2,
                   const FileName &fn_1, const FileName &fn_2,
                   const FileName &fn_out);
void radial_avg(int operand_type1, const FileName &fn_1, const FileName &fn_out);

// Operations suported
#define    OPERATE_PLUS       1
#define    OPERATE_MINUS      2
#define    MULTIPLICATION     3
#define    DIVISION           4
#define    LOG10              5
#define    SQRT               6
#define    SLICE              7
#define    COLUMN             8
#define    ROW                9
#define    RADIAL_AVG        10
#define    FORCE_POSITIVE    11
#define    OPERATE_MAX       12
#define    OPERATE_MIN       13
#define    OPERATE_WHICH_MAX 14
#define    OPERATE_WHICH_MIN 15

// types supported
#define    IMAGE     1
#define    NUMBER    2
#define    SELFILE   3

/**************************************************************************

   NAME:          main

   DESCRIPTION:

   DATE:          22-8-2001

/**************************************************************************/
int main(int argc, char **argv)
{
    FileName  fn_1, fn_2, fn_out;
    int operation = 0,
         operand_type1 = 0,
         operand_type2 = 0,
         operand_type_input2 = 0;

    // Obtain command line parameters form the program
    try
    {
        fn_out  = getParameter(argc, argv, "-oext","");
        if (fn_out=="")
            fn_out  = getParameter(argc, argv, "-o","");
        check_for_operation(argc, argv, "-i", fn_1, operand_type1);

        // Check the operations supported
        if (check_for_operation(argc, argv, "-plus", fn_2, operand_type2))
            operation = OPERATE_PLUS;
        else if (check_for_operation(argc, argv, "-minus", fn_2, operand_type2))
            operation = OPERATE_MINUS;
        else if (check_for_operation(argc, argv, "-mult", fn_2, operand_type2))
            operation = MULTIPLICATION;
        else if (check_for_operation(argc, argv, "-divide", fn_2, operand_type2))
            operation = DIVISION;
        else if (checkParameter(argc, argv, "-log10"))
            operation = LOG10;
        else if (checkParameter(argc, argv, "-sqrt"))
            operation = SQRT;
        else if (check_for_operation(argc, argv, "-slice", fn_2, operand_type2))
            operation = SLICE;
        else if (check_for_operation(argc, argv, "-column", fn_2, operand_type2))
            operation = COLUMN;
        else if (check_for_operation(argc, argv, "-row", fn_2, operand_type2))
            operation = ROW;
        else if (checkParameter(argc, argv, "-radial_avg"))
            operation = RADIAL_AVG;
        else if (checkParameter(argc, argv, "-forcePositive"))
            operation = FORCE_POSITIVE;
        else if (check_for_operation(argc, argv, "-max", fn_2, operand_type2))
            operation = OPERATE_MAX;
        else if (check_for_operation(argc, argv, "-min", fn_2, operand_type2))
            operation = OPERATE_MIN;
        else if (check_for_operation(argc, argv, "-which_max", fn_2, operand_type2))
            operation = OPERATE_WHICH_MAX;
        else if (check_for_operation(argc, argv, "-which_min", fn_2, operand_type2))
            operation = OPERATE_WHICH_MIN;
        else
            REPORT_ERROR(1, "No valid operation specified");

        compute(operation, operand_type1, operand_type2, fn_1, fn_2, fn_out);
    }
    catch (XmippError XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }
}

bool check_for_operation(int argc, char **argv, char *operation,
                         FileName &fn, int &operand_type)
{
    if (checkParameter(argc, argv, operation))
    {
        fn  = getParameter(argc, argv, operation);
        // If the file exist, tell if it is an image or a volume
        if (existsTrim(fn) && operation!="-slice" && operation!="-row"
            && operation!="-column")
        {
            if (fn.isMetaData())
                operand_type = SELFILE;
            else
                operand_type = IMAGE;

        }
        // In other case check if it's a number
        else
        {
            // check if we have a number using textToFloat
            // If a problem exist, textToFloat will throw an exception, catched by the main function
            double dummy = textToFloat(fn);
            operand_type = NUMBER;
        }

        return true;
    }

    return false;
}

// Insert new functions here and create its operations as in plus() or log10(), for example
void compute(int operation, int operand_type1, int operand_type2,
             const FileName &fn_1, const FileName &fn_2, const FileName &fn_out)
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
    case    FORCE_POSITIVE:
        forcePositive(operand_type1, fn_1, fn_out);
        break;
    case    OPERATE_MAX:
        operate_max(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    OPERATE_MIN:
        operate_min(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    OPERATE_WHICH_MAX:
        operate_which_max(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    case    OPERATE_WHICH_MIN:
        operate_which_min(operand_type1, operand_type2, fn_1, fn_2, fn_out);
        break;
    }
}

void operate_plus(int operand_type1, int operand_type2, const FileName &fn_1,
                  const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        Image<double> out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        out() = out() + number1;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_2);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        out() = out() + number2;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        Image<double> Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        out() = Op1() + out();
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());

        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);

            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            operate_plus(IMAGE, operand_type2, fn_img, fn_2, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void operate_minus(int operand_type1, int operand_type2,
                   const FileName &fn_1, const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        Image<double> out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        out() = number1 - out();
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_2);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        out() = out() - number2;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        Image<double> Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        out() = Op1() - out();
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());

        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);

            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            operate_minus(IMAGE, operand_type2, fn_img, fn_2, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void multiplication(int operand_type1, int operand_type2,
                    const FileName &fn_1, const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        Image<double> out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        out() *= number1;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_2);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        out() *= number2;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        Image<double> Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        out() *= Op1();
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());

        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            multiplication(IMAGE, operand_type2, fn_img, fn_2, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void division(int operand_type1, int operand_type2,
              const FileName &fn_1, const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        Image<double> out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = number1 / DIRECT_MULTIDIM_ELEM(out(),n);
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_2);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        out() /= number2;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        Image<double> Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = DIRECT_MULTIDIM_ELEM(Op1(),n) / DIRECT_MULTIDIM_ELEM(out(),n);
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            division(IMAGE, operand_type2, fn_img, fn_2, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void operate_max(int operand_type1, int operand_type2, const FileName &fn_1,
                 const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        Image<double> out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = XMIPP_MAX(DIRECT_MULTIDIM_ELEM(out(),n),number1);
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_2);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = XMIPP_MAX(DIRECT_MULTIDIM_ELEM(out(),n),number2);
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        Image<double> Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = XMIPP_MAX(DIRECT_MULTIDIM_ELEM(out(),n), \
                                        DIRECT_MULTIDIM_ELEM(Op1(),n));
        //        out() = Op1() + out(); //FIXME this line is badly included. Check and remove in version 2.4
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            operate_max(IMAGE, operand_type2, fn_img, fn_2, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void operate_min(int operand_type1, int operand_type2, const FileName &fn_1,
                 const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        Image<double> out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = XMIPP_MIN(DIRECT_MULTIDIM_ELEM(out(),n),number1);
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_2);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = XMIPP_MIN(DIRECT_MULTIDIM_ELEM(out(),n),number2);
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        Image<double> Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = XMIPP_MIN(DIRECT_MULTIDIM_ELEM(out(),n),DIRECT_MULTIDIM_ELEM(Op1(),n));

        //        out() = Op1() + out(); //FIXME this line is badly included. Check and remove in version 2.4
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            operate_min(IMAGE, operand_type2, fn_img, fn_2, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void operate_which_max(int operand_type1, int operand_type2, const FileName &fn_1,
                       const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        Image<double> out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = (number1>DIRECT_MULTIDIM_ELEM(out(),n))?0:1;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_2);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = (DIRECT_MULTIDIM_ELEM(out(),n)>number2)?0:1;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        Image<double> Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = (DIRECT_MULTIDIM_ELEM(Op1(),n)>DIRECT_MULTIDIM_ELEM(out(),n))?0:1;
        //        out() = Op1() + out(); FIXME the same copied error. CHECK in v 2.4
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            operate_which_max(IMAGE, operand_type2, fn_img, fn_2, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void operate_which_min(int operand_type1, int operand_type2, const FileName &fn_1,
                       const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == NUMBER && operand_type2 == IMAGE)
    {
        Image<double> out;
        out.read(fn_2);
        double number1 = textToFloat(fn_1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = (number1<DIRECT_MULTIDIM_ELEM(out(),n))?0:1;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_2);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> out;
        out.read(fn_1);
        double number2 = textToFloat(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = (DIRECT_MULTIDIM_ELEM(out(),n)<number2)?0:1;
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == IMAGE && operand_type2 == IMAGE)
    {
        Image<double> Op1, out;
        Op1.read(fn_1);
        out.read(fn_2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = (DIRECT_MULTIDIM_ELEM(Op1(),n)<DIRECT_MULTIDIM_ELEM(out(),n))?0:1;
        //        out() = Op1() + out(); // FIXME again the same problem. Check in v 2.4
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            operate_which_min(IMAGE, operand_type2, fn_img, fn_2, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void log10(int operand_type1, const FileName &fn_1, const FileName &fn_out)
{
    if (operand_type1 == IMAGE)
    {
        Image<double> out;
        out.read(fn_1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = log10(DIRECT_MULTIDIM_ELEM(out(),n));
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            log10(IMAGE, fn_img, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void forcePositive(int operand_type1, const FileName &fn_1,
                   const FileName &fn_out)
{
    if (operand_type1 == IMAGE)
    {
        Image<double> out;
        out.read(fn_1);

        bool negativeRemaining;

        if (out().zdim==1) // IMAGE
        {
            do
            {
            	negativeRemaining=false;

                FOR_ALL_ELEMENTS_IN_ARRAY2D(out())
                if (out(i, j)<=0)
                {
                    std::vector<double> neighbours;
                    for (int ii=-2; ii<=2; ii++)
                    {
                        int iii=i+ii;
                        if (iii<0 || iii>=YSIZE(out()))
                            continue;
                        for (int jj=-2; jj<=2; jj++)
                        {
                            int jjj=j+jj;
                            if (jjj<0 || jjj>=XSIZE(out()))
                                continue;
                            double val=out(iii,jjj);
                            if (val>0)
                                neighbours.push_back(val);
                        }
                    }
                    int N=neighbours.size();
                    if (N==0)
                        negativeRemaining=true;
                    else
                    {
                        std::sort(neighbours.begin(),neighbours.end());
                        if (N%2==0)
                            out(i,j)=0.5*(neighbours[N/2-1]+neighbours[N/2]);
                        else
                            out(i,j)=neighbours[N/2];
                    }
                }
            }
            while (negativeRemaining);
            out.setShifts(0., 0.);
            out.setEulerAngles(0., 0., 0.);
            if (fn_out=="")
                out.write(fn_1);
            else
                out.write(fn_out);

        }
        else // VOLUME
        {
            do
            {
            	negativeRemaining=false;

                FOR_ALL_ELEMENTS_IN_ARRAY3D(out())
                if (out(k, i, j)<=0)
                {
                    std::vector<double> neighbours;
                    for (int kk=-2; kk<=2; kk++)
                    {
                        int kkk=k+kk;
                        if (kkk<0 || kkk>=ZSIZE(out()))
                            continue;
                        for (int ii=-2; ii<=2; ii++)
                        {
                            int iii=i+ii;
                            if (iii<0 || iii>=YSIZE(out()))
                                continue;
                            for (int jj=-2; jj<=2; jj++)
                            {
                                int jjj=j+jj;
                                if (jjj<0 || jjj>=XSIZE(out()))
                                    continue;
                                double val=out(kkk,iii,jjj);
                                if (val>0)
                                    neighbours.push_back(val);
                            }
                        }
                        int N=neighbours.size();
                        if (N==0)
                            negativeRemaining=true;
                        else
                        {
                            std::sort(neighbours.begin(),neighbours.end());
                            if (N%2==0)
                                out(k,i,j)=0.5*(neighbours[N/2-1]+
                                                neighbours[N/2]);
                            else
                                out(k,i,j)=neighbours[N/2];
                        }
                    }
                }
            }
            while (negativeRemaining);
            if (fn_out=="")
                out.write(fn_1);
            else
                out.write(fn_out);
        }
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            forcePositive(IMAGE, fn_img, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void sqrt(int operand_type1, const FileName &fn_1, const FileName &fn_out)
{
    if (operand_type1 == IMAGE)
    {
        Image<double> out;
        out.read(fn_1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(out())
        DIRECT_MULTIDIM_ELEM(out(),n) = sqrt(DIRECT_MULTIDIM_ELEM(out(),n));
        out.setShifts(0., 0.);
        out.setEulerAngles(0., 0., 0.);
        if (fn_out=="")
            out.write(fn_1);
        else
            out.write(fn_out);
    }
    else if (operand_type1 == SELFILE)
    {
        MetaData SF;
        SF.read(fn_1);
        init_progress_bar(SF.size());
        int i=0;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            FileName fnl_out, fn_img;
            SF.getValue(MDL_IMAGE, fn_img);
            if (fn_out!="")
                fnl_out = fn_img.without_extension() + "." + fn_out;
            else
                fnl_out="";
            sqrt(IMAGE, fn_img, fnl_out);
            if ((i++ % 50) ==0)
                progress_bar(i);
        }
        progress_bar(SF.size());
    }
}

void extract_slice(int operand_type1, int operand_type2, const FileName &fn_1,
                   const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == IMAGE && operand_type2 == NUMBER)
    {
        Image<double> Op1;
        Image<double>  out;
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

void extract_column(int operand_type1, int operand_type2, const FileName &fn_1,
                    const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == IMAGE && operand_type2 == NUMBER) // VOLUME
    {
        Image<double> Op1;
        Image<double>  out;
        Op1.read(fn_1);
        int number2 = textToInteger(fn_2);
        // If the column requested exists
        if (number2 >= STARTINGX(Op1()) && number2 <= FINISHINGX(Op1()))
        {
            // Resize image
            out().resize(Op1().ydim, Op1().zdim); // Columns along Z axis are now along X axis
            // Copy
            for (int k = STARTINGZ(Op1());k <= FINISHINGZ(Op1());k++)
                for (int i = STARTINGY(Op1());i <= FINISHINGY(Op1());i++)
                {
                    out(i, k) = Op1(k, i, number2);
                    if (Op1().zdim==1)
                    {
                        // For 1D output: also write in text format to screen
                        std::cout <<i+FIRST_XMIPP_INDEX(Op1().colNumber())<<" "<<Op1(i, number2)<<std::endl;
                    }
                }
            // Save
            out.write(fn_out);
        }
    }
}


void extract_row(int operand_type1, int operand_type2, const FileName &fn_1,
                 const FileName &fn_2, const FileName &fn_out)
{
    if (operand_type1 == IMAGE && operand_type2 == NUMBER) // VOLUME
    {
        Image<double> Op1;
        Image<double>  out;
        Op1.read(fn_1);
        int number2 = textToInteger(fn_2);
        // If the column requested exists
        if (number2 >= STARTINGY(Op1()) && number2 <= FINISHINGY(Op1()))
        {
            // Resize image
            out().resize(Op1().zdim, Op1().xdim);
            // Copy
            for (int k = STARTINGZ(Op1());k <= FINISHINGZ(Op1());k++)
                for (int j = STARTINGX(Op1());j <= FINISHINGX(Op1());j++)
                {
                    out(k, j) = Op1(k, number2, j);
                    if (Op1().zdim==1)
                    {
                        // For 1D output: also write in text format to screen
                        std::cout <<j+FIRST_XMIPP_INDEX(Op1().ydim)<<" "<<Op1(number2,j)<<std::endl;
                    }
                }
            // Save
            out.write(fn_out);
        }
    }
}

void radial_avg(int operand_type1, const FileName &fn_1, const FileName &fn_out)
{
    if (operand_type1 == IMAGE)
    {
        Image<double> input;
        input.read(fn_1);
        input().setXmippOrigin();

        if (input().zdim==1) // IMAGE
        {
            Matrix1D<int> center(2);
            center.initZeros();
            MultidimArray<double> radial_mean;
            MultidimArray<int> radial_count;
            radialAverage(input(), center, radial_mean, radial_count);
            radial_mean.write(fn_out);


            int my_rad;
            FOR_ALL_ELEMENTS_IN_ARRAY2D(input())
            {
                my_rad = (int)floor(sqrt((double)(i * i + j * j)));
                input(i, j) = radial_mean(my_rad);
            }
        }
        else // VOLUME
        {
            Matrix1D<int> center(3);
            center.initZeros();
            MultidimArray<double> radial_mean;
            MultidimArray<int> radial_count;
            radialAverage(input(), center, radial_mean, radial_count);
            radial_mean.write(fn_out);

            int my_rad;
            FOR_ALL_ELEMENTS_IN_ARRAY3D(input())
            {
                my_rad = (int)floor(sqrt((double)(i * i + j * j + k * k)));
                input(k, i, j) = radial_mean(my_rad);
            }
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
    std::cout  << " A simple Xmipp images calculator. Binary and unary operations\n"
    << " Parameters:\n"
    << " -i xmipp image, selfile or volume. This is the input to the program. \n"
    << "                                    Only image selfiles are allowed\n"
    << "[-o <file> / -oext <extension>]     If no output is given, the input\n"
    << "                                    images are rewritten.\n"
    << "\n"
    << " CURRENTLY SUPPORTED OPERATIONS \n"
    << "================================\n"
    << " -plus <file or value>    Sums two images, volumes or adds a numerical value to an image\n"
    << " -minus <file or value>   Substracts two images, volumes or substracts a numerical value to an image\n"
    << " -mult <file or value>    Multiplies two images, volumes, or multiplies per a given number\n"
    << " -divide <file or value>  Divides two images, volumes, or divides per a given number\n"
    << " -min <file or value>     Minimum of two images, volumes, or number (pixel-wise)\n"
    << " -max <file or value>     Maximum of two images, volumes, or number (pixel-wise)\n"
    << " -which_min <file or value> Returns 0 if the first argument is the minimum of the two (pixel-wise)\n"
    << " -which_max <file or value> Returns 0 if the first argument is the maximum of the two (pixel-wise)\n"
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
