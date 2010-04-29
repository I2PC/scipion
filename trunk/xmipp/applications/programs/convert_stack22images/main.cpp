/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2007)
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
#include <data/image.h>
#include <data/selfile.h>
#include <data/volume.h>

void Usage();

int main(int argc, char **argv)
{
    FileName fn_stack, fn_sel, fn_vol;
    bool reversed, skipHeaders;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_stack = getParameter(argc,argv,"-stack","");
        fn_sel   = getParameter(argc,argv,"-sel","");
        fn_vol   = getParameter(argc,argv,"-vol","");
        reversed = checkParameter(argc,argv,"-reverse");
        skipHeaders = checkParameter(argc,argv,"-skipHeaders");

        if (fn_sel=="" && fn_vol=="" || fn_stack=="" && fn_vol=="" ||
            fn_sel=="" && fn_stack=="")
            REPORT_ERROR(1,"stack22images: Please, provide at two of -stack, -sel and -vol");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }

    // True work -----------------------------------------------------------
    try
    {
        ImageXmippStack stack;
        VolumeXmipp V;
        SelFile SF;

        // From stack to selfile or volume .................................
        if (exists(fn_stack))
        {
            stack.readFromStack(fn_stack,reversed,false,false,skipHeaders);

            // From stack to selfile
            if (fn_sel!="")
            {
                FileName fn_root=fn_sel.without_extension();
                stack.writeAsSelFile(fn_root);
            }
            // From stack to volume
            else
            {
                stack.writeAsVolume(fn_vol);
            }
        }
        // From volume to stack or selfile .................................
        else if (exists(fn_vol) && Is_VolumeXmipp(fn_vol))
        {
            V.read(fn_vol);

            // From volume to selfile
            if (fn_sel!="")
            {
                FileName fn_root=fn_sel.without_extension();
                V.writeAsSelFile(fn_root);
            }
            // From volume to stack
            else
            {
                for (int k=0; k<ZSIZE(V()); k++)
                {
                    ImageXmipp I;
                    V().getSlice(k,I());
                    stack.appendImage(I);
                }
                stack.writeAsStack(fn_stack);
            }
            // From selfile to volume or stack .................................
        }
        else
        {
            SF.read(fn_sel);
            if (fn_vol!="")
            {
                int Zdim, Ydim, Xdim;
                SF.ImgSize(Ydim,Xdim);
                Zdim=SF.ImgNo();
                V().initZeros(Zdim,Ydim,Xdim);
            }
            int k=0;
            while (!SF.eof())
            {
                FileName fn_img=SF.NextImg();
                if (fn_img=="") break;
                ImageXmipp I(fn_img);
                if (fn_vol!="")
                    V().setSlice(k++,I());
                else
                    stack.appendImage(I);
            }
            if (fn_vol!="")
                V.write(fn_vol);
            else
                stack.writeAsStack(fn_stack);
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    exit(0);
} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cerr << "Purpose:\n";
    std::cerr << "    Converts Spider stacks into images or volumes or any other\n"
              << "    combination of these three elements\n";
    std::cerr << "Usage: stack22images " << std::endl
              << "   [-stack stackFile]: Stack with the set of images\n"
              << "   [-sel selFile]    : Selfile with the set of images\n"
              << "   [-vol volume]     : Volume with the set of images\n"
              << "   [-reverse]        : Reverse endiannness for reading the stack\n"
              << "   [-skipHeaders]    : skip headers in badly formed Spider stacks\n"
              ;
}
