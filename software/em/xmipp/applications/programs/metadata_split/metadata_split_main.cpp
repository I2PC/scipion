/***************************************************************************
 *
 * Authors:     Sjors Scheres (scheres@cnb.csic.es)
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

#include <data/metadata.h>
#include <data/args.h>
#include <data/xmipp_program.h>
#include <data/numerical_tools.h>
#include <classification/ahc_classifier.h>

class ProgMetadataSplit: public XmippProgram
{
    FileName fn_in, fn_out, fn_root, fn_cc;
    String sortLabelStr, extension;
    MDLabel sortLabel;
    MetaData  mdIn;
    MetaData *mdPtr, *mdPtr2;
    std::vector<MetaData> mdVector;
    bool     dont_randomize, dont_remove_disabled, dont_sort, use_cc;
    size_t N;
    int AHCiter, AHCsubset;

protected:

    void defineParams()
    {
        addUsageLine("Split a metadata (randomly by default) in any number of equally sized output metadata.");
        addUsageLine("By default, only enabled entries in the input file will be written to the output files.");

        addParamsLine("    -i <inputSelfile>    : Input MetaData File");
        addParamsLine("  [ -n <parts=2> ] : Number of output MetaDatas");
        addParamsLine("  [ --oroot <rootname=\"\"> ] : Rootname for output MetaDatas");
        addParamsLine("                          : output will be rootname_<n>.xmd");
        addParamsLine("  [--dont_randomize ]     : Do not generate random groups");
        addParamsLine("  [--dont_sort ]          : Do not sort the outputs MetaData");
        addParamsLine("  [--dont_remove_disabled]: Do not remove disabled images from MetaData");
        addParamsLine("  [--use_correlation <fnCC=\"\"> <iter=100> <subset=16>] : Use correlation to a single volume.");
        addParamsLine("                             : This is the output of xmipp_reconstruct_significant");
        addParamsLine("                             : with a single reference volume.");
        addParamsLine("                             :+ The metadata is splitted by taking random subsets of the projection");
        addParamsLine("                             :+ directions and applying hierarchical clustering a number of");
        addParamsLine("                             :+ iterations on these random subsets");
        addParamsLine("  [-l <label=\"image\">] : sort using a label, default image");

        addSeeAlsoLine("reconstruct_significant");

        addExampleLine("Splits input.sel in two parts:", false);
        addExampleLine("  xmipp_metadata_split -i input.sel --oroot output_part");
        addExampleLine("Splits input.sel in 4 output files without randomizing input metdata:", false);
        addExampleLine("  xmipp_metadata_split -i input.sel -n 4 --dont_randomize");
        addExampleLine("Splits input.sel in two parts with sel extension:", false);
        addExampleLine("  xmipp_metadata_split -i input.sel --oroot output_part:sel");
    }

    void readParams()
    {
        fn_in = getParam("-i");
        N = getIntParam("-n");
        fn_root = getParam("--oroot");
        extension = fn_root.getFileFormat();
	    fn_root = fn_root.removeFileFormat();
        if (fn_root.empty())
            fn_root = fn_in.withoutExtension();
        if (extension.empty())
            extension = "xmd";
        dont_randomize = checkParam("--dont_randomize");
        dont_sort      = checkParam("--dont_sort");
        dont_remove_disabled = checkParam("--dont_remove_disabled");
        use_cc = checkParam("--use_correlation");
        if (use_cc)
        {
        	fn_cc = getParam("--use_correlation",0);
        	AHCiter = getIntParam("--use_correlation",1);
        	AHCsubset = getIntParam("--use_correlation",2);
        	dont_randomize=true;
        }

        sortLabelStr = "objId"; //if not sort, by default is objId
        if(!dont_sort)
            sortLabelStr = getParam("-l");

        sortLabel = MDL::str2Label(sortLabelStr);
        if (sortLabel == MDL_UNDEFINED)
            REPORT_ERROR(ERR_MD_UNDEFINED, (String)"Unrecognized label '" + sortLabelStr + "'");
    }

public:
    void run()
    {
        mdIn.read(fn_in);

        if (dont_randomize)
            mdPtr = &mdIn;
        else
        {
            mdPtr = new MetaData();
            mdPtr->randomize(mdIn);
        }

        if (!dont_remove_disabled && mdPtr->containsLabel(MDL_ENABLED)) //remove disabled entries by default
            mdPtr->removeObjects(MDValueEQ(MDL_ENABLED, -1));

        size_t Num_images = mdPtr->size();
        size_t Num_groups = XMIPP_MIN(N, Num_images);
        if (!use_cc)
        	mdPtr->split(Num_groups, mdVector); //Split MD in Num_groups
        else
        {
        	Image<double> cc;
        	cc.read(fn_cc);
        	MultidimArray<double> &mcc=cc();
        	int Ndirs=(int)XSIZE(mcc);
        	int Nvols=(int)YSIZE(mcc);
        	int Nimgs=(int)ZSIZE(mcc);
        	if (Nvols!=1)
        		REPORT_ERROR(ERR_ARG_INCORRECT,"This program is meant only for cross-correlation volumes of a single reference model");
        	Matrix2D<int> coocurrence(Nimgs,Nimgs);
        	Matrix2D<double> X(Nimgs,AHCsubset);
        	MultidimArray<int> permutation;
        	AHCClassifier ahc;

        	for (int iter=0; iter<AHCiter; ++iter)
        	{
        		// Take the correlation with random directions that will serve as probe functions
        		randomPermutation(Ndirs,permutation);
        		FOR_ALL_ELEMENTS_IN_MATRIX2D(X)
        			MAT_ELEM(X,i,j)=A3D_ELEM(mcc,i,0,A1D_ELEM(permutation,j));

        		// Classify
        		ahc.clusterData(X,(int)Num_groups);

        		// Update coocurrence
				for (int i=0; i<Nimgs-1; ++i)
					for (int j=i+1; j<Nimgs; ++j)
						if (VEC_ELEM(ahc.clusterAssigned,i)==VEC_ELEM(ahc.clusterAssigned,j))
						{
							MAT_ELEM(coocurrence,i,j)+=1;
							MAT_ELEM(coocurrence,j,i)=MAT_ELEM(coocurrence,i,j);
						}
        	}

        	// Convert the coocurrence matrix into a distance matrix
        	int maxCoocurrence=coocurrence.computeMax();
        	Matrix2D<double> D;
        	D.resizeNoCopy(coocurrence);
        	FOR_ALL_ELEMENTS_IN_MATRIX2D(coocurrence)
        		MAT_ELEM(D,i,j)=maxCoocurrence-MAT_ELEM(coocurrence,i,j);

        	// Now perform the clustering with this distance
        	ahc.clusterWithDistance(D,(int)Num_groups);

        	// Construct output metadata
        	MetaData dummy;
        	for (size_t k=0; k<Num_groups; ++k)
        		mdVector.push_back(dummy);
        	size_t i=0;
        	MDRow row;
        	FOR_ALL_OBJECTS_IN_METADATA(*mdPtr)
        	{
        		mdPtr->getRow(row,__iter.objId);
        		mdVector[VEC_ELEM(ahc.clusterAssigned,i)].addRow(row);
        		i++;
        	}
        }

        //Write splitted groups to disk
        for (size_t i = 0; i < Num_groups; ++i)
        {
            fn_out.compose(fn_root, i + 1, extension);

            if(!dont_sort)
                mdPtr->sort(mdVector[i], sortLabel);
            else
                mdPtr = &(mdVector[i]);

            mdPtr->write(fn_out);
        }
    }

}
;//end of class ProgMetadataSplit

/* MAIN -------------------------------------------------------------------- */
RUN_XMIPP_PROGRAM(ProgMetadataSplit)

