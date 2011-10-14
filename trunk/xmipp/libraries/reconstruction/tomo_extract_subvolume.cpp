/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2010)
 *
 * Unidad de Bioinformatica del Centro Nacional de Biotecnologia , CSIC
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
#include "tomo_extract_subvolume.h"

//#define DEBUG

// Read arguments ==========================================================

void ProgTomoExtractSubvolume::preProcess()
{
    produceSideInfo();

    center.resizeNoCopy(3);
    doccenter.resizeNoCopy(3);
    A.resizeNoCopy(3,3);
    R.resizeNoCopy(3,3);
    I.resizeNoCopy(3,3);
    I.initIdentity();
}
void ProgTomoExtractSubvolume::defineParams()
{
    defaultComments["-i"].addComment("Metadata file with input volumes (and rotations/shifts)");

    XmippMetadataProgram::defineParams();
    addUsageLine("Extract subvolumes of the asymmetric parts of each subtomogram");

    addUsageLine("+++This program works closely together with the ml_tomo program, which is used to align");
    addUsageLine("+++ (and classify) a series of subtomograms. Now, if each subtomogram contains some sort");
    addUsageLine("+++ of (pseudo-) symmetry, tomo_extract_subvolume can be used to extract subvolumes of the");
    addUsageLine("+++ asymmetric parts of each subtomogram. The metadata file  that is output by this program orients");
    addUsageLine("+++ each subvolume and its missing wedge in the correct orientation, so that it can be fed");
    addUsageLine("+++ directly back into ml_tomo. There, the sub-subtomograms can be further aligned and/or classified");
    addParamsLine("--oroot <root=\"out\">  : Root FileName output subvolume");
    addParamsLine("[-o <filename=\"\">]         : Name of output metadata (\"oroot\".xmd by default)");
    addParamsLine("--sym  <sym=\"c1\">     : Symmetry group");
    addParamsLine("--size     <dim>        : size output subvolumes");
    addParamsLine("--mindist  <distance>   : minimum distance between subvolume centers, usefull to avoid repetition of subvolumes place at simmetry axis");
    addParamsLine("--center  <x> <y> <z>   :  position of center of subvolume to be extracted");
    addExampleLine("Extract 12 vertices (subvolumes) in boxes of size 21x21x21 pixels from each subtomogram in the data set: ", false);
    addExampleLine("xmipp_extract_subvolume -i align/mltomo_1deg_it000001.doc -center 0 0 59 -size 21 -sym i3 -o vertices");

    addExampleLine("+++Imagine we have N subtomograms of a virus particle and we have aligned them",false);
    addExampleLine("+++ with respect to a reference structure (or we have aligned them in a reference-free manner).",false);
    addExampleLine("+++ Then, suppose we are interested in the vertices of the virus: perhaps because we suspect",false);
    addExampleLine("+++ that one of them is 'special': e.g. it is a portal. All we have to do is:",false);

    addExampleLine("+++* Pass the metadata of the aligned subtomograms of the ml_tomo prgram",false);
    addExampleLine("+++* Find the x,y,z coordinates of one of the vertices in the average (reference) map",false);
    addExampleLine("+++* Pass the (pseudo) symmetry description of the virus (icosahedral, see Symmetries)",false);

    addExampleLine("+++The program will then extract all symmetry-related vertices in each of the N aligned subtomograms.",false);
    addExampleLine("+++Each icosahedral has 12 vertices, so if the x,y,z coordinates coincided with a vertex (with an accuracy ",false);
    addExampleLine("+++given by the --mindist option, see below)",false);
    addExampleLine("+++, then we will end up with 12*N sub-subtomograms. The sub-subtomograms are not rotated ",false);
    addExampleLine("+++to avoid interpolation. Instead, the output metadata holds the correct orientations ",false);
    addExampleLine("+++and translation for each of them to subsequently superimpose them in the ml_tomo program. ",false);
    addExampleLine("+++Also the missing wedges of all sub-subtomograms will be treated correctly in this way. Then, in the ",false);
    addExampleLine("+++ml_tomo program one could attempt to classify the 12*N vertices, in an attempt to 'fish out' the N ",false);
    addExampleLine("+++supposed portal structures.",false);

    addExampleLine("+++The format of the input and output files is as described for the ml_tomo program.",false);

    addExampleLine("+++ First, lets align our set of virus subtomograms (images.sel, with missing wedges defined ",false);
    addExampleLine("+++ in images.doc and wedges.doc, see ml_tomo for details) ",false);
    addExampleLine("+++. We will use a known virus structure (myreference.vol) as reference and will use ",false);
    addExampleLine("+++ maximum cross-correlation (rather than maximum likelihood) to avoid problems with absolute ",false);
    addExampleLine("+++ greyscales and the standard deviation in the noise. ",false);

    addExampleLine("+++ ml_tomo -i images.sel -doc images.doc -ref myreference.vol --oroot align/mltomo_10deg -missing wedges.doc -iter 1 -ang 10");
    addExampleLine("+++ ml_tomo -i images.sel -doc align/mltomo_10deg_it000001.doc -ref myreference.vol --oroot align/mltomo_5deg -missing wedges.doc -iter 1 -ang 5 -ang_search 20");
    addExampleLine("+++ ml_tomo -i images.sel -doc align/mltomo_5deg_it000001.doc -ref myreference.vol --oroot align/mltomo_2deg -missing wedges.doc -iter 1 -ang 2 -ang_search 8");
    addExampleLine("+++ ml_tomo -i images.sel -doc align/mltomo_2deg_it000001.doc -ref myreference.vol --oroot align/mltomo_1deg -missing wedges.doc -iter 1 -ang 1 -ang_search 3");

    addExampleLine("+++ Now, we will identify the x,y,z coordinates of the vertex in the reference structure (which has symmetry i3): ",false);
    addExampleLine("+++ -center 0 0 59. And we will extract 12 vertices (subvolumes) in boxes of size 12x12x12 pixels ",false);
    addExampleLine("+++ from each subtomogram in the data set: ",false);

    addExampleLine("+++ xmipp_extract_subvolume -doc align/mltomo_1deg_it000001.doc -center 0 0 59 -size 21 -sym i3 -o vertices ",false);

    addExampleLine("+++ his has generated subvolumes in the same directory as the original subtomograms. For each subtomogram 12  ",false);
    addExampleLine("+++ subvolumes, which are named in the same way as the original subtomograms, but ending in _sub000001.vol etc. ",false);
    addExampleLine("+++ The orientations and translations are in a file called vertices.doc, the selection file is called vertices.sel.  ",false);
    addExampleLine("+++ These files are already in the correct format for the ml_tomo program. Therefore, classifying the vertices into  ",false);
    addExampleLine("+++ 4 classes using ML is now straightforward (see ml_tomo for details again) ",false);
    addExampleLine("+++ : ",false);

    addExampleLine("+++ ml_tomo -i vertices.sel -doc vertices.doc -dont_align -keep_angles -missing wedges.doc -nref 4 -reg0 5 -iter 25");

    addExampleLine("+++ Note that the whole functionality of ml_tomo can be employed on the sub-subtomograms: alignment, ",false);
    addExampleLine("+++ classification, local angular searches, symmetry, etc. For example, in case of the vertices, on ",false);
    addExampleLine("+++ could apply -sym c5, and/or one could allow the vertices to re-orient slightly (with respect ",false);
    addExampleLine("+++ to the i3 symmetry) by using -ang 5 ang_search 15. ",false);
}

void ProgTomoExtractSubvolume::readParams()
{
	XmippMetadataProgram::readParams();
	oroot = getParam("--oroot");
	fn_out = getParam("-o");
	if (fn_out.empty())
		fn_out = oroot.addExtension("xmd");
    // Read command line
    fn_sym  = getParam( "--sym");
    center_ref.resize(3);
    size = getIntParam("--size");
    XX(center_ref) = getIntParam( "--center",0);
    YY(center_ref) = getIntParam( "--center",1);
    ZZ(center_ref) = getIntParam( "--center",2);
    mindist = getDoubleParam("--mindist");
}

// Usage ===================================================================
void ProgTomoExtractSubvolume::show()
{
#ifdef DEBUG
    std::cerr << "start show" <<std::endl;
#endif

    if (verbose)
    {
        std::cout << " -----------------------------------------------------------------" << std::endl;
        std::cout << " | Read more about this program in the following publication:    |" << std::endl;
        std::cout << " |  Scheres ea.  (2009) Structure, 17, 1563-1572                 |" << std::endl;
        std::cout << " |                                                               |" << std::endl;
        std::cout << " |   *** Please cite it if this program is of use to you! ***    |" << std::endl;
        std::cout << " -----------------------------------------------------------------" << std::endl;

        std::cout << "Input  Metadata File " <<  fn_in       << std::endl;
        std::cout << "Output volume root   " <<  oroot      << std::endl;
        std::cout << "Output metadata file " <<  fn_out      << std::endl;
        std::cout << "Coordenate center    " <<  center_ref   << std::endl;
        std::cout << "Size subvolume       " <<  size         << std::endl;
        std::cout << "Symmetry group       " <<  fn_sym       << std::endl;
        std::cout << "Minimum distance between subvolumes " << mindist << std::endl;
    }
#ifdef DEBUG
    std::cerr << "end show" <<std::endl;
#endif
}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produce_Side_info2)
void ProgTomoExtractSubvolume::produceSideInfo()
{
	//do not write metadata

#ifdef  DEBUG
    std::cerr<<"Start produceSideInfo"<<std::endl;
#endif

    //    //Read in Docfile
    //    DF.read(fn_doc);
    //    nr_exp_images = DF.size();

    // Setup symmetry
    if (!SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
        REPORT_ERROR(ERR_PARAM_INCORRECT, "tomo_extract_subvolume: Invalid symmetry " +  fn_sym);
    SL.readSymmetryFile(fn_sym);

    // Check which symmetry operators give unique output points
    centers_subvolumes.clear();
    rotations_subvolumes.clear();
    Matrix2D<double>  L(4, 4), R(4, 4);
    Matrix2D<double>  I(3,3);
    Matrix1D<double>  newcenter(3), distcenter(3);
    double dist;
    bool is_uniq;
    I.initIdentity();

    centers_subvolumes.push_back(center_ref);
    rotations_subvolumes.push_back(I);

    for (int isym = 0; isym < SL.symsNo(); isym++)
    {
        SL.getMatrices(isym, L, R);
        L.resize(3,3);
        R.resize(3,3);
        newcenter = L * (center_ref.transpose() * R).transpose();
        is_uniq=true;
        for (int i = 0; i < centers_subvolumes.size(); i++)
        {
            distcenter = centers_subvolumes[i] - newcenter;
            dist = sqrt(distcenter.sum2());
            if (dist < mindist)
            {
                is_uniq=false;
                break;
            }
        }
        if (is_uniq)
        {
            centers_subvolumes.push_back(newcenter);
            rotations_subvolumes.push_back(R);
        }
    }
#ifdef  DEBUG
    std::cerr<<"End produceSideInfo"<<std::endl;
#endif

}
void ProgTomoExtractSubvolume::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
#ifdef  DEBUG
    std::cerr<<"Start processImages"<<std::endl;
#endif

    //    for (int imgno=imgno_start; imgno <imgno_end; imgno++)
    rowIn.getValue(MDL_ANGLEROT,rot);
    rowIn.getValue(MDL_ANGLETILT,tilt);
    rowIn.getValue(MDL_ANGLEPSI,psi);

    double auxD;
    rowIn.getValue(MDL_ORIGINX,auxD);
    XX(doccenter) = auxD;
    rowIn.getValue(MDL_ORIGINY,auxD);
    YY(doccenter) = auxD;
    rowIn.getValue(MDL_ORIGINZ,auxD);
    ZZ(doccenter) = auxD;

    // Read volume
    vol.read(fnImg);
    vol().setXmippOrigin();

    x0 = FIRST_XMIPP_INDEX(size);
    xF = LAST_XMIPP_INDEX(size);

    // Tomo_Extract each of the unique subvolumes
    size_t oId;
    for (int i = 0; i < centers_subvolumes.size(); i++)
    {
        center=centers_subvolumes[i];

        // 1. rotate center
        Euler_angles2matrix(-psi,-tilt,-rot,A,false);
        M3x3_BY_V3x1(center, A, center);
        // 2. translate center
        center -= doccenter;
        // 3. Apply possible non-integer center to volume
        //FIXME -center rob
        translate(BSPLINE3,volout(),vol(),center);

        //vol().translate(-center, volout(), DONT_WRAP);
        //4. Window operation and write subvolume to disc
        volout().selfWindow(x0,x0,x0,xF,xF,xF);
        fn_aux=fnImg.removeLastExtension();
        fn_aux+="_sub";
        fn_aux.compose(fn_aux,i+1,"vol");
        volout.write(fn_aux);


        // 5. Calculate output angles: apply symmetry rotation to rot,tilt and psi
        Euler_apply_transf(rotations_subvolumes[i], I, rot, tilt, psi, rotp, tiltp, psip);
        oId=DFout.addObject();

        DFout.setValue(MDL_IMAGE,fn_aux,oId);
        DFout.setValue(MDL_ANGLEROT,rotp,oId);
        DFout.setValue(MDL_ANGLETILT,tiltp,oId);
        DFout.setValue(MDL_ANGLEPSI,psip,oId);
    }
    // 6. Output translations will be zero because subvolumes are centered by definition
    DFout.setValueCol(MDL_ORIGINX,0.);
    DFout.setValueCol(MDL_ORIGINY,0.);
    DFout.setValueCol(MDL_ORIGINZ,0.);


#ifdef  DEBUG

    std::cerr<<"end processImages"<<std::endl;
#endif
}
void ProgTomoExtractSubvolume::postProcess()
{
    DFout.write(fn_out);
}
