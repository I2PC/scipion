#!/usr/bin/env python

import sys, os

class Tester:
    def __init__(self, fnDir):
        self.fnDir = fnDir
        self.lastProgram = ""
        self.progDict = {}
        self.addAllTests()

    def addProgram(self, program):
        self.progDict[program] = []
        self.lastProgram = program
        
    def addTest(self, test):
        self.progDict[self.lastProgram].append(test)
        
    def runProgramTests(self, program):
        tests = self.progDict[program]
        n = len(tests)
        outPath = os.path.join(self.fnDir, program)
        outDir = outPath
        testName = ""
        
        testNo = 1
        for test in tests:
            if n > 1:
                outDir = outPath + "_%02d" % testNo
                testName = "(%d of %d)" % (testNo, n)
            print "------------------------------------------------------------------------------------"
            print ">>> Running test", testName, "of", program
            print "    Output dir: "
            print "       ", outDir
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            test = test.replace("%o", outDir)
            test = test.replace("%p", program)
            test = test.replace("%d", self.fnDir)
            cmd = "%s %s > %s/stdout.txt 2> %s/stderr.txt" % (program, test, outDir, outDir)
            print "    Command: "
            print "       ", cmd
            os.system(cmd)
            testNo += 1
            
    
    def runAllTests(self):
        for program in self.progDict.keys():
            self.runProgramTests(program)
    
    def addAllTests(self):
        # Add all desired tests -------------------------------------------
    
        self.addProgram("xmipp_classify_analyze_cluster")
        self.addTest("-i input/smallStack.stk --ref 1@input/smallStack.stk -o %o/pca.xmd")

        self.addProgram("xmipp_image_align")
        self.addTest("-i input/smallStack.stk --oroot %o/aligned")

        self.addProgram("xmipp_image_convert")
        self.addTest("-i input/smallStack.stk -o %o/smallStack.mrcs -t stk")
    
        self.addProgram("xmipp_image_header")
        self.addTest("-i input/smallStack.stk --extract -o %o/header.doc")
        self.addTest("-i input/header.doc --assign -o %o/smallStack2.stk")
        
        self.addProgram("xmipp_image_statistics")
        self.addTest("-i input/smallStack.stk --image_stats %o/stats")
        
        self.addProgram("xmipp_metadata_utilities")
        self.addTest ("-i input/mD1.doc --set union input/mD2.doc  -o %o/out.doc")
        self.addTest ("-i input/mD1.doc --operate sort -o %o/out.doc")
        self.addTest ("-i input/mD1.doc --operate add_column \"shiftX shiftY\" -o %o/out.doc")
        self.addTest ("-i input/mD1.doc --fill shiftX rand_uniform 0 10 -o %o/out.doc")
        self.addTest ("-i input/mD1.doc -l \"shiftX shiftY\" constant 5 -o %o/out.doc")
        self.addTest ("-i input/mD1.doc --query select \"angleRot > 10 AND anglePsi < 0.5\" -o %o/out.doc")
        self.addTest("-i input/mD1.doc --operate modify_values \"angleRot=(angleRot*3.1416/180.)\" -o %o/out.doc")
        self.addTest("-i input/mD1.doc --operate modify_values \"image=replace(image, 'xmp','spi')\" -o %o/out.doc")
        
        self.addProgram("xmipp_ml_align2d")
        self.addTest("-i input/images_some.stk --ref input/seeds2.stk --oroot %o/ml2d --fast --mirror")
      
        self.addProgram("xmipp_phantom_project")
        self.addTest("-i input/phantomBacteriorhodopsin.vol -o %o/image.xmp --angles 0 0 0")
        self.addTest("-i input/phantomBacteriorhodopsin.vol     --oroot %o/projections --params input/clusterProjection.param")
        self.addTest("-i input/phantomBacteriorhodopsin.vol     --oroot %o/projections --params input/uniformProjection.param")
        self.addTest("-i input/Crystal/cylinder_with_axis.descr --oroot %o/MRCproj     --params input/Crystal/MRC_projection.param --crystal input/Crystal/MRC_crystal_projection.param")
    
        self.addProgram("xmipp_phantom_simulate_microscope")
        self.addTest("-i input/smallStack.stk -o %o/smallStackPlusCtf.stk --ctf input/input.ctfparam" )
    
        self.addProgram("xmipp_tomo_align_dual_tilt_series")
        self.addTest("--ref input/tomo_dual_alignment/ref.sel --dual input/tomo_dual_alignment/dual.sel --scale 1")

        self.addProgram("xmipp_tomo_align_tilt_series")
        self.addTest("-i input/tomo_dual_alignment/ref.sel --oroot %o/ref_aligned")

        self.addProgram("xmipp_tomo_align_refinement")
        self.addTest("--ref input/tomo_align_refinement/volume.vol --sel input/tomo_align_refinement/tilt_series.sel --oroot %o/refined --max_tilt_change 1 --max_rot_change 1 --tilt_step 1 --rot_step 1 --adjustGray")

        self.addProgram("xmipp_tomo_project")
        self.addTest("-i input/phantomCandida.vol -o %o/image.xmp --angles 0 90 90" )
        self.addTest("-i input/phantomCandida.vol --oroot %o/projections --params input/tomoProjection.param")
    
        self.addProgram("xmipp_transform_add_noise")
        self.addTest("-i input/cleanImage.spi --type gaussian 10 5 -o %o/noisyGaussian.spi")
    
        self.addProgram("xmipp_transform_adjust_volume_grey_levels")
        self.addTest("-i input/phantomCandida.vol -m goldStandard/xmipp_tomo_project_02/projections.sel -o %o/adjusted.vol")
    
        self.addProgram("xmipp_transform_center_image")
        self.addTest("-i input/smallStack.stk -o %o/smallStackCentered.stk")
    
        self.addProgram("xmipp_transform_geometry")
        self.addTest("-i input/header.doc --apply_transform -o %o/images.stk");
        self.addTest("-i input/phantomBacteriorhodopsin.vol --shift 10 5 -10 -o %o/volume.vol --dont_wrap");
        self.addTest("-i input/header.doc --scale factor 0.5 --oroot %o/halvedOriginal");
        self.addTest("-i input/header.doc --scale dim 64 --oroot %o/halvedDim");

        self.addProgram("xmipp_transform_normalize")
        self.addTest("-i input/smallStack.stk -o %o/smallStackNormalized.stk --method NewXmipp --background circle 32")

        self.addProgram("xmipp_transform_symmetrize")
        self.addTest("-i input/smallStack.stk -o %o/smallStackSymmetrized.stk --sym 5")
        self.addTest("-i input/phantomBacteriorhodopsin.vol -o %o/symmetrizedVolume.vol --sym C5")

        self.addProgram("xmipp_transform_window")
        self.addTest("-i input/singleImage.spi -o %o/image.xmp --size 32")
        self.addTest("-i input/singleImage.spi -o %o/image.xmp --corners -16 -16 15 15")
        self.addTest("-i input/singleImage.spi -o %o/image.xmp --corners 0 0 31 31 --physical")
        self.addTest("-i input/singleImage.spi -o %o/image.xmp --crop -10")
        self.addTest("-i input/xray_import/Images/img48949.spe -o %o/image.xmp --size 512")
    
        self.addProgram("xmipp_volume_segment")
        self.addTest("-i input/phantomBacteriorhodopsin.vol -o %o/maskOtsu.vol --method otsu")
        self.addTest("-i input/phantomBacteriorhodopsin.vol -o %o/maskVoxelMass.vol --method voxel_mass 54000")
        self.addTest("-i input/phantomBacteriorhodopsin.vol -o %o/maskProb.vol --method prob 1")

        self.addProgram("xmipp_xray_import")
        self.addTest("--data input/xray_import/Images --flat input/xray_import/Flatfields --oroot %o/stack --crop 30")
    
        self.addProgram("xmipp_xray_project")
        self.addTest("-i input/phantomCandida.vol -o %o/image.xmp --angles 0 90 90 -s 10 --psf input/xray_psf.xmd")
        self.addTest("-i input/phantomCandida.vol --oroot %o/projections --params input/tomoProjection.param -s 10 --psf input/xray_psf.xmd")
    
        self.addProgram("xmipp_xray_psf_create")
        self.addTest("-i input/xray_psf.xmd -o %o/psf.vol")    
#
# Main
#
if __name__ == '__main__':
    
    argc = len(sys.argv)
    
    if argc < 2 or argc > 3:
        print "Usage: ./batch.py <directory> [program]"
        sys.exit()

    fnDir = sys.argv[1]
    # Create tester
    tester = Tester(fnDir)
    
    if argc > 2:        
        program = sys.argv[2]
        tester.runProgramTests(program)
        if not os.path.exists(fnDir):
            os.makedirs(fnDir)
    else:
        # Remove the output directory if it is not goldStandard
        if fnDir != 'goldStandard' and fnDir != 'goldStandard/':
            if os.path.exists(fnDir):
                os.system("rm -rf " + fnDir)
            os.makedirs(fnDir)
        tester.runAllTests()
