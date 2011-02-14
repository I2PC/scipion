#!/usr/bin/env python
class Tester:
    def __init__(self, fnDir):
        self.fnDir = fnDir

    def testProgram(self, program, arguments, testNo=0):
        import os
        print "------------------------------------------------------------------------------------"
        print ">>> Testing " + program
        outDir = os.path.join(self.fnDir, program)
        if testNo!=0:
            outDir+="_%02d"%testNo
        print "   Making output directory: ", outDir
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        cmd = "%s %s > %s/stdout.txt 2> %s/stderr.txt" % (program, arguments, outDir, outDir)
        print "   Running command: ", cmd
        os.system(cmd)

#
# Main
#
import os,sys
if __name__ == '__main__':
    if not sys.argv[1:] or len(sys.argv)<=1:
        print "Usage: ./batch.py <directory>"
        sys.exit()
    args = sys.argv[1:]
    fnDir=args[0]

    # Remove the output directory if it is not goldStandard
    if fnDir!='goldStandard' and fnDir!='goldStandard/':
        if os.path.exists(fnDir):
            os.system("rm -rf "+fnDir)
        os.makedirs(fnDir)

    # Create tester
    tester=Tester(fnDir)

    # Test the programs -------------------------------------------

    program = "xmipp_convert_image"
    tester.testProgram(program, "-i input/smallStack.stk -o %s/%s/smallStack.mrcs -t stk" % (fnDir, program))

    program = "xmipp_header_extract"
    tester.testProgram(program, "-i input/smallStack.stk -o %s/%s/header.doc" % (fnDir, program))

    program = "xmipp_phantom_project"
    tester.testProgram(program, "-i input/phantomBacteriorhodopsin.vol -o %s/%s_%02d/image.xmp --angles 0 0 0" % (fnDir, program,1),1)
    tester.testProgram(program, "-i input/phantomBacteriorhodopsin.vol --oroot %s/%s_%02d/projections --params input/clusterProjection.param" % (fnDir, program,2),2)
    tester.testProgram(program, "-i input/phantomBacteriorhodopsin.vol --oroot %s/%s_%02d/projections --params input/uniformProjection.param" % (fnDir, program,3),3)

    program = "xmipp_phantom_simulate_microscope"
    tester.testProgram(program, "-i input/smallStack.stk -o %s/%s/smallStackPlusCtf.stk --ctf input/input.ctfparam" % (fnDir, program))

    program = "xmipp_tomo_project"
    tester.testProgram(program, "-i input/phantomCandida.vol -o %s/%s_%02d/image.xmp --angles 0 90 90" % (fnDir, program,1),1)
    tester.testProgram(program, "-i input/phantomCandida.vol --oroot %s/%s_%02d/projections --params input/tomoProjection.param" % (fnDir, program,2),2)

    program = "xmipp_transform_add_noise"
    tester.testProgram(program, "-i input/cleanImage.spi --type gaussian 10 5 -o %s/%s/noisyGaussian.spi" % (fnDir, program))

    program = "xmipp_transform_center_image"
    tester.testProgram(program, "-i input/smallStack.stk -o %s/%s/smallStackCentered.stk" % (fnDir, program))

    program = "xmipp_transform_window"
    tester.testProgram(program, "-i input/singleImage.spi -o %s/%s_%02d/image.xmp --size 32" % (fnDir, program,1),1)
    tester.testProgram(program, "-i input/singleImage.spi -o %s/%s_%02d/image.xmp --corners -16 -16 15 15" % (fnDir, program,2),2)
    tester.testProgram(program, "-i input/singleImage.spi -o %s/%s_%02d/image.xmp --corners 0 0 31 31 --physical" % (fnDir, program,3),3)
    tester.testProgram(program, "-i input/singleImage.spi -o %s/%s_%02d/image.xmp --crop -10" % (fnDir, program,4),4)
    tester.testProgram(program, "-i input/xray_import/Images/img48949.spe -o %s/%s_%02d/image.xmp --size 512" % (fnDir, program,5),5)

    program = "xmipp_xray_import"
    tester.testProgram(program, "--data input/xray_import/Images --flat input/xray_import/Flatfields --oroot %s/%s/stack --crop 30" % (fnDir, program))

    program = "xmipp_xray_project"
    tester.testProgram(program, "-i input/phantomCandida.vol -o %s/%s_%02d/image.xmp --angles 0 90 90 -s 10 --psf input/xray_psf.xmd" % (fnDir, program,1),1)
    tester.testProgram(program, "-i input/phantomCandida.vol --oroot %s/%s_%02d/projections --params input/tomoProjection.param -s 10 --psf input/xray_psf.xmd" % (fnDir, program,2),2)

    program = "xmipp_xray_psf_create"
    tester.testProgram(program, "-i input/xray_psf.xmd -o %s/%s/psf.vol" % (fnDir, program))
