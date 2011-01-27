#!/usr/bin/env python
class Tester:
    def __init__(self,fnDir):
        self.fnDir=fnDir

    def testProgram(self,program,arguments):
        print "Testing "+program
        if not os.path.exists(self.fnDir+"/"+program):
            os.makedirs(self.fnDir+"/"+program)
        os.system(program+" "+arguments+\
            " > "+self.fnDir+"/"+program+"/stdout.txt 2>"+\
            self.fnDir+"/"+program+"/stderr.txt")

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

    # Test the programs
    program="xmipp_xray_psf_create"
    tester.testProgram(program,
        "-i input/xray_psf.xmd -o "+fnDir+"/"+program+"/psf.vol")
