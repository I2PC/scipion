import os
import shutil
import subprocess


sizes = {32,50,64,91,128,197,256,399,512,771,1024}
for s in sizes:
	folderName = str(s)+"x"+str(s)
	if os.path.isdir(folderName):
		shutil.rmtree(folderName)
	os.makedirs(folderName)
	with open(folderName + os.sep + "axes.descr", "w") as f:
		l1 = str(s) + " " + str(s) + " " + str(s) + " 0\n"
		l2 = "cyl + 1  0 0 " + str(5/64.*s) + "  " + str(2/64.*s) + " " + str(2/64.*s) + " " + str(10/64.*s) + "  0 0 0\n"
		l3 = "cyl + 1  " + str(10/64.*s) + " 0 0  " + str(4/64.*s) + " " + str(4/64.*s) + " " + str(20/64.*s) + "  0 90 0\n"
		l4 = "cyl + 1  0 " + str(13/64.*s) + " 0  " + str(6/64.*s) + " " + str(6/64.*s) + " " + str(25/64.*s) + "  90 90 0\n"
		l5 = "sph + 1  0 0 " + str(-15/64.*s) + "  " + str(3/64.*s)
		f.writelines([l1, l2, l3, l4, l5])
	with open(folderName + os.sep + "projParams.xmd", "w") as f:
		f.write("# XMIPP_STAR_1 *\n"
				+ "data_block1\n"
				+ "_dimensions2D   '" + str(s) + " " +str(s) + "'\n"
				+ "_projRotRange    '0 360 10'\n_projRotRandomness   even\n"
				+ "_projRotNoise   '0'\n_projTiltRange    '0 180 9'\n"
				+ "_projTiltRandomness   even\n_projTiltNoise   '0'\n"
				+ "_projPsiRange    '0 360 1'\n_projPsiRandomness   even\n"
				+ "_projPsiNoise   '0'\n_noisePixelLevel   '0'\n"
				+ "_noiseCoord   '0'");
	p = subprocess.Popen(["xmipp_phantom_create",
			"-i", "axes.descr",
			"-o", "axes.vol"], cwd=folderName)
	p.wait()
	p = subprocess.Popen(["xmipp_phantom_project",
			"-i", "axes.vol",
			"-o", "projAxes.stk",
			"--params", "projParams.xmd"], cwd=folderName)
	p.wait()
	p = subprocess.Popen(["xmipp_cuda_reconstruct_fourier",
			"-i", "projAxes.xmd",
			"-o", "maxwell_" + str(s) + ".vol",
			"--bufferSize", "50"], cwd=folderName)
	p.wait()
