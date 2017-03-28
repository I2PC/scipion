 /* Authors:     Roberto Marabini (roberto@cnb.csic.es)*/

#include <stdlib.h>
#include <data/xmipp_image.h>
#include <data/xmipp_image_extension.h>
#include <iostream>
#include <gtest/gtest.h>
#include <data/phantom.h>
#include <data/xmipp_program.h>
#include <cstdio>
#include <unistd.h>

const char ico_i2[] = "# Phantom description file, (generated with phantom help)\n\
# General Volume Parameters: \n\
#      Xdim      Ydim      Zdim   Background_Density Scale \n\
     5 5 5 0 60  \n\
# Feature Parameters: \n\
#Type  +/=  Density X_Center Y_Center Z_Center \n\
sph  + 1.   1.000 0.000 1.618 .15\n\
sph  + 1.   1.000 0.000 -1.618 .15\n\
sph  + 1.   -1.000 0.000 1.618 .15\n\
sph  + 1.   -1.000 0.000 -1.618 .15\n\
sph  + 1.   0.000 1.618 1.000 .15\n\
sph  + 1.   0.000 -1.618 1.000 .15\n\
sph  + 1.   0.000 1.618 -1.000 .15\n\
sph  + 1.   0.000 -1.618 -1.000 .15\n\
sph  + 1.   1.618 1.000 0.000 .15\n\
sph  + 1.   -1.618 1.000 0.000 .15\n\
sph  + 1.   1.618 -1.000 0.000 .15\n\
sph  + 1.   -1.618 -1.000 0.000 .15\n\
sph  + 1.   0.000 -0.539 1.412 .10\n\
sph  + 1.   0.000 0.539 1.412 .10\n\
sph  + 1.   0.873 -0.873 0.873 .10\n\
sph  + 1.   0.873 0.873 0.873 .10\n\
sph  + 1.   1.412 0.000 0.539 .10\n\
sph  + 1.   0.000 0.539 -1.412 .10\n\
sph  + 1.   0.873 0.873 -0.873 .10\n\
sph  + 1.   0.000 -0.539 -1.412 .10\n\
sph  + 1.   1.412 0.000 -0.539 .10\n\
sph  + 1.   0.873 -0.873 -0.873 .10\n\
sph  + 1.   -0.873 0.873 0.873 .10\n\
sph  + 1.   -1.412 0.000 0.539 .10\n\
sph  + 1.   -0.873 -0.873 0.873 .10\n\
sph  + 1.   -0.873 0.873 -0.873 .10\n\
sph  + 1.   -1.412 0.000 -0.539 .10\n\
sph  + 1.   -0.873 -0.873 -0.873 .10\n\
sph  + 1.   -0.539 1.412 0.000 .10\n\
sph  + 1.   0.539 1.412 0.000 .10\n\
sph  + 1.   -0.539 -1.412 0.000 .10\n\
sph  + 1.   0.539 -1.412 0.000 .10\n\
sph  + 1.   0.000 0.000 1.618 .10\n\
sph  + 1.   -0.500 -0.809 1.309 .10\n\
sph  + 1.   0.500 -0.809 1.309 .10\n\
sph  + 1.   0.500 0.809 1.309 .10\n\
sph  + 1.   -0.500 0.809 1.309 .10\n\
sph  + 1.   0.000 0.000 1.618 .10\n\
sph  + 1.   0.500 -0.809 1.309 .10\n\
sph  + 1.   0.809 -1.309 0.500 .10\n\
sph  + 1.   1.309 -0.500 0.809 .10\n\
sph  + 1.   1.309 0.500 0.809 .10\n\
sph  + 1.   0.809 1.309 0.500 .10\n\
sph  + 1.   0.500 0.809 1.309 .10\n\
sph  + 1.   1.309 -0.500 0.809 .10\n\
sph  + 1.   1.618 0.000 0.000 .10\n\
sph  + 1.   1.309 0.500 0.809 .10\n\
sph  + 1.   0.000 0.000 -1.618 .10\n\
sph  + 1.   -0.500 0.809 -1.309 .10\n\
sph  + 1.   0.500 0.809 -1.309 .10\n\
sph  + 1.   0.500 0.809 -1.309 .10\n\
sph  + 1.   0.809 1.309 -0.500 .10\n\
sph  + 1.   1.309 0.500 -0.809 .10\n\
sph  + 1.   0.500 -0.809 -1.309 .10\n\
sph  + 1.   -0.500 -0.809 -1.309 .10\n\
sph  + 1.   0.000 0.000 -1.618 .10\n\
sph  + 1.   1.309 0.500 -0.809 .10\n\
sph  + 1.   1.618 0.000 0.000 .10\n\
sph  + 1.   1.309 -0.500 -0.809 .10\n\
sph  + 1.   1.309 -0.500 -0.809 .10\n\
sph  + 1.   0.809 -1.309 -0.500 .10\n\
sph  + 1.   0.500 -0.809 -1.309 .10\n\
sph  + 1.   -0.500 0.809 1.309 .10\n\
sph  + 1.   -0.809 1.309 0.500 .10\n\
sph  + 1.   -1.309 0.500 0.809 .10\n\
sph  + 1.   -1.309 0.500 0.809 .10\n\
sph  + 1.   -1.618 0.000 0.000 .10\n\
sph  + 1.   -1.309 -0.500 0.809 .10\n\
sph  + 1.   -1.309 -0.500 0.809 .10\n\
sph  + 1.   -0.809 -1.309 0.500 .10\n\
sph  + 1.   -0.500 -0.809 1.309 .10\n\
sph  + 1.   -0.500 0.809 -1.309 .10\n\
sph  + 1.   -0.809 1.309 -0.500 .10\n\
sph  + 1.   -1.309 0.500 -0.809 .10\n\
sph  + 1.   -1.309 0.500 -0.809 .10\n\
sph  + 1.   -1.618 0.000 0.000 .10\n\
sph  + 1.   -1.309 -0.500 -0.809 .10\n\
sph  + 1.   -1.309 -0.500 -0.809 .10\n\
sph  + 1.   -0.809 -1.309 -0.500 .10\n\
sph  + 1.   -0.500 -0.809 -1.309 .10\n\
sph  + 1.   0.000 1.618 0.000 .10\n\
sph  + 1.   -0.809 1.309 -0.500 .10\n\
sph  + 1.   -0.809 1.309 0.500 .10\n\
sph  + 1.   0.809 1.309 0.500 .10\n\
sph  + 1.   0.809 1.309 -0.500 .10\n\
sph  + 1.   0.000 1.618 0.000 .10\n\
sph  + 1.   0.000 -1.618 0.000 .10\n\
sph  + 1.   -0.809 -1.309 -0.500 .10\n\
sph  + 1.   -0.809 -1.309 0.500 .10\n\
sph  + 1.   0.809 -1.309 0.500 .10\n\
sph  + 1.   0.809 -1.309 -0.500 .10\n\
sph  + 1.   0.000 -1.618 0.000 .10\n";


class TransformWindowTest : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
        XMIPP_TRY
        ;
        XMIPP_CATCH
    }
};


TEST_F( TransformWindowTest, unitcell)
{
	//create phantom volume
	char filename[] = "/tmp/mytemp.XXXXXX"; // template for our file.
    int fd = mkstemp(filename);    // Creates and opens a new temp file r/w.
	                                 // Xs are replaced with a unique number.
#define DEBUG
#ifdef DEBUG
    std::cerr << "filename: " << filename << std::endl;
#endif
    if (fd == -1) {
    	std::cerr << " cannot open temporary file" << std::endl;
    	ASSERT_TRUE(false);
    	return;
    }
    write(fd, ico_i2,sizeof(ico_i2));
	Phantom           phantom;
    Image<double>     vol;
    const char  inFn[] = "/tmp/inTransformWindowTest.mrc";
    const char  outFn[] ="/tmp/outTransformWindowTest.mrc";
    XMIPP_TRY
		phantom.read(filename);
		phantom.draw_in(vol());
		vol().addNoise(0,0.1,"gaussian");
		vol.write(inFn);
	XMIPP_CATCH
	//transform window is not in a library as it should be
	//so I need to make a system call
    char commandBuff[180];

	snprintf(commandBuff, sizeof(commandBuff),
			"xmipp_transform_window -i %s -o %s --unitcell i2 80 140 .25 0",
			inFn, outFn);
    system(commandBuff);
#ifdef DEBUG
    std::cerr << "commandBuff: " << commandBuff << std::endl;
#endif
#ifndef DEBUG
	unlink(inFn);
	unlink(outFn);
	unlink(filename);
#endif
	ASSERT_TRUE(true);
#undef DEBUG
}

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
