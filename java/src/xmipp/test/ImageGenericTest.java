/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 * 				Roberto Marabini       (roberto@cnb.csic.es)
 *
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
package xmipp.test;

import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.process.ImageProcessor;

import java.awt.Rectangle;
import java.io.File;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import xmipp.ij.commons.XmippImageConverter;
import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MetaData;
import static org.junit.Assert.*;
import java.util.Arrays;

/**
 * 
 * @author Juanjo Vega
 */
public class ImageGenericTest {

	// Filenames to test.
	static String filename = XmippTest.getTestFilename("singleImage.spi");
	static String filename2 = XmippTest.getTestFilename("singleImage2.spi");
	// Indexes for arrays.
	final static int X = 0;
	final static int Y = 1;
	final static int Z = 2;
	// Predefined dimension of 'singleImage.spi' to test
	int xdim = 3;
	int ydim = 3;
	int zdim = 1;
	int ndim = 1;
	// Image pixels values
	float[] imagePixels = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f,
			9.0f };
	// Images statistics.
	// min, max, avg, std dev
	double[] imageStats = { 1.0, 9.0, 5.0, 2.5820 };
	// Image datatype
	int imageDatatype = ImageGeneric.Float;

	public ImageGenericTest() {

	}

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	@Before
	public void setUp() {
	}

	@After
	public void tearDown() {
	}

	/**
	 * Test of getFilename method, of class ImageGeneric.
	 */
	@Test
	public void testGetFilename() throws Exception {
		ImageGeneric image = new ImageGeneric(filename);
		String imageFn = image.getFilename();
		assertEquals(filename, imageFn);
	}

	/**
	 * Test of getXDim method, of class ImageGeneric.
	 */
	@Test
	public void testGetDimensions() throws Exception {
		ImageGeneric image = new ImageGeneric(filename);
		assertEquals(xdim, image.getXDim());
		assertEquals(ydim, image.getYDim());
		assertEquals(zdim, image.getZDim());
		assertEquals(ndim, image.getNDim());
	}

	/**
	 * Test of getDataType method, of class ImageGeneric.
	 */
	@Test
	public void testGetDataType() throws Exception {
		ImageGeneric image = new ImageGeneric(filename);
		assertEquals(imageDatatype, image.getDataType());
	}

	/**
	 * Test of getArrayFloat method, of class ImageGeneric.
	 */
	@Test
	public void testGetArrayFloat() throws Exception {
		ImageGeneric image = XmippTest.getImageGeneric(filename);

		float[] result = image.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.FIRST_SLICE);
		for (int i = 0; i < imagePixels.length; i++)
			assertEquals(imagePixels[i], result[i], XmippTest.EQUAL_ACCURACY);
	}

	/**
	 * Test of getDataType method, of class ImageGeneric.
	 */
	@Test
	public void testEqual() throws Exception  {
		ImageGeneric image = XmippTest.getImageGeneric(filename);
		ImageGeneric image2 = XmippTest.getImageGeneric(filename2);

		assertTrue(image.equal(image, XmippTest.EQUAL_ACCURACY));
		assertFalse(image.equal(image2, XmippTest.EQUAL_ACCURACY));
	}

	/**
	 * Test of subtract method, of class ImageGeneric.
	 */
	@Test
	public void testSubtract() throws Exception  {
		ImageGeneric image = XmippTest.getImageGeneric(filename);
		ImageGeneric image2 = XmippTest.getImageGeneric(filename2);
		ImageGeneric image3 = new ImageGeneric();
		image2.subtract(image,image3);
		float[] result = image3.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.FIRST_SLICE);
		for (int i = 0; i < imagePixels.length -1; i++)
			assertEquals(result[i],0., XmippTest.EQUAL_ACCURACY);
		assertEquals(result[8],27.0, XmippTest.EQUAL_ACCURACY);
	}

	/**
	 * Test of getStatistics method, of class ImageGeneric.
	 */
	@Test
	public void testGetStatistics() throws Exception {
		ImageGeneric image = XmippTest.getImageGeneric(filename);
		double[] stats = image.getStatistics();
		for (int i = 0; i < stats.length; i++) {
			assertEquals(imageStats[i], stats[i], XmippTest.EQUAL_ACCURACY);
		}
	}
	
	/**
	 * Test of getStatistics method, of class ImageGeneric.
	 */
	@Test
	public void testSetArrayFloat() throws Exception {
		System.out.println("TestSetArrayFloat");
		int size = 100;
		String templatesfile = "EmptyImageGeneric.stk";
		ImageGeneric templates = new ImageGeneric(ImageGeneric.Float);
		templates.resize(size, size, 1, 1);
		
		templates.write(templatesfile);
		//templates.setFilename(templatesfile);
		
		String micrographfile = XmippTest.getTestFilename("BPV_1386.mrc");
		ImagePlus mimage = new ImagePlus(micrographfile);
		Rectangle r = new Rectangle(0, 0, size, size);
		mimage.setRoi(r);
		ImageProcessor processor = mimage.getProcessor().crop();
		ImagePlus particleimg = new ImagePlus("particle", processor);
		//new FileSaver(particleimg).save();
		ImageGeneric particleig = XmippImageConverter.convertToImageGeneric(particleimg);
		float[] matrix = particleig.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.ALL_SLICES);
		templates.setArrayFloat(matrix, ImageGeneric.FIRST_IMAGE, ImageGeneric.ALL_SLICES);
		templates.write(templatesfile);
	}

	/**
	 * Test of isPSD method, of class ImageGeneric.
	 */
	// /@Test
	public void testWhatIs() throws Exception {
		ImageGeneric image = new ImageGeneric(filename);
		assertTrue(image.isSingleImage());
		assertFalse(image.isStackOrVolume());
		assertFalse(image.isPSD());
	}

}//class ImageGenericTest
