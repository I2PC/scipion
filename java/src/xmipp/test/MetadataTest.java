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

import java.io.File;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.windows.GalleryJFrame;
import static org.junit.Assert.*;

/**
 * 
 * Test metadata
 * 
 * sample metadata:
 * 
 * data_ loop_ _image _enabled _shiftX _ref 000001@tux.stk 1 0.000000 3
 * 000002@tux.stk 1 1.404060 4 000003@tux.stk 1 1.715445 1 000004@tux.stk 1
 * 1.266333 2
 */

public class MetadataTest
{
	static String mdFn = XmippTest.getTestFilename("tux.xmd");
	// Values expected from the metadata
	String[] imageValues = { "000001@tux.stk", "000002@tux.stk", "000003@tux.stk", "000004@tux.stk" };

	double[] shiftXValues = { 0.000000, 1.404060, 1.715445, 1.266333 };
	double[] shiftXSorted = { 0.000000, 1.266333, 1.404060, 1.715445 };
	int[] refValues = { 3, 4, 1, 2 };
	int[] refSorted = { 4, 3, 2, 1 };

	@BeforeClass
	public static void setUpClass() throws Exception
	{
	}

	@AfterClass
	public static void tearDownClass() throws Exception
	{
	}

	@Before
	public void setUp()
	{
	}

	@After
	public void tearDown()
	{
	}

	@Test
	public void testRead() throws Exception
	{
		MetaData md = new MetaData(mdFn);

		long[] ids = md.findObjects();
		long id;
		for (int i = 0; i < ids.length; ++i)
		{
			id = ids[i];
			assertEquals(imageValues[i], md.getValueString(MDLabel.MDL_IMAGE, id));
			assertEquals(shiftXValues[i], md.getValueDouble(MDLabel.MDL_SHIFT_X, id), XmippTest.EQUAL_ACCURACY);
			assertEquals(refValues[i], md.getValueInt(MDLabel.MDL_REF, id));
		}
		md.destroy();
	}// function testRead

	@Test
	public void testLabelType() throws Exception
	{
		int type = MetaData.getLabelType(MDLabel.MDL_IMAGE);
		assertEquals(MetaData.LABEL_STRING, type);
	}// function testLabelType

	@Test
	public void testSort() throws Exception
	{
		MetaData md = new MetaData(mdFn);
		// Sort by MDL_SHIFT_X ascending
		md.sort(MDLabel.MDL_SHIFT_X, true);
		long[] ids = md.findObjects();
		for (int i = 0; i < ids.length; ++i)
			assertEquals(shiftXSorted[i], md.getValueDouble(MDLabel.MDL_SHIFT_X, ids[i]), XmippTest.EQUAL_ACCURACY);
		// Sort by MDL_REF descending
		md.sort(MDLabel.MDL_REF, false);
		ids = md.findObjects();
		for (int i = 0; i < ids.length; ++i)
			assertEquals(refSorted[i], md.getValueInt(MDLabel.MDL_REF, ids[i]));
		md.destroy();
	}// function testSort

	@Test
	public void testMakeAbsPath() throws Exception
	{
		MetaData md = new MetaData(mdFn), md2 = new MetaData(mdFn);
		md2.makeAbsPath(MDLabel.MDL_IMAGE);
		long[] ids = md.findObjects();
		long[] ids2 = md2.findObjects();

		for (int i = 0; i < ids.length; ++i)
		{
			String path = md.getValueString(MDLabel.MDL_IMAGE, ids[i]);
			File f = new File(Filename.getFilename(path));
			path = md2.getValueString(MDLabel.MDL_IMAGE, ids2[i]);
			assertEquals(Filename.getFilename(path), f.getAbsolutePath());
		}
		md.print();
		md2.print();

		md.destroy();
		md2.destroy();

	}// function testSort

	@Test
	public void testFill() throws Exception
	{
		// Test fillConstant
		String imagePath = "resources/test/singleImage.img";
		MetaData md = new MetaData(mdFn);
		md.fillConstant(MDLabel.MDL_IMAGE1, imagePath);
		long[] ids = md.findObjects();
		for (int i = 0; i < ids.length; ++i)
			assertEquals(imagePath, md.getValueString(MDLabel.MDL_IMAGE1, ids[i]));

		// Test fillLinear
		md.fillLinear(MDLabel.MDL_SHIFT_Y, 0.0, 0.5);
		double value = 0.0;
		for (int i = 0; i < ids.length; ++i)
		{
			assertEquals(value, md.getValueDouble(MDLabel.MDL_SHIFT_Y, ids[i]), XmippTest.EQUAL_ACCURACY);
			value += 0.5;
		}
		// md.print();

		// Test fillRandom
		// md.fillRandom(MDLabel.MDL_SHIFT_Y, "uniform", 2.0, 3.0);
		// md.print();
		md.destroy();
	}// function testFill

	@Test
	public void testRemoveObject() throws Exception
	{
		double[] localShiftXValues = { shiftXValues[1], shiftXValues[3] };
		MetaData md = new MetaData(mdFn);
		long[] ids = md.findObjects();
		md.removeObject(ids[0]);
		md.removeObject(ids[2]);
		assertEquals(2, md.size());
		ids = md.findObjects();
		for (int i = 0; i < ids.length; ++i)
			assertEquals(localShiftXValues[i], md.getValueDouble(MDLabel.MDL_SHIFT_X, ids[i]), XmippTest.EQUAL_ACCURACY);
		// md.print();
		md.destroy();
	}

	@Test
	public void testReadQuotedString() throws Exception
	{
		// MetaData md = new MetaData("kk.xmd");
		// long[] ids = md.findObjects();
		// md.print();
		// System.out.println("------------------------------------");
		// for (int i=0; i < ids.length; ++i){
		// System.out.println(md.getValueString(MDLabel.MDL_IMAGE, ids[i]));
		// System.out.println(md.getValueString(MDLabel.MDL_ANGLE_DIFF,
		// ids[i]));
		// }
		// md.print();
	}

	// @Test
	// public void createBigMetadata() throws Exception {
	// MetaData md = new MetaData();
	// int n = 1000000;
	// long id;
	//
	// DEBUG.tic();
	// System.out.println("Creating big md....");
	// for (int i = 0; i < n; ++i){
	// id = md.addObject();
	// md.setValueString(MDLabel.MDL_IMAGE, String.format("image%06d.xmp", i),
	// id);
	// }
	// DEBUG.toc();
	// DEBUG.tic();
	// System.out.println("Filling big md....");
	// md.fillLinear(MDLabel.MDL_ANGLE_PSI, 0.0, 0.5);
	// md.fillRandom(MDLabel.MDL_SHIFT_X, "gaussian", 0., 3.);
	// md.fillRandom(MDLabel.MDL_SHIFT_Y, "gaussian", 0., 3.);
	// DEBUG.toc();
	// md.write("big_md.xmd");
	// }

	// public void createBigMetadata(MetaData md, int n) {
	// System.out.println("Creating big md....");
	// long id;
	// double v = 0.;
	//
	// for (int i = 0; i < n; ++i) {
	// id = md.addObject();
	// md.setValueString(MDLabel.MDL_IMAGE,
	// String.format("image%06d.xmp", i), id);
	// md.setValueDouble(MDLabel.MDL_ANGLE_PSI, v, id);
	// md.setValueDouble(MDLabel.MDL_SHIFT_X, 100.0, id);
	// md.setValueDouble(MDLabel.MDL_SHIFT_Y, 100.0, id);
	// v+=0.5;
	// }
	// }

//	@Test
//	public void readBigMetadata() throws Exception
//	{
//		try
//		{
//
//			DEBUG.tic();
//			System.out.println("Creating big md....");
//			MetaData md = new MetaData("javaMdTest/big_md.xmd");
//			// MetaData md = new MetaData();
//			// createBigMetadata(md, 1000000);
//			DEBUG.toc();
//			System.out.println("Finding objects big md....");
//			long[] ids = md.findObjects();
//			DEBUG.toc();
//			long counter = 0, total = ids.length;
//			String image;
//			double angle = 0, x = 0, y = 0;
//
//			System.out.println("Looping objects big md....");
//			for (long id : ids)
//			{
//				image = md.getValueString(MDLabel.MDL_IMAGE, id);
//				angle += md.getValueDouble(MDLabel.MDL_ANGLE_PSI, id);
//				x += md.getValueDouble(MDLabel.MDL_SHIFT_X, id);
//				y += md.getValueDouble(MDLabel.MDL_SHIFT_Y, id);
//				if (++counter % 10000 == 0)
//					System.out.format("%d of %d\n", counter, total);
//			}
//			System.out.format("sums angle: %f x: %f y: %f %n", angle, x, y);
//			md.destroy();
//		}
//		catch (Exception ex)
//		{
//			ex.printStackTrace();
//		}
//	}

	public class Worker implements Runnable
	{
		MetaData imagesmd;
		MetaData md;
		
		public Worker(MetaData md,MetaData imagesmd)
		{
			this.imagesmd = imagesmd;
			this.md = md;
		}

		public void run()
		{
			String imagepath;
			int idlabel = MDLabel.MDL_IMAGE;
			
			long id2;
			for (long id : md.findObjects())
			{
				imagepath = md.getValueString(idlabel, id, true);
				if (imagepath != null && ImageGeneric.exists(imagepath))
				{
					id2 = imagesmd.addObject();
					imagesmd.setValueString(idlabel, imagepath, id2);
				}
			}
			
			for (long id : imagesmd.findObjects())
			{
				imagepath = imagesmd.getValueString(idlabel, id, true);
				System.out.printf("%d %s\n", id, imagepath);
			}
			
			imagesmd.destroy();
		}
	}
	@Test
	public void readRuntimeMd() throws Exception
	{
		try
		{
			System.out.println("read runtime md...");
			
			MetaData md = new MetaData(XmippTest.getTestFilename("images.stk"));
			MetaData imagesmd = new MetaData();
			String imagepath;
			int idlabel = MDLabel.MDL_IMAGE;
			

			
			Worker w = new Worker(md, imagesmd);
			//w.run();
			
			Thread th = new Thread(w);
			th.start();

			md.destroy();
			
			System.out.println("read runtime md ended...");
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}
	
	@Test
	public void updateMetaDataRowWithRow() throws Exception
	{
		try
		{
			System.out.println("updateMetaDataRowWithRow");
			
			MetaData md = new MetaData(XmippTest.getTestFilename("images.stk"));
			MetaData imagesmd = new MetaData();
			String imagepath;
			int idlabel = MDLabel.MDL_IMAGE;
			MetaData mdRow = new MetaData();
			
			long id2;
			for (long id : md.findObjects())
			{
				imagepath = md.getValueString(idlabel, id, true);
				if (imagepath != null && ImageGeneric.exists(imagepath))
				{
					id2 = imagesmd.addObject();
					md.getRow(id, mdRow);
					mdRow.setValueString(idlabel, imagepath, mdRow.firstObject());
					imagesmd.setRow(mdRow, id2);
				}
			}
			
			imagesmd.print();
			
			imagesmd.destroy();
			md.destroy();
			mdRow.destroy();
			
			System.out.println("updateMetaDataRowWithRow ended...");
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

//	@Test
//	public void readTinyMetadatas() throws Exception
//	{
//		try
//		{
//			MetaData md = new MetaData("javaMdTest/all_tiny.xmd");
//			MetaData md2 = new MetaData();
//
//			long[] ids = md.findObjects();
//			long[] ids2;
//			String tinyFn;
//			double angle, x, y;
//
//			for (long id : ids)
//			{
//				tinyFn = md.getValueString(MDLabel.MDL_IMAGE, id);
//				md2.read(tinyFn);
//				System.out.println("Tiny md: " + tinyFn);
//				ids2 = md2.findObjects();
//				for (long id2 : ids2)
//				{
//					angle = md2.getValueDouble(MDLabel.MDL_ANGLE_PSI, id2);
//					x = md2.getValueDouble(MDLabel.MDL_SHIFT_X, id2);
//					y = md2.getValueDouble(MDLabel.MDL_SHIFT_Y, id2);
//					// System.out.format("    image: %s, angle: %f x: %f y: %f %n",
//					// tinyFn, angle, x, y);
//				}
//			}
//			md.destroy();
//			md2.destroy();
//		}
//		catch (Exception ex)
//		{
//			ex.printStackTrace();
//		}
//	}

	@Test
	public void testOperate() throws Exception
	{
		MetaData md = new MetaData(mdFn);
		md.operate("shiftX=2*shiftX, ref=ref+10");
		long[] ids = md.findObjects();
		long id;
		for (int i = 0; i < ids.length; ++i)
		{
			id = ids[i];
			assertEquals(shiftXValues[i] * 2, md.getValueDouble(MDLabel.MDL_SHIFT_X, id), XmippTest.EQUAL_ACCURACY);
			assertEquals(refValues[i] + 10, md.getValueInt(MDLabel.MDL_REF, id));
		}
	}// function testRead

}// class MetadataTest
