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
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;
import static org.junit.Assert.*;

/**
 *
 * Test metadata 
 * 
 * sample metadata:
 * 
data_
loop_
 _image
 _enabled
 _shiftX
 _ref
 000001@tux.stk          1     0.000000   3 
 000002@tux.stk          1     1.404060   4 
 000003@tux.stk          1     1.715445   1 
 000004@tux.stk          1     1.266333   2 
 */
 
public class MetadataTest {
	static String mdFn = XmippTest.getTestFilename("tux.xmd");
	//Values expected from the metadata
	String[] imageValues = {  "000001@tux.stk"
                             ,"000002@tux.stk"
			                 ,"000003@tux.stk"
			                 ,"000004@tux.stk"};
	
	double[] shiftXValues = {  0.000000, 1.404060, 1.715445, 1.266333};
	double[] shiftXSorted = {  0.000000, 1.266333, 1.404060, 1.715445};
	int[] refValues = { 3, 4, 1, 2};
	int[] refSorted = { 4, 3, 2, 1};		
			
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

    @Test
    public void testRead() throws Exception{
    	MetaData md = new MetaData(mdFn);
    	long[] ids = md.findObjects();
    	long id;
    	for (int i=0; i < ids.length; ++i){
    		id = ids[i];
    		assertEquals(imageValues[i], md.getValueString(MDLabel.MDL_IMAGE, id));
    		assertEquals(shiftXValues[i], md.getValueDouble(MDLabel.MDL_SHIFTX, id), XmippTest.EQUAL_ACCURACY);
    		assertEquals(refValues[i], md.getValueInt(MDLabel.MDL_REF, id));
    	}
    }//function testRead
    
    @Test
    public void testLabelType() throws Exception{    	
    	int type = MetaData.getLabelType(MDLabel.MDL_IMAGE);
    	assertEquals(MetaData.LABEL_STRING, type);
    }//function testLabelType

    @Test
    public void testSort() throws Exception{
    	MetaData md = new MetaData(mdFn);
    	//Sort by MDL_SHIFTX ascending
    	md.sort(MDLabel.MDL_SHIFTX, true);
    	long [] ids = md.findObjects();
    	for (int i=0; i < ids.length; ++i)
    		assertEquals(shiftXSorted[i], md.getValueDouble(MDLabel.MDL_SHIFTX, ids[i]), XmippTest.EQUAL_ACCURACY);
    	//Sort by MDL_REF descending
    	md.sort(MDLabel.MDL_REF, false);
    	ids = md.findObjects();
    	for (int i=0; i < ids.length; ++i)
    		assertEquals(refSorted[i], md.getValueInt(MDLabel.MDL_REF, ids[i]));
    }//function testSort
    
    @Test
    public void testMakeAbsPath() throws Exception{
    	MetaData md = new MetaData(mdFn), md2 = new MetaData(mdFn);
    	md2.makeAbsPath(MDLabel.MDL_IMAGE);
    	long [] ids = md.findObjects();
    	long [] ids2 = md2.findObjects();
    	
    	for (int i=0; i < ids.length; ++i) {
    		String path = md.getValueString(MDLabel.MDL_IMAGE, ids[i]);
    		File f = new File(Filename.getFilename(path));
    		path = md2.getValueString(MDLabel.MDL_IMAGE, ids2[i]);
    		assertEquals(Filename.getFilename(path), f.getAbsolutePath());
    	}    	
    	md.print();
    	md2.print();
    	
    }//function testSort
    
    @Test
    public void testFill() throws Exception{
    	//Test fillConstant
    	String imagePath = "resources/test/singleImage.img";
    	MetaData md = new MetaData(mdFn);
    	md.fillConstant(MDLabel.MDL_ASSOCIATED_IMAGE1, imagePath);
    	long [] ids = md.findObjects();
    	for (int i=0; i < ids.length; ++i)
    		assertEquals(imagePath, md.getValueString(MDLabel.MDL_ASSOCIATED_IMAGE1, ids[i]));
    	
    	//Test fillLinear
    	md.fillLinear(MDLabel.MDL_SHIFTY, 0.0, 0.5);
    	double value = 0.0;
    	for (int i=0; i < ids.length; ++i) {
    		assertEquals(value, md.getValueDouble(MDLabel.MDL_SHIFTY, ids[i]), XmippTest.EQUAL_ACCURACY);
    		value += 0.5;
    	}
    	md.print();
    	
    	//Test fillRandom
    	md.fillRandom(MDLabel.MDL_SHIFTY, "uniform", 2.0, 3.0);
    	md.print();
    }//function testFill
    
    @Test
    public void testRemoveObject() throws Exception {
    	double[] localShiftXValues = {  shiftXValues[1], shiftXValues[3]};
    	MetaData md = new MetaData(mdFn);
    	long[] ids = md.findObjects();
    	md.removeObject(ids[0]);
    	md.removeObject(ids[2]);
    	assertEquals(2, md.size());
    	ids = md.findObjects();
    	for (int i=0; i < ids.length; ++i)
    		assertEquals(localShiftXValues[i], md.getValueDouble(MDLabel.MDL_SHIFTX, ids[i]), XmippTest.EQUAL_ACCURACY);
    	md.print();
    }
    
    @Test
    public void testReadQuotedString() throws Exception {
    	MetaData md = new MetaData("kk.xmd");
    	long[] ids = md.findObjects(); 
    	md.print();
    	System.out.println("------------------------------------");   	
    	for (int i=0; i < ids.length; ++i){
    		System.out.println(md.getValueString(MDLabel.MDL_IMAGE, ids[i]));
    		System.out.println(md.getValueString(MDLabel.MDL_ANGLE_COMPARISON, ids[i]));
    	}
    	//md.print();
    }  
    

}//class MetadataTest
