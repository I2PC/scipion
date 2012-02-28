/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.jni;

import java.io.File;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import java.util.Arrays;

/**
 *
 * @author Juanjo Vega
 */
public class ImageGeneric_Test {

    // Xmipp dir
    static String XMIPP_HOME;
    final static String TESTS_PATH = "applications/tests/test_java/";
    // Filenames to test.
    String filenames[] = new String[]{
        "singleImage.spi"
    };
    String filename  = "singleImage.spi";
    String filename2 = "singleImage2.spi";
    // Indexes for arrays.
    final static int X = 0;
    final static int Y = 1;
    final static int Z = 2;
    // Images dimensions to test.
    int dimensions[][] = new int[][]{
        // X Y Z
        {3, 3, 1}
    };
    long N[] = new long[]{
        1
    };
    // Images statistics.
    double statistics[][] = new double[][]{
        // min, max, avg, std dev
        {1.0, 9.0, 5.0, 2.5820}
    };
    int dataType[] = new int[]{
        8
    };

    public ImageGeneric_Test() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
        findXmippHome();
        fixFilenamesPaths();
    }

    void findXmippHome() {
        XMIPP_HOME = System.getenv("XMIPP_HOME");

        if (XMIPP_HOME == null) {
            XMIPP_HOME = System.getProperty("xmipp_home");
        }

        XMIPP_HOME += File.separator + TESTS_PATH;
    }

    void fixFilenamesPaths() {
        for (int i = 0; i < filenames.length; i++) {
            filenames[i] = XMIPP_HOME + filenames[i];
        }
        filename  = XMIPP_HOME + filename;
        filename2 = XMIPP_HOME + filename2;
    }

    @After
    public void tearDown() {
        findXmippHome();
        fixFilenamesPaths();

        for (int i = 0; i < filenames.length; i++) {
            String string = filenames[i];
            System.out.println(i + ": " + string);
        }
    }

    /**
     * Test of setDimensions method, of class ImageGeneric.
     */
//     @Test
//     public void testtestResize() {
//         for (int i = 0; i < dimensions.length; i++) {
//             try {
//                 int xDim = dimensions[i][0];
//                 int yDim = dimensions[i][1];
//                 int zDim = dimensions[i][2];
// 
//                 testResize(xDim, yDim, zDim);
//             } catch (Exception ex) {
//                 fail("testGetFilename(): " + filenames[i]);
//             }
//         }
//     }
// 
//     void testResize(int xDim, int yDim, int zDim) throws Exception {
//         ImageGeneric instance = new ImageGeneric();
//         instance.resize(xDim, yDim, zDim);
// 
//         int xResult = instance.getXDim();
//         int yResult = instance.getYDim();
//         int zResult = instance.getZDim();
// 
//         assertEquals(xResult, xDim);
//         assertEquals(yResult, yDim);
//         assertEquals(zResult, zDim);
//     }

    /**
     * Test of getFilename method, of class ImageGeneric.
     */
    @Test
    public void testGetFilename() {
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetFilename(filenames[i]);
            } catch (Exception ex) {
                fail("testGetFilename(): " + filenames[i]);
            }
        }
    }

    void testGetFilename(String filename) throws Exception {
        ImageGeneric instance = new ImageGeneric(filename);
        String result = instance.getFilename();

        assertEquals(result, filename);
    }

    /**
     * Test of getXDim method, of class ImageGeneric.
     */
    @Test
    public void testGetXDim() {
        for (int i = 0; i < filenames.length; i++) {
            int result = 0;
            try {
                result = testGetXDim(filenames[i], dimensions[i][X]);
            } catch (Exception ex) {
                fail("testGetXDim(): " + filenames[i] + " // expected: " + dimensions[i][X] + " / got:" + result);
                //Logger.getLogger(ImageGenericTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    int testGetXDim(String filename, int xDim) throws Exception {
        ImageGeneric instance = new ImageGeneric(filename);
        int result = instance.getXDim();

        assertEquals(xDim, result);

        return result;
    }

    /**
     * Test of getYDim method, of class ImageGeneric.
     */
    @Test
    public void testGetYDim() {
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetYDim(filenames[i], dimensions[i][Y]);
            } catch (Exception ex) {
                fail("testGetYDim(): " + filenames[i]);
                //Logger.getLogger(ImageGenericTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    void testGetYDim(String filename, int yDim) throws Exception {
        ImageGeneric instance = new ImageGeneric(filename);
        int result = instance.getYDim();

        assertEquals(yDim, result);
    }

    /**
     * Test of getZDim method, of class ImageGeneric.
     */
    @Test
    public void testGetZDim() {
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetZDim(filenames[i], dimensions[i][Z]);
            } catch (Exception ex) {
                fail("testGetZDim(): " + filenames[i]);
                //Logger.getLogger(ImageGenericTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    void testGetZDim(String filename, int zDim) throws Exception {
        ImageGeneric instance = new ImageGeneric(filename);
        int result = instance.getZDim();

        assertEquals(zDim, result);
    }

    /**
     * Test of getNDim method, of class ImageGeneric.
     */
    @Test
    public void testGetNDim() {
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetNDim(filenames[i], N[i]);
            } catch (Exception ex) {
                fail("testGetNDim(): " + filenames[i]);
                //Logger.getLogger(ImageGenericTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    void testGetNDim(String filename, long nDim) throws Exception {
        ImageGeneric instance = new ImageGeneric(filename);
        long result = instance.getNDim();

        assertEquals(nDim, result);
    }

    /**
     * Test of getDataType method, of class ImageGeneric.
     */
    @Test
    public void testGetDataType() {
        System.out.println("testGetDataType");
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetDataType(filenames[i], dataType[i]);
            } catch (Exception ex) {
                fail("testGetDataType(): " + filenames[i]);
                //Logger.getLogger(ImageGenericTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public void testGetDataType(String filename, int dataType) throws Exception {
        ImageGeneric instance = new ImageGeneric(filename);
        int result = instance.getDataType();

        assertEquals(dataType, result);
    }

      /**
       * Test of getArrayFloat method, of class ImageGeneric.
       */
      @Test
      public void testGetArrayFloat() throws Exception {
          System.out.println("getArrayFloat");
          int slice = 1;
          ImageGeneric instance = new ImageGeneric(filename);
          instance.read(slice);
          float[] expResult = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
          float[] result = instance.getArrayFloat(slice);
          //System.out.println(Arrays.toString(result));
          //may be there is a better way but I do not know it
          for (int i = 0; i < expResult.length; i++)
        	  assertEquals(expResult[i], result[i],0.001);
      }

    /**
     * Test of getDataType method, of class ImageGeneric.
     */
    @Test
    public void testEqual() {
         try {
             System.out.println("testEqual");
             int slice = 1;
             ImageGeneric instance = new ImageGeneric(filename);
             ImageGeneric instance2 = new ImageGeneric(filename2);
             instance.read(slice);
             instance2.read(slice);
             assertTrue(instance.equal(instance));
         } catch (Exception ex) {
             fail("testGetDataType(): " + filename);
             //Logger.getLogger(ImageGenericTest.class.getName()).log(Level.SEVERE, null, ex);
         }
    }

//    /**
//     * Test of read method, of class ImageGeneric.
//     */
//    @Test
//    public void testRead_int() throws Exception {
//        System.out.println("read");
//        int slice = 0;
//        ImageGeneric instance = new ImageGeneric(filename);
//        instance.read(slice);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
// 
//    /**
//     * Test of read method, of class ImageGeneric.
//     */
//    @Test
//    public void testRead_3args_1() throws Exception {
//        System.out.println("read");
//        int width = 0;
//        int height = 0;
//        int slice = 0;
//        ImageGeneric instance = new ImageGeneric("/home/jvega/Escritorio/imgs_Roberto/DnaB_50000X1.dm3");
//        instance.read(width, height, slice);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
// //
// //    /**
// //     * Test of read method, of class ImageGeneric.
// //     */
// //    @Test
// //    public void testRead_int_long() throws Exception {
// //        System.out.println("read");
// //        int slice = 0;
// //        long image = 0L;
// //        ImageGeneric instance = new ImageGeneric(filename);
// //        instance.read(slice, image);
// //        // TODO review the generated test code and remove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }
// //
// //    /**
// //     * Test of read method, of class ImageGeneric.
// //     */
// //    ///@Test
// //    public void testRead_long() throws Exception {
// //        System.out.println("read");
// //        long image = 0L;
// //        ImageGeneric instance = new ImageGeneric();
// //        instance.read(image);
// //        // TODO review the generated test code and remlove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }
// //
// //    /**
// //     * Test of read method, of class ImageGeneric.
// //     */
// //    ///@Test
// //    public void testRead_3args_2() throws Exception {
// //        System.out.println("read");
// //        int width = 0;
// //        int height = 0;
// //        long image = 0L;
// //        ImageGeneric instance = new ImageGeneric();
// //        instance.read(width, height, image);
// //        // TODO review the generated test code and remove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }
// //
// //    /**
// //     * Test of read method, of class ImageGeneric.
// //     */
// //    @Test
// //    public void testRead_4args() throws Exception {
// //        System.out.println("read");
// //        int width = 0;
// //        int height = 0;
// //        int slice = 0;
// //        long image = 0L;
// //        ImageGeneric instance = new ImageGeneric(filename);
// //        instance.read(width, height, slice, image);
// //        int xDim = instance.getXDim();
// //        int yDim = instance.getYDim();
// //
// //
// //        // TODO review the generated test code and remove the default call to fail.
// //        ///fail("The test case is a prototype.");
// //        assertNotNull(instance);
// //        assert xDim > 0;
// //        assert yDim > 0;
// //
// //    }
// //
// //    /**
// //     * Test of getArrayByte method, of class ImageGeneric.
// //     */
// //    ///@Test
// //    public void testGetArrayByte() throws Exception {
// //        System.out.println("getArrayByte");
// //        int slice = 0;
// //        ImageGeneric instance = new ImageGeneric();
// //        byte[] expResult = null;
// //        byte[] result = instance.getArrayByte(slice);
// //        assertEquals(expResult, result);
// //        // TODO review the generated test code and remove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }
// //
// //    /**
// //     * Test of getArrayShort method, of class ImageGeneric.
// //     */
// //    ///@Test
// //    public void testGetArrayShort() throws Exception {
// //        System.out.println("getArrayShort");
// //        int slice = 0;
// //        ImageGeneric instance = new ImageGeneric();
// //        short[] expResult = null;
// //        short[] result = instance.getArrayShort(slice);
// //        assertEquals(expResult, result);
// //        // TODO review the generated test code and remove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }
// //

// //
// //    /**
// //     * Test of write method, of class ImageGeneric.
// //     */
// //    ///@Test
// //    public void testWrite() throws Exception {
// //        System.out.println("write");
// //        String filename = "";
// //        ImageGeneric instance = new ImageGeneric();
// //        instance.write(filename);
// //        // TODO review the generated test code and remove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }
// //
// //    /**
// //     * Test of setArrayByte method, of class ImageGeneric.
// //     */
// //    ///@Test
// //    public void testSetArrayByte() throws Exception {
// //        System.out.println("setArrayByte");
// //        byte[] data = null;
// //        ImageGeneric instance = new ImageGeneric();
// //        instance.setArrayByte(data);
// //        // TODO review the generated test code and remove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }
// //
// //    /**
// //     * Test of setArrayShort method, of class ImageGeneric.
// //     */
// //    ///@Test
// //    public void testSetArrayShort() throws Exception {
// //        System.out.println("setArrayShort");
// //        short[] data = null;
// //        ImageGeneric instance = new ImageGeneric();
// //        instance.setArrayShort(data);
// //        // TODO review the generated test code and remove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }
// //
// //    /**
// //     * Test of setArrayFloat method, of class ImageGeneric.
// //     */
// //    ///@Test
// //    public void testSetArrayFloat() throws Exception {
// //        System.out.println("setArrayFloat");
// //        float[] data = null;
// //        ImageGeneric instance = new ImageGeneric();
// //        instance.setArrayFloat(data);
// //        // TODO review the generated test code and remove the default call to fail.
// //        fail("The test case is a prototype.");
// //    }

    /**
     * Test of getStatistics method, of class ImageGeneric.
     */
    @Test
    public void testGetStatistics() throws Exception {
        for (int i = 0; i < filenames.length; i++) {
            testGetStatistics(filenames[i], statistics[i]);
        }
    }

    public void testGetStatistics(String filename, double statistics[]) throws Exception {
        try {
            ImageGeneric instance = new ImageGeneric(filename);

            //for (long nimage = ImageGeneric.FIRST_IMAGE; nimage <= instance.getNDim(); nimage++) {
            instance.read(ImageGeneric.FIRST_IMAGE);//nimage);

            double[] result = instance.getStatistics();
            //double[] expResult = new double[]{min, max, avg, std};

            for (int i = 0; i < statistics.length; i++) {
                assertEquals(result[i], statistics[i], 0.0002);
            }
            //}
        } catch (Exception ex) {
            fail("testGetStatistics(String filename): " + filename);
        }
    }

    /**
     * Test of isPSD method, of class ImageGeneric.
     */
    ///@Test
    public void testIsPSD() throws Exception {
        System.out.println("isPSD");
        ImageGeneric instance = new ImageGeneric();
        boolean expResult = false;
        boolean result = instance.isPSD();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isStack method, of class ImageGeneric.
     */
    ///@Test
    public void testIsStack() {
        try {
            System.out.println("isStack");
            ImageGeneric instance = new ImageGeneric();
            boolean expResult = false;
            boolean result = instance.isStack();
            assertEquals(expResult, result);
            // TODO review the generated test code and remove the default call to fail.
            fail("The test case is a prototype.");
        } catch (Exception ex) {
            fail("testIsStack()");
        }
    }

    /**
     * Test of isVolume method, of class ImageGeneric.
     */
    ///@Test
    public void testIsVolume() {
        try {
            System.out.println("isVolume");
            ImageGeneric instance = new ImageGeneric();
            boolean expResult = false;
            boolean result = instance.isVolume();
            assertEquals(expResult, result);
            // TODO review the generated test code and remove the default call to fail.
            fail("The test case is a prototype.");
        } catch (Exception ex) {
            fail("testIsVolume()");
        }
    }

    /**
     * Test of isStackOrVolume method, of class ImageGeneric.
     */
    ///@Test
    public void testIsStackOrVolume() {

        try {
            System.out.println("isStackOrVolume");
            ImageGeneric instance = new ImageGeneric();
            boolean expResult = false;
            boolean result = instance.isStackOrVolume();
            assertEquals(expResult, result);
            // TODO review the generated test code and remove the default call to fail.
            fail("The test case is a prototype.");
        } catch (Exception ex) {
            fail("testIsStackOrVolume()");
        }
    }

    /**
     * Test of isSingleImage method, of class ImageGeneric.
     */
    ///@Test
    public void testIsSingleImage() {
        try {
            System.out.println("isSingleImage");
            ImageGeneric instance = new ImageGeneric();
            boolean expResult = false;
            boolean result = instance.isSingleImage();
            assertEquals(expResult, result);
            // TODO review the generated test code and remove the default call to fail.
            fail("The test case is a prototype.");
        } catch (Exception ex) {
            fail("testIsSingleImage()");
        }
    }
}
