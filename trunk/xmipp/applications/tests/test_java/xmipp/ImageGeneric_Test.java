/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp;

import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

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
     * Test of getFilename method, of class ImageGeneric_.
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
        ImageGeneric_ instance = new ImageGeneric_(filename);
        String result = instance.getFilename();

        assertEquals(result, filename);
    }

    /**
     * Test of getXDim method, of class ImageGeneric_.
     */
    @Test
    public void testGetXDim() {
        for (int i = 0; i < filenames.length; i++) {
            int result = 0;
            try {
                result = testGetXDim(filenames[i], dimensions[i][X]);
            } catch (Exception ex) {
                fail("testGetXDim(): " + filenames[i] + " // expected: " + dimensions[i][X] + " / got:" + result);
                //Logger.getLogger(ImageGeneric_Test.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    int testGetXDim(String filename, int xDim) throws Exception {
        ImageGeneric_ instance = new ImageGeneric_(filename);
        int result = instance.getXDim();

        assertEquals(xDim, result);

        return result;
    }

    /**
     * Test of getYDim method, of class ImageGeneric_.
     */
    @Test
    public void testGetYDim() {
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetYDim(filenames[i], dimensions[i][Y]);
            } catch (Exception ex) {
                fail("testGetYDim(): " + filenames[i]);
                //Logger.getLogger(ImageGeneric_Test.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    void testGetYDim(String filename, int yDim) throws Exception {
        ImageGeneric_ instance = new ImageGeneric_(filename);
        int result = instance.getYDim();

        assertEquals(yDim, result);
    }

    /**
     * Test of getZDim method, of class ImageGeneric_.
     */
    @Test
    public void testGetZDim() {
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetZDim(filenames[i], dimensions[i][Z]);
            } catch (Exception ex) {
                fail("testGetZDim(): " + filenames[i]);
                //Logger.getLogger(ImageGeneric_Test.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    void testGetZDim(String filename, int zDim) throws Exception {
        ImageGeneric_ instance = new ImageGeneric_(filename);
        int result = instance.getZDim();

        assertEquals(zDim, result);
    }

    /**
     * Test of getNDim method, of class ImageGeneric_.
     */
    @Test
    public void testGetNDim() {
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetNDim(filenames[i], N[i]);
            } catch (Exception ex) {
                fail("testGetNDim(): " + filenames[i]);
                //Logger.getLogger(ImageGeneric_Test.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    void testGetNDim(String filename, long nDim) throws Exception {
        ImageGeneric_ instance = new ImageGeneric_(filename);
        long result = instance.getNDim();

        assertEquals(nDim, result);
    }

    /**
     * Test of getDataType method, of class ImageGeneric_.
     */
    @Test
    public void testGetDataType() {
        for (int i = 0; i < filenames.length; i++) {
            try {
                testGetDataType(filenames[i], dataType[i]);
            } catch (Exception ex) {
                fail("testGetDataType(): " + filenames[i]);
                //Logger.getLogger(ImageGeneric_Test.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public void testGetDataType(String filename, int dataType) throws Exception {
        ImageGeneric_ instance = new ImageGeneric_(filename);
        int result = instance.getDataType();

        assertEquals(dataType, result);
    }
//
//    /**
//     * Test of read method, of class ImageGeneric_.
//     */
//    @Test
//    public void testRead_int() throws Exception {
//        System.out.println("read");
//        int slice = 0;
//        ImageGeneric_ instance = new ImageGeneric_(filename);
//        instance.read(slice);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of read method, of class ImageGeneric_.
//     */
//    @Test
//    public void testRead_3args_1() throws Exception {
//        System.out.println("read");
//        int width = 0;
//        int height = 0;
//        int slice = 0;
//        ImageGeneric_ instance = new ImageGeneric_("/home/jvega/Escritorio/imgs_Roberto/DnaB_50000X1.dm3");
//        instance.read(width, height, slice);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of read method, of class ImageGeneric_.
//     */
//    @Test
//    public void testRead_int_long() throws Exception {
//        System.out.println("read");
//        int slice = 0;
//        long image = 0L;
//        ImageGeneric_ instance = new ImageGeneric_(filename);
//        instance.read(slice, image);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of read method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testRead_long() throws Exception {
//        System.out.println("read");
//        long image = 0L;
//        ImageGeneric_ instance = new ImageGeneric_();
//        instance.read(image);
//        // TODO review the generated test code and remlove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of read method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testRead_3args_2() throws Exception {
//        System.out.println("read");
//        int width = 0;
//        int height = 0;
//        long image = 0L;
//        ImageGeneric_ instance = new ImageGeneric_();
//        instance.read(width, height, image);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of read method, of class ImageGeneric_.
//     */
//    @Test
//    public void testRead_4args() throws Exception {
//        System.out.println("read");
//        int width = 0;
//        int height = 0;
//        int slice = 0;
//        long image = 0L;
//        ImageGeneric_ instance = new ImageGeneric_(filename);
//        instance.read(width, height, slice, image);
//        int xDim = instance.getXDim();
//        int yDim = instance.getYDim();
//
//
//        // TODO review the generated test code and remove the default call to fail.
//        ///fail("The test case is a prototype.");
//        assertNotNull(instance);
//        assert xDim > 0;
//        assert yDim > 0;
//
//    }
//
//    /**
//     * Test of getArrayByte method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testGetArrayByte() throws Exception {
//        System.out.println("getArrayByte");
//        int slice = 0;
//        ImageGeneric_ instance = new ImageGeneric_();
//        byte[] expResult = null;
//        byte[] result = instance.getArrayByte(slice);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getArrayShort method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testGetArrayShort() throws Exception {
//        System.out.println("getArrayShort");
//        int slice = 0;
//        ImageGeneric_ instance = new ImageGeneric_();
//        short[] expResult = null;
//        short[] result = instance.getArrayShort(slice);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getArrayFloat method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testGetArrayFloat() throws Exception {
//        System.out.println("getArrayFloat");
//        int slice = 0;
//        ImageGeneric_ instance = new ImageGeneric_();
//        float[] expResult = null;
//        float[] result = instance.getArrayFloat(slice);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of write method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testWrite() throws Exception {
//        System.out.println("write");
//        String filename = "";
//        ImageGeneric_ instance = new ImageGeneric_();
//        instance.write(filename);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setArrayByte method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testSetArrayByte() throws Exception {
//        System.out.println("setArrayByte");
//        byte[] data = null;
//        ImageGeneric_ instance = new ImageGeneric_();
//        instance.setArrayByte(data);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setArrayShort method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testSetArrayShort() throws Exception {
//        System.out.println("setArrayShort");
//        short[] data = null;
//        ImageGeneric_ instance = new ImageGeneric_();
//        instance.setArrayShort(data);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of setArrayFloat method, of class ImageGeneric_.
//     */
//    ///@Test
//    public void testSetArrayFloat() throws Exception {
//        System.out.println("setArrayFloat");
//        float[] data = null;
//        ImageGeneric_ instance = new ImageGeneric_();
//        instance.setArrayFloat(data);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }

    /**
     * Test of getStatistics method, of class ImageGeneric_.
     */
    @Test
    public void testGetStatistics() throws Exception {
        for (int i = 0; i < filenames.length; i++) {
            testGetStatistics(filenames[i], statistics[i]);
        }
    }

    public void testGetStatistics(String filename, double statistics[]) throws Exception {
        try {
            ImageGeneric_ instance = new ImageGeneric_(filename);

            //for (long nimage = ImageGeneric_.FIRST_IMAGE; nimage <= instance.getNDim(); nimage++) {
            instance.read(ImageGeneric_.FIRST_IMAGE);//nimage);

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
     * Test of isPSD method, of class ImageGeneric_.
     */
    ///@Test
    public void testIsPSD() {
        System.out.println("isPSD");
        ImageGeneric_ instance = new ImageGeneric_();
        boolean expResult = false;
        boolean result = instance.isPSD();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isStack method, of class ImageGeneric_.
     */
    ///@Test
    public void testIsStack() {
        System.out.println("isStack");
        ImageGeneric_ instance = new ImageGeneric_();
        boolean expResult = false;
        boolean result = instance.isStack();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isVolume method, of class ImageGeneric_.
     */
    ///@Test
    public void testIsVolume() {
        System.out.println("isVolume");
        ImageGeneric_ instance = new ImageGeneric_();
        boolean expResult = false;
        boolean result = instance.isVolume();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isStackOrVolume method, of class ImageGeneric_.
     */
    ///@Test
    public void testIsStackOrVolume() {
        System.out.println("isStackOrVolume");
        ImageGeneric_ instance = new ImageGeneric_();
        boolean expResult = false;
        boolean result = instance.isStackOrVolume();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isSingleImage method, of class ImageGeneric_.
     */
    ///@Test
    public void testIsSingleImage() {
        System.out.println("isSingleImage");
        ImageGeneric_ instance = new ImageGeneric_();
        boolean expResult = false;
        boolean result = instance.isSingleImage();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}
