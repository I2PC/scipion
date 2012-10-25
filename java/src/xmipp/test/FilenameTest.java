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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import org.junit.*;
import xmipp.jni.Filename;

public class FilenameTest {

    // Xmipp dir
    static String TESTS_PATH;
    
    String filename  = "singleImage.spi";
    String filename2 = "singleImage2.spi";


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
    public void testGetXmippPath() {
            try {
                String xmipp_home = Filename.getXmippPath();
                String xmipp_home2 = System.getenv("XMIPP_HOME");
                assertEquals(xmipp_home, xmipp_home2);
            } catch (Exception ex) {
                fail("testGetFilename");
            }
    }//function testGetXmippPath

}//class FilenameTest
