
import ij.IJ;
import ij.ImageJ;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Main {

    static {
        System.setProperty("plugins.dir", "/home/juanjo/Desktop/ImageJ/plugins");
    }

    public static void main(String args[]) {
        new ImageJ();

        String filename = "class_000002@/data2/MicrographPreprocessing/2D/CL2D/run_001/results_level_00_classes.xmd";

        try {
            MetaDataReader xr = new MetaDataReader();
            xr.run(filename);
            xr.show();

            String output = "/home/juanjo/Desktop/stackfull.stk";
            ImageWriter iw = new ImageWriter();
            iw.write(xr, output);

            // Reload.
            ImageReader ir = new ImageReader();
            ir.run(output);
            ir.show();
        } catch (Exception ex) {
            ex.printStackTrace();
            IJ.error(ex.getMessage());
        }
    }
//
//    public static void main(String args[]) {
//        new ImageJ();
//
//        String filename = "/home/juanjo/Desktop/imgs_Roberto/kk.vol";
//
//        try {
//            ImagePlus imp = IJ.openImage(filename);
//            imp.show();
//
//            String output = "/home/juanjo/Desktop/full.vol";
//            IJ.run(imp, "Xmipp writer", "save=" + output);
//
//            ImagePlus imp2 = IJ.openImage(output);
//            imp2.show();
//        } catch (Exception ex) {
//            IJ.error(ex.getMessage());
//        }
//    }
}
