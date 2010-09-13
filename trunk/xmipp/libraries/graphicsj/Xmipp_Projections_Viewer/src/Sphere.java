
import constants.LABELS;
import ij.IJ;
import ij.ImagePlus;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class Sphere {

    public Sphere() {
    }

    public ImagePlus getProjection(MultidimArrayd xmippVolume, double rot, double tilt, int w, int h) {
        ImageDouble id = new ImageDouble(w, h);

        Projection p = new Projection();
        p.reset(h, w);

        IJ.showStatus(LABELS.MESSAGE_RETRIEVING_PROJECTION);
        XmippData.project_Volume(xmippVolume, p, h, w, rot, tilt, 0);

        ImagePlus projection = ImageConverter.xmipp2imagej(p, w, h);   // To imageJ.

        p.delete();

        return projection;
    }
}
