/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package table;

import ij.IJ;
import ij.ImagePlus;
import java.awt.Image;
import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class ScoreItem {

    public String fileName;
    public double score = 0.0;
    public Image image;
    public String label;
    public int w, h;

    public ScoreItem(String fileName, double score) {
        this.fileName = fileName;
        this.score = score;

        label = (new File(fileName)).getName();
    }

    public Image getImage() {
        if (image == null) {
            ImagePlus ip = IJ.openImage(fileName);
            image = ip.getImage();

            w = ip.getWidth();
            h = ip.getHeight();
        }

        return image;
    }

    public String getLabel() {
        return label;
    }

    public String getTooltip() {
        return "score = " + score;
    }
}
