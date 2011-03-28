/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package table;

import ij.ImagePlus;
import java.io.File;
import sphere.ImageConverter;
import xmipp.ImageDouble;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class ScoreItem implements Comparable<ScoreItem> {

    protected MetaData md;
    protected long id;
    protected String fileName;
    protected double score = 0.0;
    protected boolean good;
    protected ImagePlus ip;
    protected String label;

    public ScoreItem(MetaData md, long id) {
        this.md = md;
        this.id = id;
        fileName = md.getValueString(MDLabel.MDL_IMAGE, id);
        score = md.getValueDouble(MDLabel.MDL_ZSCORE, id);
        good = md.getValueInt(MDLabel.MDL_ENABLED, id) != 0;

        label = (new File(fileName)).getName();
    }

    public ImagePlus getImagePlus() {
        try {
            if (ip == null) {
                ImageDouble image = new ImageDouble(md, id);

                ip = ImageConverter.convertToImagej(image, fileName);
            }
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        } finally {
            return ip;
        }
    }

    public int getWidth() {
        return getImagePlus().getWidth();
    }

    public int getHeight() {
        return getImagePlus().getHeight();
    }

    public String getLabel() {
        return label;
    }

    public String getTooltip() {
        return "score = " + score;
    }

    public int compareTo(ScoreItem item) {
        if (good != item.good) {
            return good ? -1 : 1;
        } else {
            if (score != item.score) {
                return score < item.score ? -1 : 1;
            } else {
                return fileName.compareTo(item.fileName);
            }
        }
    }
}
