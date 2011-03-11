/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import browser.ICONS_MANAGER;
import browser.LABELS;
import browser.imageitems.ImageDimension;
import browser.table.ImageOperations;
import ij.ImagePlus;
import java.io.File;
import java.text.DecimalFormat;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public abstract class AbstractImageItem extends FileItem {

    //public int width, height;
    protected ImageDimension dimension;
    protected Cache cache;
    protected double min, max, mean, stdDev;

    public AbstractImageItem(File file, Cache cache) {
        super(file);

        this.cache = cache;

        loadImageData();
    }

    public abstract ImagePlus getImagePlus();

    public ImagePlus getPreview(int w, int h) {
        ImagePlus preview;

        if (dimension.getWidth() > 0 && dimension.getHeight() > 0) {
            // Tries to load from cache.
            preview = (ImagePlus) cache.get(getKey());

            // If not in cache.
            if (preview == null) {
                System.out.println("Reading from disk: " + getKey());
                preview = loadPreview(w, h);

                if (preview != null) {
                    cache.put(getKey(), preview);

                    // Stores image info.
                    min = preview.getProcessor().getMin();
                    max = preview.getProcessor().getMax();
                    mean = ImageOperations.mean(preview);
                    stdDev = ImageOperations.std_dev(preview);
                }
            }
        } else {    // Null preview.
            preview = new ImagePlus("", ICONS_MANAGER.MISSING_ITEM.getImage());
        }

        return preview;
    }

    public String getKey() {
        return file.getAbsolutePath() + "_" + dimension.getWidth() + "_" + dimension.getHeight();
    }

    public int getWidth() {
        return dimension.getWidth();
    }

    public int getHeight() {
        return dimension.getHeight();
    }

    public int getDepth() {
        return dimension.getDepth();
    }

    public long getNImages() {
        return dimension.getNimages();
    }

    public boolean isSingleImage() {
        return !isStack() & !isVolume();
    }

    public boolean isStack() {
        return dimension.getNimages() > 1;
    }

    public boolean isVolume() {
        return dimension.getDepth() > 1;
    }

    //protected abstract void loadImageData();
    protected void loadImageData() {
        try {
            ImageDouble image = new ImageDouble();
            image.readHeader(file.getAbsolutePath());

            dimension = new ImageDimension(image);
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }

    //public abstract String getImageInfo();
    public String getImageInfo() {
        loadImageData();

        DecimalFormat myFormatter = new DecimalFormat("#.###");
        String strMin = myFormatter.format(min);
        String strMax = myFormatter.format(max);
        String strMean = myFormatter.format(mean);
        String strStdDev = myFormatter.format(stdDev);

        return "<html>"
                + LABELS.LABEL_WIDTH + dimension.getWidth() + "<br>"
                + LABELS.LABEL_HEIGHT + dimension.getHeight() + "<br>"
                + (isVolume() ? LABELS.LABEL_DEPTH + dimension.getDepth() + "<br>" : "")
                + (isStack() ? LABELS.LABEL_NIMAGES + dimension.getNimages() : "")
                + "<br>" + "<br>"
                + "Min=" + strMin + "<br>"
                + "Max=" + strMax + "<br>"
                + "Mean=" + strMean + "<br>"
                + "Std. dev.=" + strStdDev + "<br>"
                + "</p>"
                + "</html>";
    }

    protected abstract ImagePlus loadPreview(int w, int h);
}
