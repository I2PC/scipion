/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table.coss;

import browser.Cache;
import browser.imageitems.ImageConverter;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import java.io.File;
import xmipp.FileName;
import xmipp.ImageDouble;
import xmipp.MultidimArrayd;
import xmipp.XmippData;
import xmipp.io.ij.FileOpener;

/**
 *
 * @author Juanjo Vega
 */
public class TableImageItemCOSS {

    public static final int CELL_WIDTH = 128;
    public static final int CELL_HEIGHT = 128;
    protected boolean enabled = true;
    protected boolean selected = false;
    protected Cache cache;
    protected String parentDir, fileName, absolutePath;
    protected FileName fn;
    protected ImageDouble image;
    protected TableModelCOSS imagesTableModel;
    protected int slice = FileOpener.MID_SLICE;

    public TableImageItemCOSS(String parentDir, String fileName, Cache cache, TableModelCOSS imagesTableModel) {
        this(parentDir, fileName, cache, FileOpener.MID_SLICE, imagesTableModel);
    }

    public TableImageItemCOSS(String parentDir, String fileName, Cache cache, int slice, TableModelCOSS imagesTableModel) {
        this.parentDir = parentDir;
        this.fileName = fileName;
        absolutePath = parentDir + File.separator + fileName;
        this.fn = new FileName(absolutePath);
        this.cache = cache;
        this.slice = slice;
        this.imagesTableModel = imagesTableModel;

        image = XmippData.readImageHeader(fileName);
    }

    @Override
    public String toString() {
        return getLabel();
    }

    public String getKey() {
        return absolutePath;//item.getKey() + (slice != FileOpener.MID_SLICE ? "-" + slice : "");
    }

    public String getLabel() {
        return fileName;//item.getLabel() + (slice != FileOpener.MID_SLICE ? "[" + slice + "]" : "");
    }

    public String getFileName() {
        return fileName.toString();
    }

    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    public void setSelected(boolean selected) {
        this.selected = selected;
    }

    public boolean isEnabled() {
        return enabled;
    }

    public boolean isSelected() {
        return selected;
    }

    protected void removeFromCache() {
        cache.remove(getKey());
    }

    public int getWidth() {
        return image.getData().getXdim();
    }

    public int getHeight() {
        return image.getData().getYdim();
    }

    public ImagePlus getFullZiseImagePlus() {
        ImagePlus ip = readImage(absolutePath);
        ip.setTitle(fileName);

        return ip;
    }

    public ImagePlus getImagePlus() {
        ImagePlus img = (ImagePlus) cache.get(getKey());

        // If not in cache...
        if (img == null) {
            // ...loads it from disk.
            img = readImage(absolutePath);

            // Resize...
            ImageProcessor ipr = img.getProcessor();
            img = new ImagePlus(img.getTitle(), ipr.resize(CELL_WIDTH, CELL_HEIGHT));

            cache.put(getKey(), img);   // Stores it
        }

        return img;
    }

    public static ImagePlus readImage(String path) {
        ImageDouble image = XmippData.readFullImage(path);
        FileName fn = new FileName(path);

        // Xmipp(MultidimArray) -> ImageJ
        ImagePlus ip = null;

        if (isPSD(fn)) {
            // ImageJ -> Xmipp(MultidimArray)
            MultidimArrayd mda_in = image.getData();
            MultidimArrayd mda_out = new MultidimArrayd();

            XmippData.xmipp2PSD(mda_in, mda_out);

            // Xmipp(MultidimArray) -> ImageJ
            ip = ImageConverter.xmipp2imagej(mda_out);
        } else {
            ip = ImageConverter.xmipp2Imagej(image, path);
        }

        return ip;
    }

    private static boolean isPSD(FileName fn) {
        if (fn.getExtension().endsWith("psd")) {
            return true;
        }

        return false;
    }
}
