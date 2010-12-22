/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import ij.ImagePlus;
import ij.io.FileInfo;
import java.io.File;
import xmipp.io.Spider_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class SpiderItem extends AbstractImageItem {

    protected static Spider_Reader spider_reader;

    public SpiderItem(File file, Cache cache) {
        super(file, cache);
    }

    public SpiderItem(File file, Cache cache, int slice) {
        super(file, cache, slice);
    }

    @Override
    public String getKey() {
        String key = super.getKey();
        if (slice >= 0) {
            key += "-" + slice;
        }

        return key;
    }

    public boolean loadImageData() {
        try {
            if (spider_reader == null) {
                spider_reader = new Spider_Reader();
            }

            spider_reader.parseSPIDER(getDirectory(), getFileName());

            FileInfo fi = spider_reader.getFileInfo();

            // Stores image info.
            width = fi.width;
            height = fi.height;
            nslices = fi.nImages;
        } catch (Exception e) {
            return false;
        }

        return true;
    }

    @Override
    public ImagePlus loadPreview(int w, int h, int slice) {
        ImagePlus ip = null;

        try {
            loadImageData();

            ip = spider_reader.loadThumbnail(getDirectory(), getFileName(), w, h, slice);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return ip;
    }

    public ImagePlus[] loadPreviews(int w, int h) {
        ImagePlus previews[] = null;

        try {
            loadImageData();

            previews = spider_reader.loadVolumeSlices(file.getParent(), getFileName(), w, h);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return previews;
    }
}
