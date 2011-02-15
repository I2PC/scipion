/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import browser.ICONS_MANAGER;
import browser.imageitems.ImageDimension;
import ij.ImagePlus;
import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public abstract class AbstractImageItem extends FileItem {

    //public int width, height;
    public ImageDimension dimension;
    protected Cache cache;

    public AbstractImageItem(File file, Cache cache) {
        super(file);

        this.cache = cache;

        loadImageData();
    }

    public abstract ImagePlus getImagePlus();

    public ImagePlus getPreview(int w, int h) {
        ImagePlus preview;

        if (dimension.width > 0 && dimension.height > 0) {
            // Tries to load from cache.
            preview = (ImagePlus) cache.get(getKey());

            // If not in cache.
            if (preview == null) {
                preview = loadPreview(w, h);

                if (preview != null) {
                    cache.put(getKey(), preview);
                }
            }
        } else {    // Null preview.
            preview = new ImagePlus("", ICONS_MANAGER.MISSING_ITEM.getImage());
        }

        return preview;
    }

    public String getKey() {
        return file.getAbsolutePath() + "-" + dimension.width + "-" + dimension.height;
    }

    public boolean isSingleImage() {
        return !isStack() & !isVolume();
    }

    public boolean isStack() {
        return dimension.nimages > 1;
    }

    public boolean isVolume() {
        return dimension.depth > 1;
    }

    protected abstract void loadImageData();

    public abstract String getImageInfo();

    protected abstract ImagePlus loadPreview(int w, int h);
}
