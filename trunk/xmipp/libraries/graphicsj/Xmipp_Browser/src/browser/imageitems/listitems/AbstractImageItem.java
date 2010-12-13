/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import browser.ICONS_MANAGER;
import browser.LABELS;
import ij.ImagePlus;
import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public abstract class AbstractImageItem extends FileItem {

    public final static int MID_SLICE = -1;
    public int width, height, nslices = 1;
    public int slice;
    protected Cache cache;

    public AbstractImageItem(File file, Cache cache) {
        this(file, cache, MID_SLICE);
    }

    public AbstractImageItem(File file, Cache cache, int slice) {
        super(file);

        this.cache = cache;
        this.slice = slice;

        loadImageData();    // Loads size at start up.
    }

    public abstract boolean loadImageData();

    public String getKey() {
        return file.getAbsolutePath();
    }

    @Override
    public ImagePlus getPreview(int w, int h) {
        ImagePlus preview;

        if (width > 0 && height > 0) {
            // Tries to load from cache.
            preview = (ImagePlus) cache.get(getKey());

            // If not in cache.
            if (preview == null) {
                preview = loadPreview(w, h, slice);

                if (preview != null) {
                    cache.put(getKey(), preview);
                }
            }
        } else {    // Null preview.
            preview = new ImagePlus("", ICONS_MANAGER.MISSING_ITEM.getImage());
        }

        return preview;
    }

    public abstract ImagePlus loadPreview(int w, int h, int slice);

    @Override
    public String getImageInfo() {
        loadImageData();

        return "<html>" + LABELS.LABEL_WIDTH + width
                + "<br>" + LABELS.LABEL_HEIGHT + height
                + (nslices > 1 ? "<br>" + LABELS.LABEL_NSLICES + nslices : "") + "</html>";
    }
}
