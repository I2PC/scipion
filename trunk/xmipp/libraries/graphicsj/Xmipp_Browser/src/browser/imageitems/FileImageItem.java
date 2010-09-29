/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import browser.Cache;
import browser.LABELS;
import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class FileImageItem extends ImageItem {

    public FileImageItem(File file, Cache cache) {
        super(file, cache);

        this.cache = cache;
    }

    @Override
    public String getFileInfo() {
        loadImageData();

        return "<html>" + LABELS.LABEL_WIDTH + width
                + "<br>" + LABELS.LABEL_HEIGHT + height
                + (nslices > 1 ? "<br>" + LABELS.LABEL_NSLICES + nslices : "") + "</html>";
    }

    @Override
    protected String getKey() {
        return getFile().getAbsolutePath();
    }
}
