/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import ij.ImagePlus;

/**
 *
 * @author Juanjo Vega
 */
public class MetadataImageItem extends AbstractImageItem {

    String absoluteFilename;
    String originalFilename;

    public MetadataImageItem(String absoluteFilename, String originalFilename, Cache cache) {
        super();

        this.absoluteFilename = absoluteFilename;
        this.originalFilename = originalFilename;
    }

    public String getValue() {
        return absoluteFilename;
    }

    @Override
    public String toString() {
        return originalFilename;
    }

    @Override
    public String getAbsoluteFileName() {
        return absoluteFilename;
    }

    @Override
    public String getFileName() {
        return absoluteFilename;
    }

    @Override
    public ImagePlus getImagePlus() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected ImagePlus loadPreview(int w, int h) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
