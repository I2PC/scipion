/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileInfo;
import java.io.File;
import xmipp.io.TIA_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class SERFileItem extends AbstractImageItem {

    protected TIA_Reader tia_reader;

    public SERFileItem(File file, Cache cache) {
        super(file, cache);
    }

    @Override
    public boolean loadImageData() {
        try {
            if (tia_reader == null) {
                tia_reader = new TIA_Reader();
            }

            String path = getFile().getCanonicalPath();

            tia_reader.parseTIA(path);

            FileInfo fi = tia_reader.getTIAFileInfo(path);

            // Stores image info.
            width = fi.width;
            height = fi.height;
            nslices = fi.nImages;
        } catch (Exception e) {
            IJ.write(getDirectory());
            IJ.write(getFileName());
            e.printStackTrace();
            return false;
        }

        return true;
    }

    @Override
    public ImagePlus loadPreview(int w, int h, int slice) {
        ImagePlus ip = null;

        try {
            loadImageData();

            ip = tia_reader.loadThumbnail(getDirectory(), getFileName(), w, h, slice);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return ip;
    }
}
