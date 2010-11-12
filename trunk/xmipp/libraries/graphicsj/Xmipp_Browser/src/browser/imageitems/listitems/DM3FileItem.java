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
import xmipp.io.DM3_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class DM3FileItem extends AbstractImageItem {

    protected DM3_Reader dm3_reader;

    public DM3FileItem(File file, Cache cache) {
        super(file, cache);
    }

    @Override
    public boolean loadImageData() {
        try {
            if (dm3_reader == null) {
                dm3_reader = new DM3_Reader();
            }

            dm3_reader.parseDM3(getDirectory(), getFileName());

            FileInfo fi = dm3_reader.getDM3FileInfo(getDirectory(), getFileName());

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

            ip = dm3_reader.loadThumbnail(getDirectory(), getFileName(), w, h, slice);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return ip;
    }
}
