/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import ij.ImagePlus;
import ij.io.FileInfo;
import java.io.File;
import xmipp.io.OpenSPE_;

/**
 *
 * @author Juanjo Vega
 */
public class SPEFileItem extends AbstractImageItem {

    protected OpenSPE_ spe_reader;

    public SPEFileItem(File file, Cache cache) {
        super(file, cache);
    }

    @Override
    public boolean loadImageData() {
        try {
            if (spe_reader == null) {
                spe_reader = new OpenSPE_();
            }

            spe_reader.parseSPE(getDirectory(), getFileName());

            FileInfo fi = OpenSPE_.getSPEFileInfo(getDirectory(), getFileName());

            // Stores image info.
            width = fi.width;
            height = fi.height;
            nslices = fi.nImages;
        } catch (Exception e) {
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

            ip = spe_reader.loadThumbnail(getDirectory(), getFileName(), w, h, slice);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return ip;
    }
}
