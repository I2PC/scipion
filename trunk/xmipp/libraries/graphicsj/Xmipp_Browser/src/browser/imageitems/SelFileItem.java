/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import browser.imageitems.listitems.FileImageItem;
import browser.Cache;
import ij.ImagePlus;
import java.io.File;
import xmipp.Sel_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class SelFileItem extends FileImageItem {

    protected File previewFile;

    public SelFileItem(File file, Cache cache) {
        super(file, cache);

        nslices = loadPreviewFile();

        this.cache = cache;
    }

    protected int loadPreviewFile() {
        String files[] = Sel_Reader.loadFileNames(file.getParent(), file.getName());

        // Loads preview (skipping while reference is not satisfied).
        int index = 0;
        do {
            // Sets file to force super class to load it.
            previewFile = new File(files[index]);

            index++;
        } while (!previewFile.exists() && index < files.length);

        return files.length;
    }

    @Override
    protected ImagePlus loadXmippPreview(int w, int h, int slice) {
        ImagePlus ip = null;

        if (previewFile == null) {
            nslices = loadPreviewFile();
        }

        File originalFile = file;
        file = previewFile;

        if (file != null) {
            ip = super.loadXmippPreview(w, h, slice);
        }

        file = originalFile;

        return ip;
    }

    @Override
    protected boolean loadXmippImageData() {
        boolean result = false;

        if (previewFile == null) {
            nslices = loadPreviewFile();
        }

        File originalFile = file;
        file = previewFile;

        if (file != null) {
            result = super.loadXmippImageData();
        }

        file = originalFile;

        return result;
    }
}
