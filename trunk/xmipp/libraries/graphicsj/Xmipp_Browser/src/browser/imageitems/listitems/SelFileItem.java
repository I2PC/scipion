/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import ij.ImagePlus;
import java.io.File;
import xmipp.io.Sel_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class SelFileItem extends SpiderItem {

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
    public ImagePlus loadPreview(int w, int h, int slice) {
        ImagePlus ip = null;

        if (previewFile == null) {
            loadPreviewFile();
        }

        File originalFile = file;
        file = previewFile;

        if (file != null) {
            ip = super.loadPreview(w, h, slice);
        }

        file = originalFile;

        return ip;
    }

    @Override
    public boolean loadImageData() {
        boolean result = false;
        int slices = nslices;


        if (previewFile == null) {
            slices = loadPreviewFile();
        }

        File originalFile = file;
        file = previewFile;

        if (file != null) {
            result = super.loadImageData();
            nslices = slices;
        }

        file = originalFile;

        return result;
    }
}
