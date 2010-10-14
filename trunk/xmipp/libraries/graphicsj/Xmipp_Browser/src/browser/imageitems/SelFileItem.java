/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems;

import browser.Cache;
import ij.IJ;
import java.io.File;
import xmipp.Sel_Reader;

/**
 *
 * @author Juanjo Vega
 */
public class SelFileItem extends FileImageItem {

    protected File originalFile;

    public SelFileItem(File file, Cache cache) {
        super(file, cache);

        originalFile = file;
        this.cache = cache;
    }

    @Override
    public File getFile() {
        return originalFile;
    }

    @Override
    public String getDirectory() {
        try {
            return originalFile.getParentFile().getCanonicalPath();
        } catch (Exception ex) {
            return null;
        }
    }

    @Override
    public String getLabel() {
        return originalFile.getName();
    }

    @Override
    protected void loadImageData() {
        if (file.getName().endsWith(".sel")) {  // Only for the first time.
            try {
                String files[] = Sel_Reader.loadFileNames(file.getParent(), file.getName());

                // Loads preview (skipping while reference is not satisfied).
                int index = 0;
                do {
                    // Sets file to force super class to load it.
                    file = new File(files[index]);

                    index++;
                } while (!file.exists() && index < files.length);

                if (file.exists()) {
                    super.loadImageData();
                }

                nslices = files.length; // Store # slices after, to overwrite super class result.
            } catch (Exception e) {
                IJ.error("SelFileItem:: File " + getFileName() + " can't be read.");
            }
        }
    }
}
