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
                /**
                int dist = 0;//, way = 1;

                while(!carga && dist>size/2[o tambien: i-dist>0])
                if(!load(i+dist)){
                if (!load(i-dist)){
                dist++;
                }
                }
                 */
                String files[] = Sel_Reader.loadFileNames(file.getParent(), file.getName());

                nslices = files.length; // Store # slices.

                // Loads preview (skipping while reference is not satisfied).
                int index = -1;//nslices / 2;
                do {
                    index++;
                    // Sets file to force super class to load it.
                    file = new File(files[index]);
                } while (!file.exists() && index < files.length);

                if (file.exists()) {
                    System.out.println(" >> " + file.getName());
                    super.loadImageData();
                } else {
                    // Not valid reference, load a "X" preview
                    width = 0;
                    height = 0;
                }
            } catch (Exception e) {
                IJ.error("File " + getFileName() + " can't be read.");
            }
        }
    }
}
