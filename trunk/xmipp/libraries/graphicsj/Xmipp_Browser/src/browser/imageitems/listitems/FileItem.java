/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * FileItem.java
 *
 * Created on 18-ene-2010, 11:31:30
 */
package browser.imageitems.listitems;

import browser.files.FileBrowser;
import ij.ImagePlus;
import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class FileItem {

    protected File file;

    public FileItem(File file) {
        this.file = file;
    }

    @Override
    public String toString() {
        return file.toString();
    }

    public String getLabel() {
        return file.getName();
    }

    public File getFile() {
        return file;
    }

    public String getFileName() {
        return file.getName();
    }

    public String getDirectory() {
        try {
            return file.getParentFile().getCanonicalPath();
        } catch (Exception ex) {
            return null;
        }
    }

    public String getDescription() {
        return FileBrowser.getFileSizeString(getFile().length());
    }

    // By default there is no preview. It will be added later just for images.
    public ImagePlus getPreview(int w, int h) {
        return null;
    }

    public String getImageInfo() {
        return null;
    }
}
