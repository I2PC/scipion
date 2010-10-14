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
import browser.ICONS_MANAGER;
import ij.ImagePlus;
import java.io.File;
import java.util.Iterator;
import javax.imageio.ImageIO;
import javax.imageio.stream.ImageInputStream;
import javax.swing.ImageIcon;
import xmipp.ij.io.Opener;

/**
 *
 * @author Juanjo Vega
 */
public class FileItem {

    protected File file;
    protected ImageIcon icon;

    public FileItem(File file) {
        this.file = file;
        icon = getIconForFile(file);
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

    public ImageIcon getIcon() {
        return icon;
    }

    // By default there is no preview. It will be added later just for images.
    public ImagePlus getPreview(int w, int h) {
        return null;
    }

    public String getFileInfo() {
        return null;
    }

    protected ImageIcon getIconForFile(File file) {
        if (file.isDirectory()) {
            return ICONS_MANAGER.DIRECTORY_FILE_TYPE;
        }

        // Tries to get a reader for the given file...
        try {
            final ImageInputStream iis = ImageIO.createImageInputStream(file);
            final Iterator readers = ImageIO.getImageReaders(iis);

            if (readers.hasNext()) {
                return ICONS_MANAGER.IMAGEJ_IMAGE_FILE_TYPE;   // ...so it is a image file.
            } else {    // No image so...
                int type = (new Opener()).getFileType(file.getCanonicalPath());

                if (file.getName().toLowerCase().endsWith(".sel")) {
                    return ICONS_MANAGER.XMIPP_SELFILE_TYPE;
                }

                if (type != Opener.UNKNOWN) {   // ...let's see if it's a supported file type...
                    return ICONS_MANAGER.IMAGEJ_FILE_TYPE;
                }
            }
        } catch (Exception ex) {
        }

        return ICONS_MANAGER.UNKNOWN_FILE_TYPE;
    }
}
