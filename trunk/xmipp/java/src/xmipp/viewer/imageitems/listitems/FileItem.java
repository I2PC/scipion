/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * FileItem.java
 *
 * Created on 18-ene-2010, 11:31:30
 */
package xmipp.viewer.imageitems.listitems;

import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class FileItem {

    protected File file;

    public FileItem() {
    }

    public FileItem(File file) {
        this.file = file;
    }

    public boolean exists() {
        return file.exists();
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

    public String getAbsoluteFileName() {
        return file.getAbsolutePath();
    }

    public boolean isDirectory() {
        return file.isDirectory();
    }

    public String getDirectory() {
        try {
            return file.getParentFile().getCanonicalPath();
        } catch (Exception ex) {
            return null;
        }
    }
}
