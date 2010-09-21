/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import java.io.File;

/**
 *
 * @author Juanjo Vega
 */
public class FileItemFolder extends FileItem {

    protected int items;

    public FileItemFolder(File file) {
        super(file);

        // Not accesible folders will return null listFiles()
        items = file.listFiles() != null ? file.listFiles().length : -1;
    }

    @Override
    public String getDescription() {
        return items >= 0 ? items + " items" : "-";
    }
}
