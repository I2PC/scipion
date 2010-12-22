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
public class FolderFileItem extends FileItem {

    protected int items;

    public FolderFileItem(File file) {
        super(file);

        // Not accessible folders will return null listFiles()
        items = file.listFiles() != null ? file.listFiles().length : -1;
    }

    @Override
    public String getDescription() {
        return items >= 0 ? items + " items" : "[not accessible]";
    }
}
