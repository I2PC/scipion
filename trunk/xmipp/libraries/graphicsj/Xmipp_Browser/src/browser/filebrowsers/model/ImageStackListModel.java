/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.filebrowsers.model;

import browser.Cache;
import browser.imageitems.listitems.MetadataImageItem;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import java.util.ArrayList;
import javax.swing.AbstractListModel;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class ImageStackListModel extends AbstractListModel {

    ArrayList<MetadataImageItem> list;
    protected static Cache<String, ImagePlus> cache = new Cache<String, ImagePlus>();

    public ImageStackListModel(String metadata) {
        super();

        list = loadMetadata(metadata);
    }

    final ArrayList<MetadataImageItem> loadMetadata(String metadata) {
        try {
            ArrayList<MetadataImageItem> items = new ArrayList<MetadataImageItem>();
            MetaData md = new MetaData(metadata);

            long ids[] = md.findObjects();
            for (int i = 0; i < ids.length; i++) {
                String original = md.getValueString(MDLabel.MDL_IMAGE, ids[i]);
                String absolute = md.getValueString(MDLabel.MDL_IMAGE, ids[i], true);

                File file = new File(absolute);
                MetadataImageItem item = new MetadataImageItem(file, original, cache);

                items.add(item);
            }

            return items;
        } catch (Exception e) {
            IJ.error(e.getMessage());
        }

        return null;
    }

    public int getSize() {
        return list.size();
    }

    public Object getElementAt(int i) {
        return list.get(i);
    }
}
