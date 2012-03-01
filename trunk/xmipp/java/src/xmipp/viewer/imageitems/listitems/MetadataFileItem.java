/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package xmipp.viewer.imageitems.listitems;

import xmipp.utils.Cache;
import xmipp.utils.DEBUG;
import xmipp.utils.XmippResource;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.ij.XmippImageConverter;

/**
 *
 * @author Juanjo Vega
 */
public class MetadataFileItem extends XmippImageItem {

    protected String previewFile;
    protected MetaData md;

    public MetadataFileItem(File file, Cache cache) {
        super(file, cache);
    }

    @Override
    void loadImageData() {
        try {
            String path = file.getAbsolutePath();
            md = new MetaData();
            md.read(path);

            previewFile = getPreviewFile(md);
            if (previewFile != null) {
                loadPreviewFileData();

                dimension.setNImages(md.size());
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
            IJ.error(ex.getMessage());
        }
    }

    // Gets first available filename if any.
    protected String getPreviewFile(MetaData md) {
        try {
            if (md.containsLabel(MDLabel.MDL_IMAGE)) {
                File f;

                long objs[] = md.findObjects();

                for (long id : objs) {
                    String field = md.getValueString(MDLabel.MDL_IMAGE, id, true);

                    field = Filename.getFilename(field);      // Avoids image@filename format ;)
                    nimage = Filename.getNimage(field);

                    f = new File(field);

                    if (f.exists()) {
//                    System.err.println(" *** preview File: " + f.getAbsolutePath());
                        return f.getAbsolutePath();
                    }
                }
            }
        } catch (Exception ex) {
            DEBUG.printException(ex);
        }

        return null;
    }

    protected void loadPreviewFileData() {
        // Tricky way to get preview file's data.
        File originalFile = file;

        file = new File(previewFile);

        super.loadImageData();

        file = originalFile;
    }

    @Override
    public ImagePlus getPreview(int w, int h) {
        ImagePlus ip = XmippResource.MISSING_ITEM;

        if (previewFile != null) {
            // Tricky way to get preview.
            File originalFile = file;
            file = new File(Filename.getFilename(previewFile)); // skips image@filename format

            ip = super.getPreview(w, h);

            file = originalFile;
        }

        return ip;
    }

    @Override
    public ImagePlus getImagePlus() {
        try {
            return XmippImageConverter.readMetadataToImagePlus(md);
        } catch (Exception ex) {
            return null;
        }
    }
}
