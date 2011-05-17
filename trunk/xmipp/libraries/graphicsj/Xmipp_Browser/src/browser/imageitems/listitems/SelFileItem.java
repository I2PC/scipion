/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.imageitems.listitems;

import browser.Cache;
import browser.imageitems.ImageConverter;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import java.io.File;
import xmipp.ImageDouble;
import xmipp.MDLabel;
import xmipp.MetaData;

/**
 *
 * @author Juanjo Vega
 */
public class SelFileItem extends XmippImageItem {

    protected String previewFile;
    protected MetaData md;

    public SelFileItem(File file, Cache cache) {
        super(file, cache);
    }
    /*        loadImageData();
    }*/

    @Override
    protected void loadImageData() {
        try {
            String path = file.getAbsolutePath();
            md = new MetaData();
            md.read(path);

            previewFile = getPreviewFile(md);
            if (previewFile != null) {
                loadPreviewFileData();
            }

            dimension.setNimages(md.size());
        } catch (Exception ex) {
            ex.printStackTrace();
            IJ.error(ex.getMessage());
            throw new RuntimeException(ex);
        }
    }

    // Gets first available filename if any.
    protected String getPreviewFile(MetaData md) {
        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            File f;

            long objs[] = md.findObjects();

            for (long id : objs) {
                String field = md.getValueString(MDLabel.MDL_IMAGE, id);
                field = MetaData.getFilename(field);      // Avoids image@filename format ;)
                nimage = MetaData.getNimage(field);

                f = new File(field);

                if (f.exists()) {
                    return f.getAbsolutePath();
                }
            }
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
        // Tricky way to get preview.
        File originalFile = file;
        file = new File(MetaData.getFilename(previewFile)); // skips image@filename format

        ImagePlus ip = super.getPreview(w, h);

        file = originalFile;

        return ip;
    }

    @Override
    public ImagePlus getImagePlus(int image) {
        ImagePlus ip = null;

        try {
            ImageStack is = null;
            String images[] = getFileNames(md, file.getParent() + File.separator);

            for (int i = 0; i < images.length; i++) {
                ImageDouble img = new ImageDouble();
                img.read(images[i]);

                ImagePlus slice_ = ImageConverter.convertToImagej(img, images[i]);

                if (is == null) {
                    is = new ImageStack(slice_.getWidth(), slice_.getHeight());
                }

                is.addSlice(images[i], slice_.getProcessor());

                ip = new ImagePlus(file.getAbsolutePath(), is);
            }
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            throw new RuntimeException(ex);
        }

        return ip;
    }

    public String[] getFileNames() {
        return getFileNames(md, file.getParent() + File.separator);
    }

    protected static String[] getFileNames(MetaData md, String parent) {
        String filenames[] = null;

        // Skips if there are no images.
        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            filenames = new String[md.size()];
            File f;
            //String parent = file.getParent() + File.separator;

            int i = 0;
            //md.iteratorBegin();
            long objs[] = md.findObjects();

            //do {
            for (long id : objs) {
                //md.getStrFromValue(MDLabel.MDL_IMAGE, field);
                String field = md.getValueString(MDLabel.MDL_IMAGE, id);

                f = new File(field);
                if (!f.isAbsolute()) {
                    f = new File(parent + field);
                }

                filenames[i++] = f.getAbsolutePath();

                //md.iteratorNext();
            }
            //while (!md.iteratorEnd());
        }

        return filenames;
    }
}
