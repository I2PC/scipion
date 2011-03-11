/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package browser.windows;

import browser.files.FilterFilesModel;
import browser.imageitems.ImageConverter;
import browser.imageitems.TableImageItem;
import browser.imageitems.listitems.FileItem;
import browser.imageitems.listitems.XmippImageItem;
import browser.imageitems.listitems.SelFileItem;
import browser.table.JFrameImagesTable;
import browser.table.JFrameVolumeTable;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import java.io.File;
import java.util.Vector;
import xmipp.ImageDouble;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesWindowFactory {

    private final static String TEMPDIR_PATH = System.getProperty("java.io.tmpdir");

    public static void openImageFiles(String files[]) {
        Vector<XmippImageItem> xmippItems = new Vector<XmippImageItem>();

        for (int i = 0; i < files.length; i++) {
            FileItem item = FilterFilesModel.createSuitableFileItem(new File(files[i]));
            if (item instanceof XmippImageItem) {
                xmippItems.add((XmippImageItem) item);
            }
        }

        for (int i = 0; i < xmippItems.size(); i++) {
            openImage(xmippItems.elementAt(i));
        }
    }

    public static void openVolumeFile(String filename) {
        FileItem item = FilterFilesModel.createSuitableFileItem(new File(filename));
        if (item instanceof SelFileItem) {
            openTable((SelFileItem) item);
        } else if (item instanceof XmippImageItem) {
            openTable((XmippImageItem) item);
        }
    }

    public static void openImage(TableImageItem item) {
        ImageWindow iw = openImage(item.getImagePlus());
        iw.setTitle(item.getLabel());
    }

    public static void openImage(SelFileItem item) {
        ImageWindow iw = openImage(item.getImagePlus());
        iw.setTitle(item.getLabel());
    }

    public static void openImage(XmippImageItem item) {
        openImage(item, 0);
    }

    public static void openImage(XmippImageItem item, int n) {
        ImageWindow iw = openImage(item.getImagePlus(n));
        iw.setTitle(item.getLabel());
    }

    public static ImageWindow openImage(ImagePlus ip) {
        return openImage(ip, false);
    }

    public static ImageWindow openImage(ImagePlus imp, boolean poll) {
        ImageWindow iw = null;

        if (imp != null) {
            if (imp.getStackSize() > 1) {
                iw = new StackWindowOperations(imp, poll);
            } else {
                iw = new ImageWindowOperations(imp, poll);
            }
        } else {
            IJ.error("Trying to open a null image.");
        }

        return iw;
    }

    public static void openImage(Vector<TableImageItem> items, String title) {
        openImage(ImageConverter.toImagePlus(items, title));
    }

    public static void openTable(ImagePlus item) {
        try {
            // Saves a temporary file...
            String filename = TEMPDIR_PATH + File.separator + item.getFileInfo().fileName + System.currentTimeMillis() + ".xmp";

            System.err.println(" >>> Saving: " + filename);
            ImageDouble img = ImageConverter.convertToXmipp(item);
            img.write(filename);

            openVolumeFile(filename);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static void openTable(Object items[]) {
        Vector<XmippImageItem> images = new Vector<XmippImageItem>();

        for (int i = 0; i < items.length; i++) {
            Object item = items[i];

            if (item instanceof XmippImageItem) {
                XmippImageItem imageItem = (XmippImageItem) item;

                // Volumes are opened in a independent table for each one
                if (!imageItem.isSingleImage()) {
                    if (item instanceof SelFileItem) {
                        openTable((SelFileItem) imageItem);
                    } else {
                        openTable(imageItem);
                    }
                } else {    // Images will be at the same table.
                    images.add(imageItem);
                }
            } else {
                IJ.error(((FileItem) item).getFile().getName() + " is not an xmipp file.");
            }
        }

        // If there was any image file.
        if (images.size() > 0) {
            openTable(images);
        }
    }

    public static void openTable(Vector<XmippImageItem> images) {
        openTable(images, -1, -1);
    }

    public static void openTable(Vector<XmippImageItem> images, int w, int h) {
        JFrameImagesTable imagesTable = new JFrameImagesTable(h, w);

        for (int i = 0; i < images.size(); i++) {
            imagesTable.addImageItem(images.elementAt(i));
        }

        imagesTable.setVisible(true);
    }

    public static void openTable(XmippImageItem item) {
        JFrameVolumeTable volumeTable = new JFrameVolumeTable();

        volumeTable.addImageItem(item);
        volumeTable.setVisible(true);

        //@TODO Normalize
        //volumeTable.setNormalizedAuto();    // Volumes are normalized at startup.
    }

    public static void openTable(SelFileItem item) {
        JFrameVolumeTable volumeTable = new JFrameVolumeTable();

        volumeTable.addImageItem(item);
        volumeTable.setVisible(true);
        volumeTable.setNormalizedAuto();    // Volumes are normalized at startup.
    }
}
