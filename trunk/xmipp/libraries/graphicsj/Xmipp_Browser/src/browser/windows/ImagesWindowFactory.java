/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package browser.windows;

import browser.imageitems.ImageItem;
import browser.imageitems.SelFileItem;
import browser.imageitems.listitems.FileImageItem;
import browser.table.JFrameImagesTable;
import browser.table.JFrameVolumeTable;
import ij.IJ;
import ij.ImagePlus;
import java.io.File;
import java.util.Vector;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesWindowFactory {

    public static ImagePlus openImageWindow(File file) {
        return openImageWindow(file, -1);
    }

    public static ImagePlus openImageWindow(String file, boolean poll) {
        File f = new File(file);
        ImagePlus ip = IJ.openImage(f.getAbsolutePath());
        openImageWindow(ip, poll);

        return ip;
    }

    public static ImagePlus openImageWindow(File file, int slice) {
        ImagePlus ip = null;

        if (file != null) {
            ip = IJ.openImage(file.getAbsolutePath());

            if (ip != null) {
                ImagePlus ip2 = ip;

                if (slice >= 0) {
                    ip2 = new ImagePlus();
                    ip2.getProcessor().setPixels(ip.getStack().getPixels(slice));
                }

                openImageWindow(ip2, false);
            } else {
                // Not an image, so let ImageJ openTableFileImageItem it as usually.
                IJ.open(file.getAbsolutePath());
            }
        } else {
            IJ.showMessage("Trying to open a null file.");
        }

        return ip;
    }

    public static void openImageWindow(ImagePlus ip, boolean poll) {
        if (ip != null) {
            if (ip.getStackSize() > 1) {
                new StackWindowOperations(ip, poll);
            } else {
                new ImageWindowOperations(ip, poll);
            }
        } else {
            IJ.showMessage("Trying to open a null image.");
        }
    }

    public static void openTableImages(Vector<ImageItem> images) {
        JFrameImagesTable imagesTable = new JFrameImagesTable();

        for (int i = 0; i < images.size(); i++) {
            imagesTable.addImage(images.elementAt(i));
        }

        imagesTable.setVisible(true);
    }

    public static void openTableVolumeFile(String file) {
        JFrameVolumeTable volumeTable = new JFrameVolumeTable();

        volumeTable.addVolumeFile(new File(file));

        volumeTable.setVisible(true);
    }

    public static void openTableSelFile(String file) {
        JFrameImagesTable selFileTable = new JFrameImagesTable();

        File selfile = new File(file);

        String name = selfile.getName();
        String parent = selfile.getParent();
        if (parent == null) {
            parent = System.getProperty("user.dir");
        }

        selFileTable.addSelFile(parent, name);

        selFileTable.setVisible(true);
    }

    public static void openTableSelFileItem(SelFileItem selfile) {
        JFrameImagesTable selFileTable = new JFrameImagesTable();

        selFileTable.addVolume(selfile.getDirectory(), selfile.getFile().getName(), selfile.nslices);

        selFileTable.setVisible(true);
    }

    public static void openTableFileImageItem(FileImageItem imageItem) {
        JFrameVolumeTable volumeTable = new JFrameVolumeTable();

        volumeTable.addVolume(imageItem);

        volumeTable.setVisible(true);

        volumeTable.normalizeAuto();    // Volumes are normalized at startup.
    }
}
