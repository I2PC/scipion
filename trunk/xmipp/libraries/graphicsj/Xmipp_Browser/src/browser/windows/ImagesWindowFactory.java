/*
 * To change this template, choose Tools | Templates
 * and openTableFileImageItem the template in the editor.
 */
package browser.windows;

import browser.imageitems.listitems.AbstractImageItem;
import browser.imageitems.listitems.SelFileItem;
import browser.imageitems.TableImageItem;
import browser.table.JFrameImagesTable;
import browser.table.JFrameVolumeTable;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
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

    public static void openImageWindow(TableImageItem item) {
        if (item != null) {
            ImageWindowOperations iwo = new ImageWindowOperations(item.getImagePlus(), false);
            iwo.setTitle(item.getLabel());
        }
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

    public static void openTableImages(Vector<AbstractImageItem> images) {
        JFrameImagesTable imagesTable = new JFrameImagesTable(-1, -1);

        for (int i = 0; i < images.size(); i++) {
            imagesTable.addImage(images.elementAt(i));
        }

        imagesTable.setVisible(true);
    }

    public static JFrameVolumeTable openTableVolumeFile(String fileName, int w, int h) {
        JFrameVolumeTable volumeTable = new JFrameVolumeTable(w, h);

        volumeTable.addVolumeFile(new File(fileName));

        volumeTable.setVisible(true);

        volumeTable.setNormalizedAuto();

        return volumeTable;
    }

    public static JFrameImagesTable openTableSelFile(String fileName, int w, int h) {
        JFrameImagesTable selFileTable = new JFrameImagesTable(w, h);

        selFileTable.addSelFile(fileName);

        selFileTable.setVisible(true);

        return selFileTable;
    }

    public static void openTableSelFileItem(SelFileItem selfile) {
        JFrameImagesTable selFileTable = new JFrameImagesTable(-1, -1);

        selFileTable.addVolume(selfile.getDirectory(), selfile.getFile().getName(), selfile.nslices);

        selFileTable.setVisible(true);
    }

    public static void openTableFileImageItem(AbstractImageItem imageItem) {
        JFrameVolumeTable volumeTable = new JFrameVolumeTable();

        volumeTable.addVolume(imageItem);

        volumeTable.setVisible(true);

        volumeTable.setNormalizedAuto();    // Volumes are normalized at startup.
    }

    public static void openTableAsStack(Vector<TableImageItem> items) {
        int width = items.firstElement().getWidth();
        int height = items.firstElement().getHeight();
        String fileName = items.firstElement().getFileName();

        ImageStack stack = new ImageStack(width, height);

        for (int i = 0; i < items.size(); i++) {
            stack.addSlice(
                    String.valueOf(i),
                    items.elementAt(i).getImagePlus().getProcessor());
        }

        openImageWindow(new ImagePlus(fileName, stack), false);
    }

    //openStackAsTable(imp);
    public static void openStackAsTable(ImagePlus ip, boolean autoAdjust) {
        JFrameVolumeTable volumeTable = new JFrameVolumeTable();
        System.out.println("FN: " + ip.getFileInfo().fileName);
        for (int i = 0; i < ip.getStackSize(); i++) {
//            volumeTable.addVolume(new ImageItem(, null));
        }
    }
}
