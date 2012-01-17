package xmipp.viewer;

import java.util.LinkedList;

import ij.IJ;
import xmipp.jni.Filename;
import xmipp.utils.Param;
import xmipp.viewer.windows.ImagesWindowFactory;


public class Viewer {
	public static void main(String args[]) {
		Param parameters = new Param(args);
		try{
	        if (parameters.files != null) 
	            openFiles(parameters);	        
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

    static void openFiles(Param parameters) {
        if (parameters.mode.trim().toLowerCase().compareTo(Param.OPENING_MODE_DEFAULT) == 0) {
            ImagesWindowFactory.openFilesAsDefault(parameters.files, parameters);
        } else if (parameters.mode.trim().toLowerCase().compareTo(Param.OPENING_MODE_IMAGE) == 0) {
            ImagesWindowFactory.openFilesAsImages(parameters.files, parameters);
        } else if (parameters.mode.trim().toLowerCase().compareTo(Param.OPENING_MODE_GALLERY) == 0) {
            // Separates single images (index 0) from the rest of files (index 1).
            String filesList[][] = getSeparatedFiles(parameters.files);

            // Single images are opened in the same table...
            if (filesList[0] != null) {
                ImagesWindowFactory.openFilesAsGallery(filesList[0], true, parameters);
            }

            // ...while rest of files are opened in separated ones.
            if (filesList[1] != null) {
                ImagesWindowFactory.openFilesAsGallery(filesList[1], false, parameters);
            }
        } else if (parameters.mode.trim().toLowerCase().compareTo(Param.OPENING_MODE_METADATA) == 0) {
            ImagesWindowFactory.openFilesAsMetadata(parameters.files, parameters);
        }
    }

    static String[][] getSeparatedFiles(String items[]) {
        LinkedList<String> images = new LinkedList<String>();
        LinkedList<String> files = new LinkedList<String>();
        LinkedList<String> storeTo;

        for (int i = 0; i < items.length; i++) {
            String filename = items[i];

            try {
                storeTo = Filename.isSingleImage(filename) ? images : files;
                storeTo.add(filename);
            } catch (Exception e) {
                IJ.error(e.getMessage());
            }
        }

        String result[][] = new String[2][];

        if (!images.isEmpty()) {
            result[0] = images.toArray(new String[images.size()]);
        }

        if (!files.isEmpty()) {
            result[1] = files.toArray(new String[files.size()]);
        }

        return result;
    }
}