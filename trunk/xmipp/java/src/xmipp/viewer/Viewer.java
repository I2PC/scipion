package xmipp.viewer;

import java.util.LinkedList;

import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import ij.IJ;
import xmipp.jni.Filename;
import xmipp.utils.DEBUG;
import xmipp.utils.Param;
import xmipp.viewer.windows.ImagesWindowFactory;


public class Viewer {
	public static void main(String args[]) {
		//Schedule a job for the event dispatch thread:
        //creating and showing this application's GUI.
//        SwingUtilities.invokeLater(new Runnable() {
//            public void run() {
//                //Turn off metal's use of bold fonts
//	        UIManager.put("swing.boldMetal", Boolean.FALSE);
//	       	ViewerTest.createAndShowGUI();
//            }
//        });
		Param parameters = new Param(args);
		try{
			if (parameters.debug)
	    		DEBUG.enableDebug(true);
			
	        if (parameters.files != null) {
	        	//DEBUG.enableDebug(true);
	            openFiles(parameters);	    
	        }
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

    static void openFiles(Param parameters) throws Exception {
    	DEBUG.printMessage(String.format("Viewer.openFiles"));
    	ImagesWindowFactory.openFilesAsDefault(parameters.files, parameters);
//        if (parameters.mode.trim().toLowerCase().compareTo(Param.OPENING_MODE_DEFAULT) == 0) {
//            ImagesWindowFactory.openFilesAsDefault(parameters.files, parameters);
//        } else if (parameters.mode.trim().toLowerCase().compareTo(Param.OPENING_MODE_IMAGE) == 0) {
//            ImagesWindowFactory.openFilesAsImages(parameters.files, parameters);
//        } else if (parameters.mode.trim().toLowerCase().compareTo(Param.OPENING_MODE_GALLERY) == 0) {
//            // Separates single images (index 0) from the rest of files (index 1).
//            String filesList[][] = getSeparatedFiles(parameters.files);
//
//            // Single images are opened in the same table...
//            if (filesList[0] != null) {
//                ImagesWindowFactory.openFilesAsGallery(filesList[0], true, parameters);
//            }
//
//            // ...while rest of files are opened in separated ones.
//            if (filesList[1] != null) {
//                ImagesWindowFactory.openFilesAsGallery(filesList[1], false, parameters);
//            }
//        } else if (parameters.mode.trim().toLowerCase().compareTo(Param.OPENING_MODE_METADATA) == 0) {
//            ImagesWindowFactory.openFilesAsMetadata(parameters.files, parameters);
//        }
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