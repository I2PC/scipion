package tests;

import xmipp.MDLabel;
import xmipp.MetaData;

public class MetaDataTest {

    public static void main(String args[]) {
        if (args.length < 1) {
            System.out.println("Usage: java MetaDataTest <xmipp_metadata_file>");
            System.exit(0);
        }

        try {
            String file = args[0];

            // Load metadata file.
            MetaData md = new MetaData(file);

            // For all items...
            long ids[] = md.findObjects();

            for (int i = 0; i < ids.length; i++) {
                // Shows image field.
                if (md.containsLabel(MDLabel.MDL_IMAGE)) {
                    String imagename = md.getValueString(MDLabel.MDL_IMAGE, ids[i]);
                    System.out.println(ids[i] + ": " + imagename);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
