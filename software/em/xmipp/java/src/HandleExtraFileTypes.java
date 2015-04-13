
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import xmipp.utils.Params;
import xmipp.viewer.windows.ImagesWindowFactory;

// Plugin to handle file types which are not implemented
// directly in ImageJ through io.Opener
// nb since there is no _ in the name it will not appear in Plugins menu
// -----
// Can be user modified so that your own specialised file types
// can be opened through File>Open
// or by drag and drop onto the ImageJ main window
// or by drag and drop onto the ImageJ icon
// or by double clicking on the file icon (requires that
// the file extension be associated with ImageJ).
// -----
// Go to the point marked MODIFY HERE and modify to
// recognise and load your own file type
// I have implemented two file types as examples:
// 	Biorad PIC and Gatan DM3
// -----
// Gregory Jefferis - 030629
// jefferis@stanford.edu
//
// Note that the File>Import>Image Sequence command and
// IJ.openImage() method will not work as expected unless the
// reader plugin extends the ImagePlus class.
/**
 *  Description of the Class
 *
 *@author     cedric
 *@created    26 juin 2006
 */
public class HandleExtraFileTypes extends ImagePlus implements PlugIn {

    final static int IMAGE_OPENED = -1;
    final static int PLUGIN_NOT_FOUND = -1;

    // Called from io/Opener.java
    /**
     *  Main processing method for the HandleExtraFileTypes object
     *
     *@param  path  Description of the Parameter
     */
    public void run(String path) {
        if (IJ.versionLessThan("1.30u")) {
            return;
        }
        if (path.equals("")) {
            return;
        }
        File theFile = new File(path);
        String directory = theFile.getParent();
        String fileName = theFile.getName();
        if (directory == null) {
            directory = "";
        }

        // Try and recognise file type and load the file if recognised
        ImagePlus imp = openImage(directory, fileName);
        if (imp == null) {
            // failed to load file or plugin has opened and displayed it
            IJ.showStatus("");
            return;
            // failed to load file or plugin has opened and displayed it
        }
        ImageStack stack = imp.getStack();
        // Set the stack of this HandleExtraFileTypes object
        // to that attached to the ImagePlus object returned by openImage()
        setStack(fileName, stack);
        // Copy over the calibration info since it doesn't come with the ImageProcessor
        setCalibration(imp.getCalibration());
        // Also copy the Show Info field over if it exists
        if (imp.getProperty("Info") != null) {
            setProperty("Info", imp.getProperty("Info"));
        }
        // Copy over the FileInfo
        setFileInfo(imp.getOriginalFileInfo());
    }

    /**
     *  Description of the Method
     *
     *@param  directory  Description of the Parameter
     *@param  name       Description of the Parameter
     *@return            Description of the Return Value
     */
    private ImagePlus openImage(String directory, String name) {
        ImagePlus imp;

        // Set out file name and path
        if (directory.length() > 0 && !directory.endsWith(Prefs.separator)) {
            directory += Prefs.separator;
        }
        String path = directory + name;

        // set up a stream to read in 132 bytes from the file header
        // These can be checked for "magic" values which are diagnostic
        // of some image types
        InputStream is;
        byte[] buf = new byte[132];
        try {
            is = new FileInputStream(path);
            is.read(buf, 0, 132);
            is.close();
        } catch (IOException e) {
            // Couldn't open the file for reading
            return null;
        }
        name = name.toLowerCase(); 

        // OK now we get to the interesting bit

        // CTR: try opening the file with LOCI Bio-Formats
        // http://www.loci.wisc.edu/ome/formats.html
		/*
        Object loci = IJ.runPlugIn("LociPlugin", path);
        if (loci != null) {
        plugin exists and was launched
        try {
        check whether plugin was successful
        boolean ok = loci.getClass().getField("success").getBoolean(loci);
        if (ok) {
        width = IMAGE_OPENED;
        return null;
        }
        }
        catch (Exception exc) { }
        }
         */
        // GJ: added Biorad PIC confocal file handler
        // Note that the Biorad_Reader plugin extends the ImagePlus class,
        // which is why the IJ.runPlugIn() call below returns an ImagePlus object.
        // ------------------------------------------
        // These make 12345 if you read them as the right kind of short
        // and should have this value in every Biorad PIC file
        if (buf[54] == 57 && buf[55] == 48) {
            // Ok we've identified the file type
            // Now load it using the relevant plugin
            imp = (ImagePlus) IJ.runPlugIn("Biorad_Reader", path);
            if (imp == null) {
                width = PLUGIN_NOT_FOUND;
            }
            if (imp != null && imp.getWidth() == 0) {
                imp = null;
            }
            return imp;
        }

        // Analyze format (.img/.hdr) handler
        // Note that the Analyze_Reader plugin opens and displays the
        // image and does not implement the ImagePlus class.
        if (name.endsWith(".img") || name.endsWith(".hdr")) {
            // Open Analyze image and display it
            IJ.runPlugIn("Analyze_Reader", path);
            // Set flag so Opener.openImage() does not display error
            width = IMAGE_OPENED;
            return null;
        }

        // IPLab file handler
        // Note that the IPLab_Reader plugin extends the ImagePlus class.
        // Little-endian IPLab files start with "iiii" or "mmmm".
        if ((buf[0] == 105 && buf[1] == 105 && buf[2] == 105 && buf[3] == 105)
                || (buf[0] == 109 && buf[1] == 109 && buf[2] == 109 && buf[3] == 109)) {
            imp = (ImagePlus) IJ.runPlugIn("IPLab_Reader", path);
            if (imp == null) {
                width = PLUGIN_NOT_FOUND;
            }
            if (imp != null && imp.getWidth() == 0) {
                imp = null;
            }
            return imp;
        }

        // Packard InstantImager format (.img) handler -> check HERE before Analyze check below!
        // Note that the InstantImager_Reader plugin extends the ImagePlus class.
        // Check extension and signature bytes KAJ_
        if (name.endsWith(".img") && buf[0] == 75 && buf[1] == 65 && buf[2] == 74 && buf[3] == 0) {
            imp = (ImagePlus) IJ.runPlugIn("InstantImager_Reader", path);
            if (imp == null) {
                width = PLUGIN_NOT_FOUND;
            }
            if (imp != null && imp.getWidth() == 0) {
                imp = null;
            }
            return imp;
        }

        // Analyze format (.img/.hdr) handler
        // Note that the Analyze_Reader plugin opens and displays the
        // image and does not implement the ImagePlus class.
        if (name.endsWith(".img") || name.endsWith(".hdr")) {
            // Open Analyze image and display it
            IJ.runPlugIn("Analyze_Reader", path);
            // Set flag so Opener.openImage() does not display error
            width = IMAGE_OPENED;
            return null;
        }

        // Image Cytometry Standard (.ics) handler
        // http://valelab.ucsf.edu/~nico/IJplugins/Ics_Opener.html
        if (name.endsWith(".ics")) {
            // Open ICS image and display it
            IJ.runPlugIn("Ics_Opener", path);
            // Set flag so Opener.openImage() does not display error
            width = IMAGE_OPENED;
            return null;
        }

        //  Zeiss Confocal LSM 510 image file (.lsm) handler
        //  http://rsb.info.nih.gov/ij/plugins/lsm-reader.html
        if (name.endsWith(".lsm")) {
            IJ.runPlugIn("LSM_Reader", path);
            width = IMAGE_OPENED;
            return null;
        }

        // BM: added Bruker  file handler 29.07.04
        if (name.equals("ser") || name.equals("fid") || name.equals("2rr") || name.equals("2ii") || name.equals("3rrr")
                || name.equals("3iii") || name.equals("2dseq")) {
            IJ.showStatus("Opening Bruker " + name + " File");
            IJ.runPlugIn("BrukerOpener", name + "|" + path);
            width = IMAGE_OPENED;
            return null;
        }

        // AVI: open AVI files using AVI_Reader plugin
        if (name.endsWith(".avi")) {
            IJ.runPlugIn("AVI_Reader", path);
            width = IMAGE_OPENED;
            return null;
        }

        // QuickTime: open .mov files using QT_Movie_Opener plugin
        if (name.endsWith(".mov")) {
            IJ.runPlugIn("QT_Movie_Opener", path);
            width = IMAGE_OPENED;
            return null;
        }

        // ZVI file handler
        // Little-endian ZVI and Thumbs.db files start with d0 cf 11 e0
        // so we can only look at the extension.
        if (name.endsWith(".zvi")) {
            IJ.runPlugIn("ZVI_Reader", path);
            width = IMAGE_OPENED;
            return null;
        }

        //University of North Carolina (UNC) file format handler
        // 'magic' numbers are (int) offsets to data structures and may change in future releases.
        if (name.endsWith(".unc") || (buf[3] == 117 && buf[7] == -127 && buf[11] == 36 && buf[14] == 32 && buf[15] == -127)) {
            IJ.runPlugIn("UNC_Reader", path);
            width = IMAGE_OPENED;
            return null;
        }

        //  Leica SP confocal .lei file handler
        if (name.endsWith(".lei")) {
            int dotIndex = name.lastIndexOf(".");
            if (dotIndex >= 0) {
                name = name.substring(0, dotIndex);
            }
            path = directory + name + ".txt";
            File f = new File(path);
            if (!f.exists()) {
                IJ.error("Cannot find the Leica information file: " + path);
                return null;
            }
            IJ.runPlugIn("Leica_TIFF_sequence", path);
            width = IMAGE_OPENED;
            return null;
        }

        // Amira file handler
        // http://wbgn013.biozentrum.uni-wuerzburg.de/ImageJ/amira-io.html
        if (buf[0] == 0x23 && buf[1] == 0x20 && buf[2] == 0x41
                && buf[3] == 0x6d && buf[4] == 0x69 && buf[5] == 0x72
                && buf[6] == 0x61 && buf[7] == 0x4d && buf[8] == 0x65
                && buf[9] == 0x73 && buf[10] == 0x68 && buf[11] == 0x20) {
            ImagePlus i = (ImagePlus) IJ.runPlugIn("AmiraMeshReader_", path);
            width = (i == null) ? PLUGIN_NOT_FOUND : IMAGE_OPENED;
            return null;
        }
        /* // Deltavision file handler
        // Open DV files generated on Applied Precision DeltaVision systems
        if (name.endsWith(".dv") || name.endsWith(".r3d")) {
        return tryPlugIn("Deltavision_Opener", path);
        }
        // Albert Cardona: read .dat files from the EMMENU software
        if (name.endsWith(".dat") && 1 == buf[1] && 0 == buf[2]) { // 'new format' only
        return tryPlugIn("Open_DAT_EMMENU", path);
        }        */

        if (name.toLowerCase().endsWith(".sdf")) {
            imp = (ImagePlus) IJ.runPlugIn("SDF_reader", path);
            if (imp == null) {
                width = PLUGIN_NOT_FOUND;
            }
            if (imp != null && imp.getWidth() == 0) {
                imp = null;
            }
            return imp;
        }
        // ****************** MODIFY HERE **	****************
        // Do what ever you have to do to recognise your own file type
        // and then call appropriate plugin
        // using the above as models
        // eg:

        
        // Obviouslly, IJ.Opener needs a valid ImagePlus. If for some reason we don't return it,
        // an alert will apear complaining something like "File is not supported format, a reader plugin
        // is not available..."
        try{
            System.err.println(path);

            Params param = new Params();

        	imp=ImagesWindowFactory.openFileAsImagePlus(path, param);
            if (imp == null) {
                width = PLUGIN_NOT_FOUND;
            }
            
            if (imp != null){
            	if(imp.getWidth() == 0) 
            		imp = null;
            	return imp;
            }
        	
        }catch (Exception ex){
        	ex.printStackTrace();
        }


        /*if (Filename.isMetadata(name)) {
            imp = (ImagePlus) tryPlugIn(xmipp.io.readers.MetaDataReader.class.getCanonicalName(), path);


        } else if (Filename.isXmippType(name)) {
            imp = (ImagePlus) tryPlugIn(xmipp.io.readers.ImageReader.class.getCanonicalName(), path);

            if (imp == null) {
                width = PLUGIN_NOT_FOUND;
            }
            if (imp != null && imp.getWidth() == 0) {
                imp = null;
            }

            return imp;
        }*/

        /*
        A. Dent: Added XYZ handler
        ----------------------------------------------
        Check if the file ends in .XYZ or .xyz
        if (name.toUpperCase().endsWith(".XYZ")) {
        Bytes 0 and 1 must equal 42 for this file type
        if(buf[0]==42 && buf[1]==42) {
        Ok we've identified the file type - now load it
        imp = (ImagePlus)IJ.runPlugIn("XYZ_Reader", path);
        if (imp==null) width = PLUGIN_NOT_FOUND;
        if (imp!=null&&imp.getWidth()==0) imp = null;
        return imp;
        }
        }
         */
        // If we got this far we didn't recognise the file type
        return null;
    }

    /**
     * Attempts to open the specified path with the given plugin. If the
     * plugin extends the ImagePlus class (e.g., BioRad_Reader), set
     * extendsImagePlus to true, otherwise (e.g., LSM_Reader) set it to false.
     *
     * @return A reference to the plugin, if it was successful.
     */
    private Object tryPlugIn(String className, String path) {
        Object o = IJ.runPlugIn(className, path);
        if (o instanceof ImagePlus) {
            // plugin extends ImagePlus class
            ImagePlus imp = (ImagePlus) o;
            if (imp.getWidth() == 0) {
                o = null; // invalid image
            } else {
                width = IMAGE_OPENED; // success
            }
        } else {
            // plugin does not extend ImagePlus; assume success
            width = IMAGE_OPENED;
        }
        return o;
    }
    
    
}
