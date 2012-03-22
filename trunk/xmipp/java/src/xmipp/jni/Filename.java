package xmipp.jni;

import java.io.File;
import java.io.FilenameFilter;
import java.net.URI;

public class Filename {

    public final static String PROJECT_FILE = ".project.sqlite";
    public final static String SEPARATOR = "@";
    public final static String EXT_XMP = ".xmp";
    public final static String EXT_IMG = ".img";
    public final static String EXT_HED = ".hed";
    public final static String EXT_PSD = ".psd";
    public final static String EXT_SER = ".ser";
    public final static String EXT_DM3 = ".dm3";
    public final static String EXT_RAW = ".raw";
    public final static String EXT_INF = ".inf";
    public final static String EXT_SPE = ".spe";
    public final static String EXT_MRC = ".mrc";
    public final static String EXT_MRCS = ".mrcs";
    public final static String EXT_MRCS2 = ".st";
    public final static String EXT_STK = ".stk";
    //metadata extensions
    public final static String EXT_XMD = ".xmd";
    public final static String EXT_SEL = ".sel";
    public final static String EXT_DOC = ".doc";
    public final static String EXT_CTFPARAM = ".ctfparam";
    public final static String EXT_CTFDAT = ".ctfdat";
    public final static String EXT_POS = ".pos";
    public final static String EXT_VOL = ".vol";
    public final static String EXT_SPI = ".spi";
    public final static String EXT_TIF = ".tif";
    public final static String EXT_LOG = ".log";
    public final static String EXT_OUT = ".out";
    public final static String EXT_ERR = ".err";
    public final static String EXT_PY = ".py";
    public final static String EXT_TXT = ".txt";
    
    
    // Initialize library.
    static {
        System.loadLibrary("XmippJNI");
        //storeIds();
    }

    public final static String[] SINGLE_IMAGES = new String[]{
        EXT_XMP,
        EXT_IMG,
        EXT_HED,
        EXT_PSD,
        EXT_SER,
        EXT_DM3,
        EXT_RAW,
        EXT_INF,
        EXT_SPE,
        EXT_SPI,
        EXT_TIF
    };
    public final static String[] VOLUMES = new String[]{
        EXT_MRC,
        EXT_VOL
    };
    public final static String[] STACKS = new String[]{
        EXT_MRCS,
        EXT_MRCS2,
        EXT_STK
    };
    
    public final static String[] METADATAS = new String[]{
    	EXT_XMD,
        EXT_SEL,
        EXT_DOC,
        EXT_CTFPARAM,
        EXT_CTFDAT,
        EXT_POS
    };
    
    public final static String[] SPIDER = new String[] {
    	EXT_SPI,
    	EXT_VOL
    };
    
    public final static String[] TEXT = new String[]{
    	EXT_TXT, EXT_LOG, EXT_ERR, EXT_OUT
    };

    public static boolean isPSD(String filename) {
        return filename != null && filename.endsWith(EXT_PSD);
    }
    
    public static boolean isSpiderVolume(String filename){
    	return filename != null && isFileType(filename, SPIDER);
    }
//
//    public static boolean isSingleImage(String filename) {
//        return filename != null && (filename.contains(SEPARATOR) || isFileType(filename, SINGLE_IMAGES));
//    }
//
//    public static boolean isVolume(String filename) {
//        return filename != null && isFileType(filename, VOLUMES);
//    }
//
//    public static boolean isStack(String filename) {
//        return filename != null && isFileType(filename, STACKS);
//    }
//
//    public static boolean isStackOrVolume(String filename) {
//        return filename != null && (isStack(filename) || isVolume(filename));
//    }
//

    public static native boolean hasStackExtension(String filename) throws Exception;

    public static native boolean hasVolumeExtension(String filename) throws Exception;
    
    public static native String compose(int slice,String path) throws Exception;
    
    public static native String getXmippPath() throws Exception;
    
    public static String getXmippPath(String relpath){
    	try {
			return getXmippPath() + File.separator + relpath;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	return null;
    }

    public static boolean isSingleImage(String filename) throws Exception {
        try {
            return (new ImageGeneric(filename)).isSingleImage();
        } catch (Exception ex) {
            return filename != null && isFileType(filename, SINGLE_IMAGES);
        }
    }

    public static boolean isVolume(String filename) {
        try {
            return (new ImageGeneric(filename)).isVolume();
        } catch (Exception ex) {
            return filename != null && isFileType(filename, VOLUMES);
        }
    }

    public static boolean isStack(String filename) throws Exception {
        try {
            return (new ImageGeneric(filename)).isStack();
        } catch (Exception ex) {
            return filename != null && isFileType(filename, STACKS);
        }
    }

    public static boolean isMetadata(String filename) {
        return filename != null && isFileType(filename, METADATAS);
    }
    
    public static boolean isTextfile(String filename){
    	return isMetadata(filename) || isFileType(filename, TEXT);
    }

    private static boolean isFileType(String filename, String filetypes[]) {
        for (int i = 0; i < filetypes.length; i++) {
            if (filename.endsWith(filetypes[i])) {
                return true;
            }
        }

        return false;
    }

    // Auxiliary methods.
    public static String fixPath(String filename, String MDdir, boolean shouldExist) {
        MDdir += !MDdir.endsWith(File.separator) ? File.separator : "";
        String fixed = filename;

        if (!filename.startsWith(File.separator)) { // Absolute path?
            String name = Filename.getFilename(filename);
            String strprefix = "";

            if (filename.contains(Filename.SEPARATOR)) { // Has #image?
                String prefix = Filename.getPrefix(filename);
                strprefix = prefix + Filename.SEPARATOR;
            }

            // Checks if path is absolute...
            if (!name.startsWith(File.separator)) {
                // ...if not: tries to build the absolute path:
                // 1st case: Relative to metadata file (metadata_path + file)
                String aux = URI.create(MDdir + name).normalize().getPath();
                File f = new File(aux);
                if (shouldExist && !f.exists()) {
                    // 2nd case: Relative to current dir.
                    aux = URI.create(System.getProperty("user.dir") + File.separatorChar + name).normalize().getPath();
                    f = new File(aux);
                    if (!f.exists()) {
                        // 3rd case: find "project dir" (the one containing a file called ".project.sqlite")
                        String projectdir = findProjectDir(MDdir, Filename.PROJECT_FILE);
                        if (projectdir != null) {
                            aux = URI.create(projectdir + name).normalize().getPath();
                        }
                    }
                }

                fixed = strprefix + aux;
            }
        }

        return fixed;
    }
    
    public static boolean exists(String path){
    	File f = new File(Filename.getFilename(path));
    	return f.exists();
    }

    public String findProjectDir(String metadata) {
        File f = new File(metadata);
        String startingdir = f.isDirectory() ? metadata : f.getParent();

        return findProjectDir(startingdir, PROJECT_FILE);
    }

    private static String findProjectDir(String current, final String PROJECT_FILE) {
        FilenameFilter filter = new FilenameFilter() {

            @Override
            public boolean accept(File dir, String name) {
                return PROJECT_FILE.compareTo(name) == 0;
            }
        };

        File dir = new File(current);
        String files[] = dir.list(filter);

        if (files == null || files.length == 0) {
            String parentdir = dir.getParent();

            if (parentdir != null) {
                return findProjectDir(dir.getParent(), PROJECT_FILE);
            } else {
                return null;
            }
        }

        return dir.toURI().normalize().getPath();
    }

    public static String getFilename(String filename) {
        if (filename.contains(SEPARATOR)) {
            return filename.split(SEPARATOR)[1];
        }

        return filename;
    }
    
    public static String getPath(String baseDir,String fileName,int slice) throws Exception{
    	return compose(slice,baseDir+File.separatorChar+fileName);
    }

    public static long getNimage(String filename) {
        String prefix = getPrefix(filename);

        return prefix != null ? Long.valueOf(prefix).longValue() : ImageGeneric.ALL_IMAGES;
    }

    public static boolean hasPrefix(String filename) {
        return filename.contains(SEPARATOR);
    }

    public static String getPrefix(String filename) {
        if (hasPrefix(filename)) {
            String prefix = "";
            String str = filename.split(SEPARATOR)[0];
            if (!str.isEmpty()) {
                // str may have a string prefix before the number, so
                // grab the leftmost part
                int i = str.length() - 1;
                while (i >= 0) {
                    if (str.charAt(i) == File.separatorChar) {
                        break;
                    }
                    i--;
                }

                prefix = str.substring(i + 1, str.length());
            }
            return prefix;
        }

        return null;
    }
    
    /** Return the name of the block associated with
     * class 'ref' in a results_classes.xmd metadata
     */
    public static String getClassBlockName(int ref){
    	return String.format("class%06d_images", ref);
    }
}
