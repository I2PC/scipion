package xmipp.jni;

import java.io.File;
import java.io.FilenameFilter;
import java.net.URI;

public class Filename {

    public final static String PROJECT_FILE = ".project.sqlite";
    public final static String SEPARATOR = "@";
    private final static String EXTENSION_XMP = ".xmp";
    private final static String EXTENSION_IMG = ".img";
    private final static String EXTENSION_HED = ".hed";
    private final static String EXTENSION_PSD = ".psd";
    private final static String EXTENSION_SER = ".ser";
    private final static String EXTENSION_DM3 = ".dm3";
    private final static String EXTENSION_RAW = ".raw";
    private final static String EXTENSION_INF = ".inf";
    private final static String EXTENSION_SPE = ".spe";
    private final static String EXTENSION_MRC = ".mrc";
    private final static String EXTENSION_MRCS = ".mrcs";
    private final static String EXTENSION_STK = ".stk";
    //metadata extensions
    private final static String EXTENSION_XMD = ".xmd";
    private final static String EXTENSION_SEL = ".sel";
    private final static String EXTENSION_DOC = ".doc";
    private final static String EXTENSION_CTFPARAM = ".ctfparam";
    private final static String EXTENSION_CTFDAT = ".ctfdat";
    private final static String EXTENSION_POS = ".pos";
    private final static String EXTENSION_VOL = ".vol";
    private final static String EXTENSION_SPI = ".spi";
    private final static String EXTENSION_TIF = ".tif";
    // Types of images contained by each file type.
//    private final static String[] XMIPP_TYPES = new String[]{
//        EXTENSION_XMP,
//        EXTENSION_IMG,
//        EXTENSION_HED,
//        EXTENSION_PSD,
//        EXTENSION_SER,
//        EXTENSION_DM3,
//        EXTENSION_RAW,
//        EXTENSION_INF,
//        EXTENSION_SPE,
//        EXTENSION_MRC,
//        EXTENSION_MRCS,
//        EXTENSION_STK,
//        EXTENSION_SEL,
//        EXTENSION_VOL,
//        EXTENSION_XMD,
//        EXTENSION_SPI,
//        EXTENSION_TIF
//    };
    public final static String[] SINGLE_IMAGES = new String[]{
        EXTENSION_XMP,
        EXTENSION_IMG,
        EXTENSION_HED,
        EXTENSION_PSD,
        EXTENSION_SER,
        EXTENSION_DM3,
        EXTENSION_RAW,
        EXTENSION_INF,
        EXTENSION_SPE,
        EXTENSION_SPI,
        EXTENSION_TIF
    };
    public final static String[] VOLUMES = new String[]{
        EXTENSION_MRC,
        EXTENSION_VOL
    };
    public final static String[] STACKS = new String[]{
        EXTENSION_MRCS,
        EXTENSION_STK
    };
    
    public final static String[] METADATAS = new String[]{
    	EXTENSION_XMD,
        EXTENSION_SEL,
        EXTENSION_DOC,
        EXTENSION_CTFPARAM,
        EXTENSION_CTFDAT,
        EXTENSION_POS
    };
    
    public final static String[] SPIDER = new String[] {
    	EXTENSION_SPI,
    	EXTENSION_VOL
    };

    public static boolean isPSD(String filename) {
        return filename != null && filename.endsWith(EXTENSION_PSD);
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
//
//    public static boolean isXmippType(String filename) {
//        return isFileType(filename, XMIPP_TYPES);
//    }

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
}
