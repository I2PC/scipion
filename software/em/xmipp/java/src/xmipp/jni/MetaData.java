/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

package xmipp.jni;

import ij.IJ;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.IllegalFormatCodePointException;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * Binding class for accessing C++ MetaData implementation.
 */
public class MetaData {
    
	/** Enum values with Labels possible types */
	public static final int LABEL_NOTYPE = -1;
	public static final int LABEL_INT = 0;
	public static final int LABEL_BOOL = 1;
	public static final int LABEL_DOUBLE = 2;
	public static final int LABEL_STRING = 3;
	public static final int LABEL_VECTOR_DOUBLE = 4;
	public static final int LABEL_SIZET = 5;
	public static final int LABEL_VECTOR_SIZET = 6;

	public static final String FILL_CONSTANT = "constant";
	public static final String FILL_LINEAR = "linear";
	public static final String FILL_RAND_UNIFORM = "random uniform";
	public static final String FILL_RAND_GAUSSIAN = "random gaussian";
        
       

	//
	// // Fields whose content is a path. They will be "fixed" conveniently.
	// private final static int PATHS_FIELDS[] = {
	// MDLabel.MDL_IMAGE1,
	// MDLabel.MDL_IMAGE2,
	// MDLabel.MDL_IMAGE3,
	// MDLabel.MDL_IMAGE,
	// MDLabel.MDL_PSD,
	// MDLabel.MDL_CTF_MODEL
	// };
	//
	// static {
	// // Sorts it to use binary search later.
	// // (It's executed just for the first time)
	// Arrays.sort(PATHS_FIELDS);
	// }

	public final static int GEOMETRY_LABELS[] = { MDLabel.MDL_FLIP,
			MDLabel.MDL_ANGLE_PSI, MDLabel.MDL_SHIFT_X, MDLabel.MDL_SHIFT_Y};
	
	public final static int GEOMETRY_RELION_LABELS[] = { MDLabel.RLN_ORIENT_ORIGIN_X,  MDLabel.RLN_ORIENT_ORIGIN_Y, MDLabel.RLN_ORIENT_PSI, MDLabel.RLN_ORIENT_FLIP};

	public final static int MICROGRAPH_BASIC_LABELS[] = {
			MDLabel.MDL_MICROGRAPH, MDLabel.MDL_PSD, MDLabel.MDL_CTF_MODEL };

	public static final int MD_OVERWRITE = 0;
	public static final int MD_APPEND = 1;
	protected String filename;
	// hold pointer to Image class in C++ space
	private long peer;

	// keep labels for avoid read all the time
	int[] activeLabels;
        
        protected HashMap<Long, EllipseCTF> ctfs = new HashMap<Long, EllipseCTF>();
        

	static {
		System.loadLibrary("XmippJNI");
		// storeIds();
	}

        	/** Create empty metadata */
	public MetaData() {
		//System.out.format("Java: Creating metadata\n");
		//Exception ex = new Exception();
    		//ex.printStackTrace();
		create();
	}

	
	// caching some ids
	// private static native void storeIds();
	// functions to create images
	private native void create();

	// destructor
	public synchronized native void destroy();
	
	// clear metadata
	public native void clear();

	// reading
	public native void read_(String filename);

	public final void read(String filename) {
		this.filename = filename;
		read_(filename);
		activeLabels = getActiveLabels();
	}

	public native int size();

	public native void setColumnFormat(boolean format);

	public native void write(String filename);

	public native void writeBlock(String filename);

	public native void print();

	public native boolean isColumnFormat();

	public native boolean containsLabel(int label);

	public native boolean removeLabel(int label);

	public boolean containsGeometryInfo(String type) {
		
		
		if (type.equals("2D")) {
			return containsGeometryInfo2D();
		
		} else {
		
			return containsGeometryInfo3D();
		}
		
	}
	
	
	private boolean containsGeometryInfo2D(){

		try {
			boolean contains = containsLabel(GEOMETRY_LABELS) && !containsLabel(MDLabel.MDL_ANGLE_ROT);
			if(!contains)
				contains = containsLabel(GEOMETRY_RELION_LABELS) &&
				((!containsLabel(MDLabel.RLN_ORIENT_ROT)) || ( getColumnMax(MDLabel.RLN_ORIENT_ROT) == 0 )) ;
			return contains;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;	
		
	}
	
	private boolean containsGeometryInfo3D(){

		try {
			boolean contains = containsLabel(GEOMETRY_LABELS) && containsLabel(MDLabel.MDL_ANGLE_ROT);
			if(!contains)
				contains = containsLabel(GEOMETRY_RELION_LABELS) &&
				((containsLabel(MDLabel.RLN_ORIENT_ROT)) && ( getColumnMax(MDLabel.RLN_ORIENT_ROT) != 0 )) ;
			return contains;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;	
		
	}

	
	
	public boolean containsLabel(int[] labels) {
		try {
			for (int i = 0; i < labels.length; i++)
				if (containsLabel(labels[i]))
					return true;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}
	
	
	public boolean containsMicrographParticles() {
		try {
			int[] labels = new int[]{MDLabel.MDL_MICROGRAPH, MDLabel.MDL_XCOOR, MDLabel.MDL_YCOOR, MDLabel.MDL_ENABLED};
			for (int i = 0; i < labels.length; i++)
				if (!containsLabel(labels[i]))
					return false;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return true;
	}

	public boolean containsMicrographsInfo() {
		try {
			// Should contain some micrographs labels
			for (int i = 0; i < MICROGRAPH_BASIC_LABELS.length; i++)
				if (!containsLabel(MICROGRAPH_BASIC_LABELS[i]))
					return false;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return true;
	}
	
	public static boolean containsBlock(String file, String block)
	{
		try
		{
			return Arrays.asList(MetaData.getBlocksInMetaDataFile(file)).contains(block);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
	}// function containsBlock

	public static native String label2Str(int label);

	public static native int str2Label(String labelName);

	/** Same of before but handling the exception */
	public static String getLabelName(int label) {
		try {
			return label2Str(label);
		} catch (Exception e) {
			return null;
		}
	}

	public static native String[] getBlocksInMetaDataFile(String filename);

	public native int[] getActiveLabels();

	public static native int getLabelType(int label);
        
        public static Class getClassForType(int type) {
            switch (type) {
		case LABEL_INT:
			return Integer.class;
		case LABEL_BOOL:
			return Boolean.class;
		case LABEL_DOUBLE:
			return Double.class;
		case LABEL_SIZET:
			return Long.class;
		case LABEL_STRING:
		case LABEL_VECTOR_DOUBLE:
		case LABEL_VECTOR_SIZET:
			return String.class;

		}
		return null;
        }

	public static Class getLabelClass(int label) {
		int type = getLabelType(label);
		return getClassForType(type);
	}

	/** Return an String representing the label type */
	public static String getLabelTypeString(int labelType) {
		switch (labelType) {
		case LABEL_STRING:
			return "STRING";
		case LABEL_DOUBLE:
			return "DOUBLE";
		case LABEL_INT:
			return "INT";
		case LABEL_BOOL:
			return "BOOL";
		case LABEL_VECTOR_DOUBLE:
			return "VECTOR(DOUBLE)";
		case LABEL_SIZET:
			return "SIZE_T";
		case LABEL_VECTOR_SIZET:
			return "VECTOR(SIZE_T)";
		}
		return "UNKNOWN";
	}// function getLabelTypeString

	public static native String getLabelComment(int label);

	public static native boolean isTextFile(int label);

	public static native boolean isMetadata(int label);

	public static native boolean isCtfParam(int label);

	public static native boolean isImage(int label);

	public static native boolean isStack(int label);

	public static native boolean isMicrograph(int label);

	public static native boolean isPSD(int label);

	public native boolean isMetadataFile();

	public native void makeAbsPath(int label);

	// get values from metadata
	public native int getValueInt(int label, long objId);

	public native long getValueLong(int label, long objId);

	public native double getValueDouble(int label, long objId);

	/** Return the value of some label as String */
	public native String getValueString(int label, long objId);

	public String getValueString(int label, long objId, boolean fixPaths) {
		String value = getValueString(label, objId);

//		if(value == null)
//			DEBUG.printFormat("label: %d, id: %d", label, objId);
		// Try to fix paths.
		if (fixPaths && filename != null && isPathField(label)) {
			value = fixPath(value);
		}

		return value;
	}

	/**
	 * Return all values of the row as String[]
	 * 
	 * @
	 */
	public String[] getRowValues(long objId) {
		String[] values = new String[activeLabels.length];
		for (int i = 0; i < activeLabels.length; ++i)
			values[i] = getValueString(activeLabels[i], objId);
		return values;
	}
	
//	/** Return a new metadata containing the values of one row */
//	public void getRow(long objId, MetaData mdRow){
//		mdRow.importObjects(this, new long[]{objId});
//	}
//	
//	public void setRow(MetaData mdRow, long objId){
//		int[] labels = mdRow.getActiveLabels();
//		String value;
//		long rowId = mdRow.firstObject();
//		for (int l : labels){
//			value = mdRow.getValueString(l, rowId);
//			setValueString(l, value, objId);
//		}
//	}
	
	public native void getRow(MDRow mdRow, long objId);
	
	public native void setRow(MDRow mdRow, long objId);

	/**
	 * Create a metadata row from another metadata Adding the activeLabels and
	 * adding an object
	 */
	public MetaData getMetaDataRow() {
		MetaData md = null;
		try {
			md = new MetaData();
			md.setColumnFormat(false);
			int[] labels = getActiveLabels();
			for (int l : labels)
				md.addLabel(l);
			md.addObject();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return md;
	}

	public static boolean isPathField(int label) {
		return isTextFile(label) || isMetadata(label) || isCtfParam(label)
				|| isImage(label) || isStack(label) || isMicrograph(label)
				|| isPSD(label);
	}

	// public static boolean isPathField(int label) {
	// return Arrays.binarySearch(PATHS_FIELDS, label) >= 0;
	// }
	public String fixPath(String value) {
		return Filename.findImagePath(value, getBaseDir(), true);
	}

	public String fixPath(String value, String baseDir) {
		return Filename.findImagePath(value, baseDir, false);
	}

	public native boolean getValueBoolean(int label, long objId);

	public String getFilename() {
		return filename != null ? Filename.getFilename(filename) : "";
	}

	public String getBlock() {
		return Filename.getPrefix(filename);
	}

	public String getPath() {
		return filename;
	}

	public String getBaseDir() {
		File f = new File(getFilename());
		f = new File(f.toURI().normalize().getPath());

		return f.getParent();
	}

	public native double[] getStatistics(boolean applyGeo, int label);

	public native double[] getColumnValues(int label);

	// set functions connection with MetaData class in C++
	public boolean setEnabled(boolean value, long objId) {
		
		return setValueInt(MDLabel.MDL_ENABLED, value ? 1 : -1, objId);
		
	}

	public boolean getEnabled(long objId) {
            if(!containsLabel(MDLabel.MDL_ENABLED))
                return true;
            return getValueInt(MDLabel.MDL_ENABLED, objId) > -1;
	}

	public native boolean setValueInt(int label, int value, long objId);

	public native boolean setValueLong(int label, long value, long objId);

	public native boolean setValueDouble(int label, double value, long objId);

        //NOTE: if value contains spaces it should be quoted, otherwise string will be cut!!!
	public native boolean setValueString(int label, String value, long objId);

	public native boolean setValueBoolean(int label, boolean value, long objId);

	/** Obtain all the objects ids in the MetaData */
	public native long[] findObjects();

	/**
	 * Order the metadata by some label. You can order ASCending or DESCending.
	 */
	public native void sort(int sortLabel, boolean ascending);

	/** Import objects from other MetaData */
	public native void importObjects(MetaData from, long ids[]);

	/** Return the id of the first object */
	public native long firstObject();

	/** Add a new object entry and return new id */
	public native long addObject();

	/** Remove an existing object from metadata */
	public native boolean removeObject(long objId);

	/** Remove disabled objects */
	public native void removeDisabled();

	/** Add new column to MetaData */
	public native void addLabel(int label);

	/** Get the average and std images, result is left on input image */
	public native void getStatsImages(ImageGeneric imageAvg,
			ImageGeneric imageStd, boolean applyGeo, boolean wrap, int label);

	public native void writeMdToStack(String filename, boolean applyGeo, boolean wrap, int label);

        
   	public native void getPCAbasis(ImageGeneric basis, int label);

	public native void computeFourierStatistics(MetaData mdIn, int label);

	/**
	 * Union of all elements in two Metadata, duplicating common elements.
	 * Result in calling metadata object, repetition are allowed
	 */
	public native void unionAll(MetaData mdIn);

	/**
	 * Fill all values in a column(label). Different types of fill are: Constant
	 * Linear Random uniform Random gaussian
	 * 
	 * @
	 */
	public native void fillConstant(int label, String value);

	/**
	 * Fill column in a random way
	 * 
	 * @param label
	 *            column to be filled
	 * @param mode
	 *            can be "uniform" or "gaussian"
	 * @param op1
	 *            min for uniform and mean for gaussian
	 * @param op2
	 *            max for uniform and std for gaussian
	 */
	public native void fillRandom(int label, String mode, double op1, double op2);

	/** Fill label starting at some value and with some step */
	public native void fillLinear(int label, double start, double step);

	/** Just for debugging purposes */
	public native void enableDebug();

	/*********** Non-native functions ********************/

	/** Create a metadata and read data from filename */
	public MetaData(String filename) {
		//DEBUG.printFormat("Java: Creating metadata from filename: %s\n", filename);
		//DEBUG.printStackTrace();
            
		create();
		read(filename);
	}

	/**
	 * Will be called by GarbageCollector before destroying. Needed to free C++
	 * object memory.
	 */
	@Override
	protected void finalize() throws Throwable {
		super.finalize();
		destroy();
	}

	/**
	 * Read a metadata from plain textfile.
	 * 
	 * @param file filename from where to read
	 * @param columns expected columns(labels) in file @
	 */
	public native void readPlain(String file, String columns);


        
        
	
	/**
	 * Operate directly in metadata
	 * @param operateStr sql syntanx string of the operation 
	 */
	public native void operate(String operateStr);
	
	/** Returns Max and Min values from a column in metadata
    * These functions can only be used for labels of type double
	*/
	public native double getColumnMin(int column);
	
	public native double getColumnMax(int column);


       

	
	public boolean isCTFMd() {
		try {
			return containsLabel(MDLabel.MDL_PSD_ENHANCED) && containsLabel(MDLabel.MDL_PSD) && containsLabel(MDLabel.MDL_CTF_MODEL);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return false;
	}
        
    public String getCTFFile(long id)
    {
        return getValueString(MDLabel.MDL_CTF_MODEL, id);
    }
    
    public EllipseCTF getEllipseCTF(long id) {
        return getEllipseCTF(id, -1);
    }
    
    public EllipseCTF getEllipseCTF(long id, int D) {
        MetaData md = new MetaData(getCTFFile(id));
        try {
            double Q0, Cs, Ts, kV, downsampleFactor, defU, defV, defAngle;

            Q0 = md.getValueDouble(MDLabel.MDL_CTF_Q0, id);
            Cs = md.getValueDouble(MDLabel.MDL_CTF_CS, id);
            downsampleFactor = md.getValueDouble(MDLabel.MDL_CTF_DOWNSAMPLE_PERFORMED, id);
            Ts = md.getValueDouble(MDLabel.MDL_CTF_SAMPLING_RATE, id) * downsampleFactor;
            kV = md.getValueDouble(MDLabel.MDL_CTF_VOLTAGE, id);

            defU = md.getValueDouble(MDLabel.MDL_CTF_DEFOCUSU, id);
            defV = md.getValueDouble(MDLabel.MDL_CTF_DEFOCUSV, id);
            defAngle = md.getValueDouble(MDLabel.MDL_CTF_DEFOCUS_ANGLE, id);

            return new EllipseCTF(id, Q0, Cs, downsampleFactor, Ts, kV, defU, defV, defAngle, D);
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
            throw new IllegalArgumentException(ex);
        }
        
    }
        
        
    public String getPSDFile(long id) {
        return getValueString(MDLabel.MDL_PSD, id);
    }
    
    public String getPSDEnhanced(long id) {
        return getValueString(MDLabel.MDL_PSD_ENHANCED, id);
    }
    
    public CTFDescription getCTFDescription(long id)
     {
            try {
                return new CTFDescription(getCTFFile(id));
            } catch (Exception ex) {
                Logger.getLogger(MetaData.class.getName()).log(Level.SEVERE, null, ex);
                throw new IllegalArgumentException(ex);
            }
     }
     
	
	
    public static boolean isPlainPos(String file)
    {
        try
        {
            InputStream in = new FileInputStream(file);
            Reader reader = new InputStreamReader(in);
         // buffer for efficiency
            BufferedReader buffer = new BufferedReader(reader);
            String line = buffer.readLine();//dismiss first line
            if(line == null)
                return true;//Empty file
            int count = 0;//small statistics to detect if pos file
            while ((line = buffer.readLine()) != null) {
                for(char c: line.toCharArray())
                    if(!(Character.isDigit(c) || Character.isSpaceChar(c)))
                        return false;


                if(count == 100)
                        break;
                count ++;
            }
        } catch (IOException ex) {
            return false;
        }
        return true;
    }

   
    public boolean hasRecalculateCTF() {
        return !ctfs.isEmpty();
    }

    public void removeCTF(long id) {
        ctfs.remove(id);
    }

    public boolean containsCTF(long id) {
        return ctfs.containsKey(id);
    }

    public void putCTF(long id, EllipseCTF ellipseCTF) {
        ctfs.put(id, ellipseCTF);
    }

    private long getFileTimeStamp(String filename){
        File file = new File(filename);
        return getFileTimeStamp(file);
    }
    private long getFileTimeStamp(File file){
        return file.lastModified();
    }

    public void writeWithException(String xmippFilename) throws IOException{

        String OSFilename = Filename.getFilename(xmippFilename);

        File fileObj = new File(OSFilename);

        boolean exist = fileObj.exists();
        long fileTs = 0;
        if (exist) {

            fileTs = getFileTimeStamp(fileObj);

            if (!fileObj.canWrite()){
                throw new IOException("Can't  write on " + OSFilename + ". Please check you have permissions to write.");
            }
        }

        this.write(xmippFilename);

        if (!exist && !fileObj.exists()) {

            throw new IOException("File: " + OSFilename + " couldn't be created. Check you have permission on the folder.");

        } else if (fileTs == getFileTimeStamp(fileObj)){

            throw new IOException("File: " + OSFilename + " couldn't be written. Something went wrong.");
        }


    }
    
       
       
       
}
