package xmipptomo;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.io.FileOpener;
import ij.io.ImageReader;
import ij.measure.Calibration;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.image.ColorModel;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.Properties;

/**
 * @author jcuenca
 * Collection of methods for handling I/O
 * - make all methods static?
 */
public class TiltSeriesOpener {

	public static int MRC_HEADER_SIZE = 1024;
	
	private class ByteHeader{
		private byte [] bytes;
		private boolean littleEndian=false;
	
		
		public ByteHeader(int size){
			bytes = new byte[size];
		}
		
		public void read(String path) throws java.io.IOException{
			RandomAccessFile f = new RandomAccessFile(path, "r");

	        for (int i = 0; i < bytes.length; i++) {
	            bytes[i] = f.readByte();
	        }
	        f.close();
		}
		
	    public int getInt(int index) {
	        byte b1 = bytes[index];
	        byte b2 = bytes[index + 1];
	        byte b3 = bytes[index + 2];
	        byte b4 = bytes[index + 3];
	        if (isLittleEndian()) {
	            return ((((b4 & 0xff) << 24) | ((b3 & 0xff) << 16) | ((b2 & 0xff) << 8) | (b1 & 0xff)));
	        }
	        return ((((b1 & 0xff) << 24) | ((b2 & 0xff) << 16) | ((b3 & 0xff) << 8) | (b4 & 0xff)));
	    }

		/**
		 * @return the littleEndian
		 */
		public boolean isLittleEndian() {
			return littleEndian;
		}

		/**
		 * @param littleEndian the littleEndian to set
		 */
		public void setLittleEndian(boolean littleEndian) {
			this.littleEndian = littleEndian;
		}
	}
	
	public ImagePlus read(String path,TomoData model) throws java.io.IOException{
		ImagePlus imp=null;
		if ((path == null) || (path.equals(""))) 
			throw new IOException("Empty path");
		
			if (path.endsWith(".mrc")) {       
				// imp = (ImagePlus) IJ.runPlugIn("MRC_Reader", path);
				// if (imp != null && imp.getWidth() != 0) {
					// mrc_reader does not save Fileinfo inside the ImagePlus
				readMRC(path,model);
					// debug(String.valueOf(ti.elementAt(3).getTilt()));
				// }
			}else if (path.toLowerCase().endsWith(".spi") || path.toLowerCase().endsWith(".xmp")|| path.toLowerCase().endsWith(".vol")) {
				imp = (ImagePlus) IJ.runPlugIn("Spider_Reader", path);
				if (imp == null) {
					//width = PLUGIN_NOT_FOUND;
				}
				if (imp != null && imp.getWidth() == 0) {
					imp = null;
				}
			}else if (path.endsWith(".sel")) {       
				imp = (ImagePlus) IJ.runPlugIn("Sel_Reader", path);
				if (imp == null) {
					//width = PLUGIN_NOT_FOUND;
				}
				if (imp != null && imp.getWidth() == 0) {
					imp = null;
				}
			}
			
			return imp;
	}
	
	public void readMRC(String path,TomoData model) throws java.io.IOException{

		// Parse header
		ByteHeader header=new ByteHeader(MRC_HEADER_SIZE);
		header.read(path);        
        FileInfo fi=prepareFileInfo(path, header);
        
        
		// read each slice, resize/scale and add it to the stack
        FileOpener fo=new FileOpener(fi);
        
        ColorModel cm=fo.createColorModel(fi);
        
		ImageStack stack = new ImageStack(fi.width, fi.height, cm),stackResized=new ImageStack(model.getDefaultWidth(), model.getDefaultHeight());
		long skip = fi.longOffset>0?fi.longOffset:fi.offset;
		Object pixels;
		try {
			ImageReader reader = new ImageReader(fi);
			InputStream is =  fo.createInputStream(fi);
			if (is==null) throw new IOException("readMRC - Null input stream") ;
			// IJ.resetEscape();
			model.setNumberOfProjections(fi.nImages);

			// Save original size to model
			model.setHeight(fi.height);
			model.setWidth(fi.width);
			
			// Adjust file info to new size
			fi.width= model.getDefaultWidth();
			fi.height = model.getDefaultHeight();
			
			// load images
			for (int i=1; i<=fi.nImages; i++) {
				/* if (!silentMode)
					IJ.showStatus("Reading: " + i + "/" + fi.nImages);
				if (IJ.escapePressed()) {
					IJ.beep();
					IJ.showProgress(1.0);
					silentMode = false;
					return null;
				} */
	
				pixels = reader.readPixels(is, skip);
				if (pixels==null) break;
				
				// resize/scale - use an aux image processor for all projections
		        // create appropiate image processor
				ImageProcessor ip=null,ipresized=null;
				
				// create the proper ImageProcessor initialized with the pixels
		        switch (fi.fileType){
		        	case FileInfo.GRAY8:
		        		ip = new ByteProcessor(fi.width, fi.height, (byte[])pixels, cm);
		        		break;
		        	case FileInfo.GRAY16_UNSIGNED:
		        		ip = new ShortProcessor(fi.width, fi.height, (short[])pixels, cm);
		        		break;
		        	case FileInfo.GRAY32_FLOAT:	
		        		ip = new FloatProcessor(fi.width, fi.height, (float[])pixels, cm);
		        		break;		        		
		        }
				
		        ipresized = ip.resize(model.getDefaultWidth(), model.getDefaultHeight());
		        
				stack.addSlice(null, ip);
				stackResized.addSlice(null,ipresized);
				skip = fi.gapBetweenImages;
				
				
				/* if (!silentMode)
					IJ.showProgress(i, fi.nImages); */
			}
			is.close();
		}catch(OutOfMemoryError e) {
			IJ.outOfMemory(fi.fileName);
			stack.trim();
		}
		
		// normal or resized?
		// ImagePlus imp = new ImagePlus(fi.fileName, stack);
		ImagePlus imp = new ImagePlus(fi.fileName, stackResized);
		if (fi.info!=null)
			imp.setProperty("Info", fi.info);
		
		imp.setFileInfo(fi);
		
		/* if (!silentMode) IJ.showProgress(1.0); */
		if (stack.getSize()==0)
			throw new IOException("Empty stack");
		if (fi.sliceLabels!=null && fi.sliceLabels.length<=stack.getSize()) {
			for (int i=0; i<fi.sliceLabels.length; i++)
				stack.setSliceLabel(fi.sliceLabels[i], i+1);
		}
		
		setCalibration(imp,fi,fo);
		ImageProcessor ip = imp.getProcessor();
		if (ip.getMin()==ip.getMax())  // find stack min and max if first slice is blank
			setStackDisplayRange(imp);

		// read tilt angles
		String tltFilePath=path.replace(".mrc", ".tlt");
		readTlt(tltFilePath,model);
		
		model.setImage(imp);
		model.setFile(path);
	}
	
	public FileInfo prepareFileInfo(String path,ByteHeader header) throws java.io.IOException {
		// get MRC DM3 specific file information into FileInfo
		FileInfo fi = new FileInfo();
        fi.fileFormat = FileInfo.RAW;
        File dest = new File(path);
        fi.fileName =dest.getName();
        fi.directory = dest.getParent();

        int mode = header.getInt(3 * 4);
        if (mode == 0) {
            fi.fileType = FileInfo.GRAY8;
            fi.width = header.getInt(0);
            fi.height = header.getInt(1 * 4);
            fi.nImages =header.getInt(2 * 4);
            if (fi.width < 0 || fi.height < 0 || fi.nImages < 0) {
                header.setLittleEndian(true);
            }
        } else if (mode == 1) {
            fi.fileType = FileInfo.GRAY16_UNSIGNED;
        } else if (mode == 2) {
            fi.fileType = FileInfo.GRAY32_FLOAT;
        } else {
            //it is not a big endian data
            header.setLittleEndian(true);
            mode = header.getInt(3 * 4);
            if (mode == 1) {
                fi.fileType = FileInfo.GRAY16_UNSIGNED;
            } else if (mode == 2) {
                fi.fileType = FileInfo.GRAY32_FLOAT;
            } else {
                throw new IOException("Unimplemented ImageData mode=" + mode + " in MRC file.");
            }
        }

        fi.intelByteOrder = header.isLittleEndian();
        // duplicated when mode = 0 ?
        fi.width = header.getInt(0);
        fi.height = header.getInt(1 * 4);
        fi.nImages = header.getInt(2 * 4);
        
        int extrabytes = header.getInt(23 * 4);
        fi.offset = 1024 + extrabytes;
        /* if (extrabytes == 131072) {
            IJ.write("FEI");
           
            FEI = true;
        }*/
        return fi;
	}
	
	// .tlt syntax: one float angle per line, stored as text
	public void readTlt(String tltFilePath,TomoData model) throws java.io.IOException{
		// String mrcFilePath=imp.getFileInfo().directory + imp.getFileInfo().fileName;

		// IJ.write(tltFilePath);
		BufferedReader brin=new BufferedReader(new FileReader(tltFilePath));
		String line=null;
		while ( (line=brin.readLine()) != null){
			model.addTiltAngle(new Float(line));
		}
	}
	
	
	// code from FileOpener - due to access modifiers restrictions...
	void setCalibration(ImagePlus imp,FileInfo fi,FileOpener fo) {
		
		
		if (fi.fileType==FileInfo.GRAY16_SIGNED) {
			if (IJ.debugMode) IJ.log("16-bit signed");
			double[] coeff = new double[2];
			coeff[0] = -32768.0;
			coeff[1] = 1.0;
 			imp.getLocalCalibration().setFunction(Calibration.STRAIGHT_LINE, coeff, "gray value");
		}
		
		Properties props = fo.decodeDescriptionString(fi);
		Calibration cal = imp.getCalibration();
		boolean calibrated = false;
		if (fi.pixelWidth>0.0 && fi.unit!=null) {
			cal.pixelWidth = fi.pixelWidth;
			cal.pixelHeight = fi.pixelHeight;
			cal.pixelDepth = fi.pixelDepth;
			cal.setUnit(fi.unit);
			calibrated = true;
		}
		
		if (fi.valueUnit!=null) {
			int f = fi.calibrationFunction;
			if ((f>=Calibration.STRAIGHT_LINE && f<=Calibration.RODBARD2 && fi.coefficients!=null)
			|| f==Calibration.UNCALIBRATED_OD) {
				boolean zeroClip = props!=null && props.getProperty("zeroclip", "false").equals("true");	
				cal.setFunction(f, fi.coefficients, fi.valueUnit, zeroClip);
				calibrated = true;
			}
		}
		
		/* if (calibrated)
			checkForCalibrationConflict(imp, cal); */
		
		if (fi.frameInterval!=0.0)
			cal.frameInterval = fi.frameInterval;
		
		if (props==null)
			return;
					
		cal.xOrigin = getDouble(props,"xorigin");
		cal.yOrigin = getDouble(props,"yorigin");
		cal.zOrigin = getDouble(props,"zorigin");
		cal.info = props.getProperty("info");		
				
		cal.fps = getDouble(props,"fps");
		cal.loop = getBoolean(props, "loop");
		cal.frameInterval = getDouble(props,"finterval");
		cal.setTimeUnit(props.getProperty("tunit", "sec"));

		double displayMin = getDouble(props,"min");
		double displayMax = getDouble(props,"max");
		if (!(displayMin==0.0&&displayMax==0.0)) {
			int type = imp.getType();
			ImageProcessor ip = imp.getProcessor();
			if (type==ImagePlus.GRAY8 || type==ImagePlus.COLOR_256)
				ip.setMinAndMax(displayMin, displayMax);
			else if (type==ImagePlus.GRAY16 || type==ImagePlus.GRAY32) {
				if (ip.getMin()!=displayMin || ip.getMax()!=displayMax)
					ip.setMinAndMax(displayMin, displayMax);
			}
		}
		
		int stackSize = imp.getStackSize();
		if (stackSize>1) {
			int channels = (int)getDouble(props,"channels");
			int slices = (int)getDouble(props,"slices");
			int frames = (int)getDouble(props,"frames");
			if (channels==0) channels = 1;
			if (slices==0) slices = 1;
			if (frames==0) frames = 1;
			//IJ.log("setCalibration: "+channels+"  "+slices+"  "+frames);
			if (channels*slices*frames==stackSize) {
				imp.setDimensions(channels, slices, frames);
				if (getBoolean(props, "hyperstack"))
					imp.setOpenAsHyperStack(true);
			}
		}
	}
	
	private Double getNumber(Properties props, String key) {
		String s = props.getProperty(key);
		if (s!=null) {
			try {
				return Double.valueOf(s);
			} catch (NumberFormatException e) {}
		}	
		return null;
	}
	
	private double getDouble(Properties props, String key) {
		Double n = getNumber(props, key);
		return n!=null?n.doubleValue():0.0;
	}
	
	private boolean getBoolean(Properties props, String key) {
		String s = props.getProperty(key);
		return s!=null&&s.equals("true")?true:false;
	}
	
	void setStackDisplayRange(ImagePlus imp) {
		ImageStack stack = imp.getStack();
		double min = Double.MAX_VALUE;
		double max = -Double.MAX_VALUE;
		int n = stack.getSize();
		for (int i=1; i<=n; i++) {
			/* if (!silentMode)
				IJ.showStatus("Calculating stack min and max: "+i+"/"+n); */
			ImageProcessor ip = stack.getProcessor(i);
			ip.resetMinAndMax();
			if (ip.getMin()<min)
				min = ip.getMin();
			if (ip.getMax()>max)
				max = ip.getMax();
		}
		imp.getProcessor().setMinAndMax(min, max);
		imp.updateAndDraw();
	}
}
