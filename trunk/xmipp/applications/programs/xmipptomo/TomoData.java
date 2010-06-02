package xmipptomo;
import ij.IJ;
import ij.ImagePlus;
import ij.io.OpenDialog;

import java.util.Collections;

import javax.swing.text.BadLocationException;
import javax.swing.text.Document;

/**
 * 
 */

/**
 * @author jcuenca
 *
 */
public class TomoData {
	
	// maybe it's better to move currentSlice to the viewer - in case different views can show different slices
	private int currentSlice=0;

    private ImagePlus imp=null;
	private java.util.List <TomoInfo> ti=Collections.emptyList();
	
	private Document tiltTextModel=null;
	
	public void setTiltModel(Document model){
		if(tiltTextModel==null){
			// initialize model with default value
			try{
				model.insertString(0, String.valueOf(getInitialTilt()), null);
			}catch (Exception ex){}
		}
		
		tiltTextModel=model;
	}
	
	public Document getTiltModel(){
		return tiltTextModel;
	}
	
	public int getCurrentSlice() {
		return currentSlice;
	}

	public void setCurrentSlice(int currentSlice) {
		this.currentSlice = currentSlice;
		// update all things depending on current slice
		getImage().setSlice(getCurrentSlice());
		setTiltText(getCurrentTilt());
	}
	
	/**
	 * @return the ti
	 */
	public java.util.List<TomoInfo> getTomoInfoList() {
		return ti;
	}

	/**
	 * @param ti the ti to set
	 */
	public void setTomoInfoList(java.util.List<TomoInfo> ti) {
		this.ti = ti;
	}

	public ImagePlus getImage(){
		return imp;
	}
	
	private void setImage(ImagePlus i){
		imp=i;
	}
	
	public TomoInfo getCurrentTomoInfo(){
		// if(getTomoInfoList().size()>0)
		return getTomoInfoList().get(getCurrentSlice());
	}
	
	public float getCurrentTilt(){
		return getCurrentTomoInfo().getTilt();
	}
	
	public void setCurrentTilt(float t){
		getCurrentTomoInfo().setTilt(t);
	}
	
	// there may be no tilts defined - change tilt to Float (so it can handle nulls)
	private float getInitialTilt(){
		return getCurrentTilt();
	}
	
	public int getNSlices(){
		return getImage().getNSlices();
	}
	
	// right now return only grayscale value as double
	// more complete method at ImagePlus.getValueAsString()
	public double getPixelValue(int x,int y){
		int[] v = getImage().getPixel(x, y);
		// 32bit images
		return Float.intBitsToFloat(v[0]);
		//return getImage().getCalibration().getCValue(v[0]);
	}
	
	public void import_data(String path) throws java.io.IOException{
			// Get ImagePlus with the help of the appropiate I/O plugin
			if (path.endsWith(".mrc")) {       
				imp = (ImagePlus) IJ.runPlugIn("MRC_Reader", path);
				if (imp != null && imp.getWidth() != 0) {
					// mrc_reader does not save Fileinfo inside the ImagePlus
					setTomoInfoList(TomoInfo.readMRC(path));
					// debug(String.valueOf(ti.elementAt(3).getTilt()));
				}
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
			
			setImage(imp);

	}
	
	private void setDocumentText(Document d,String t){
		try{
			d.remove(0,d.getLength());
			d.insertString(0, t, null);
		}catch (BadLocationException e){}
	}
	
	private void setTiltText(float tilt){
		setDocumentText(getTiltModel(),String.valueOf(tilt));	
	}
}
