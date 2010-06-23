package xmipptomo;
import ij.IJ;
import java.util.*;
import java.io.*;

/**
 * Specific features of each projection
 */

/**
 * @author jcuenca
 *
 */
public class TomoInfo{
	private float tilt;

	// Reader plugins don't read tilt info, so custom methods are required
		
	// MRC: tilt info is in .tlt file
	// .tlt syntax: one float angle per line, stored as text
	public static List <TomoInfo> readMRC(String mrcFilePath) throws java.io.IOException{
		// by now use a vector to store the info
		List <TomoInfo> res=new Vector<TomoInfo> (); //imp.getStackSize());
		
		// String mrcFilePath=imp.getFileInfo().directory + imp.getFileInfo().fileName;
		String tltFilePath=mrcFilePath.replace(".mrc", ".tlt");
		IJ.write(tltFilePath);
		BufferedReader brin=new BufferedReader(new FileReader(tltFilePath));
		String line=null;
		while ( (line=brin.readLine()) != null){
			TomoInfo e=new TomoInfo();
			e.setTilt(new Float(line).floatValue());
			res.add(e);
		}
			
		
		return res;
	}
	
	// Spider: tilt info is in the header
	// Xmipp: tilt info is in Metadata
	
	public float getTilt(){ return tilt;}
	public void setTilt(float t){tilt=t;}
	
	public static float getMinTilt(List <TomoInfo> ti){
		return ti.get(0).getTilt();
	}
	
	public static float getMaxTilt(List <TomoInfo> ti){
		return ti.get(ti.size()).getTilt();
	}
	
}
