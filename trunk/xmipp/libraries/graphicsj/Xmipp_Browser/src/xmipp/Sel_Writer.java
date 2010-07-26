package xmipp;

import ij.ImagePlus;
import ij.IJ;
import ij.plugin.filter.PlugInFilter;
import ij.io.SaveDialog;
import ij.process.ImageProcessor;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
 *
 * @author MESSAOUDI Cedric
 *         Date: 5 juin 2009
 *         Time: 12:49:06
 */
public class Sel_Writer implements PlugInFilter {
	ImagePlus myimp;
	String directory;
	String name;

	/**
	 *  Main processing method for the SelWriter_ object
	 *
	 *@param  arg  Description of the Parameter
	 *@param  imp  Description of the Parameter
	 *@return      Description of the Return Value
	 */
	public int setup(String arg, ImagePlus imp) {
		myimp = imp;
		if(!(arg==null||arg.equals(""))){
			File dest = new File(arg);
			directory = dest.getParent()+ System.getProperty("file.separator");
			name = dest.getName();
			System.out.println("saving as asked "+arg);
			System.out.println("saving as sel "+directory+name);
		}else{
			name=null;
		}
		return DOES_32 + NO_CHANGES;
	}


	/**
	 *  Main processing method for the SelWriter_ object
	 *
	 *@param  ip  Description of the Parameter
	 */
	public void run(ImageProcessor ip) {
		if(name==null){
			SaveDialog sd = new SaveDialog("save as...", myimp.getTitle(), ".sel");
			directory=sd.getDirectory();
			name=sd.getFileName();
		}
		if (name == null) {
			return;
		}
		IJ.showStatus("saving: " + directory + name);
		//System.out.println("saving: " + directory + name);
		try{
			BufferedWriter selfile=new BufferedWriter(new FileWriter(directory+name));
			for(int i=1;i<=myimp.getNSlices();i++){
				myimp.setSlice(i);
				ImagePlus tmp=new ImagePlus(""+i,myimp.getProcessor());
				String filename=name.substring(0,name.length()-4)+i+".xmp";
				//System.out.println("saving: " + directory + filename);
				IJ.runPlugIn(tmp,"Spider_Writer",directory+filename);
				selfile.write(filename+" 1\n");
			}
			selfile.close();


		}catch (Exception e){
			System.out.println("error with sel file : "+e);
		}


	}
}
