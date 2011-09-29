/***************************************************************************
 *
 * @author: Jesus Cuenca (jcuenca@cnb.csic.es)
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

/**
 * XmippTomo package-level documentation 
 * 
 * XmippTomo is an ImageJ plugin that simplifies 3D-EM & X-Ray tomography tasks with a simple yet powerful workflow
 * 
 * Dependences: Xmipp Java bindings, which means...
 * - XmippJavaInterface.jar inside ImageJ plugins directory
 * - LD_LIBRARY_PATH must include Xmipp's lib directory (since we also need libXmippJavaInterface.so)
 * 
 * - Why "package" is not used?
 * 
 * - Why "package" may be used?
 * Because then the plugin is shown on its own submenu in ImageJ
 */
// TODO: start the plugin from the command line as xmipptomo
// TODO: Project management. Nodes of graph are the results (Series, alignment Parameters, Landmarks, Volumes...)
// Transitions are actions (crop, align...)
// Now all the results are in disk, so it goes like this:
// - open the first (or current) image of the series
// - apply the action to this first image and show the result (like a preview)
// - if the user agrees, then apply the action to all the series in disk
// - show the results in the canvas as each image is saved
// TODO: organize a directory tree that replicates the graph structure (each subdirectory stores the results of 1 action)
// TODO: visualization - add a "Thumbnail" checkbox. When enabled, display the scaled series from memory (faster).
// When disabled, display the original from disk (slower)
// TODO: Project management - action to delete the results file of a node (to save disk space)
// Objects of same type (for example, Landmarks) share the same Color.
// Compact view: display nodes as Type-Sequence # (action). For example:
// - S0
// - S1 (crop S0)
// TODO: remove deprecated methods (once they are not needed anymore)
import ij.*;
import java.io.*;

import ij.plugin.Options;
import ij.plugin.PlugIn;

/**
 * implements PlugIn (ImageJ plugins base interface)
 * 
 * Main class, responsible for plugin initialization.
 * It serves as entry point to the plugin.
 */

// underscore in the name is required for automatic installation in the plugins menu...
public class Xmipp_Tomo implements PlugIn{

	/**
	 * Required because of inherited "serializationability"...
	 */
	private static final long serialVersionUID = -4063711977454855701L;
	
	// enable debuging tests - release versions should have this set to 0
	final static int TESTING = 1;

	// UI elements
	private static TomoWindow tw=null;
	private final static String COMMAND_OPTION_INPUT = "i";

	
	private int lastWindowId=0;
	
	/* entry point for PlugIn implementors
	 */
	public void run(String arg){
		createTomoWindow();
		// debug("library path: " + System.getProperty("java.library.path"));
		
		// IJ.debugMode = true;
		
        if (IJ.isMacro() && Macro.getOptions() != null && !Macro.getOptions().trim().isEmpty()) { // From macro.
            // "string" is used when called from another plugin or installed command.
            // "Macro.getOptions()" used when called from a run("command", arg) macro function.
            String file = processArgs(Macro.getOptions().trim());
            if (file != null)
            	tw.getController().loadEM(file);
        }
		
		if(tw != null)
			tw.display();
	}
	
/**
 * Basic parsing of args: "-i file"
 * @param args
 * @return
 */
	 private String processArgs(String args) {
	        String argsList[] = args.split(" ");

	        String ret=null;
	        
	        if("-i".equals(argsList[0]))
	        	ret=argsList[1];
	        return ret;
	    }

	

	/** Run cmdline in a shell and show standard output & error via debug()
	 * @param cmdline
	 * @return the exit value of cmdline (@see Process.waitFor())
	 */
	public static ExitValue exec(String cmdline,boolean readStderr){
		// execution details may change with each OS...
		//String osName = System.getProperty("os.name" );
		
		ExitValue exitValue=ExitValue.OK;
		Process proc=null;
		
		Runtime rt = Runtime.getRuntime();
		
		// Xmipp_Tomo.debug(System.getenv("LD_LIBRARY_PATH"));
		Logger.debug(cmdline);
		
		try{
			proc = rt.exec(cmdline);
		}catch (IOException ex){
			Logger.debug(ex.toString());
			// one improvement would be to extract the error code from the exception and return the exitvalue accordingly
			return ExitValue.PROGRAM_NOT_FOUND;
		}
		
		// Prepare buffered readers from inputstreams (stderr, stdout) ...
		InputStream stderr=proc.getErrorStream();
		InputStreamReader stderr_isr = new InputStreamReader(stderr);
        BufferedReader stderr_br = new BufferedReader(stderr_isr,10);

		InputStream stdout=proc.getInputStream();
		InputStreamReader stdout_isr = new InputStreamReader(stdout);
		// try reading small chunks of output
        BufferedReader stdout_br = new BufferedReader(stdout_isr,10);
		
		String line=null;
		
		try{
			while ( (line = stdout_br.readLine()) != null)
	            Logger.debug(line);    
	        
			if(readStderr)
				while ( (line = stderr_br.readLine()) != null)
		            Logger.debug(line);    
			
        } catch (IOException ex){
            ex.printStackTrace();  
        }
        
        try{
        	int ev = proc.waitFor();
        	// convert from int ev to ExitValue
        	switch(ev){
        		case 0:
        			exitValue = ExitValue.OK;
        			break;
        		case 127:
        			// recompile the program
        			exitValue = ExitValue.EXTERNAL_PROGRAM_BROKEN;
        			break;
        		default:
        			Logger.debug(String.valueOf(ev));
        			exitValue = ExitValue.ERROR;
        	}
        	
        }catch (java.lang.InterruptedException ex){
        	exitValue=ExitValue.ERROR;
        }
        return exitValue;
	} // exec end
	


    /********************* Trace methods for debugging ****************************/
    
    
	private int getNextWindowId(){
		lastWindowId++;
		return lastWindowId;
	}
	
	private void createTomoWindow(){
			tw= new TomoWindow(getNextWindowId());

	}
		
	
	// ij.IJ class offers a lot of useful methods...
	private void ijExamples(){
	
		IJ.showStatus("Xmipp_Tomo started.");
		IJ.showProgress(0.0);
		String name = IJ.getString("Please enter your name: ","I.J. User");
		IJ.showProgress(0.5);
		IJ.write("Starting sample plugin Red And Blue ... ");
		IJ.runPlugIn("Red_And_Blue","");
		IJ.showProgress(1.0);
		IJ.showMessage("Finished.",name+", thank you for running this plugin");
		
	}


	
}

