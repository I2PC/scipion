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
 * - Why "package" is not used?
 * 
 * - Why "package" may be used?
 * Because then the plugin is shown on its own submenu in ImageJ
 */

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
// TODO: Project management - workflow import/export to disk (SQLite format?)
// TODO: Project management - action to delete the results file of a node (to save disk space)
// Objects of same type (for example, Landmarks) share the same Color.
// Compact view: display nodes as Type-Sequence # (action). For example:
// - S0
// - S1 (crop S0)
// TODO: load ImageDoubles into a cache on demand (cache miss)
/**
 * TODO: changes in MVC design
 * - WindowModel: not needed (no relevant window data other than menus)
 * - StackModel: knows how to get the current projection to display,
 * metadata (tilt), # of projections... No need to store data in the model
 * - StackView: display projection image, tilt, and player controls
 * - WorkflowModel: stores all the user operations
 * - WorkflowOperation: input and output files, command, parameters, status (progress)...
 * - WorkflowView: displays the workflow + buttons (load,save,discard)
 * - ImageDoubleCache: get("003@F2")
 */

import ij.*;
import java.io.*;
import java.util.Calendar;
import java.util.Enumeration;
import java.util.LinkedList;
import java.util.List;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeNode;

import ij.plugin.PlugIn;

/**
 * implements PlugIn (ImageJ plugins base interface)
 * 
 * Main class, responsible for plugin initialization
 */

// underscore in the name is required for automatic installation in the plugins menu...
public class Xmipp_Tomo implements PlugIn{

	/**
	 * Required because of inherited "serializationability"...
	 */
	private static final long serialVersionUID = -4063711977454855701L;
	
	public static enum ExitValues {
		OK(0),ERROR(1),YES(2),NO(3),CANCEL(4),RUNTIME_ERROR(5),PROGRAM_NOT_FOUND(6),
		EXTERNAL_PROGRAM_BROKEN(7);
		private final int value;
		ExitValues(int err) { value=err;}
		public int value(){return value;}
	};
	
	// enable debuging tests - release versions should have this set to 0
	final static int TESTING = 1;

	// UI elements
	private static TomoWindow tw=null;
	
	private static DefaultMutableTreeNode workflow;
	
	private int lastWindowId=0;
	
	/* entry point for PlugIn implementors
	 */
	public void run(String arg){
		setWorkflow(null);
		createTomoWindow();
		// debug("library path: " + System.getProperty("java.library.path"));
		
		// IJ.debugMode = true;
		
		if(tw != null)
			tw.display();
	}
	
	public static DefaultMutableTreeNode addUserAction(DefaultMutableTreeNode last, UserAction newAction){
		if(last == null){
			Xmipp_Tomo.debug("Xmipp_Tomo.addUserAction - last is null");
		}else
			last.add(new DefaultMutableTreeNode(newAction));
		return last;
	}
	
	/** 
	 * @param last
	 * @return the section of the workflow which last action is the one passed as a parameter
	 */
	public static List<UserAction> getWorkflow(DefaultMutableTreeNode last){
		LinkedList <UserAction> res=new LinkedList<UserAction>();
		res.push((UserAction)last.getUserObject());
		DefaultMutableTreeNode current_node=last;
		while(current_node != getWorkflow().getRoot()){
			DefaultMutableTreeNode parent_node= (DefaultMutableTreeNode)current_node.getParent();
			if (parent_node == null)
				debug("Null parent");
			else{
				current_node=parent_node;
				res.push((UserAction)current_node.getUserObject());
			}
		}
		
		return res;
	}
	
	// workflow is a singleton
	public static DefaultMutableTreeNode getWorkflow(){
		if(workflow == null){
			UserAction workflowRoot = new UserAction(UserAction.ROOT_WINDOWID,"Project");
			setWorkflow(new DefaultMutableTreeNode(workflowRoot));
		}
		return workflow;
	}
	
	private static UserAction getWorkflowRoot(){
		return (UserAction)((DefaultMutableTreeNode)getWorkflow().getRoot()).getUserObject();
	}
	
	public static void printWorkflow(){
		Enumeration e = getWorkflow().breadthFirstEnumeration();
		while(e.hasMoreElements())
			debug((e.nextElement()).toString());
	}

	
	/** Run cmdline in a shell and show standard output & error via debug()
	 * @param cmdline
	 * @return the exit value of cmdline (@see Process.waitFor())
	 */
	public static ExitValues exec(String cmdline,boolean readStderr){
		// execution details may change with each OS...
		//String osName = System.getProperty("os.name" );
		
		ExitValues exitValue=ExitValues.OK;
		Process proc=null;
		
		Runtime rt = Runtime.getRuntime();
		
		// Xmipp_Tomo.debug(System.getenv("LD_LIBRARY_PATH"));
		Xmipp_Tomo.debug(cmdline);
		
		try{
			proc = rt.exec(cmdline);
		}catch (IOException ex){
			debug(ex.toString());
			// one improvement would be to extract the error code from the exception and return the exitvalue accordingly
			return ExitValues.PROGRAM_NOT_FOUND;
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
	            debug(line);    
	        
			if(readStderr)
				while ( (line = stderr_br.readLine()) != null)
		            debug(line);    
			
        } catch (IOException ex){
            ex.printStackTrace();  
        }
        
        try{
        	int ev = proc.waitFor();
        	// convert from int ev to ExitValue
        	switch(ev){
        		case 0:
        			exitValue = ExitValues.OK;
        			break;
        		case 127:
        			// recompile the program
        			exitValue = ExitValues.EXTERNAL_PROGRAM_BROKEN;
        			break;
        		default:
        			debug(String.valueOf(ev));
        			exitValue = ExitValues.ERROR;
        	}
        	
        }catch (java.lang.InterruptedException ex){
        	exitValue=ExitValues.ERROR;
        }
        return exitValue;
	} // exec end
	


    /********************* Trace methods for debugging ****************************/
    
    
	/** print s with a timestamp, using IJ Results window (@see IJ.write)
	 * @param s
	 */
	public static void debug(String s){
		Calendar calendar = Calendar.getInstance();
		java.util.Date now = calendar.getTime();
		String output="" + now.getTime() + " > " + s;
		if(TESTING == 1){
			/*if(tw != null)
				IJ.write(output);
			else */
				System.err.println(output);
		}
	}

	/** same as Debug, plus ex stack trace
	 * @param s
	 * @param ex Exception
	 */
	public static void debug(String s, Exception ex){
		debug(s + ex.toString());
		ex.printStackTrace();
	}
	
	private int getNextWindowId(){
		lastWindowId++;
		return lastWindowId;
	}
	
	private void createTomoWindow(){
			tw= new TomoWindow(getNextWindowId());
			UserAction newWindow= new UserAction(UserAction.ROOT_WINDOWID,"New Window");
			addUserAction(getWorkflow(),newWindow);
			tw.setFirstAction(getWorkflow());
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
	
	public void testWorkflow(){
		/* // create tree
		DefaultMutableTreeNode root=new DefaultMutableTreeNode("1"),n2=new DefaultMutableTreeNode("2"),
		n3=new DefaultMutableTreeNode("3"),n4=new DefaultMutableTreeNode("4"),n5=new DefaultMutableTreeNode("5");
		// get one lineage
		addUserAction(root, n2);
		addUserAction(root, n4);
		addUserAction(n2, n3);
		addUserAction(n4, n5);
		printWorkflow();
		
		List <UserAction> l=getWorkflow(n5);
		// should print root - n4 - n5
		for(UserAction ua:l)
			debug(ua.toString());*/
	}
	
	public static void main(String[] args) {
		Xmipp_Tomo xt= new Xmipp_Tomo();
		xt.testWorkflow();
	}

	/**
	 * @param workflow the workflow to set
	 */
	private static void setWorkflow(DefaultMutableTreeNode workflow) {
		Xmipp_Tomo.workflow = workflow;
	}
	
}

