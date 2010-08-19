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

// some recommend that plugin classes must be in the default package...
// here we use it has its own submenu
package xmipptomo;


import ij.*;
import ij.io.*;

import java.awt.Dimension;
import java.io.*;
import java.security.acl.LastOwnerException;
import java.util.Calendar;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import ij.plugin.PlugIn;

/**
 * @author jcuenca
 * implements PlugIn (ImageJ plugins base interface)
 * 
 * Main class, responsible for plugin initialization
 */

// underscore in the name is required for automatic installation in the plugins menu...
// Requirements: http://u759.curie.u-psud.fr/compteur/download.php?Fichier=software/update/20090928/U759_InputOutput.jar
public class Xmipp_Tomo implements PlugIn{

	/**
	 * Required because of inherited "serializationability"...
	 */
	private static final long serialVersionUID = -4063711977454855701L;
	
	// if image is bigger than this threshold, resize it to this size
	public static Dimension resizeThreshold = new Dimension(400,400);
	
	private static enum CmdExitValues {
		OK(0),ERROR(1);
		private final int value;
		CmdExitValues(int err) { value=err;}
		public int value(){return value;}
	};
	
	// enable debuging tests - release versions should have this set to 0
	final static int TESTING = 1;

	// UI elements
	private static TomoWindow tw=null;
	
	private static Tree<UserAction> workflow;
	
	private int lastWindowId=0;
	
	/*

	public Xmipp_Tomo() {
		// super("Xmipp_Tomo");
		
		// dataModel=new TomoData();
		
		/*setLayout(new FlowLayout());
		panel = new Panel();
		panel.setLayout(new GridLayout(4, 4, 5, 5));
		panel.add(new Label("Xmipp Tomo"));
		addButton(ButtonLabels.CMD_INFO.label());
		addButton(ButtonLabels.CMD_PREPROC.label());
		addButton(ButtonLabels.CMD_ALIGN.label());
		
		add(panel);
		
		pack();
		GUI.center(this);
		// show is deprecated...
		show(); 
	} */
	
	/* entry point for PlugIn implementors
	 */
	public void run(String arg){
		createTomoWindow();
		
		// IJ.debugMode = true;
		
		if(tw != null)
			tw.display();
	}
	
	public static void addUserAction(UserAction last, UserAction newAction){
		if(last == null){
			Xmipp_Tomo.debug("Xmipp_Tomo.addUserAction - last is null");
			return;
		}
		getWorkflow().addLeaf(last,newAction);
		
	}
	
	public static List<UserAction> getWorkflow(UserAction last){
		LinkedList <UserAction> res=new LinkedList<UserAction>();
		res.push(last);
		UserAction current=last;
		Tree<UserAction> current_node=getWorkflow().getTree(last);
		while(current != getWorkflow().getHead()){
			Tree<UserAction> parent_node= current_node.getParent();
			if (parent_node == null)
				debug("Null parent");
			else{
				current=parent_node.getHead();
				current_node = parent_node;
				res.push(current);
			}
		}
		
		return res;
	}
	
	// workflow is a singleton
	private static Tree<UserAction> getWorkflow(){
		if(workflow == null){
			UserAction workflowRoot = new UserAction(UserAction.ROOT_WINDOWID,"Plugin start");
			workflow = new Tree<UserAction>(workflowRoot);
		}
		return workflow;
	}
	
	private static UserAction getWorkflowRoot(){
		return getWorkflow().getHead();
	}
	
	public static void printWorkflow(){
		debug(getWorkflow().toString());
	}
	

	
	/** Add button to main panel
	 * @param label
	 *
	void addButton(String label) {
		/* Button b = new Button(label);
		b.addActionListener(this);
		b.addKeyListener(IJ.getInstance());
		panel.add(b); 
	}*/

	
	/** Run cmdline in a shell and show standard output & error via debug()
	 * @param cmdline
	 * @return the exit value of cmdline (@see Process.waitFor())
	 */
	int exec(String cmdline){
		// execution details may change with each OS...
		//String osName = System.getProperty("os.name" );
		
		int exitValue=CmdExitValues.OK.value();
		Process proc=null;
		
		Runtime rt = Runtime.getRuntime();
		
		try{
			proc = rt.exec(cmdline);
		}catch (IOException ex){
			return -1;
		}
		
		// Prepare buffered readers from inputstreams (stderr, stdout) ...
		InputStream stderr=proc.getErrorStream();
		InputStreamReader stderr_isr = new InputStreamReader(stderr);
        BufferedReader stderr_br = new BufferedReader(stderr_isr);

		InputStream stdout=proc.getInputStream();
		InputStreamReader stdout_isr = new InputStreamReader(stdout);
        BufferedReader stdout_br = new BufferedReader(stdout_isr);
		
		String line=null;
		
		try{
			// print stdout
			while ( (line = stdout_br.readLine()) != null)
	            debug(line);    
	        
			// print stderr
			while ( (line = stderr_br.readLine()) != null)
	            debug(line);    
			
        } catch (IOException ex){
            ex.printStackTrace();  
        }
        
        try{
        	exitValue = proc.waitFor();
        }catch (java.lang.InterruptedException ex){
        	exitValue=CmdExitValues.ERROR.value();;
        }
        return exitValue;
	} // exec end
	
	
	/** Show a file browser and... 
	 * @return the path of the file chosen by the user
	 *
	public String browseFile(){
		OpenDialog od = new OpenDialog("Import file",null);
        String directory = od.getDirectory();
		String fileName = od.getFileName();
		String path= directory + fileName;
		return path;
	}*/
	
	
	/* handle button/keyboard pressing, from both this plugin and the windows it opens
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 *
	public void actionPerformed(ActionEvent e) {
	
		String label = e.getActionCommand();
		/*
		// select proper action method based on button's label
		if (label==null)
			return;
		else if (label.equals(ButtonLabels.CMD_INFO.label())){
			this.infoAction();
		}else if (label.equals(ButtonLabels.CMD_PREPROC.label())){
			exec("date");
		}else if (label.equals(ButtonLabels.CMD_ALIGN.label())){
			String directory=IJ.getDirectory("");
			exec("ls "+directory);
		}
	} // actionPerformed end */


    /********************* Trace methods for debugging ****************************/
    
    
	/** print s with a timestamp, using IJ Results window (@see IJ.write)
	 * @param s
	 */
	public static void debug(String s){
		Calendar calendar = Calendar.getInstance();
		java.util.Date now = calendar.getTime();
		String output="" + now.getTime() + " > " + s;
		if(tw != null)
			IJ.write(output);
		else
			System.err.println(output);
	}

	/** same as Debug, plus ex stack trace
	 * @param s
	 * @param ex Exception
	 */
	public static void debug(String s, Exception ex){
		debug(s);
		ex.printStackTrace();
	}
	
	private int getNextWindowId(){
		lastWindowId++;
		return lastWindowId;
	}
	
	private void createTomoWindow(){
			if(tw==null){
				tw= new TomoWindow(getNextWindowId());
			}
			UserAction newWindow= new UserAction(UserAction.ROOT_WINDOWID,"New Window");
			addUserAction(getWorkflowRoot(),newWindow);
			tw.setFirstAction(newWindow);
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
		// create tree
		UserAction root=new UserAction(1),n2=new UserAction(2),n3=new UserAction(3),n4=new UserAction(4),n5=new UserAction(5);
		// get one lineage
		addUserAction(root, n2);
		addUserAction(root, n4);
		addUserAction(n2, n3);
		addUserAction(n4, n5);
		printWorkflow();
		
		List <UserAction> l=getWorkflow(n5);
		// should print root - n4 - n5
		for(UserAction ua:l)
			debug(ua.toString());
	}
	
	public static void main(String[] args) {
		Xmipp_Tomo xt= new Xmipp_Tomo();
		xt.testWorkflow();
	}
	
}

