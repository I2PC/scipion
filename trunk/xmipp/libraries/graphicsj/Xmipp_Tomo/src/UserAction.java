import java.io.File;
import java.util.Date;


import javax.swing.text.BadLocationException;
import javax.swing.text.PlainDocument;
import javax.swing.tree.TreePath;

import xmipp.MDLabel;
import xmipp.MetaData;

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
 * - Why?
 * A simple data structure to save user action's info (for the user action workflow)
 */
/* TODO: map user actions to a SQLite database (using Java code)
 * (reuse ideas from xmipp protocols python code)
 * Define the parameters each action may use, and when they may or may not use them
 * (for instance, some operations produce volumes, other need alignedstacks as input...).
 * In short, define the rules of (in)compatibility between parameters and actions/steps
 *  Ana action can have one or many inputs. By now an action has only one output
 * Then, once all the posibilities are defined, a user can build a specific workflow, choosing among the steps and setting 
 * the parameters
 * The database may also hold the user interface definition
 * 
 * STATIC TEMPLATES (possibilities definition, always the same for every project)
 * Action: id, name, label, method, iconName, enabled (at startup), panel (where to place the button), panelOrder
 * Parameter: id, name, type {Volume, Stack, AlignedStack, Float, String}
 * ActionParameters: actionId, paramId, defaultParamValue, paramType{input/output}
 * 
 * DYNAMIC STEPS (workflow, changes from project to project)
 * StepParameters: stepId, paramId, paramValue
 * Step: id, order, stepId, workingDir (each project has its working dir where the workflow database
 * is stored as .project.sqlite, and all the files are stored too there), comments 
 */
public class UserAction {
	public static int ROOT_WINDOWID=0;
	private static String WORKFLOW_ROOT_NAME="Start";
	
	private int windowId, actionId;
	

	private String command=null,parameters=null,name=null;
	String progress=null;
	private Plugin plugin;
	private boolean neededForFile=false;
	private PlainDocument comments;
	private UserActionIO output,parentOutput=null;

	public UserActionIO getParentOutput() {
		return parentOutput;
	}

	private void setParentOutput(UserActionIO parentOutput) {
		this.parentOutput = parentOutput;
	}
	
	public void setParentOutput(UserAction parent) {
		if(parent != null)
			setParentOutput(parent.getOutput());
	}

	private UserActionIO getOutput() {
		if(output == null)
			output = new UserActionIO();
		return output;
	}

	/**
	 * Specify the output path manually (typically when loading a new file)
	 * @param path
	 */
	public void setOutputFilePath(String path){
		File output = new File(path);
		String workingDir=output.getParent();
		if (workingDir == null)
			workingDir = "./";
		getOutput().setWorkingDir(workingDir,false);
		getOutput().setFileName(output.getName());
	}
	
	/*public String getInputFilePath(){
		return getIoDetails().getInputFilePath();
	}*/
	
	public String getOutputFilePath(int projection){
		return getOutput().getFilePath(projection);
	}
	
	public String getOutputFileName(){
		return getOutput().getFileName();
	}
	
	
	public boolean isEnabled(long id) {
		return getOutput().isEnabled(id);
	}
	
	public void setEnabled(long id,int enabled){
		getOutput().setEnabled(id, enabled);
	}
	
	public String getInputFilePath(int projection){
		return getParentOutput().getFilePath(projection);
	}
	
	public int getNumberOfProjections(){
		return getOutput().getNumberOfProjections();
	}
	
	public long getProjectionId(int projection){
		return getOutput().getProjectionId(projection);
	}
	
	public double getTiltAngle(long id){
		return getOutput().getTiltAngle(id);
	}
	
	public String getInfo(int projectionNumber){
		return getOutput().getInfo(projectionNumber);
	}
	

	/**
	 * Call this method once the working dir has been set (when the action has been added to the workflow)
	 * @param name
	 */
	public void setOutputFileName(String name){
		getOutput().setFileName(name);
	}
	
	/**
	 * If you need to re-set the working dir, use setFilePath
	 * @param workflowWorkingDir
	 * @throws NullPointerException if the parent output has not been set before
	 */
	// TODO: maybe it's better to set the working dir in the constructor...
	public void createWorkingDir() throws NullPointerException {
		String uaWorkingSubdir = getWorkingSubdirName();
		String parentWorkingDir = getParentOutput().getWorkingDir(); 
		String workingDir =	parentWorkingDir + "/" + uaWorkingSubdir;
		getOutput().setWorkingDir(workingDir,true);
	}
	
	// TODO: -current- use four digits for Id
	private String getWorkingSubdirName(){
		//TODO: it may be a good idea to also strip all punctuation signs with replaceAll
		return String.valueOf(getId()) + "-" + getName().replaceAll("\\s+", "");
	}

	public String getComments() {
		String result="";
		try{
			result= comments.getText(0, comments.getLength());
		}catch (BadLocationException ex){}
		return result;
	}
	
	public PlainDocument getCommentsDocument(){
		if(comments == null)
			comments=new PlainDocument();
		return comments;
	}

	public void setComments(String text) {
		try{
			comments.replace(0,comments.getLength(),text,null);
		}catch (BadLocationException ex){}
	}

	public UserAction(int windowId){
		setWindowId(windowId);
		
	}
	
	public UserAction(int windowId,String cmd){
		this(windowId);
		setCommand(cmd);
	}
	
	public UserAction(int windowId,String cmd,String params){
		this(windowId,cmd);
		setParameters(params);
	}
	
	public UserAction(int windowId,String name,String cmd,String params){
		this(windowId,cmd,params);
		setName(name);
	}
	
	public UserAction(int windowId,String cmd,Plugin plugin){
		this(windowId,cmd);
		setPlugin(plugin);
		setNeededForFile(true); 
	}
	
	public UserAction(int windowId,String name, String cmd,Plugin plugin){
		this(windowId,cmd,plugin);
		setName(name);
	}
	
	public static UserAction start(){
		return new UserAction(UserAction.ROOT_WINDOWID,WORKFLOW_ROOT_NAME,"","");
		
	}
	
	public int getWindowId() {
		return windowId;
	}

	public void setWindowId(int windowId) {
		this.windowId = windowId;
	}

	public String getCommand() {
		return command;
	}

	public void setCommand(String command) {
		this.command = command;
	}

	public String getParameters() {
		return parameters;
	}

	public void setParameters(String parameters) {
		this.parameters = parameters;
	}
	
	public String toString(){
		String actionName = getName(), actionText="";
		if(actionName != null)
		  actionText = actionText + actionName + ".";

		if(getProgress() != null)
			actionText = actionText + " :: " + getProgress();
		
		return actionText;
	}

	public String getCommandDetails(){
		String details="";

		String actionCommand = getCommand();
		if(actionCommand != null){
		  details = details + actionCommand;
		  String actionParameters = getParameters();
		  if(actionParameters!= null)
			  details = details + " ["+actionParameters+"]";
		  if(getPlugin()!= null)
			  details = details + " ["+getPlugin().getOptions()+"]"; 
		  
		}
		return details;
	}
	
	public Plugin getPlugin() {
		return plugin;
	}

	public void setPlugin(Plugin plugin) {
		this.plugin = plugin;
	}

	public boolean isNeededForFile() {
		return neededForFile;
	}

	public void setNeededForFile(boolean neededForFile) {
		this.neededForFile = neededForFile;
	}

	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public String getProgress(){
		if(WORKFLOW_ROOT_NAME.equals(getName()))
			return null;
		String ret = "100%";
		if (ret == null){
			int p = (getCommand().length() % 4) * 25;
			ret = String.valueOf(p) + "%";
		}
		return ret;
	}
	
	public void setProgress(String progress){
		this.progress = progress;
	}
	
	public int getId() {
		return actionId;
	}

	public void setId(int actionId) {
		this.actionId = actionId;
	}
	
	public void applySelFile() {
		getOutput().applySelFile(null);
	}
	
	/**
	 * 
	 * @throws NullPointerException if parent output has not been set
	 */
	public void createSelFile() throws NullPointerException{
		getOutput().applySelFile(getParentOutput());
	}
	

}
