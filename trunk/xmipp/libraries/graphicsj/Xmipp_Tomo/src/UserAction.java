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
public class UserAction {
	public static int ROOT_WINDOWID=0;
	
	private int windowId;
	private String command=null,parameters=null,name="";
	private Plugin plugin;
	private boolean neededForFile=false;
	
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
	
	public UserAction(int windowId,String cmd,Plugin plugin){
		this(windowId,cmd);
		setPlugin(plugin);
		setNeededForFile(true); 
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
		String res= getName() + ". " + getCommand();
		return res;
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
	
}
