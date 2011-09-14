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
 * - Why a Command class?
 * To provide a simple declaration (using this class constructors) of the program's commands,
 * independent of whether the program implements them as buttons, menus... 
 * @see XmippTomoCommands
 */

import javax.swing.ImageIcon;

public class Command {

	/* Complex commands may transit through some states.
	 * For example, loading a file takes time - the command may be canceled.
	 * Furthermore, the results of other commands may change depending upon the load command is still running,
	 * or it has finished.
	 */
	public static enum State{IDLE,LOADING, LOADED, RELOADING,CANCELED;};
	
	private String id;
	private String label;
	private String method;
	private String iconName;
	private boolean enabled = true;

	/**
	 * 
	 * @param id unique identifier string, using a dotted syntax, like "align.correlation"
	 * @param label what the user will read
	 * @param method to be called when this command takes place
	 * @param enabled at startup?
	 * @param iconName (optional) URL like ""resources/icon-play.png", or NULL
	 */
	Command(String id,String label, String method, boolean enabled, String iconName) {
		this.id = id;
		this.method = method;
		this.label = label;
		this.enabled = enabled;
		this.iconName = iconName;
	}

	public String getLabel() {
		return label;
	}

	public boolean isEnabled() {
		return enabled;
	}

	public String getMethod() {
		return method;
	}

	public ImageIcon getIcon(){
		if(iconName == null)
			return null;
		java.net.URL imageURL = Command.class.getResource(iconName);
		if (imageURL != null) {
			return new ImageIcon(imageURL);
		}else
			return null;
	}

	public String getId() {
		return id;
	}

}
