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

import javax.swing.ImageIcon;

/**
 * Data structure that relates an action (method) and its parameters iwith a
 * label, an icon...
 * @author jcuenca
 */

public class Command {
	
	public static Command PLAY= new Command("controls.play","", "playPause", false, TomoWindow.PLAY_ICON),
	PLAY_LOOP= new Command("controls.play_loop","Loop", "changePlayMode", true, null),
	LOAD=new Command("file.emload","Load","loadEM",true,null),
	XRAY=new Command("file.xrayload","Import X-Ray","loadXray",true,null),
	SAVE=new Command("file.save","Save","save",false,null),
	NORMALIZE_SERIES=new Command("file.normalize","Normalize series","normalize",false,null),
	DEFINE_TILT = new Command("file.tilt","Set tilt angles","setTilt",false,null),
	GAUSSIAN = new Command("proc.gaussian","Gaussian","gaussian", false,null),
	MEDIAN = new Command("proc.median","Median","median", false,null),
	SUB_BACKGROUND = new Command("proc.sub_background","Substract Background","subBackground",false,null),
	ENHANCE_CONTRAST = new Command("proc.enhance_contrast","Enhance Contrast","enhanceContrast",false,null),
	GAMMA_CORRECTION = new Command("proc.gamma_correction","Gamma Correction","gammaCorrection",false,null),
	MEASURE = new Command("file.measure","Measure","measure",false,null),
	APPLY = new Command("proc.apply","Apply","apply",false,null),
	ALIGN_AUTO = new Command("align.auto","Auto","alignAuto",false,null),
	ALIGN_MANUAL = new Command("align.manual","Manual","alignManual",false,null),
	ALIGN_CORRELATION = new Command("align.correlation","Quick","alignCorrelation",false,null),
	PRINT_WORKFLOW = new Command("debug.print_workflow","Print workflow","printWorkflow",false,null);
	
	public static enum State{IDLE,LOADING, LOADED, RELOADING,CANCELED;};
	
	private String id;
	private String label;
	private String method;
	private String iconName;
	private boolean enabled = true;

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

	/**
	 * @return the method
	 */
	public String getMethod() {
		return method;
	}

	public ImageIcon getIcon() throws Exception {
		if(iconName == null)
			return null;
		java.net.URL imageURL = Xmipp_Tomo.class.getResource(iconName);
		if (imageURL != null) {
			return new ImageIcon(imageURL);
		}else
			throw new Exception("Command.getIcon() - icon not found");
	}

	/**
	 * @return the id
	 */
	public String getId() {
		return id;
	}

}
