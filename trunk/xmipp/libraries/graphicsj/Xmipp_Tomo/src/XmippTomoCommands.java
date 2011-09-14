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
 * - Why XmippTomoCommand?
 * To provide a central repository of all the app commands 
 * - Alternatives:
 * A configuration file (or database) parsed at runtime before GUI construction
 */


public class XmippTomoCommands {
	public static Command 
	PLAY= new Command("controls.play","", "playPause", false, TomoWindow.PLAY_ICON),
	PLAY_LOOP= new Command("controls.play_loop","Loop", "changePlayMode", true, null),
	LOAD=new Command("file.emload","Load","loadEM",true,null),
	XRAY=new Command("file.xrayload","Import X-Ray","loadXray",true,null),
	CONVERT=new Command("file.convert","Convert","convert",false,null),
	ADJUSTBC=new Command("controls.adjustbc","Adjust Brightness/Contrast","adjustbc",false,null),
	NORMALIZE_SERIES=new Command("preproc.normalize","Normalize series","normalize",false,null),
	DEFINE_TILT = new Command("file.tilt","Set tilt angles","setTilt",false,null),
	DISCARD_PROJECTION = new Command("file.discard_projection","Discard Projection","discardProjection",false,null),
	// defined here only for the label
	UNDO_DISCARD_PROJECTION = new Command("file.undo_discard_projection","Undo Discard Proj.","discardProjection",false,null),
	HOTSPOT_REMOVAL = new Command("proc.hotspotremoval","Hotspot Removal","hotspotRemoval", false,null),
	GAUSSIAN = new Command("proc.gaussian","Gaussian","gaussian", false,null),
	MEDIAN = new Command("proc.median","Median","median", false,null),
	SUB_BACKGROUND = new Command("proc.sub_background","Substract Background","subBackground",false,null),
	ENHANCE_CONTRAST = new Command("proc.enhance_contrast","Enhance Contrast","enhanceContrast",false,null),
	GAMMA_CORRECTION = new Command("proc.gamma_correction","Gamma Correction","gammaCorrection",false,null),
	HISTOGRAM_EQUALIZATION = new Command("proc.histogram_equalization","Histogram Equalization","histogramEqualization",false,null),
	CROP = new Command("proc.crop","Crop","crop",false,null),
	BANDPASS = new Command("proc.bandpass","Bandpass Filter","bandpass",false,null),
	MEASURE = new Command("controls.measure","Measure","measure",false,null),
	// Apply now does not make sense, every step is written to disk
	// APPLY = new Command("proc.apply","Apply&Save","apply",false,null),
	ALIGN_AUTO = new Command("align.auto","Auto","alignAuto",false,null),
	ALIGN_MANUAL = new Command("align.manual","Manual","alignManual",false,null),
	ALIGN_CORRELATION = new Command("align.correlation","Quick","alignCorrelation",false,null),
	PRINT_WORKFLOW = new Command("debug.print_workflow","Print workflow","printWorkflow",false,null), 
	CURRENT_PROJECTION_INFO = new Command("debug.current_projection_info","Current Projection Info","currentProjectionInfo",false,null);;
	
}
