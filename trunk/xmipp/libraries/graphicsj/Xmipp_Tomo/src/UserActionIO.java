import java.io.File;

import xmipp.Filename;
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
 * Why? To isolate useractions from low-level IO related details (like metadatas)
 *
 */

public class UserActionIO {
	private String inputFileName=null;
	MetaData imageMetadata;
	private String workingDir=null;
	
	/**
	 * Precondition: the working dir must have been assigned...
	 * @return
	 */
	public String getInputFilePath() {
		return getWorkingDir() + "/" + inputFileName;
	}

	// TODO: import more code details from Tomodata.readMetadata()
	public void setWorkingDir(String workingDir){
		// TODO: delete old working dir in case it already existed?
		this.workingDir = workingDir;
		new File(this.workingDir).mkdirs();
	}
	
	private MetaData getMetadata(){
		// metadata is initialized on metadata.read()
		if (getInputFilePath() == null)
			return null;
		
		if(imageMetadata==null){
			imageMetadata = new MetaData();
			imageMetadata.enableDebug();
			try{
				imageMetadata.read(getInputFilePath());
			}catch (Exception ex){
				imageMetadata = null;
			}
		}
		return imageMetadata;
	}
	
	/**
	 * 
	 * @param projection 1..N
	 * @return
	 */
	public String getFilePath(int projection){
		String filename= null;
		if(getMetadata() != null)
			filename = getMetadata().getValueString(MDLabel.MDL_IMAGE,getProjectionId(projection));
		return filename;
	}
	
	public long getProjectionId(int projection){
		return getStackIds()[projection-1];
	}
	
	public String getInputFileName(){
		return inputFileName;
	}
	
	public void setInputFileName(String name){
		inputFileName = name;
	}
	
	public String getWorkingDir(){
		return workingDir;
	}
	
	public int getNumberOfProjections(){
		int ret=0;
		if(getMetadata() != null)
			ret = getMetadata().size();
		return ret;
	}
	
	
	private long[] getStackIds(){
		long [] result=new long[1];
		if(getMetadata() != null)
			result= getMetadata().findObjects();
		return result;
	}

	public boolean isEnabled(long id) {
		boolean enabled = true;
		if (getMetadata().getValueInt(MDLabel.MDL_ENABLED, id) != 1)
			enabled = false;
		return enabled;
	}
	
	public void setEnabled(long id,int enabled){
		getMetadata().setValueInt(MDLabel.MDL_ENABLED, enabled, id);
	}
	
	public String getInfo(int projectionNumber){
			long id=getProjectionId(projectionNumber);
			return "(" + getMetadata().getValueString(MDLabel.MDL_IMAGE, id) + "). Enabled: " + getMetadata().getValueInt(MDLabel.MDL_ENABLED, id);
	}
}
