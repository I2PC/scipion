import java.io.File;

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
	private String inputFilePath;
	MetaData imageMetadata;
	
	public String getInputFilePath() {
		return inputFilePath;
	}

	// TODO: import more code details from Tomodata.readMetadata()
	public void setInputFilePath(String f) {
		this.inputFilePath = f;
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
		String ret = null;
		if(getInputFilePath() != null)
			ret = new File(getInputFilePath()).getName();
		return ret;
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
}
