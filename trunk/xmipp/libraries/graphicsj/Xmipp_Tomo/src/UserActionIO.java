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
	private String inputFile,outputFile;
	MetaData imageMetadata;
	
	public String getInputFile() {
		return inputFile;
	}

	// TODO: import more code details from Tomodata.readMetadata()
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
		try{
			getMetadata().read(inputFile);
		}catch (Exception ex){
			Logger.debug("readMetadata - problem with "+inputFile+". Please check that the file exists and its permissions", ex);
		}
	}
	
	public MetaData getMetadata(){
		// metadata is initialized on metadata.read()
		
		if(imageMetadata==null){
			imageMetadata = new MetaData();
			imageMetadata.enableDebug();
		}
		return imageMetadata;
	}
	
	/**
	 * 
	 * @param projection 1..N
	 * @return
	 */
	public String getFilePath(int projection){
		if(!initialized())
			return null;
		long projectionId = getStackIds()[projection-1];
		String filename= getMetadata().getValueString(MDLabel.MDL_IMAGE,projectionId);
		return filename;
	}
	
	// TODO: improve the test to check if metadata has been read
	private boolean initialized(){
		return inputFile != null;
	}
	
	private long[] getStackIds(){
		long [] result=new long[1];
		if(getMetadata() != null)
			result= getMetadata().findObjects();
		return result;
	}
}
