import java.io.File;

import xmipp.jni.Filename;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;

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
 * Why? To isolate useractions from low-level IO related details (like
 * metadatas)
 * 
 */

public class UserActionIO {
	private String fileName = null;
	MetaData imageMetadata;
	private String workingDir = null;

	public void applySelFile(UserActionIO parentIo) {
		String outputFilePath = getFilePath();
		if (isSelFile(outputFilePath)) {
			try {
				if(parentIo != null){
					MetaData metadata= parentIo.getMetadata();
					if(metadata != null)
						metadata.write(outputFilePath);
					else
						getMetadata().write(outputFilePath);
				}else
					getMetadata().write(outputFilePath);
			} catch (Exception ex) {
				Logger.debug("can not write metadata",ex);
			}
		}
	}

	/**
	 * 
	 * @return the path to the output file, which can be a .sel file or a stack file
	 */
	public String getFilePath() {
		return getWorkingDir() + "/" + fileName;
	}

	public static boolean isSelFile(String path) {
		// right now, identify file type by its extension
		// should check also uppercase like .SEL
		return path.endsWith(".sel");
	}

	// TODO: import more code details from Tomodata.readMetadata()
	public void setWorkingDir(String workingDir,boolean createWorkingDir) {
		// TODO: delete old working dir in case it already existed?
		this.workingDir = workingDir;
		if(createWorkingDir)
			new File(this.workingDir).mkdirs();
	}

	private MetaData getMetadata() {
		// metadata is initialized on metadata.read()
		if (getFilePath() == null)
			return null;

		if (imageMetadata == null) {
			imageMetadata = new MetaData();
			imageMetadata.enableDebug();
			try {
				imageMetadata.read(getFilePath());
			} catch (Exception ex) {
				Logger.debug("getmetadata",ex);
				imageMetadata = null;
			}
		}
		return imageMetadata;
	}

	/**
	 * 
	 * @param projection
	 *            1..N
	 * @return
	 */
	public String getFilePath(int projection) {
		String filename = null;
		if (getMetadata() != null){
			filename = getMetadata().getValueString(MDLabel.MDL_IMAGE,getProjectionId(projection));
			filename = getMetadata().fixPath(filename,getWorkingDir()); 
		}
		return filename;
	}
	
	/**
	 * Precondition: the working dir must have been assigned...
	 * @deprecated
	 * @return the path to a projection of the image stack 
	 */
	public String getOutputFilePath(int projection) {
		String ret = null;
		if (isSelFile(fileName))
			ret = getFilePath(projection);
		else
			ret = String.valueOf(getProjectionId(projection)) + "@" + getWorkingDir() + "/" + getFileName();
		return ret;
	}

	public long getProjectionId(int projection) {
		return getStackIds()[projection - 1];
	}

	public String getFileName() {
		return fileName;
	}

	public void setFileName(String name) {
		fileName = name;
	}

	public String getWorkingDir() {
		return workingDir;
	}

	public int getNumberOfProjections() {
		int ret = 0;
		if (getMetadata() != null)
			ret = getMetadata().size();
		return ret;
	}

	private long[] getStackIds() {
		long[] result = new long[1];
		if (getMetadata() != null)
			result = getMetadata().findObjects();
		return result;
	}

	public boolean isEnabled(long id) {
		boolean enabled = true;
		if (getMetadata().getValueInt(MDLabel.MDL_ENABLED, id) != 1)
			enabled = false;
		return enabled;
	}

	public void setEnabled(long id, int enabled) {
		getMetadata().setValueInt(MDLabel.MDL_ENABLED, enabled, id);
	}

	public String getInfo(int projectionNumber) {
		long id = getProjectionId(projectionNumber);
		return "(" + getMetadata().getValueString(MDLabel.MDL_IMAGE, id)
				+ "). Enabled: "
				+ getMetadata().getValueInt(MDLabel.MDL_ENABLED, id);
	}

	// Helper methods to manage tilt angles "list"
	// TODO: -low- get angles from .tlt files too (and save them like sel files)
	// - @see TiltSeriesIO.readTiltAngles
	public double getTiltAngle(long id) {
		return getMetadata().getValueDouble(MDLabel.MDL_ANGLETILT, id);
	}

}
