/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
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

package xmipp.jni;

/**
 * Binding class for accessing C++ MDRow implementation.
 */
public class MDRow {
	// hold pointer to Image class in C++ space
	private long peer;

	static {
		System.loadLibrary("XmippJNI");
		// storeIds();
	}

	// caching some ids
	// private static native void storeIds();
	// functions to create images
	private native void create();

	// destructor
	public synchronized native void destroy();
	
	// clear metadata
	public native void clear();

	public native int size();

	public native boolean containsLabel(int label);

	public native boolean removeLabel(int label);

	// get values from metadata
	public native int getValueInt(int label);

	public native long getValueLong(int label);

	public native double getValueDouble(int label);

	/** Return the value of some label as String */
	public native String getValueString(int label);
	

	// set functions conection with MDRow class in C++
	public boolean setEnabled(boolean value) {
		return setValueInt(MDLabel.MDL_ENABLED, value ? 1 : -1);
	}

	public boolean getEnabled(long objId) {
		return getValueInt(MDLabel.MDL_ENABLED) > 0;
	}

	public native boolean setValueInt(int label, int value);

	public native boolean setValueLong(int label, long value);

	public native boolean setValueDouble(int label, double value);

	public native boolean setValueString(int label, String value);

	public native boolean setValueBoolean(int label, boolean value);

	public native void addLabel(int label);

	/*********** Non-native functions ********************/
	/** Create empty metadata */
	public MDRow() {
		create();
	}

	/**
	 * Will be called by GarbageCollector before destroying. Needed to free C++
	 * object memory.
	 */
	@Override
	protected void finalize() throws Throwable {
		super.finalize();
		destroy();
	}

}
