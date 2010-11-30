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

import ij.Prefs;

import java.awt.Component;
import java.io.File;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * @author jcuenca Swing file dialog
 */
public class FileDialog {
	private JFileChooser fc;
	private Component parent;
	private static File defaultDirectory;

	public void setupOpenDialog(String title, String path, final String fileName) {
		fc = new JFileChooser();
		fc.setDialogTitle(title);
		fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		fc.setCurrentDirectory(getDefaultDirectory());

		File fdir = null;
		if (path != null)
			fdir = new File(path);
		if (fdir != null)
			fc.setCurrentDirectory(fdir);
		if (fileName != null)
			fc.setSelectedFile(new File(fileName));
	}

	public static String openDialog(String title, Component parent) {
		FileDialog fd = new FileDialog();
		fd.setParent(parent);
		fd.setupOpenDialog(title, null, null);
		FileNameExtensionFilter filter = new FileNameExtensionFilter(
				"Tomo images", "sel", "mrcs");
		fd.setFilter(filter);
		int status = fd.showOpenDialog();
		if (status != JFileChooser.APPROVE_OPTION)
			return "";
		return fd.getPath();
	}

	public int showOpenDialog() {
		return fc.showOpenDialog(getParent());
	}

	public String getPath() {
		File file = fc.getSelectedFile();
		if (file == null)
			return "";
		String name = file.getName();
		String dir = fc.getCurrentDirectory().getPath() + File.separator;

		String returnPath = "";
		if (name != null)
			returnPath = dir + name;
		return returnPath;
	}

	public void setFilter(FileNameExtensionFilter filter) {
		fc.setFileFilter(filter);
	}

	/**
	 * @return the parent
	 */
	private Component getParent() {
		return parent;
	}

	/**
	 * @param parent
	 *            the parent to set
	 */
	private void setParent(Component parent) {
		this.parent = parent;
	}

	/**
	 * Returns the current working directory, which may be null. The returned
	 * string always ends with the separator character ("/" or "\").
	 */
	public static File getDefaultDirectory() {
		if (defaultDirectory == null)
			defaultDirectory = new File(Prefs.getString(Prefs.DIR_IMAGE));
		return defaultDirectory;
	}

}
