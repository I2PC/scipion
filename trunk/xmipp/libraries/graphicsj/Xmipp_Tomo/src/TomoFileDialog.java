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
 * Why a FileDialog?
 * ImageJ default file dialogs are simple AWT. This class uses the more beautiful and powerful Swing dialogs.
 * 
 * - Features:
 * FileNameExtensionFilters
 * Remember last path
 */
import ij.Prefs;
import java.awt.Component;
import java.io.File;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileNameExtensionFilter;

// TODO: -low-means low priority task
// TODO: -low- reuse FileDialog code by extending it
public class TomoFileDialog {
	private JFileChooser fileChooser;
	private static File lastChosenDirectory;
	private static String prefsPrefix="xmipptomo.";

	private void setupDialog(String title, String path, final String fileName,int type) {
		fileChooser = new JFileChooser();
		fileChooser.setDialogTitle(title);
		fileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		File lastChosenDirectory = getLastChosenDirectory(title);
		if (lastChosenDirectory == null)
			lastChosenDirectory=new File(".");
		fileChooser.setCurrentDirectory(TomoFileDialog.lastChosenDirectory);
		fileChooser.setDialogType(type);

		File fdir = null;
		if (path != null)
			fdir = new File(path);
		if (fdir != null)
			fileChooser.setCurrentDirectory(fdir);
		if (fileName != null)
			fileChooser.setSelectedFile(new File(fileName));
	}

	public static String openDialog(String title, Component parent) {
		TomoFileDialog fd = new TomoFileDialog();
		fd.setupDialog(title, null, null,JFileChooser.OPEN_DIALOG);
		fd.addTomoImagesFileExtensionFilters();
		int status = fd.showOpenDialog(parent);
		if (status != JFileChooser.APPROVE_OPTION)
			return "";
		String path=fd.getPath();
		setLastChosenDirectory(title, path);
		return path;
	}
	
	public static String saveDialog(String title, Component parent) {
		TomoFileDialog fd = new TomoFileDialog();
		fd.setupDialog(title, null, null,JFileChooser.SAVE_DIALOG);
		fd.addTomoImagesFileExtensionFilters();
		int status = fd.showSaveDialog(parent);
		if (status != JFileChooser.APPROVE_OPTION)
			return "";
		String path=fd.getPath();
		setLastChosenDirectory(title, path);
		return path;
	}
	
	public void addTomoImagesFileExtensionFilters(){
		addFilter(new FileNameExtensionFilter("Spider", "spi"));
		addFilter(new FileNameExtensionFilter("Sel file", "sel"));
		addFilter(new FileNameExtensionFilter("MRC", "mrcs"));
		addFilter(new FileNameExtensionFilter("VOL", "vol"));
		addFilter(new FileNameExtensionFilter("Xmipp Stack", "stk"));
	}

	/**
	 * @depends on setupOpenDialog
	 * @return Cancel, Approve or Error - @see fileChooser.showOpenDialog
	 */
	private int showOpenDialog(Component parent) {
		return fileChooser.showOpenDialog(parent);
	}
	
	/**
	 * @depends on setupOpenDialog
	 * @return Cancel, Approve or Error - @see fileChooser.showOpenDialog
	 */
	private int showSaveDialog(Component parent) {
		return fileChooser.showSaveDialog(parent);
	}

	public String getPath() {
		boolean removeFile=false;
		File selectedFile = fileChooser.getSelectedFile();
		if (selectedFile == null)
			return null;

		if (fileChooser.getDialogType() == JFileChooser.SAVE_DIALOG)
			if ((selectedFile != null) && selectedFile.exists()) {
				int response = JOptionPane.showConfirmDialog(fileChooser,"The file " + selectedFile.getName() + " already exists. Do you want to replace it?",
						"Ovewrite file", JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
				if (response != JOptionPane.YES_OPTION)
					return null;
				// remove file to prevent overwriting problems
				removeFile=true;
			}


		String name = selectedFile.getName();
		String dir = fileChooser.getCurrentDirectory().getPath() + File.separator;

		String returnPath = "";
		if (name != null){
			returnPath = dir + name;
			if(removeFile)
				new File(returnPath).delete();
		}
		return returnPath;
	}

	public void addFilter(FileNameExtensionFilter filter) {
		// fileChooser.setFileFilter(filter);
		fileChooser.addChoosableFileFilter(filter);
	}

	/**
	 * @return path, which may be null.
	 * This path always ends with the separator character ("/" or "\").
	 */
	private static File getLastChosenDirectory(String title) {
		if (lastChosenDirectory == null){
			String dirImage = Prefs.get(prefsPrefix+title,".");
			lastChosenDirectory = new File(dirImage);
		}
		return lastChosenDirectory;
	}

	private static void setLastChosenDirectory(String title,String filePath){
		String directory=new File(filePath).getParent();
		Prefs.set(prefsPrefix+title,directory);
	}
}
