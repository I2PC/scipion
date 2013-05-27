

package xmipp.utils;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;
import java.util.Map;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComponent;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.KeyStroke;

/** This a menu creator to provide some useful functions to define menu structures
 * It will allow creation of JMenuBar and JPopupMenu
 * */
public abstract class XmippMenuCreator implements ActionListener {

	protected Map<String, JMenuItem> items;
	/** This will be used to group Radiobuttons under a Menu
	 * A new menu defines a new group and also a separator 
	 * in the same menu defines a new group;
	 */
	protected ButtonGroup group; 
	protected JComponent menu;
	
	/** Constructor */
	public XmippMenuCreator(JComponent menu){
		try {
			this.menu = menu;
			items = new Hashtable<String, JMenuItem>();
			group = null;
			createItems();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**
	 * The idea behind this method is simplify the creation of menu options. The
	 * name of the menu item should be Parent.Item to refer to parent, and also	
	 * the prefix cb_ and rb_ shoulb be used to create JCheckBoxMenuItem and
	 * JRadioButtonMenuItem respectively Example: addItem("File", "File");
	 * addItem("File.SaveAs", "Save as...", "save_as.gif")
	 * 
	 * @return the created MenuItem, null if an error occurs
	 * @throws Exception
	 */
	public JMenuItem addItem(String key, String text, String... values){
		String name = key;
		JMenuItem item = null;
		JMenu parent = null;
		boolean isMenu = false;

		try {
		if (key.contains(".")) {
			int pos = key.lastIndexOf(".");
			name = key.substring(pos + 1);
			parent = (JMenu)getItem(key.substring(0, pos));
		}

		// Create the right type of menu item
		if (name.endsWith("cb"))
			item = new JCheckBoxMenuItem();
		else if (name.endsWith("rb")){
			if (group == null)
				group = new ButtonGroup();
			item = new JRadioButtonMenuItem();
			group.add(item);
		}
		else if (name.endsWith("mi"))
			item = new JMenuItem();
		else {
			item = new JMenu();
			isMenu = true;
			group = null; //define a new group
		}
		// Setup the item
		item.setText(text);
		item.setActionCommand(key);
		item.addActionListener(this);
		// Store item and add to menu
		items.put(key, item);

		if (!isMenu) {
			int n = values.length;
			switch (n) {
			case 2:
				// The second argument could be the keystroke, be careful
				item.setAccelerator(KeyStroke.getKeyStroke(values[1]));
			case 1:
				if (values[0] != null)
					item.setIcon(XmippResource.getIcon(values[0]));
				break;
			case 0:
				break; // ignore if no arguments
			default:
				throw new Exception(
						"No more than 2 arguments to create menu item!!!");
			}
		}

		if (parent == null)
			menu.add(item);
		else
			parent.add(item);
		} catch (Exception e){
			e.printStackTrace();
		}
		return item;
	}
	
	public JMenuItem getItem(String key) {
		return items.get(key);
	}

	public boolean getItemEnabled(String key){
		return getItem(key).isEnabled();
	}
	
	public void setItemEnabled(String key, boolean value) {
		getItem(key).setEnabled(value);
	}
	
	public boolean getItemSelected(String key){
		// Create the right type of menu item
		if (key.endsWith("rb"))
			return ((JRadioButtonMenuItem)getItem(key)).isSelected();
		if (key.endsWith("cb"))
			return ((JCheckBoxMenuItem)getItem(key)).isSelected();
		return false;
	}
	
	public void setItemSelected(String key, boolean value){
		//DEBUG.printMessage("setItemSelected -> key:" + key);
		// Create the right type of menu item
		if (key.endsWith("rb"))
			((JRadioButtonMenuItem)getItem(key)).setSelected(value);
		else if (key.endsWith("cb"))
			((JCheckBoxMenuItem)getItem(key)).setSelected(value);		
	}
	
	public void setItemVisible(String key, boolean value){
		getItem(key).setVisible(value);
	}

	@Override
	public void actionPerformed(ActionEvent evt) {
		handleActionPerformed(evt);
	}
	
	/** Abstract methods */
	abstract protected void handleActionPerformed(ActionEvent evt);

	/** This abstract method will be called from constructor
	 * to create the items 
	 */
	abstract protected void createItems() throws Exception;
	
	/** Method to initialize items state */
	public void initItems() { }
	
	
	/** Public menu items constants */
	public final String FILE = "File";
	public final String FILE_OPEN = "File.Open_mi";
	public final String FILE_OPENWITH_IJ = "File.OpenWithIJ_mi";
	public final String FILE_OPENMICROGRAPHS = "File.OpenMicrographs_mi";
	public final String FILE_OPENWITH_CHIMERA = "File.OpenWithChimera_mi";
	public final String FILE_INFO = "File.FileInfo_mi";
	public final String FILE_SAVE = "File.Save_mi";
	public final String FILE_SAVEAS = "File.SaveAs_mi";
	public final String FILE_REFRESH = "File.Refresh_mi";
	public final String FILE_EXIT = "File.Exit_mi";
	public final String DISPLAY = "Display";
	public final String DISPLAY_NORMALIZE = "Display.Normalize_cb";
	public final String DISPLAY_SHOWLABELS = "Display.ShowLabels_cb";
	public final String DISPLAY_RENDERIMAGES = "Display.RenderImages_cb";
	
	public final String DISPLAY_RENDERIMAGECOLUMN = "Display.RenderImagesColumn";
	
	
	public final String DISPLAY_APPLYGEO = "Display.ApplyGeo_cb";
	public final String DISPLAY_WRAP = "Display.Wrap_cb";
	public final String DISPLAY_COLUMNS = "Display.Columns_mi";
	public final String DISPLAY_RESLICE = "Display.Reslice";
	public final String DISPLAY_RESLICE_VIEWS[] = { 
			"Display.Reslice.ZNeg_rb",
			"Display.Reslice.YNeg_rb",
			"Display.Reslice.XNeg_rb",
			"Display.Reslice.YPos_rb",
			"Display.Reslice.XPos_rb"};

	public final String METADATA = "Metadata";
	public final String MD_CLASSES = "Metadata.Classes_mi";
	public final String MD_EDIT_COLS = "Metadata.EditCols_mi";
	public final String MD_ADD_OBJECT = "Metadata.AddObject_mi";
	public final String MD_REMOVE_SELECTION = "Metadata.RemoveSelection_mi";
	public final String MD_SAVE_SELECTION = "Metadata.SaveSelection_mi";
	public final String MD_FIND_REPLACE = "Metadata.FindReplace_mi";
	public final String MD_REMOVE_DISABLED = "Metadata.RemoveDisabled_mi";
	public final String STATS = "Metadata.Stats";
	public final String STATS_AVGSTD = "Metadata.Stats.AvgStd_mi";
	public final String STATS_PCA = "Metadata.Stats.Pca_mi";
	public final String STATS_FSC = "Metadata.Stats.Fsc_mi";
	public final String MD_PLOT = "Metadata.Plot_mi";
	
	public final String HELP = "Help";
	public final String HELP_ONLINE = "Help.Online_mi";
	public final String KEY_ASSIST = "Help.Key_mi";
	//public final String HELP_ABOUT = "Help.About_mi";
	

	public final static String ENABLED = "Enabled_mi";
	public final static String DISABLED = "Disabled_mi";
	public final static String REFRESH = "Refresh_mi";
	public final static String OPEN = "Open_mi";
	public final static String OPEN_ASTEXT = "OpenAsText_mi";
	public final static String CTF_PROFILE = "CTFProfile_mi";
	public final static String CTF_RECALCULATE = "CTFRecalculate_mi"; 
	public final static String SELECT = "Select";
	public final static String SELECT_ALL = "Select.All_mi";
	public final static String SELECT_TOHERE = "Select.ToHere_mi";
	public final static String SELECT_FROMHERE = "Select.FromHere_mi";
	public final static String SET_CLASS = "SetClass_mi";
	public final static String OPEN_IMAGES = "OpenImages_mi";

}
