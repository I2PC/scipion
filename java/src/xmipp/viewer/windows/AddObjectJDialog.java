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

package xmipp.viewer.windows;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.LookAndFeel;

import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippWindowUtil;
import xmipp.viewer.RowHeaderRenderer;
import xmipp.viewer.models.GalleryRowHeaderModel;
import xmipp.viewer.models.ImageGalleryTableModel;
import xmipp.viewer.models.MetadataRowTableModel;

public class AddObjectJDialog extends XmippDialog {
	private static final long serialVersionUID = 1L;
	private JTable table;
	private MetadataRowTableModel model;
	// This will be used for check for results from the dialog
	public MetaData md;
	GridBagConstraints gbc = new GridBagConstraints();
	ImageGalleryTableModel gallery;
	
	public AddObjectJDialog(GalleryJFrame parent) {
		super(parent, "Classes", true);
		this.gallery = parent.gallery;
		this.md = gallery.data.md.getMetaDataRow();
		initComponents();
	}// constructor AddObjectJDialog
	
	@Override
	protected void createContent(JPanel panel) {
		setResizable(false);
		panel.setLayout(new GridBagLayout());
		gbc.anchor = GridBagConstraints.EAST;

		JPanel groupstbpn = new JPanel();
		JScrollPane sp = new JScrollPane();
		//groupstbpn.setBorder(BorderFactory.createTitledBorder("Classes"));
		groupstbpn.add(sp);
		sp.setOpaque(true);
		try {
			model = new MetadataRowTableModel(this, md);
			JList rowHeader = new JList();
			rowHeader.setModel(new GalleryRowHeaderModel(model.data));
			LookAndFeel.installColorsAndFont(rowHeader, "TableHeader.background",
					"TableHeader.foreground", "TableHeader.font");
			rowHeader.setCellRenderer(new RowHeaderRenderer());
			
			sp.setRowHeaderView(rowHeader);
			table = new JTable(model);
			model.setupTable(table);
			table.setPreferredScrollableViewportSize(new Dimension(350, 200));
			sp.setViewportView(table);
			int h = model.getCellSize().height;
			rowHeader.setFixedCellHeight(h);
			table.setRowHeight(h);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		panel.add(groupstbpn, XmippWindowUtil.getConstraints(gbc, 0, 1, 2));
		
		// listen to selection changes (only one row selected)
		table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
	}// function initComponents

}// class AddObjectJDialog
