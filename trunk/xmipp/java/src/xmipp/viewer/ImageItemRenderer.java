package xmipp.viewer;

import ij.ImagePlus;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.border.Border;
import javax.swing.table.DefaultTableCellRenderer;

import xmipp.utils.XmippResource;

/**
 * 
 * @author Juanjo Vega
 */
public class ImageItemRenderer extends DefaultTableCellRenderer {
	private static final long serialVersionUID = 1L;

	Border BORDER_SELECTED = new StrokeBorder(Color.RED, 3);
	Border BORDER_FOCUSED = BorderFactory.createLineBorder(Color.RED, 3);
	boolean hackBorders = true;

	public ImageItemRenderer() {
		super();
		setOpaque(true);
		setHorizontalAlignment(JLabel.CENTER);
		setHorizontalTextPosition(JLabel.CENTER);
		setVerticalTextPosition(JLabel.BOTTOM);		
	}

	public ImageItemRenderer(boolean hackBorders) {
		this();
		this.hackBorders = false;
	}

	@Override
	public Component getTableCellRendererComponent(JTable table, Object object,
			boolean isSelected, boolean hasFocus, int row, int column) {
		ImageItem item = (ImageItem) object;
		// DEBUG.printMessage("*** 1. Rendering: " + item.getLabel());
		// Calls super class so foreground, background, borders and rest of
		// stuff is set.
		super.getTableCellRendererComponent(table, null, item != null
				&& item.isSelected, item != null && hasFocus, row, column);

		if (item != null) {
			// DEBUG.printMessage("*** 2. Rendering: " + item.getLabel());
			// AbstractXmippTableModel tableModel = (AbstractXmippTableModel)
			// table.getModel();

			setPreferredSize(item.cellDim);

			// ... and sets it.
			setEnabled(item.isEnabled);

			// Loads image...
			ImagePlus imp = item.getImage();
			if (imp != null)
				setIcon(new ImageIcon(imp.getImage()));
			else
				setIcon(XmippResource.MISSING_ICON);

			// Tooltip.
			setToolTipText(item.getLabel());

			// (Shows label only when required).
			if (item.showLabel) {
				String label = cutString(item.getLabel(), table
						.getColumnModel().getColumn(column).getWidth());
				setText(label);
			} else {
				setText(null);
			}

			if (hackBorders) {
				// Hacking borders to enhance the default one.
				if (item.isSelected) {
					setBorder(BORDER_SELECTED);
				}

				if (hasFocus) {
					setBorder(BORDER_FOCUSED);
				}
			}
		} else {
			setIcon(null);
			setText(null);
			setToolTipText(null);
		}

		return this;
	}

	protected String cutString(String string, int columnWidth) {
		StringBuilder sb = new StringBuilder(string);
		String str = sb.toString();

		Font font = getFont();
		int w = getFontMetrics(font).stringWidth(str);

		int i = 0;
		while (w > columnWidth) {
			str = "..." + sb.substring(i++);
			w = getFontMetrics(font).stringWidth(str);
		}

		return str;
	}
}
