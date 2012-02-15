package xmipp.viewer;

import xmipp.utils.DEBUG;
import xmipp.viewer.ImageItem;
import ij.ImagePlus;
import java.awt.Component;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import xmipp.viewer.gallery.models.AbstractXmippTableModel;
import xmipp.viewer.gallery.renderers.StrokeBorder;

import java.awt.Color;
import java.awt.Font;
import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.table.DefaultTableCellRenderer;

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
			// Loads image...
			ImagePlus imp = item.getImage();

			// ... and sets it.
			setEnabled(item.isEnabled);

			setIcon(new ImageIcon(imp.getImage()));

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

	private void normalize(ImagePlus image, AbstractXmippTableModel tableModel) {
		if (tableModel.isNormalizing()) {
			image.getProcessor().setMinAndMax(tableModel.getNormalizeMin(),
					tableModel.getNormalizeMax());
		} else {
			image.getProcessor().resetMinAndMax();
		}

		image.updateImage(); // Repaint
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
