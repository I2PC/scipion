package particlepicker.training.gui;

import java.awt.Color;
import java.awt.Component;

import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

/**
 * Used by {@link CLegendJDialog} to display colors on legend.
 * 
 * @author Airen
 *
 */
public class ColorRenderer  extends JLabel implements TableCellRenderer
{
	private Color color;

		 public ColorRenderer() {
		        setOpaque(true); //MUST do this for background to show up.
		 }

		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			color = (Color)value;
			setIcon(new ColorIcon(color));
			return this;
		}

}
