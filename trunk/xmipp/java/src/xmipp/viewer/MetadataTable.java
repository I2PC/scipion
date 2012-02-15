package xmipp.viewer;

import java.awt.Component;

import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;

public class MetadataTable extends MetadataGallery {

	public MetadataTable(String fn, int zoom) throws Exception {
		super(fn, zoom);
		cols = visibleLabels.size();
	}

	@Override
	public String getColumnName(int column) {
		return visibleLabels.get(column).getLabelName();
	}

	@Override
	public Class getColumnClass(int column) {
		ColumnInfo col = visibleLabels.get(column);
		if (col.render)
			return ImageItem.class;
		return Object.class;
	}

	@Override
	public int getIndex(int row, int col) {
		return row;
	}

	@Override
	public Object getValueAt(int row, int column) {
		try {
			ColumnInfo col = visibleLabels.get(column);
			if (col.render) {
				String key = getItemKey(row, col.getLabel());
				return createImageItem(row, col.getLabel(), col.getLabel(), key);
				// return super.getValueAt(row, col);
			}
			return md.getValueString(col.getLabel(), ids[row]);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return null;
	}

	@Override
	public void setColumns(int cols) {

	}

	@Override
	public void adjustColumn(int width) {

	}
	
	@Override
	/** Whether to display the labels */
	public void setRenderImages(boolean value) {
		if (renderLabels != value) {
			renderLabels = value;
			calculateCellSize();
			fireTableDataChanged();
		}
	}

	@Override
	protected void calculateCellSize() {
		if (renderLabels)
			super.calculateCellSize();
		else {
			int font_height;
			font_height = renderer.getFontMetrics(renderer.getFont())
					.getHeight();
			font_height += renderer.getIconTextGap(); // Adds the extra gap.
			cellDim.setSize(cellDim.width, font_height + 5);
		}
	}

	@Override
	public void setupTable(JTable table) {
		table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		table.setDefaultRenderer(ImageItem.class, renderer);
		table.setDefaultRenderer(Object.class, new TestRenderer());
	}

	@Override
	/** Return the column model to be used with this table model */
	public GalleryColumnModel createColumnModel() {
		return new MetadataColumnModel();
	}

	public class MetadataColumnModel extends GalleryColumnModel {
		public MetadataColumnModel() {
			super(0);
		}

		@Override
		public void adjustColumnsWidth(JTable table) {
			try {
				calculateCellSize();
				// String[] row = md.getRowValues(ids[0]);
				//table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
				int width = 0;
				for (int i = 0; i < visibleLabels.size(); ++i) {
					ColumnInfo col = visibleLabels.get(i);
					if (col.render) {
						width = cellDim.width;
						DEBUG.printMessage("is render -->");
					}
					else {
						TableCellRenderer rend = table.getCellRenderer(3, i);
						Component comp = rend.getTableCellRendererComponent(
								table, getValueAt(3, i), false, false, 0, 0);
						width = comp.getPreferredSize().width + 10;
					}
					getColumn(i).setPreferredWidth(width);
					DEBUG.printMessage(String.format("col: %d, width: %d", i, width));
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
}
