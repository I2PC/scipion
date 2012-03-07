package xmipp.viewer;

import java.awt.Component;
import java.util.ArrayList;

import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.DEBUG;

public class MetadataTable extends MetadataGallery {
	private static final long serialVersionUID = 1L;

	public MetadataTable(GalleryData data) throws Exception {
		super(data);
		cols = visibleLabels.size();
	}

	@Override
	public String getColumnName(int column) {
		return visibleLabels.get(column).getLabelName();
	}

	@Override
	public Class getColumnClass(int column) {
		ColumnInfo ci = visibleLabels.get(column);
		if (ci.render)
			return ImageItem.class;
		else if (ci.getLabel() == MDLabel.MDL_ENABLED)
			return Boolean.class;
		try {
			return MetaData.getLabelClass(ci.getLabel());
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public int getIndex(int row, int col) {
		return row;
	}

	@Override
	public Object getValueAt(int row, int column) {
		// DEBUG.printMessage(String.format("MetadataTable.getValueAt(%d, %d)",
		// row, column));
		try {
			ColumnInfo ci = visibleLabels.get(column);
			if (ci.render) {
				String key = getItemKey(row, ci.getLabel());
				return createImageItem(row, ci.getLabel(), ci.getLabel(), key);
				// return super.getValueAt(row, col);
			}
			int label = ci.getLabel();
			long id = data.ids[row];
			int type = MetaData.getLabelType(label);
			MetaData md = data.md;
			switch (type) {
			case MetaData.LABEL_INT:
				int value = md.getValueInt(label, id);
				// treat special case of MDL_ENABLED
				if (label == MDLabel.MDL_ENABLED)
					return (value > 0);
				return value;
			case MetaData.LABEL_BOOL:
				return md.getValueBoolean(label, id);
			case MetaData.LABEL_FLOAT:
			case MetaData.LABEL_DOUBLE:
				return md.getValueDouble(label, id);
			case MetaData.LABEL_LONG:
				return md.getValueLong(label, id);
			case MetaData.LABEL_STRING:
			case MetaData.LABEL_VECTOR:
			case MetaData.LABEL_VECTOR_LONG:
				return md.getValueString(ci.getLabel(), data.ids[row]);

			}
			return null;

		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;
	}

	@Override
	public void setValueAt(Object value, int row, int column) {
		try {
			ColumnInfo ci = visibleLabels.get(column);
			if (!ci.render) {
				int label = ci.getLabel();
				long id = data.ids[row];
				int type = MetaData.getLabelType(label);
				MetaData md = data.md;
				switch (type) {
				case MetaData.LABEL_BOOL:
					md.setValueBoolean(label, (Boolean) value, id);
					break;
				case MetaData.LABEL_INT:
					if (label == MDLabel.MDL_ENABLED) {
						md.setEnabled((Boolean)value, id);
					} else
						md.setValueInt(label, ((Integer) value).intValue(), id);
					break;
				case MetaData.LABEL_FLOAT:
				case MetaData.LABEL_DOUBLE:
					md.setValueDouble(label, ((Double) value).doubleValue(), id);
					break;
				case MetaData.LABEL_LONG:
					md.setValueInt(label, ((Integer) value).intValue(), id);
					break;
				case MetaData.LABEL_STRING:
				case MetaData.LABEL_VECTOR:
				case MetaData.LABEL_VECTOR_LONG:
					// TODO: Implement a better editor for vectors
					md.setValueString(label, value.toString(), id);
					break;
				}
				fireTableRowsUpdated(row, row);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

	}// function setValueAt

	@Override
	public boolean isCellEditable(int row, int column) {
		ColumnInfo ci = visibleLabels.get(column);
		return ci.allowEdit && !ci.render;
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
		boolean changed = false;
		for (ColumnInfo ci : data.labels)
			if (ci.allowRender && ci.render != value) {
				ci.render = value;
				changed = true;
			}
		if (changed) {
			data.globalRender = value;
			calculateCellSize();
			fireTableDataChanged();
		}
	}

	@Override
	protected void calculateCellSize() {
		DEBUG.printMessage(String.format("MetadataTable:calculateSize"));
		if (data.globalRender) {
			super.calculateCellSize();
			DEBUG.printMessage(String.format(
					"MetadataTable:calculateSize w:%d, h:%d", cellDim.width,
					cellDim.height));

		} else {
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
		table.setDefaultRenderer(Double.class, new TestRenderer());
		table.setAutoCreateRowSorter(true);
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
				// String[] row = md.getRowValues(data.ids[0]);
				// table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
				int width = 0;
				for (int i = 0; i < visibleLabels.size(); ++i) {
					ColumnInfo col = visibleLabels.get(i);
					if (col.render) {
						width = cellDim.width;
						DEBUG.printMessage("is render -->");
					} else {
						TableCellRenderer rend = table.getCellRenderer(0, i);
						Component comp = rend.getTableCellRendererComponent(
								table, getValueAt(0, i), false, false, 0, 0);
						width = comp.getPreferredSize().width + 10;
					}
					getColumn(i).setPreferredWidth(width);
//					DEBUG.printMessage(String.format("col: %d, width: %d", i,
//							width));
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
}
