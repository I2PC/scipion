package gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import model.Micrograph;

public class MicrographsTableModel extends AbstractTableModel {
	
	
	
	private List<Micrograph> micrographs;
	private String[] columns = new String[]{"Name", "Particles"};

	public MicrographsTableModel(List<Micrograph> micrographs)
	{
		this.micrographs = micrographs;
	}
	
	@Override
	public int getColumnCount() {
		return columns.length;
	}
	@Override
	public String getColumnName(int c)
	{
		return columns[c];
	}

	@Override
	public int getRowCount() {
		return micrographs.size();
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		Micrograph m = micrographs.get(rowIndex);
		if(columnIndex == 0)
			return m.getName();
		if(columnIndex == 1)
			return m.getParticles().size();
		return null;
	}

}
