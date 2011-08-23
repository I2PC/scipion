package gui;

import java.awt.Frame;
import java.util.List;

import javax.swing.table.AbstractTableModel;

import model.Family;
import model.Micrograph;

public class MicrographsTableModel extends AbstractTableModel {
	
	
	
	private List<Micrograph> micrographs;
	private String[] columns = new String[]{"Name", "Particles"};
	private ParticlePickerJFrame frame;

	public MicrographsTableModel(ParticlePickerJFrame frame)
	{
		this.micrographs = frame.getParticlePicker().getMicrographs();
		this.frame = frame;
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
			return m.getFamilyData(frame.getFamily()).getParticles().size();
		return null;
	}

}
