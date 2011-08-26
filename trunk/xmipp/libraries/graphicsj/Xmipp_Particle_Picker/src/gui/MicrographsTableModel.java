package gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import model.Micrograph;
import model.MicrographFamilyData;

public class MicrographsTableModel extends AbstractTableModel {
	
	
	
	private List<Micrograph> micrographs;
	private String[] columns = new String[]{"Name", "Particles", "State"};
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
		if(c == 1)
			return frame.getFamily().getName();
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
		MicrographFamilyData mfd = m.getFamilyData(frame.getFamily()); 
		if(columnIndex == 1)
			return mfd.getManualParticles().size();
		if(columnIndex == 2)
			return mfd.getState();
		return null;
	}

}
