package gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import model.Micrograph;
import model.MicrographFamilyData;
import model.FamilyState;

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
		{
			if(mfd.getStep() == FamilyState.Manual)
				return Integer.toString(mfd.getManualParticles().size());
			if(mfd.getStep() == FamilyState.Available)
				return "0";
			if(mfd.getStep() == FamilyState.Supervised)
				return String.format("%s + %s", mfd.getManualParticles().size(), mfd.getAutomaticParticlesCount());
		}
		if(columnIndex == 2)
			return mfd.getState();
		return null;
	}

}
