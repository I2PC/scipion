package xmipp.viewer.particlepicker.training.gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.SingleParticlePickerMicrograph;


public class MicrographsTableModel extends AbstractTableModel {
	
	
	
	private List<SingleParticlePickerMicrograph> micrographs;
	private String[] columns = new String[]{"", "Name", "Particles", "State"};
	private SingleParticlePickerJFrame frame;

	public MicrographsTableModel(SingleParticlePickerJFrame frame)
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
		SingleParticlePickerMicrograph m = micrographs.get(rowIndex);
		if(columnIndex == 0)
			return rowIndex + 1;
		if(columnIndex == 1)
			return m.getName();
		if(columnIndex == 2)
		{
			if(m.getStep() == Mode.Manual)
				return Integer.toString(m.getManualParticles().size());
			if(m.getStep() == Mode.Available)
				return "0";
			if(m.getStep() == Mode.Supervised)
				return String.format("%s + %s", m.getManualParticles().size(), m.getAutomaticParticlesNumber(m.getThreshold()));
		}
		if(columnIndex == 3)
			return m.getState();

		return null;
	}
	
	

}
