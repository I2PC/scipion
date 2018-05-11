package xmipp.viewer.particlepicker.training.gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import xmipp.viewer.particlepicker.CtfInfo;
import xmipp.viewer.particlepicker.training.model.Mode;
import xmipp.viewer.particlepicker.training.model.SupervisedPickerMicrograph;


public class MicrographsTableModel extends AbstractTableModel {
	
	
	
	private List<SupervisedPickerMicrograph> micrographs;
	private String[] columns;
	private SupervisedPickerJFrame frame;

	public MicrographsTableModel(SupervisedPickerJFrame frame)
	{
		this.micrographs = frame.getParticlePicker().getMicrographs();
		CtfInfo ctfInfo = micrographs.get(0).getCtfInfo();
		if (ctfInfo == null || ctfInfo.defocusU == null)
			columns = new String[]{"", "Name", "Particles", "State"};
		else
			columns = new String[]{"", "Name", "Particles", "State", "DefocusU"};
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
		SupervisedPickerMicrograph m = micrographs.get(rowIndex);
		switch (columnIndex) {
			case 0: return rowIndex + 1;
			case 1: return m.getName();
			case 2:
				if(m.getStep() == Mode.Manual)
					return Integer.toString(m.getManualParticles().size());
				if(m.getStep() == Mode.Available)
					return "0";
				if(m.getStep() == Mode.Supervised)
					return String.format("%s + %s", m.getManualParticles().size(),
							             m.getAutomaticParticlesNumber(m.getThreshold()));
			case 3: return m.getState();
			case 4:
			    CtfInfo ctfInfo = m.getCtfInfo();
			    return ctfInfo.defocusU.toString();
			default:
				return null;
		}
	}
	
}
