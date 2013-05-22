package xmipp.viewer.particlepicker.training.gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import xmipp.viewer.particlepicker.training.model.FamilyState;
import xmipp.viewer.particlepicker.training.model.MicrographFamilyData;
import xmipp.viewer.particlepicker.training.model.TrainingMicrograph;


public class MicrographsTableModel extends AbstractTableModel {
	
	
	
	private List<TrainingMicrograph> micrographs;
	private String[] columns = new String[]{"", "Name", "Particles", "State"};
	private TrainingPickerJFrame frame;

	public MicrographsTableModel(TrainingPickerJFrame frame)
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
		if(getParticlesPosition() == c)
			return frame.getFamily().getName();
		return columns[c];
	}

	@Override
	public int getRowCount() {
		return micrographs.size();
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		TrainingMicrograph m = micrographs.get(rowIndex);
		if(columnIndex == 0)
			return rowIndex + 1;
		if(columnIndex == 1)
			return m.getName();
		MicrographFamilyData mfd = m.getFamilyData(frame.getFamily()); 
		if(columnIndex == 2)
		{
			if(mfd.getStep() == FamilyState.Manual)
				return Integer.toString(mfd.getManualParticles().size());
			if(mfd.getStep() == FamilyState.Available)
				return "0";
			if(mfd.getStep() == FamilyState.Supervised)
				return String.format("%s + %s", mfd.getManualParticles().size(), mfd.getAutomaticParticlesCount(frame.getThreshold()));
		}
		if(columnIndex == 3)
			return mfd.getState();
		return null;
	}
	
	public int getParticlesPosition()
	{
		return 2;
	}

}
