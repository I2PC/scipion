package tiltpairpicker.gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;


import tiltpairpicker.model.UntiltedMicrograph;
import trainingpicker.model.FamilyState;
import trainingpicker.model.TrainingMicrograph;
import trainingpicker.model.MicrographFamilyData;


public class MicrographPairsTableModel extends AbstractTableModel {
	
	
	
	private List<UntiltedMicrograph> micrographs;
	private String[] columns = new String[]{"", "Name", "Pair Name", "Particles"};
	private TiltPairPickerJFrame frame;

	public MicrographPairsTableModel(TiltPairPickerJFrame frame)
	{
		this.micrographs = frame.getParticlePairPicker().getMicrographs();
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
		UntiltedMicrograph m = micrographs.get(rowIndex);
		if(columnIndex == 0)
			return rowIndex + 1;
		if(columnIndex == 1)
			return m.getName();
		if(columnIndex == 2)
			return m.getTiltedMicrograph().getName();
		return null;
	}
	
	public int getParticlesPosition()
	{
		return 2;
	}

}
