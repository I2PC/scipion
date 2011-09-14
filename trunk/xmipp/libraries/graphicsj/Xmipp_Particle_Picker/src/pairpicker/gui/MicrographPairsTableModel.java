package pairpicker.gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import pairpicker.model.MicrographPair;
import picker.model.FamilyState;
import picker.model.Micrograph;
import picker.model.MicrographFamilyData;


public class MicrographPairsTableModel extends AbstractTableModel {
	
	
	
	private List<MicrographPair> micrographs;
	private String[] columns = new String[]{"", "Name", "Pair Name", "Particles"};
	private ParticlePairPickerJFrame frame;

	public MicrographPairsTableModel(ParticlePairPickerJFrame frame)
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
		MicrographPair m = micrographs.get(rowIndex);
		if(columnIndex == 0)
			return rowIndex + 1;
		if(columnIndex == 1)
			return m.getName();
		if(columnIndex == 2)
		{
			return m.getTiltedName();
		}
			
		
		return null;
	}
	
	public int getParticlesPosition()
	{
		return 2;
	}

}
