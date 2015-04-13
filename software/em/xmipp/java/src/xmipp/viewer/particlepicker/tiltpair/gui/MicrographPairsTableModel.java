package xmipp.viewer.particlepicker.tiltpair.gui;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import xmipp.viewer.particlepicker.tiltpair.model.UntiltedMicrograph;


public class MicrographPairsTableModel extends AbstractTableModel {
	
	
	
	private List<UntiltedMicrograph> micrographs;
	private String[] columns = new String[]{"", "Name", "Pair Name", "Particles", "Tilt Angle"};

	public MicrographPairsTableModel(TiltPairPickerJFrame frame)
	{
		this.micrographs = frame.getParticlePicker().getMicrographs();
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
		if(columnIndex == 3)
			return m.getParticles().size();
		if(columnIndex == 4)
			return String.format("%.2f", m.getTiltAngle());
		return null;
	}
	
	public int getParticlesPosition()
	{
		return 2;
	}

}
