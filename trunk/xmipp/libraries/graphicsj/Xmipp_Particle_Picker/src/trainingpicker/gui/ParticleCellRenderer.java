package trainingpicker.gui;

import java.awt.Component;

import javax.swing.JList;
import javax.swing.JTable;
import javax.swing.ListCellRenderer;
import javax.swing.table.TableCellRenderer;

import trainingpicker.model.TrainingParticle;



public class ParticleCellRenderer implements TableCellRenderer {

	@Override
	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
	{
		ParticleCanvas c = (ParticleCanvas)value;
		return c;
	}
	

	

}
