package trainingpicker.gui;

import java.awt.Component;

import javax.swing.JList;
import javax.swing.ListCellRenderer;

import trainingpicker.model.TrainingParticle;



public class ParticleCellRenderer implements ListCellRenderer {
	
	@Override
	public Component getListCellRendererComponent(JList list, Object o,
			int arg2, boolean arg3, boolean arg4) {
		ParticleCanvas c = (ParticleCanvas)o;
		return c;
	}
	

}
