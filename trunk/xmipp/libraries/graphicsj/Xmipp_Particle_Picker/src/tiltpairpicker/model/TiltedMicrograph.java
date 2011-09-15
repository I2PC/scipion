package tiltpairpicker.model;

import ij.ImagePlus;

import java.io.File;
import java.util.List;

import javax.swing.Icon;

import trainingpicker.model.Constants;
import trainingpicker.model.Family;
import trainingpicker.model.FamilyState;
import trainingpicker.model.Micrograph;
import trainingpicker.model.MicrographFamilyData;
import trainingpicker.model.MicrographFamilyState;
import trainingpicker.model.Particle;
import trainingpicker.model.TrainingParticle;

public class TiltedMicrograph extends Micrograph{
	
	
	public TiltedMicrograph(String file) {
		super(file);

	}

	public TiltedParticle getParticle(int x, int y, int size) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean hasData() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void reset() {
		// TODO Auto-generated method stub
		
	}

}
