package xmipp.viewer.particlepicker;

import java.util.List;

import xmipp.viewer.particlepicker.training.model.Mode;

public class SingleParticlePicker extends ParticlePicker
{

	public SingleParticlePicker(String block, String selfile, String outputdir, String fname, Mode mode)
	{
		super(block, selfile, outputdir, fname, mode);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void loadEmptyMicrographs()
	{
		// TODO Auto-generated method stub

	}

	@Override
	public List<? extends Micrograph> getMicrographs()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void saveData(Micrograph m)
	{
		// TODO Auto-generated method stub

	}

	@Override
	public void saveConfig()
	{
		// TODO Auto-generated method stub

	}

	@Override
	public Micrograph getMicrograph()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setMicrograph(Micrograph m)
	{
		// TODO Auto-generated method stub

	}

	@Override
	public boolean isValidSize(int size)
	{
		// TODO Auto-generated method stub
		return false;
	}

}
