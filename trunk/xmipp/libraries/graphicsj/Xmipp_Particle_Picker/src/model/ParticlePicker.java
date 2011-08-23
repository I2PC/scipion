package model;

import gui.ParticlePickerJFrame;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import javax.swing.JOptionPane;

import xmipp.MDLabel;
import xmipp.MetaData;
import xmipp.Program;

public class ParticlePicker {

	private List<Family> families;
	private List<Micrograph> micrographs;
	private static ParticlePicker ppdata;

	private ParticlePicker() {
		this.families = new ArrayList<Family>();
		this.micrographs = new ArrayList<Micrograph>();
		loadFamilyData();
		loadMicrographsData();
	}

	public void train() {

		String args;
		for (Micrograph micrograph : micrographs) {
			for (Family f : families) {
				args = String
						.format("-i %s --particleSize %s --model %s --outputRoot %s --mode train %s",
								micrograph.getFilename(),// -i
								f.getSize(), // --particleSize
								ExecutionEnvironment.getOutputPath(f.getName()),// --model
								ExecutionEnvironment.getOutputPath(micrograph
										.getName()), // --outputRoot
								f.getName() + "@" + micrograph.getOutputFName());// train
																					// parameter
				if (ExecutionEnvironment.isFastMode())
					args += " --fast";
				if (ExecutionEnvironment.isIncore())
					args += " --in_core";
				executeProgram("xmipp_micrograph_automatic_picking", args);
			}
		}
	}

	private void executeProgram(String name, String args) {
		System.out.println(args);
		try {
			Program.runByName(name, args);
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
					e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public int getParticlesNumber() {
		int count = 0;
		for (Micrograph m : micrographs)
			count += m.getParticles().size();
		return count;
	}

	public void classify() {
		String args;
		for (Micrograph micrograph : micrographs) {
			for (Family f : families) {
				args = String
						.format("-i %s --particleSize %s --model %s --outputRoot %s --mode try --thr %s",
								micrograph.getFilename(),// -i
								f.getSize(), // --particleSize
								ExecutionEnvironment.getOutputPath(f.getName()),// --model
								ExecutionEnvironment.getOutputPath(micrograph
										.getName()),// --outputRoot
								ExecutionEnvironment.getThreads()// --thr
						);

				if (ExecutionEnvironment.isFastMode())
					args += " --fast";
				if (ExecutionEnvironment.isIncore())
					args += " --in_core";
				executeProgram("xmipp_micrograph_automatic_picking", args);
			}
		}
	}

	public static ParticlePicker getInstance() {
		if (ppdata == null)
			ppdata = new ParticlePicker();
		return ppdata;
	}

	public List<Family> getFamilies() {
		return families;
	}

	public List<Micrograph> getMicrographs() {
		return micrographs;
	}

	public void saveFamilyData() {
		long id;
		String filename = Family.getOFilename();
		try {
			MetaData md = new MetaData();
			for (Family f : families) {
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, f.getName(), id);
				md.setValueInt(MDLabel.MDL_PICKING_COLOR,
						f.getColor().getRGB(), id);
				md.setValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, f.getSize(),
						id);
			}
			md.write(filename);
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
					e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void loadFamilyData() {
		families.clear();
		String filename = Family.getOFilename();
		if (!new File(filename).exists()) {
			families.add(Family.getDefaultFamily());
			return;
		}

		Family family;
		int rgb, size;
		String gname;
		try {
			MetaData md = new MetaData(filename);
			long[] ids = md.findObjects();
			for (long id : ids) {
				gname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				rgb = md.getValueInt(MDLabel.MDL_PICKING_COLOR, id);
				size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
				family = new Family(gname, new Color(rgb), size);
				families.add(family);
			}
			if (families.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No families specified on %s", filename));
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
					e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}

	public Family getFamily(String name) {
		for (Family f : getFamilies())
			if (f.getName().equalsIgnoreCase(name))
				return f;
		return null;
	}

	public boolean existsFamilyName(String name) {
		return getFamily(name) != null;
	}

	public void loadMicrographsData() {
		String xmd = Micrograph.getIFilename();
		micrographs.clear();
		Micrograph micrograph;
		String ctf, filename;
		try {
			MetaData md = new MetaData(xmd);
			int count = 1;
			long[] ids = md.findObjects();
			for (long id : ids) {

				filename = ExecutionEnvironment.getMicrographPath(md
						.getValueString(MDLabel.MDL_IMAGE, id));
				ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new Micrograph(filename, ctf);
				loadParticles(micrograph);
				micrographs.add(micrograph);
				count++;
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No micrographs specified on %s", xmd));

		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
					e);
		}

	}

	public void saveParticles() {
		long id;
		try {
			int count;
			MetaData md;
			String block = null;
			for (Micrograph m : micrographs) {
				count = 0;
				for (Family f : families) {
					md = new MetaData();
					for (Particle p : m.getParticles()) {
						if (f.equals(p.getFamily())) {
							id = md.addObject();
							md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
							md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
						}
					}
					block = f.getName() + "@" + m.getOutputFName();
					if (count == 0)
						md.write(block);
					else
						md.writeBlock(block);
					count++;
				}
			}

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void loadParticles(Micrograph micrograph) {
		try {
			int x, y;
			Particle particle;
			MetaData md;
			if (!new File(micrograph.getOutputFName()).exists())
				return;
			for (Family f : families) {
				md = new MetaData(f.getName() + "@"
						+ micrograph.getOutputFName());
				long[] ids = md.findObjects();
				for (long id : ids) {

					x = md.getValueInt(MDLabel.MDL_XINT, id);
					y = md.getValueInt(MDLabel.MDL_YINT, id);
					particle = new Particle(x, y, f, micrograph);
					micrograph.addParticle(particle);
				}
			}
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
					e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

	public void loadAutomaticParticles(Micrograph micrograph) {
		try {
			int x, y;
			double cost;
			Particle particle;
			MetaData md;
			Family f;
			String filename = micrograph.getAutoOutputFName();
			if (!new File(filename).exists())
				return;

			String[] blocks = MetaData.getBlocksInMetaDataFile(filename);
			for (String blockname : blocks) {
				md = new MetaData(blockname + "@" + filename);
				long[] ids = md.findObjects();
				for (long id : ids) {

					x = md.getValueInt(MDLabel.MDL_XINT, id);
					y = md.getValueInt(MDLabel.MDL_YINT, id);
					cost = md.getValueDouble(MDLabel.MDL_COST, id);
					f = getFamily(blockname);
					if (f == null)
						throw new IllegalArgumentException("Unknown family "
								+ blockname);
					particle = new Particle(x, y, f, micrograph, true, cost);
					micrograph.addParticle(particle);
				}
			}
		} catch (Exception e) {
			ExecutionEnvironment.getLogger().log(Level.SEVERE, e.getMessage(),
					e);
			throw new IllegalArgumentException(e.getMessage());
		}

	}

}
