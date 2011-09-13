package model;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;

import xmipp.MDLabel;
import xmipp.MetaData;


public class ReviewParticlePicker extends ParticlePicker {
	
	private String reviewfile;
	private Family reviewfamily;
	
	public String getReviewFile() {
		return reviewfile;
	}
	

	public ReviewParticlePicker(String selfile, String outputdir, String reviewfile) {
		super(selfile, outputdir, FamilyState.Review);
		if (!new File(reviewfile).exists())
			throw new IllegalArgumentException(
					Constants.getNoSuchFieldValueMsg("review file", reviewfile));
		this.reviewfile = reviewfile;
		String[] parts = reviewfile.split(File.separator);
		String familyname = parts[parts.length - 1].split("_")[0];
		this.reviewfamily = getFamily(familyname);
		if (reviewfamily == null)
			throw new IllegalArgumentException(
					Constants.getNoSuchFieldValueMsg("family", familyname));
		loadMicrographs();
	}
	
	@Override
	public void persistMicrographs() {
		try {
			MetaData md;
			MicrographFamilyData mfd;
			long id;
			new File(reviewfile).delete();// to ensure I clean file and can
											// append later
			for (Micrograph m : micrographs) {

				mfd = m.getFamilyData(reviewfamily);
				if (!mfd.isEmpty()) {
					md = new MetaData();
					for (Particle p : mfd.getParticles()) {
						id = md.addObject();
						md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
						md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
					}
					md.writeBlock("mic_" + m.getName() + "@" + reviewfile);
				}
			}
		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}
	}

	@Override
	public void loadMicrographs() {
		micrographs.clear();
		Micrograph micrograph;
		String ctf = null, filename;
		try {
			MetaData md = new MetaData(getMicrographsSelFile());
			boolean existsctf = md.containsLabel(MDLabel.MDL_PSD_ENHANCED);
			long[] ids = md.findObjects();
			for (long id : ids) {

				filename = getMicrographPath(md.getValueString(
						MDLabel.MDL_IMAGE, id));
				if (existsctf)
					ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new Micrograph(filename, ctf, families, getMode());
				micrographs.add(micrograph);
			}
			if (micrographs.size() == 0)
				throw new IllegalArgumentException(String.format(
						"No micrographs specified on %s", getMicrographsSelFile()));
			int x, y;
			Particle particle;
			List<String> blocks = Arrays.asList(MetaData
					.getBlocksInMetaDataFile(reviewfile));
			String block;
			for (Micrograph m : micrographs) {
				block = "mic_" + m.getName();
				if (blocks.contains(block)) {
					md = new MetaData(block + "@" + reviewfile);

					ids = md.findObjects();
					for (long id : ids) {

						x = md.getValueInt(MDLabel.MDL_XINT, id);
						y = md.getValueInt(MDLabel.MDL_YINT, id);
						particle = new Particle(x, y, reviewfamily, m);
						m.addManualParticle(particle);
					}
				}
			}

		} catch (Exception e) {
			getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e);
		}

	}
	



}
