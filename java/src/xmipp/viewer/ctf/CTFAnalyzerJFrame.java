package xmipp.viewer.ctf;

import ij.ImagePlus;
import ij.gui.Line;
import ij.gui.ProfilePlot;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTabbedPane;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import xmipp.jni.CTFDescription;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.XmippDialog;
import xmipp.utils.XmippFileChooser;
import xmipp.utils.XmippLabel;
import xmipp.utils.XmippWindowUtil;

public class CTFAnalyzerJFrame extends JFrame implements ActionListener
{

	private ImagePlus imp;
	private String ctffile;
	private JPanel graphicpn;
	private JPanel imagepn;
	private CTFDescription ctfmodel;
	private int samples;
	private double[] xvalues;
	private XYSeriesCollection collection;
	private ChartPanel radialchartpn;
	private JPanel actionspn;
	private JRadioButton ctfrb;
	private JRadioButton psdrb;
	private GridBagConstraints constraints;
	private JCheckBox theorethicalpsdchb;
	private JCheckBox envelopepsdchb;
	private JCheckBox noisechb;
	private CTFAnalyzerImagePane imageprofilepn;
	private JTabbedPane tabspn;
	private String psdfile;
	private ChartPanel avgchartpn;
	double[] ctf_avgplot, theorethicalpsd_avgplot, envelope_avgplot, bgnoise_avgplot, psdprofile_avgplot, difference_avgplot;
	double[] ctfplot, theorethicalpsdplot, envelopeplot, bgnoiseplot, psdprofileplot, difference;
	private JButton exportbt;
	private JButton exportavgbt;
	private XmippFileChooser fc;
	private JCheckBox differencechb;

	public CTFAnalyzerJFrame(ImagePlus imp, String ctffile, String psdfile)
	{
		try
		{
			setTitle("CTF Analyzer");
			this.imp = imp;
			this.ctffile = ctffile;
			this.psdfile = psdfile;

			ctfmodel = new CTFDescription(ctffile);

			samples = imp.getWidth() / 2;
			xvalues = getXValues(samples, getSamplingRate(ctffile));
			fc = new XmippFileChooser();

			constraints = new GridBagConstraints();
			constraints.insets = new Insets(0, 5, 0, 5);
			constraints.anchor = GridBagConstraints.NORTHWEST;
			// constraints.fill = GridBagConstraints.HORIZONTAL;
			setLayout(new GridBagLayout());
			JPanel contentpn = new JPanel(new GridBagLayout());
			contentpn.setBorder(BorderFactory.createEtchedBorder());
			add(contentpn, XmippWindowUtil.getConstraints(constraints, 0, 0));
			imagepn = new JPanel(new GridBagLayout());
			imagepn.setBorder(BorderFactory.createTitledBorder("PSD Image"));
			imageprofilepn = new CTFAnalyzerImagePane(imp, this);
			imagepn.add(imageprofilepn, XmippWindowUtil.getConstraints(constraints, 0, 0));
			contentpn.add(imagepn, XmippWindowUtil.getConstraints(constraints, 0, 0, 1));

			initGraphicPanes();
			contentpn.add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 1, 1, 1));
			contentpn.add(graphicpn, XmippWindowUtil.getConstraints(constraints, 1, 0, 1, 3));

			JPanel actionspn = new JPanel();
			exportbt = XmippWindowUtil.getTextButton("Export Graphics", this);
			actionspn.add(exportbt, XmippWindowUtil.getConstraints(constraints, 0, 1));
//			exportavgbt = XmippWindowUtil.getTextButton("Export Average Graphics", this);
//			actionspn.add(exportavgbt, XmippWindowUtil.getConstraints(constraints, 1, 1));
			add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 1));
			enableDisplay();
			pack();
			setVisible(true);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
	}

	private void initGraphicPanes()
	{
		actionspn = new JPanel(new GridBagLayout());
		actionspn.setBorder(BorderFactory.createTitledBorder("Display"));
		ctfrb = new JRadioButton(getCTFLabel());
		ctfrb.addActionListener(new FillGraphicsActionListener());
		ctfrb.addActionListener(new EnableDisplayActionListener());

		psdrb = new JRadioButton(getPSDProfileLabel());
		psdrb.addActionListener(new FillGraphicsActionListener());
		psdrb.addActionListener(new EnableDisplayActionListener());
		ButtonGroup drawbg = new ButtonGroup();
		drawbg.add(ctfrb);
		drawbg.add(psdrb);
		ctfrb.setSelected(true);
		actionspn.add(ctfrb, XmippWindowUtil.getConstraints(constraints, 0, 0));
		actionspn.add(psdrb, XmippWindowUtil.getConstraints(constraints, 0, 1));

		constraints.insets = new Insets(0, 20, 0, 5);
		theorethicalpsdchb = new JCheckBox(getTheorethicalPSDLabel());
		theorethicalpsdchb.addActionListener(new FillGraphicsActionListener());
		envelopepsdchb = new JCheckBox(getEnvelopeLabel());
		envelopepsdchb.addActionListener(new FillGraphicsActionListener());
		noisechb = new JCheckBox(getBGNoiseLabel());
		noisechb.addActionListener(new FillGraphicsActionListener());
		differencechb = new JCheckBox(getDifferenceLabel());
		differencechb.addActionListener(new FillGraphicsActionListener());

		actionspn.add(theorethicalpsdchb, XmippWindowUtil.getConstraints(constraints, 0, 2));
		actionspn.add(envelopepsdchb, XmippWindowUtil.getConstraints(constraints, 0, 3));
		actionspn.add(noisechb, XmippWindowUtil.getConstraints(constraints, 0, 4));
		actionspn.add(differencechb, XmippWindowUtil.getConstraints(constraints, 0, 5));
		
		graphicpn = new JPanel();
		tabspn = new JTabbedPane();
		JPanel tab1pn = new JPanel();
		JPanel tab2pn = new JPanel();
		tabspn.addTab("Radial Profile", tab1pn);
		tabspn.addTab("Radial Average", tab2pn);
		graphicpn.add(tabspn);
		JFreeChart chart = ChartFactory
				.createXYLineChart("", XmippLabel.LABEL_SAMPLING, showCTF() ? getCTFLabel() : getPSDProfileLabel(), collection, PlotOrientation.VERTICAL, true, true, false);

		radialchartpn = createChartPanel(chart);
		JFreeChart avgchart = ChartFactory
				.createXYLineChart("", XmippLabel.LABEL_SAMPLING, showCTF() ? getCTFLabel() : getPSDProfileLabel(), collection, PlotOrientation.VERTICAL, true, true, false);

		avgchartpn = createChartPanel(avgchart);
		fillGraphics();
		tab1pn.add(radialchartpn);
		tab2pn.add(avgchartpn);
	}

	

	private static ChartPanel createChartPanel(JFreeChart chart)
	{
		XYPlot plot = (XYPlot) chart.getPlot();
		plot.setDomainPannable(true);
		plot.setRangePannable(true);

		List<Integer> list = Arrays.asList(new Integer[] { new Integer(0), new Integer(1) });
		plot.mapDatasetToDomainAxes(0, list);
		plot.mapDatasetToRangeAxes(0, list);
		ChartUtilities.applyCurrentTheme(chart);

		return new ChartPanel(chart);
	}

	public boolean isRadialProfile()
	{
		return tabspn.getSelectedIndex() == 0;
	}

	void fillGraphics()
	{
		fillRadialGraphics();
		fillAvgGraphics();
	}

	private void fillAvgGraphics()
	{
		if (psdprofile_avgplot == null)// avg data does not change, calculated
										// only once
		{
			int x0 = imageprofilepn.getX0();
			int y0 = imageprofilepn.getY0();
			// // Get profile.
			Line line;
			double[] avgprofileplot = new double[samples];
			double angle;
			int x1, y1;
			for (int i = 0; i < 360; i++)
			{
				angle = Math.toRadians(-i);

				// Make line roi.
				x1 = x0 + (int) (samples * Math.cos(angle));
				y1 = y0 + (int) (samples * Math.sin(angle));

				line = new Line(x0, y0, x1, y1);
				imp.setRoi(line);

				// Get profile.
				ProfilePlot profilePlot = new ProfilePlot(imp);
				psdprofile_avgplot = profilePlot.getProfile();

				// Total summatory.
				sum(avgprofileplot, psdprofile_avgplot);
			}

			// Get profile average.
			for (int i = 0; i < avgprofileplot.length; i++)
				avgprofileplot[i] /= 360;
			psdprofile_avgplot = avgprofileplot;

			// Get rest of profiles average.
			ctfmodel.CTFAverageProfiles(samples);

			bgnoise_avgplot = ctfmodel.avgprofiles[CTFDescription.BACKGROUND_NOISE];
			envelope_avgplot = ctfmodel.avgprofiles[CTFDescription.ENVELOPE];
			theorethicalpsd_avgplot = ctfmodel.avgprofiles[CTFDescription.PSD];
			ctf_avgplot = ctfmodel.avgprofiles[CTFDescription.CTF];
			difference_avgplot = new double[bgnoiseplot.length];
			for (int i = 0; i < xvalues.length; i++)
				difference_avgplot[i] = psdprofileplot[i] - bgnoiseplot[i];
		}
		XYSeriesCollection collection = getXYSeriesCollection(psdprofile_avgplot, ctf_avgplot, theorethicalpsd_avgplot, envelope_avgplot, bgnoise_avgplot);
		XYPlot plot = ((XYPlot) avgchartpn.getChart().getPlot());
		plot.setDataset(0, collection);
		plot.getRangeAxis().setLabel(showCTF() ? getCTFLabel() : getPSDProfileLabel());
	}

	private void fillRadialGraphics()
	{
		

		int x0 = imageprofilepn.getX0();
		int y0 = imageprofilepn.getY0();
		// // Get profile.
		Line line;

		line = new Line(x0, y0, imageprofilepn.getX1(), imageprofilepn.getY1());
		imp.setRoi(line);
		psdprofileplot = new ProfilePlot(imp).getProfile();

		ctfmodel.CTFProfile(imageprofilepn.getProfileangle(), samples);
		bgnoiseplot = ctfmodel.profiles[CTFDescription.BACKGROUND_NOISE];
		envelopeplot = ctfmodel.profiles[CTFDescription.ENVELOPE];
		theorethicalpsdplot = ctfmodel.profiles[CTFDescription.PSD];
		ctfplot = ctfmodel.profiles[CTFDescription.CTF];
		difference = new double[bgnoiseplot.length];
		for (int i = 0; i < xvalues.length; i++)
			difference[i] = psdprofileplot[i] - bgnoiseplot[i];

		XYSeriesCollection collection = getXYSeriesCollection(psdprofileplot, ctfplot, theorethicalpsdplot, envelopeplot, bgnoiseplot);

		XYPlot plot = ((XYPlot) radialchartpn.getChart().getPlot());
		plot.setDataset(0, collection);
		plot.getRangeAxis().setLabel(showCTF() ? getCTFLabel(): getPSDProfileLabel());
	}

	private static void sum(double a[], double b[])
	{
		for (int i = 0; i < a.length; i++)
		{
			a[i] += b[i];
		}
	}

	private XYSeriesCollection getXYSeriesCollection(double[] psdprofileplot, double[] ctfplot, double[] theorethicalpsdplot, double[] envelopeplot,
			double[] bgnoiseplot)
	{
		double max = Double.MAX_VALUE;
		double min = -1;
		collection = new XYSeriesCollection();

		if (showCTF())
			collection.addSeries(createSeries(XmippLabel.CB_PLOT_CTF, xvalues, ctfplot, max, min));
		else if (showPSD())
		{
			collection.addSeries(createSeries(getPSDProfileLabel(), xvalues, psdprofileplot, max, min));

			/////////////////some max value is established and values displayed are filtered accordingly////////////////////////////////
			// min = max;
			max = -max;
			for (int i = 0; i < xvalues.length; i++)
			{
				max = Math.max(max, psdprofileplot[i]);
				// min = Math.min(min, profile[i]);
			}
			max += Math.abs(max) * 0.1;
			// min -= Math.abs(min) * 0.1;
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			if (showBGNoise())
			{
				collection.addSeries(createSeries(getBGNoiseLabel(), xvalues, bgnoiseplot, max, min));
				
			}
			if(showDifference())
			{
				
				collection.addSeries(createSeries(getDifferenceLabel(), xvalues, difference, max, min));
			}
			if (showEnvelope())
				collection.addSeries(createSeries(getEnvelopeLabel(), xvalues, envelopeplot, max, min));

			if (showTheorethicalPSD())
				collection.addSeries(createSeries(getTheorethicalPSDLabel(), xvalues, theorethicalpsdplot, max, min));
		}
		return collection;

	}

	private static XYSeries createSeries(String name, double[] xs, double[] values, double max, double min)
	{
		XYSeries series = new XYSeries(name);
		for (int i = 0; i < xs.length; i++)
		{
			if (values[i] > min && values[i] < max)
				series.add(xs[i], values[i]);
		}

		return series;
	}

	private boolean showPSD()
	{
		return psdrb.isSelected();
	}

	private boolean showTheorethicalPSD()
	{
		return theorethicalpsdchb.isSelected();
	}

	private boolean showEnvelope()
	{
		return envelopepsdchb.isSelected();
	}

	private boolean showBGNoise()
	{
		return noisechb.isSelected();
	}

	private boolean showCTF()
	{
		return ctfrb.isSelected();
	}
	
	private boolean showDifference()
	{
		return differencechb.isSelected();
	}

	private static double[] getXValues(int length, double Tm)
	{
		double values[] = new double[length];
		double factor = 1 / (length * 2 * Tm);

		for (int i = 0; i < length; i++)
		{
			values[i] = (double) i * factor;
		}

		return values;
	}

	private double getSamplingRate(String CTFFilename)
	{
		double samplingRate = 1.0;
		double downsampling = 1.0;

		try
		{
			MetaData md = new MetaData(CTFFilename);
			samplingRate = md.getValueDouble(MDLabel.MDL_CTF_SAMPLING_RATE, md.firstObject());
			if (md.containsLabel(MDLabel.MDL_CTF_DOWNSAMPLE_PERFORMED))
				downsampling = md.getValueDouble(MDLabel.MDL_CTF_DOWNSAMPLE_PERFORMED, md.firstObject());
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}

		return samplingRate * downsampling;
	}

	class FillGraphicsActionListener implements ActionListener
	{

		@Override
		public void actionPerformed(ActionEvent e)
		{

			fillGraphics();
			envelopepsdchb.setEnabled(psdrb.isSelected());

		}

	}
	
	class EnableDisplayActionListener implements ActionListener
	{

		@Override
		public void actionPerformed(ActionEvent e)
		{
			enableDisplay();
		}

	}
	
	private void enableDisplay()
	{
		theorethicalpsdchb.setEnabled(showPSD());
		envelopepsdchb.setEnabled(showPSD());
		noisechb.setEnabled(showPSD());
		differencechb.setEnabled(showPSD());
	
	}

	@Override
	public void actionPerformed(ActionEvent e)
	{
		int result = fc.showSaveDialog(this);
		boolean export = result == XmippFileChooser.APPROVE_OPTION;
		if (export)
		{
			File file = fc.getSelectedFile();
			if (file != null)
			{
				
				boolean saved = export(file.getAbsolutePath());
				if (saved)
					XmippDialog.showInfo(this, "File saved sucesfully.");
			}
		}

	}
	
    private boolean export(String path) {


        try {
            FileWriter fstream = new FileWriter(path);
            BufferedWriter out = new BufferedWriter(fstream);

            String hformat = "%40s%36s AVG";
            String vformat = "%40.4f%40.4f";
            // Header.
            out.write(String.format("%40s", "X"));
            if(showCTF())
            	out.write(String.format(hformat, getCTFLabel(), getCTFLabel()));
            if(showPSD())
            	out.write(String.format(hformat, getPSDProfileLabel(), getPSDProfileLabel()));
            if(showTheorethicalPSD())
            	out.write(String.format(hformat, getTheorethicalPSDLabel(), getTheorethicalPSDLabel()));
            if(showEnvelope())
            	out.write(String.format(hformat, getEnvelopeLabel(), getEnvelopeLabel()));
            if(showBGNoise())
            	out.write(String.format(hformat, getBGNoiseLabel(), getBGNoiseLabel()));
            if(showDifference())
            	out.write(String.format(hformat, getDifferenceLabel(), getDifferenceLabel()));
            out.newLine();
            
            for (int i = 0; i < ctfplot.length; i++) {
            	out.write(String.format("%40.4f", xvalues[i]));
            	if(showCTF())
                	out.write(String.format(vformat, ctfplot[i], ctf_avgplot[i]));
                if(showPSD())
                	out.write(String.format(vformat, psdprofileplot[i], psdprofile_avgplot[i]));
                if(showTheorethicalPSD())
                	out.write(String.format(vformat, theorethicalpsdplot[i], theorethicalpsd_avgplot[i]));
                if(showEnvelope())
                	out.write(String.format(vformat, envelopeplot[i], envelope_avgplot[i]));
                if(showBGNoise())
                	out.write(String.format(vformat, bgnoiseplot[i], bgnoise_avgplot[i]));
                if(showDifference())
                	out.write(String.format(vformat, difference[i], difference_avgplot[i]));

            	out.newLine();
            }
            
            
            out.close();

            return true;
        } catch (Exception ex) {
            ex.printStackTrace();
            XmippDialog.showError(this, ex.getMessage());
        }

        return false;
    }
    
    public String getDifferenceLabel()
    {
    	return XmippLabel.CB_PLOT_BGNOISE + " Corrected " + XmippLabel.LABEL_PSD + " Profile";
    }
    
    private String getTheorethicalPSDLabel()
	{
		// TODO Auto-generated method stub
		return "Theorethical " + XmippLabel.LABEL_PSD;
	}
    
    public String getPSDProfileLabel()
    {
    	return XmippLabel.LABEL_PSD + " Profile";
    }
    
    public String getEnvelopeLabel()
    {
    	return XmippLabel.CB_PLOT_ENVELOPE;
    }
    
    public String getCTFLabel()
    {
    	return XmippLabel.LABEL_CTF;
    }
    
    
    public String getBGNoiseLabel()
    {
    	return XmippLabel.CB_PLOT_BGNOISE;
    }
    
    
}
