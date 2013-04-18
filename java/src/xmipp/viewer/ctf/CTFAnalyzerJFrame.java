package xmipp.viewer.ctf;

import ij.ImagePlus;
import ij.gui.Line;
import ij.gui.ProfilePlot;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
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
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import xmipp.jni.CTFDescription;
import xmipp.jni.MDLabel;
import xmipp.jni.MetaData;
import xmipp.utils.XmippLabel;
import xmipp.utils.XmippWindowUtil;

public class CTFAnalyzerJFrame extends JFrame
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
	double[] ctf_avgplot, theorethicalpsd_avgplot, envelope_avgplot, bgnoise_avgplot, psdprofile_avgplot;

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

			constraints = new GridBagConstraints();
			constraints.insets = new Insets(0, 5, 0, 5);
			constraints.anchor = GridBagConstraints.NORTHWEST;
			// constraints.fill = GridBagConstraints.HORIZONTAL;
			setLayout(new GridBagLayout());

			initImagePane();
			add(imagepn, XmippWindowUtil.getConstraints(constraints, 0, 0, 1));

			initGraphicPane();
			add(actionspn, XmippWindowUtil.getConstraints(constraints, 0, 1, 1, 1));
			add(graphicpn, XmippWindowUtil.getConstraints(constraints, 1, 0, 1, 2));
			pack();
			setVisible(true);
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
	}

	private void initImagePane()
	{
		imagepn = new JPanel(new GridBagLayout());
		imagepn.setBorder(BorderFactory.createTitledBorder("PSD Image"));
		imageprofilepn = new CTFAnalyzerImagePane(imp, this);
		imagepn.add(imageprofilepn, XmippWindowUtil.getConstraints(constraints, 0, 0));

	}

	private void initGraphicPane()
	{
		actionspn = new JPanel(new GridBagLayout());
		actionspn.setBorder(BorderFactory.createTitledBorder("Display"));
		ctfrb = new JRadioButton(XmippLabel.LABEL_CTF);
		ctfrb.addActionListener(new FillGraphicsActionListener());

		psdrb = new JRadioButton(XmippLabel.LABEL_PSD + " Profile");
		psdrb.addActionListener(new FillGraphicsActionListener());
		ButtonGroup drawbg = new ButtonGroup();
		drawbg.add(ctfrb);
		drawbg.add(psdrb);
		ctfrb.setSelected(true);
		actionspn.add(ctfrb, XmippWindowUtil.getConstraints(constraints, 0, 0));
		actionspn.add(psdrb, XmippWindowUtil.getConstraints(constraints, 0, 1));

		constraints.insets = new Insets(0, 20, 0, 5);
		theorethicalpsdchb = new JCheckBox("Theorethical PSD");
		theorethicalpsdchb.addActionListener(new FillGraphicsActionListener());
		envelopepsdchb = new JCheckBox("Envelope");
		envelopepsdchb.addActionListener(new FillGraphicsActionListener());
		noisechb = new JCheckBox("Noise");
		noisechb.addActionListener(new FillGraphicsActionListener());

		actionspn.add(theorethicalpsdchb, XmippWindowUtil.getConstraints(constraints, 0, 2));
		actionspn.add(envelopepsdchb, XmippWindowUtil.getConstraints(constraints, 0, 3));
		actionspn.add(noisechb, XmippWindowUtil.getConstraints(constraints, 0, 4));

		graphicpn = new JPanel();
		tabspn = new JTabbedPane();
		JPanel tab1pn = new JPanel();
		JPanel tab2pn = new JPanel();
		tabspn.addTab("Radial Profile", tab1pn);
		tabspn.addTab("Radial Average", tab2pn);
		graphicpn.add(tabspn);
		JFreeChart chart = ChartFactory
				.createXYLineChart("", XmippLabel.LABEL_SAMPLING, showCTF() ? XmippLabel.LABEL_CTF : XmippLabel.LABEL_PSD, collection, PlotOrientation.VERTICAL, true, true, false);

		radialchartpn = createChartPanel(chart);
		JFreeChart avgchart = ChartFactory
				.createXYLineChart("", XmippLabel.LABEL_SAMPLING, showCTF() ? XmippLabel.LABEL_CTF : XmippLabel.LABEL_PSD, collection, PlotOrientation.VERTICAL, true, true, false);

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
		if (psdprofile_avgplot == null)//avg data does not change, calculated only once
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
		}
		XYSeriesCollection collection = getXYSeriesCollection(psdprofile_avgplot, ctf_avgplot, theorethicalpsd_avgplot, envelope_avgplot, bgnoise_avgplot);
		XYPlot plot = ((XYPlot) avgchartpn.getChart().getPlot());
		plot.setDataset(0, collection);
		plot.getRangeAxis().setLabel(showCTF() ? XmippLabel.LABEL_CTF : XmippLabel.LABEL_PSD);
	}

	private void fillRadialGraphics()
	{
		double[] ctfplot, theorethicalpsdplot, envelopeplot, bgnoiseplot, psdprofileplot;
		
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

		XYSeriesCollection collection = getXYSeriesCollection(psdprofileplot, ctfplot, theorethicalpsdplot, envelopeplot, bgnoiseplot);

		XYPlot plot = ((XYPlot) radialchartpn.getChart().getPlot());
		plot.setDataset(0, collection);
		plot.getRangeAxis().setLabel(showCTF() ? XmippLabel.LABEL_CTF : XmippLabel.LABEL_PSD);
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
		else
		{
			collection.addSeries(createSeries(XmippLabel.LABEL_PSD + " " + XmippLabel.CB_PLOT_PROFILE, xvalues, psdprofileplot, max, min));

			// min = max;
			max = -max;
			for (int i = 0; i < xvalues.length; i++)
			{
				max = Math.max(max, psdprofileplot[i]);
				// min = Math.min(min, profile[i]);
			}
			max += Math.abs(max) * 0.1;
			// min -= Math.abs(min) * 0.1;

			if (showBGNoise())
			{
				collection.addSeries(createSeries(XmippLabel.CB_PLOT_BGNOISE, xvalues, bgnoiseplot, max, min));
				double[] difference = new double[bgnoiseplot.length];
				for (int i = 0; i < xvalues.length; i++)
					difference[i] = psdprofileplot[i] - bgnoiseplot[i];
				collection.addSeries(createSeries(XmippLabel.CB_PLOT_DIFFERENCE, xvalues, difference, max, min));
			}
			if (showEnvelope())
				collection.addSeries(createSeries(XmippLabel.CB_PLOT_ENVELOPE, xvalues, envelopeplot, max, min));

			if (showTheorethicalPSD())
				collection.addSeries(createSeries("Theorethical " + XmippLabel.CB_PLOT_PSD, xvalues, theorethicalpsdplot, max, min));
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
}
