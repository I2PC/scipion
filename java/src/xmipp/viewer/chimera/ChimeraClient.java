package xmipp.viewer.chimera;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;

import xmipp.jni.Filename;
import xmipp.jni.ImageGeneric;

public class ChimeraClient
{
	
	public static void main(String[] args)
	{
		Socket client;
		BufferedReader in;
		PrintWriter out;
		try
		{
			String var = "XMIPP_CHIMERA_PORT";
			String value = "6000";
			
			int port = Integer.parseInt(value);
			String volfile = String.format("%1$shome%1$sairen%1$sxprojects%1$sshowj%1$shand.vol", File.separator);
			System.setProperty(var, value);
			Runtime.getRuntime().exec("chimera " + Filename.getXmippPath("libraries/bindings/chimera/xmipp_chimera_server.py"));
			Thread.sleep(3000);
			
			client = new Socket("", port);
			out = new PrintWriter(client.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(client.getInputStream()));
			ImageGeneric ig = new ImageGeneric(volfile);
			ig.read(ImageGeneric.ALL_SLICES);
			//ig.setDataType(ImageGeneric.Double);
			float[] data = ig.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.ALL_SLICES);
			out.write("open_volume");
			//out.wri(data);
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

}
