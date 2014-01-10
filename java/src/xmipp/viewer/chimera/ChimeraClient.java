package xmipp.viewer.chimera;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;

import xmipp.ij.commons.XmippImageConverter;
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
			
			String volfile = String.format("%1$shome%1$sairen%1$sxprojects%1$sBPV_Project%1$sBPV_scale_filtered_windowed.vol", File.separator);
			System.out.println(volfile);
			String serverpy = Filename.getXmippPath("libraries/bindings/chimera/xmipp_chimera_server.py");
			int port = 6000;
//			String command = String.format("chimera %s %s ", serverpy, port);
//			String command = String.format("xmipp_python %s", "/home/airen/socketserver.py");
			String command = String.format("chimera %s", serverpy);
			System.out.println(command);
			Runtime.getRuntime().exec(command);
			Thread.sleep(3000);
			client = new Socket("", port);
			out = new PrintWriter(client.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(client.getInputStream()));
            
//            ImageGeneric ig = new ImageGeneric(volfile);
//			ig.setDataType(ImageGeneric.Double);
//			ig.read(volfile, false);
//			
//			float[] data = ig.getArrayFloat(ImageGeneric.ALL_IMAGES, ImageGeneric.ALL_SLICES);
//			XmippImageConverter.convertToImagePlus(ig).show();
            String msg = "open_volume" ;
            sendMessage(out, msg);
			String cube = "[";
			for(int i = 0; i < 2; i ++)
			{
				cube += "[";
				for(int j = 0; j < 2; j ++)
				{
					cube += "[";
					for(int k = 0; k < 2; k ++)
						cube += " " + k;
					cube += "]";
				}
				cube += "]";
			}
			cube += "]";
			System.out.println(cube);			
			//sendMessage(out, cube);
            while ((msg = in.readLine()) != null) {
                System.out.println("Server: " + msg);
                sendMessage(out, cube);
                break;
               

               
            }
			
			//out.write(data);
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
	
	
	public static void sendMessage(PrintWriter out, String msg)
	{
		out.print(msg);
		System.out.println(msg);
	}

}
