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
			
			String volfile = String.format("%1$shome%1$sairen%1$sxprojects%1$sshowj%1$shand.vol", File.separator);

			int port = 6000;
//			String command = String.format("chimera %s %s ", Filename.getXmippPath("libraries/bindings/chimera/xmipp_chimera_server.py"), port);
			String command = String.format("chimera %s", Filename.getXmippPath("libraries/bindings/chimera/xmipp_chimera_server.py"));
			System.out.println(command);
			Runtime.getRuntime().exec(command);
			Thread.sleep(3000);
			client = new Socket("", port);
			out = new PrintWriter(client.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(client.getInputStream()));
            
//            ImageGeneric ig = new ImageGeneric(volfile);
//			ig.setDataType(ImageGeneric.Double);
//			ig.read(ImageGeneric.ALL_SLICES);
//			
//			float[] data = ig.getArrayFloat(ImageGeneric.FIRST_IMAGE, ImageGeneric.ALL_SLICES);
            String msg = "open_volume" ;
			
            sendMessage(out, msg);
            while ((msg = in.readLine()) != null) {
                System.out.println("Server: " + msg);
               

               
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
		out.write(msg);
		System.out.println(msg);
	}

}
