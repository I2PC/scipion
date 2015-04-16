package xmipp.viewer.chimera;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.ServerSocket;
import java.net.Socket;

public class ChimeraClient
{
	private Socket client;
        private BufferedReader in;
        private PrintWriter out;
        private String volfile;
        private int port;
        
        public ChimeraClient(String volfile)
        {
            this.volfile = volfile;
            
            
        }
        
        public void connect()
        {
            try
            {
                port = findFreePort();
                System.out.println(port);
                String command = String.format("xmipp_chimera_client -i %s --mode projector --client %s", this.volfile, port);
                Runtime.getRuntime().exec(command);
                Thread.sleep(3000);
                client = new Socket("", this.port);
                out = new PrintWriter(client.getOutputStream(), true);
                in = new BufferedReader(new InputStreamReader(client.getInputStream()));
                sendMessage("test");
            }
            catch(Exception e)
            {
                e.printStackTrace();
            }
        }
        
	public static void main(String[] args)
	{
		
                //String volfile = String.format("/home/airen/pyworkflow-code/data/tests/xmipp_tutorial/volumes/BPV_scale_filtered_windowed_64.vol", File.separator);
                String volfile = String.format("/home/airen/pyworkflow-code/data/tests/xmipp_tutorial/volumes/BPV_scale_filtered_windowed_64.vol", File.separator);
                ChimeraClient client = new ChimeraClient(volfile);
		client.connect();	
		
	}
	
        public static int findFreePort() 
             throws IOException {
            ServerSocket server = new ServerSocket(0);
            int port = server.getLocalPort();
            server.close();
            return port;
        }
	
	public void sendMessage(String msg) throws Exception
	{
		out.println(msg);
		System.out.println("client:" + msg);
                System.out.println("server:" + in.readLine());
	}

}
