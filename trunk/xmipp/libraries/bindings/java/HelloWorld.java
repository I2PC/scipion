import xmipp.*;

class HelloWorld {
	public static void main(String[] args) {
		HelloWorld h = new HelloWorld();

			String path = args[0];
			String s = path + "osoImages/oso.xmp";
		    
		    ImageDouble image = new ImageDouble();
		    image.read(s);
		    image.printShape();
		    
		    ImageDouble image2 = image;
		    image2.printShape();
		    int[] size = image2.getSize();
		    System.out.println("x: " + size[0] + " y: " + size[1] + " z: " + size[2]);
		    double [] data = image2.getData();
		    s = path + "oso.sel";
		    MetaData md = new MetaData(s);
		    md.print();

		    if (md.containsLabel(MDLabel.MDL_IMAGE))
		    	System.out.println("contains image");
		    
		    int label = MDLabel.MDL_SHIFTX;
		    boolean hasLabel = md.containsLabel(label);
		    if (hasLabel)
		    	System.out.println("contains rot");
		    long[] ids = md.findObjects();
		    System.out.println(ids.length);
		    for (long id: ids)
		    {
		    	System.out.print(id);
		    	System.out.print("-->");
		    	if (hasLabel)
		    		System.out.print(md.getValueDouble(label, id));
		    	System.out.println();
		    		
		    }
		
	}
}