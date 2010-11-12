package xmipp.io;



public class SpeHeader {

	public static final int headerSize = 4100;

	private static final int DIODEN = 0;
	private static final int AVGEXP = 2;
	private static final int EXPOSURE = 4;
	private static final int MODE = 8;
	private static final int ASYSEQ = 16;
	private static final int DATE = 20;
	private static final int NOSCAN = 34;
	private static final int FACCOUNT = 42;
	private static final int DATATYPE = 108;
	private static final int CALIBNAN = 110;
	private static final int LFLOAT = 620;
	private static final int BKGDFILE = 638;
	private static final int STRIPE = 656;
	private static final int SCRAMBLE = 658;
	private static final int LEXPOS = 660;
	private static final int LNOSCAN = 664;
	private static final int LAVGEXP = 668;
	private static final int VERSION = 688;
	private static final int TYPE = 704;
	private static final int NumberofFrames = 1446;
	private static final int WINVIEW_ID = 2996;

	public static final int RAW_OFFSET = 4100;
	public static final int FLOAT = 0;
	public static final int LONG = 1;
	public static final int INT = 2;
	public static final int UNINT = 3;

	private byte[] header;

	public SpeHeader() {
		this(0, 0, 0, 0);
	}

	public SpeHeader(int datatype, int width, int height, int stackSize) {
		header = new byte[headerSize];
//		setShort(DIODEN, width);
		setShort(AVGEXP, 1);
		setShort(EXPOSURE, -1);
		setShort(MODE, 0);
		setShort(ASYSEQ, 1); // ?
		setString(DATE, "12/31/99\0\0");
//		setShort(NOSCAN, -1);
//		setShort(FACCOUNT, width);
//		setShort(DATATYPE, datatype);
		setByte(CALIBNAN + 6, 3); // ?
		setShort(LFLOAT, 0); // ?
		setString(BKGDFILE, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0");
//		setShort(STRIPE, height);
		setShort(SCRAMBLE, 1); // ?
		setInt(LEXPOS, 1000); // exposure time : 1000 msec
//		setInt(LNOSCAN, height * stackSize);
		setInt(LAVGEXP, 1);
		setString(VERSION, "01.000 02/01/90\0");
		setShort(TYPE, 0); // 0 : ST1000
		setInt(WINVIEW_ID, 0x01234567);

		setProperties(datatype, width, height, stackSize);
	}

	public SpeHeader(byte[] header) {
		this.header = header;
	}

	public byte[] getHeader() {
		return header;
	}

	public void setHeader(byte[] header) {
		this.header = header;
	}

	private int getByte(int index) {
		int value = header[index];
		if (value < 0) value += 256;
		return value;
	}

	private int getShort(int index) {
		int b1 = getByte(index);
		int b2 = getByte(index + 1);
		return ((b2 << 8) | b1);
	}

	private int getInt(int index) {
		int b1 = getByte(index);
		int b2 = getByte(index + 1);
		int b3 = getByte(index + 2);
		int b4 = getByte(index + 3);
		return ((b4 << 24) | (b3 << 16) | (b2 << 8) | b1);
	}

	private String getString(int index, int length) {
		StringBuffer buffer = new StringBuffer(length);
		for (int i=0; i<length; i++) {
			buffer.setCharAt(i, (char)header[index + i]);
		}
		return buffer.toString();
	}

	private boolean setByte(int index, int value) {
		if (value >= 0 && value < 128) {
			header[index] = (byte)value;
			return true;
		} else if (value >= 128 && value <256) {
			header[index] = (byte)(value - 256);
			return true;
		}
		return false;
	}

	private boolean setShort(int index, int value) {
		if (value >= -32768 && value < 32768) {
			int b1 = value & 0xff;
			int b2 = (value >> 8) & 0xff;
			setByte(index, b1);
			setByte(index + 1, b2);
			return true;
		} else {
			return false;
		}
	}

	private void setInt(int index, int value) {
		int b1 = value & 0xff;
		int b2 = (value >> 8) & 0xff;
		int b3 = (value >> 16) & 0xff;
		int b4 = (value >> 24) & 0xff;
		setByte(index, b1);
		setByte(index + 1, b2);
		setByte(index + 2, b3);
		setByte(index + 3, b4);
	}

	private void setString(int index, String string) {
		for (int i=0; i<string.length(); i++) {
			setByte(index + i, (int)string.charAt(i));
		}
	}

	public int getDatatype() {
		return getShort(DATATYPE);
	}

	public int getWidth() {
		return getShort(FACCOUNT);
	}

	public int getHeight() {
		return getShort(STRIPE);
	}

	public int getStackSize() {
		int stripe = getShort(STRIPE);
		if (stripe == 0) return -1;
		int noscan = getShort(NOSCAN);
		int NumFrames = getInt(NumberofFrames);
		if (noscan == 65535) {
			int lnoscan = getInt(LNOSCAN);
			if (lnoscan == -1) {
		 		return NumFrames;} else {
			return lnoscan / stripe;}
		} else {return noscan / stripe;}
	}

	public void setProperties(int datatype, int width, int height, int stackSize) {
		setDatatype(datatype);
		setWidth(width);
		setHeight(height);
		setStackSize(stackSize);
	}

	public void setDatatype(int datatype) {
		setShort(DATATYPE, datatype);
	}

	public void setWidth(int width) {
		setShort(DIODEN, width);
		setShort(FACCOUNT, width);
	}

	public void setHeight(int height) {
		int stackSize = getStackSize();
		if (stackSize < 0) stackSize = 0;
		setShort(STRIPE, height);
		setStackSize(stackSize);
	}

	public void setStackSize(int stackSize) {
		setShort(NOSCAN, -1);
		setInt(LNOSCAN, getHeight() * stackSize);
	}

	public int getPixelSize() {
		return getPixelSize(getDatatype());
	}

	public int getPixelSize(int datatype) {
		switch (datatype) {
			case FLOAT:
				return 4;
			case LONG:
				return 4;
			case INT:
				return 2;
			case UNINT:
				return 2;
		}
		return -1;
	}

}

