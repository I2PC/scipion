
public enum ExitValue {
	OK(0),ERROR(1),YES(2),NO(3),
	CANCEL(4),RUNTIME_ERROR(5),
	PROGRAM_NOT_FOUND(6),EXTERNAL_PROGRAM_BROKEN(7);
	private int value;
	
	public int getValue() {
		return value;
	}

	public void setValue(int value) {
		this.value = value;
	}

	// TODO: check v is valid
	ExitValue(int v){
		value=v;
	}
	
}
