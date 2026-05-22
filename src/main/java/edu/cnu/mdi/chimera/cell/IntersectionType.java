package edu.cnu.mdi.chimera.cell;

public enum IntersectionType {

	CORNERIN("Corner In"),
	CORNEROUT("Corner Out"),
	DOUBLECORNERIN("Double Corner In"),
	DOUBLECORNEROUT("Double Corner Out"),
	FACECUT("Face Cut"),
	CORNERPULL("Corner Pull"),
	CORNERPUSH("Corner Push"),
	SKEWCUT("Skew Cut"),
	KISS("Kiss"),
	UNKNOWN("Unknown");
	

    private final String label;

	/**
	 * Constructor for the enum constants.
	 *
	 * @param label the label to be associated with the enum constant
	 */
	IntersectionType(String label) {
		this.label = label;
	}
	
    @Override
    public String toString() {
        return label;
    }

}
