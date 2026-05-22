package edu.cnu.mdi.chimera.model;

/**
 * Length units supported by this app.
 * <p>
 * The default and preferred unit for the GSM application is Earth
 * radii. Thus a coordinate value of {@code x = 10.0} means ten Earth radii
 * in the positive GSM x direction.
 * </p>
 */
public enum LengthUnit {

    /**
     * Earth radii.
     */
    EARTH_RADII("Earth radii", "R\u2091");

    /** Human-readable unit name. */
    private final String displayName;

    /** Short symbol. */
    private final String symbol;

    /**
     * Creates a length unit descriptor.
     *
     * @param displayName the human-readable unit name
     * @param symbol the short symbol
     */
    LengthUnit(String displayName, String symbol) {
        this.displayName = displayName;
        this.symbol = symbol;
    }

    /**
     * Gets the human-readable unit name.
     *
     * @return the display name
     */
    public String getDisplayName() {
        return displayName;
    }

    /**
     * Gets the unit symbol.
     *
     * @return the unit symbol
     */
    public String getSymbol() {
        return symbol;
    }

    @Override
    public String toString() {
        return displayName + " (" + symbol + ")";
    }
}