package edu.cnu.mdi.chimera.model;

/**
 * Identifies the physical coordinate system used by a grid.
 * <p>
 * The implementation uses Geocentric Solar Magnetospheric
 * coordinates, abbreviated GSM. In GSM coordinates, the origin is at the
 * center of Earth, the positive x-axis points from Earth toward the Sun, the
 * y-axis lies perpendicular to Earth's magnetic dipole axis in the plane
 * containing the Sun-Earth line, and the z-axis completes the right-handed
 * system, tilted toward the north magnetic pole.
 * </p>
 */
public enum CoordinateSystem {

    /**
     * Geocentric Solar Magnetospheric coordinates.
     */
    GSM("Geocentric Solar Magnetospheric", "GSM");

    /** Long descriptive name. */
    private final String displayName;

    /** Short abbreviation. */
    private final String abbreviation;

    /**
     * Creates a coordinate system descriptor.
     *
     * @param displayName the human-readable name
     * @param abbreviation the short abbreviation
     */
    CoordinateSystem(String displayName, String abbreviation) {
        this.displayName = displayName;
        this.abbreviation = abbreviation;
    }

    /**
     * Gets the human-readable display name.
     *
     * @return the display name
     */
    public String getDisplayName() {
        return displayName;
    }

    /**
     * Gets the short abbreviation.
     *
     * @return the abbreviation
     */
    public String getAbbreviation() {
        return abbreviation;
    }

    @Override
    public String toString() {
        return abbreviation + " (" + displayName + ")";
    }
}