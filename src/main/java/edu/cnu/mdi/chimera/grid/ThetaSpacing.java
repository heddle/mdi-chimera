package edu.cnu.mdi.chimera.grid;

/**
 * Describes how theta grid points are generated for a spherical grid.
 */
public enum ThetaSpacing {

    /**
     * Equally spaced values of theta from 0 to pi.
     */
    UNIFORM_THETA("Uniform theta"),

    /**
     * Equally spaced values of cos(theta), which gives more nearly equal-area
     * latitude bands on the sphere.
     */
    UNIFORM_COS_THETA("Uniform cos(theta)");

    /** Human-readable label. */
    private final String label;

    /**
     * Creates a theta spacing option.
     *
     * @param label the display label
     */
    ThetaSpacing(String label) {
        this.label = label;
    }

    @Override
    public String toString() {
        return label;
    }
}