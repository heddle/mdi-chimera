package edu.cnu.mdi.chimera.map;

import java.awt.Color;

import edu.cnu.mdi.ui.colors.ScientificColorMap;

/**
 * Heuristic color assignment for Monte Carlo patches.
 * <p>
 * Monte Carlo points are colored without consulting the analytic patch objects
 * or an adjacency graph. This class therefore assigns colors directly from the
 * five indices:
 * </p>
 *
 * <pre>
 * (nx, ny, nz, ntheta, nphi)
 * </pre>
 *
 * <p>
 * The color class is computed by a modular affine hash. With a prime modulus
 * and nonzero coefficients, any two bins that differ by exactly one of the five
 * indices are guaranteed to receive different color classes. This is not a
 * full graph coloring of the final patches, but it is a good deterministic
 * approximation for Monte Carlo visualization.
 * </p>
 */
public final class ChimeraPatchColor {

    /**
     * Number of categorical color classes.
     * <p>
     * Use a prime number. Eleven gives enough variety while still keeping colors
     * visually distinguishable when sampled from a perceptual color map.
     * </p>
     */
    private static final int NUM_CLASSES = 11;

    /**
     * Contrast-friendly order for sampling a continuous color map.
     * <p>
     * These are intentionally not monotonic. A monotonic Viridis sequence can
     * make nearby classes look too similar. This order jumps around the color
     * map to improve contrast.
     * </p>
     */
    private static final double[] COLOR_POSITIONS = {
            0.06, 0.86, 0.24, 0.98, 0.42, 0.70,
            0.14, 0.78, 0.34, 0.58, 0.92
    };

    /**
     * Hidden constructor for utility class.
     */
    private ChimeraPatchColor() {
    }

    /**
     * Gets a deterministic color for a five-index patch ID.
     *
     * @param nx the Cartesian x-cell index
     * @param ny the Cartesian y-cell index
     * @param nz the Cartesian z-cell index
     * @param ntheta the spherical theta-cell index
     * @param nphi the spherical phi-cell index
     * @param colorMap the scientific color map to sample; if {@code null},
     *        {@link ScientificColorMap#VIRIDIS} is used
     * @return the assigned color
     */
    public static Color getColor(int nx, int ny, int nz,
            int ntheta, int nphi, ScientificColorMap colorMap) {

        ScientificColorMap cmap = (colorMap == null)
                ? ScientificColorMap.VIRIDIS
                : colorMap;

        int colorClass = colorClass(nx, ny, nz, ntheta, nphi);
        return cmap.colorAt(COLOR_POSITIONS[colorClass]);
    }

    /**
     * Gets a deterministic color with an alpha value.
     *
     * @param nx the Cartesian x-cell index
     * @param ny the Cartesian y-cell index
     * @param nz the Cartesian z-cell index
     * @param ntheta the spherical theta-cell index
     * @param nphi the spherical phi-cell index
     * @param colorMap the color map; if {@code null}, Viridis is used
     * @param alpha alpha in [0, 255]
     * @return the assigned color with alpha
     */
    public static Color getColor(int nx, int ny, int nz,
            int ntheta, int nphi, ScientificColorMap colorMap, int alpha) {

        Color base = getColor(nx, ny, nz, ntheta, nphi, colorMap);
        int a = Math.max(0, Math.min(255, alpha));
        return new Color(base.getRed(), base.getGreen(), base.getBlue(), a);
    }

    /**
     * Computes the discrete color class for the five-index patch ID.
     * <p>
     * The coefficients are all nonzero modulo {@link #NUM_CLASSES}. Therefore,
     * changing any one index by one changes the color class.
     * </p>
     *
     * @param nx the Cartesian x-cell index
     * @param ny the Cartesian y-cell index
     * @param nz the Cartesian z-cell index
     * @param ntheta the spherical theta-cell index
     * @param nphi the spherical phi-cell index
     * @return a color class in [0, {@link #NUM_CLASSES} - 1]
     */
    public static int colorClass(int nx, int ny, int nz, int ntheta, int nphi) {
        int value = nx
                + 2 * ny
                + 3 * nz
                + 5 * ntheta
                + 7 * nphi;

        return Math.floorMod(value, NUM_CLASSES);
    }

    /**
     * Gets the number of color classes used by this heuristic.
     *
     * @return the number of color classes
     */
    public static int getNumClasses() {
        return NUM_CLASSES;
    }
    
    /**
     * Gets the color-map sample value for a five-index patch ID.
     * <p>
     * This is useful when storing large Monte Carlo point clouds. The model can
     * store one {@code double} color value per point instead of storing a full
     * {@link java.awt.Color} object. The active color map can then be applied later
     * at draw time.
     * </p>
     *
     * @param nx the Cartesian x-cell index
     * @param ny the Cartesian y-cell index
     * @param nz the Cartesian z-cell index
     * @param ntheta the spherical theta-cell index
     * @param nphi the spherical phi-cell index
     * @return a color-map sample value in {@code [0, 1]}
     */
    public static double colorValue(int nx, int ny, int nz, int ntheta, int nphi) {
        return COLOR_POSITIONS[colorClass(nx, ny, nz, ntheta, nphi)];
    }
}
