package edu.cnu.mdi.chimera.mc;

/**
 * One Monte Carlo sample point on the spherical shell.
 * <p>
 * The point is stored by its GSM spherical coordinates and by the five
 * cell indices that identify the Monte Carlo patch bin:
 * </p>
 *
 * <pre>
 * (nx, ny, nz, ntheta, nphi)
 * </pre>
 *
 * <p>
 * The {@code colorValue} is a value in {@code [0, 1]} used to sample the active
 * scientific color map at draw time. Storing the color value rather than a
 * {@code Color} object keeps the point independent of the selected color map.
 * </p>
 *
 * @param theta GSM polar angle measured from +GSM Z, in radians
 * @param phi GSM azimuth measured from +GSM X toward +GSM Y, in radians
 * @param nx Cartesian x-cell index
 * @param ny Cartesian y-cell index
 * @param nz Cartesian z-cell index
 * @param ntheta spherical theta-cell index
 * @param nphi spherical phi-cell index
 * @param colorValue color-map sample value in {@code [0, 1]}
 */
public record MonteCarloPoint(
        double theta,
        double phi,
        int nx,
        int ny,
        int nz,
        int ntheta,
        int nphi,
        double colorValue) {
}