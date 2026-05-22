package edu.cnu.mdi.chimera.mc;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.IntConsumer;

import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.map.ChimeraPatchColor;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;


/**
 * Generates uniformly distributed Monte Carlo points on the spherical
 * shell.
 * <p>
 * Uniformity on the sphere is achieved by drawing {@code phi} uniformly in
 * {@code [-pi, pi)} and {@code cos(theta)} uniformly in {@code [-1, 1]}.
 * </p>
 */
public final class MonteCarloGenerator {

    /**
     * Hidden constructor for utility class.
     */
    private MonteCarloGenerator() {
    }

    /**
     * Generates Monte Carlo points for the supplied grid specification.
     *
     * @param gridSpec active grid specification
     * @param count number of points to generate
     * @param random random number generator
     * @param progressConsumer optional progress callback receiving integer
     *        percentages from 0 to 100; may be {@code null}
     * @param cancelCheck optional cancellation check; may be {@code null}
     * @return generated points
     */
    public static List<MonteCarloPoint> generate(
            ChimeraGridSpec gridSpec,
            int count,
            Random random,
            IntConsumer progressConsumer,
            CancelCheck cancelCheck) {

        if (gridSpec == null) {
            throw new IllegalArgumentException("gridSpec must not be null.");
        }
        if (count < 0) {
            throw new IllegalArgumentException("count must be nonnegative.");
        }
        if (random == null) {
            random = new Random();
        }

        CartesianGrid cartesianGrid = gridSpec.getCartesianGrid();
        SphericalGrid sphericalGrid = gridSpec.getSphericalGrid();

        Grid1D xGrid = cartesianGrid.getXGrid();
        Grid1D yGrid = cartesianGrid.getYGrid();
        Grid1D zGrid = cartesianGrid.getZGrid();
        Grid1D thetaGrid = sphericalGrid.getThetaGrid();
        Grid1D phiGrid = sphericalGrid.getPhiGrid();

        double radius = sphericalGrid.getRadius();

        ArrayList<MonteCarloPoint> points = new ArrayList<>(count);

        int lastProgress = -1;

        for (int i = 0; i < count; i++) {
            if (cancelCheck != null && cancelCheck.isCancelled()) {
                break;
            }

            double phi = -Math.PI + 2.0 * Math.PI * random.nextDouble();

            /*
             * Uniform sphere sampling: u = cos(theta) is uniform in [-1, 1].
             */
            double u = -1.0 + 2.0 * random.nextDouble();
            double theta = Math.acos(u);

            double sinTheta = Math.sin(theta);
            double x = radius * sinTheta * Math.cos(phi);
            double y = radius * sinTheta * Math.sin(phi);
            double z = radius * Math.cos(theta);

            int nx = xGrid.cellIndex(x);
            int ny = yGrid.cellIndex(y);
            int nz = zGrid.cellIndex(z);
            int ntheta = thetaGrid.cellIndex(theta);
            int nphi = phiGrid.cellIndex(phi);

            /*
             * These should normally all be valid if the Cartesian grid encloses
             * the spherical shell and the angular grids span theta=[0,pi],
             * phi=[-pi,pi]. The guard keeps pathological boundary cases from
             * poisoning the visualization.
             */
            if (nx >= 0 && ny >= 0 && nz >= 0 && ntheta >= 0 && nphi >= 0) {
                double colorValue = ChimeraPatchColor.colorValue(nx, ny, nz, ntheta, nphi);
                points.add(new MonteCarloPoint(theta, phi, nx, ny, nz, ntheta, nphi, colorValue));
            }

            if (progressConsumer != null) {
                int progress = (count == 0) ? 100 : (int) Math.round(100.0 * (i + 1) / count);
                if (progress != lastProgress) {
                    progressConsumer.accept(progress);
                    lastProgress = progress;
                }
            }
        }

        if (progressConsumer != null) {
            progressConsumer.accept(100);
        }

        return points;
    }

    /**
     * Small cancellation callback used by the generator.
     */
    @FunctionalInterface
    public interface CancelCheck {

        /**
         * Checks whether generation has been cancelled.
         *
         * @return {@code true} if generation should stop
         */
        boolean isCancelled();
    }
}