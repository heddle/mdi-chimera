package edu.cnu.mdi.chimera.model;

import java.util.Arrays;

import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.grid.ThetaSpacing;


/**
 * Factory methods for common grid specifications.
 */
public final class ChimeraGridPresets {

    /**
     * Hidden constructor for utility class.
     */
    private ChimeraGridPresets() {
    }

    /**
     * Creates the paper/test grid from the original prototype.
     *
     * @return the paper test grid specification
     */
    public static ChimeraGridSpec paperTestGrid() {
        double[] xgrid = {
                -5.5000, -5.1667, -4.8333, -4.5000, -4.1667, -3.8333, -3.5000,
                -3.1667, -2.8333, -2.5000, -2.1667, -1.8333, -1.5000, -1.1667,
                -0.8333, -0.5000, -0.1667, 0.1667, 0.5000, 0.8333, 1.1667,
                1.5000, 1.8333, 2.1667, 2.5000, 2.8333, 3.1667, 3.5000,
                3.8333, 4.1667, 4.5000, 4.8333, 5.1667, 5.5000 };

        double[] ygrid = {
                -5.6667, -5.3333, -5.0000, -4.6667, -4.3333, -4.0000, -3.6667,
                -3.3333, -3.0000, -2.6667, -2.3333, -2.0000, -1.6667, -1.3333,
                -1.0000, -0.6667, -0.3333, 0.00001, 0.3333, 0.6667, 1.0000,
                1.3333, 1.6667, 2.0000, 2.3333, 2.6667, 3.0000, 3.3333,
                3.6667, 4.0000, 4.3333, 4.6667, 5.0000, 5.3333, 5.6803 };

        double[] zgrid = Arrays.copyOf(xgrid, xgrid.length);

        double[] thetagrid = {
                0.0, 0.0262, 0.0785, 0.1287, 0.1767, 0.2226, 0.2666,
                0.3087, 0.3491, 0.3894, 0.4297, 0.4701, 0.5104, 0.5507,
                0.5910, 0.6377, 0.6917, 0.7541, 0.8264, 0.9099, 1.0066,
                1.1185, 1.2479, 1.3976, 1.5708, 1.7440, 1.8937, 2.0231,
                2.1350, 2.2317, 2.3152, 2.3875, 2.4499, 2.5039, 2.5506,
                2.5909, 2.6312, 2.6715, 2.7119, 2.7522, 2.7925, 2.8329,
                2.8750, 2.9190, 2.9649, 3.0129, 3.0631, 3.1154, Math.PI };

        int numPhiPoints = 33;
        double[] phigrid = new double[numPhiPoints];
        for (int i = 0; i < numPhiPoints; i++) {
            phigrid[i] = -Math.PI + 2.0 * Math.PI * i / (numPhiPoints - 1);
        }

        CartesianGrid cartesianGrid = new CartesianGrid(xgrid, ygrid, zgrid, 0.0, 0.0, 0.0);
        SphericalGrid sphericalGrid = new SphericalGrid(thetagrid, phigrid, 4.60993);

        return new ChimeraGridSpec("Paper Test Grid", CoordinateSystem.GSM,
                LengthUnit.EARTH_RADII, cartesianGrid, sphericalGrid);
    }

    /**
     * Creates a small, symmetric debug grid.
     *
     * @return a small debug grid
     */
    public static ChimeraGridSpec smallDebugGrid() {
        CartesianGrid cartesianGrid = CartesianGrid.uniform(
                -6.0, 6.0, 12,
                -6.0, 6.0, 12,
                -6.0, 6.0, 12);

        SphericalGrid sphericalGrid = SphericalGrid.generated(
                4.0, 12, 16, ThetaSpacing.UNIFORM_THETA);

        return new ChimeraGridSpec("Small Debug Grid", CoordinateSystem.GSM,
                LengthUnit.EARTH_RADII, cartesianGrid, sphericalGrid);
    }

    /**
     * Creates a coarse GSM grid suitable for early visualization testing.
     *
     * @return a coarse GSM grid
     */
    public static ChimeraGridSpec coarseGsmGrid() {
        CartesianGrid cartesianGrid = CartesianGrid.uniform(
                -10.0, 10.0, 30,
                -10.0, 10.0, 30,
                -10.0, 10.0, 30);

        SphericalGrid sphericalGrid = SphericalGrid.generated(
                6.0, 24, 32, ThetaSpacing.UNIFORM_COS_THETA);

        return new ChimeraGridSpec("Coarse GSM Grid", CoordinateSystem.GSM,
                LengthUnit.EARTH_RADII, cartesianGrid, sphericalGrid);
    }
}