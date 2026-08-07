package edu.cnu.mdi.chimera.model;

import java.util.Objects;

import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphericalGrid;

/**
 * Provides the grid specification used by geometry objects during an algorithm
 * run.
 *
 * <p>The original implementation reached through the Swing application
 * singleton whenever a curve or patch needed a grid. Besides coupling the
 * geometry layer to the desktop UI, that could make
 * {@link edu.cnu.mdi.chimera.alg.ChimeraAlgorithm#run(ChimeraGridSpec)} use a
 * different grid from the method argument. This context is activated from the
 * explicit algorithm input and therefore also permits deterministic headless
 * tests.</p>
 *
 * <p>The current algorithm is synchronous, so a single active immutable grid
 * specification is sufficient. If runs become concurrent, this class should be
 * replaced by explicit constructor/method dependencies rather than expanded
 * into mutable per-thread state.</p>
 */
public final class ChimeraGridContext {

    private static volatile ChimeraGridSpec activeSpec;

    private ChimeraGridContext() {
    }

    /**
     * Activates the grids used by subsequently constructed curves and patches.
     *
     * @param gridSpec immutable specification supplied to the algorithm
     * @throws NullPointerException if {@code gridSpec} is {@code null}
     */
    public static void activate(ChimeraGridSpec gridSpec) {
        activeSpec = Objects.requireNonNull(gridSpec, "gridSpec must not be null");
    }

    /** @return a defensive copy of the active Cartesian grid */
    public static CartesianGrid cartesianGrid() {
        return activeSpec().getCartesianGrid();
    }

    /** @return a defensive copy of the active spherical grid */
    public static SphericalGrid sphericalGrid() {
        return activeSpec().getSphericalGrid();
    }

    /** @return the active sphere radius */
    public static double radius() {
        return activeSpec().getSphericalGrid().getRadius();
    }

    /**
     * Returns the active specification or fails explicitly when geometry is used
     * outside an algorithm run. A silent default would recreate the original bug
     * in which geometry could be evaluated against an unrelated grid.
     */
    private static ChimeraGridSpec activeSpec() {
        ChimeraGridSpec result = activeSpec;
        if (result == null) {
            throw new IllegalStateException(
                    "No Chimera grid is active; call ChimeraAlgorithm.run first.");
        }
        return result;
    }
}
