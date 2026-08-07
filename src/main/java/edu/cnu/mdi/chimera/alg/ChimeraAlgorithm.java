package edu.cnu.mdi.chimera.alg;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import edu.cnu.mdi.chimera.cell.Cell;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphereIntersectionScanner;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.model.ChimeraGridContext;
import edu.cnu.mdi.chimera.patch.Patch;
import edu.cnu.mdi.chimera.patch.PrePatch;
import edu.cnu.mdi.chimera.patch.ThetaPatch;

public class ChimeraAlgorithm {

    /**
     * Computes the analytic intersection of a Cartesian grid and a spherical
     * surface grid.
     *
     * <p>The pipeline classifies sphere-intersecting Cartesian cells, constructs
     * one directed closed prepatch boundary for every supported non-kiss cell,
     * splices those boundaries at theta grid lines, and finally splices each
     * theta patch at phi grid lines. Returned final patches carry complete
     * {@code (nx, ny, nz, nTheta, nPhi)} indices and closed directed boundaries.</p>
     *
     * <p>Kiss intersections are retained in the result for diagnostics but do
     * not yet generate patches, matching the limitation documented in the
     * thesis.</p>
     *
     * @param gridSpec Cartesian and spherical grids to intersect
     * @return cells and every intermediate and final patch stage
     * @throws NullPointerException if {@code gridSpec} is {@code null}
     * @throws IllegalStateException if any supported cell cannot form a closed
     *         patch or a splice cannot preserve directed-loop topology
     */
    public static ChimeraAlgorithmResult run(ChimeraGridSpec gridSpec) {
        Objects.requireNonNull(gridSpec, "gridSpec must not be null.");
        ChimeraGridContext.activate(gridSpec);
        ChimeraAlgorithmResult result = ChimeraAlgorithmResult.empty();

        CartesianGrid cartGrid = gridSpec.getCartesianGrid();
        SphericalGrid sphGrid  = gridSpec.getSphericalGrid();

        // Step 1: find intersecting cells.
        List<Cell> intersectingCells = findIntersectingCells(cartGrid, sphGrid);
        System.out.println("Found " + intersectingCells.size() + " intersecting cells.");
        result.setIntersectingCells(intersectingCells);

        // Step 2: build prepatches from all non-Kiss cells.
        List<PrePatch> prePatches = PrePatch.buildAll(intersectingCells);
        result.setPrePatches(prePatches);

        // Step 3: theta splice — cut each prepatch by the spherical theta grid.
        List<ThetaPatch> thetaPatches = new ArrayList<>();
        for (PrePatch pre : prePatches) {
            thetaPatches.addAll(ThetaPatch.splice(pre));
        }
        result.setThetaPatches(thetaPatches);
        
        // Step 4: phi splice — cut each theta patch by the spherical phi grid.
        System.out.println("Theta splice produced " + thetaPatches.size() + " theta patches.");
        List<Patch> patches = new ArrayList<>();
        for (ThetaPatch theta : thetaPatches) {
			patches.addAll(Patch.splice(theta));
		}
        result.setPatches(patches);
        System.out.println("Phi splice produced " + patches.size() + " final patches.");
        
        return result;
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    private static List<Cell> findIntersectingCells(CartesianGrid cartGrid,
                                                     SphericalGrid sphGrid) {
        double radius = sphGrid.getRadius();
        SphereIntersectionScanner scanner =
                new SphereIntersectionScanner(cartGrid, radius);
        return scanner.scan();
    }
}
