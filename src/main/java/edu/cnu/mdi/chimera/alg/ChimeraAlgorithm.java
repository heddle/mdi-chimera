package edu.cnu.mdi.chimera.alg;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import edu.cnu.mdi.chimera.cell.Cell;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphereIntersectionScanner;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.patch.Patch;
import edu.cnu.mdi.chimera.patch.PrePatch;
import edu.cnu.mdi.chimera.patch.ThetaPatch;

public class ChimeraAlgorithm {

    /**
     * The main entry point for the algorithm. It takes a ChimeraGridSpec as
     * input, which contains both the Cartesian grid and the Spherical grid
     * specifications. It returns a ChimeraAlgorithmResult containing the list
     * of intersecting cells and the prepatches built from them.
     */
    public static ChimeraAlgorithmResult run(ChimeraGridSpec gridSpec) {
        Objects.requireNonNull(gridSpec, "gridSpec must not be null.");
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
        
        // TODO: Step 4 — phi splice
        // TODO: Step 5 — area and perimeter

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
