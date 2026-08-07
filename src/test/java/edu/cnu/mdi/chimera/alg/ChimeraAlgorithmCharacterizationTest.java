package edu.cnu.mdi.chimera.alg;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import java.util.HashSet;
import java.util.Set;

import org.junit.jupiter.api.Test;

import edu.cnu.mdi.chimera.model.ChimeraGridPresets;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.patch.BasePatch;
import edu.cnu.mdi.chimera.patch.Patch;
import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.cell.Cell;
import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.util.MathUtil;

/**
 * End-to-end characterization tests for the patch construction pipeline.
 *
 * <p>These tests assert topological invariants instead of exact patch counts.
 * Exact counts are sensitive to grid choice and numerical tolerances, whereas a
 * successful splice must preserve closed directed loops, assign both spherical
 * indices, and preserve spherical area.</p>
 */
class ChimeraAlgorithmCharacterizationTest {

    /**
     * Verifies that the small deterministic grid reaches the final-patch stage
     * and that every returned boundary is a closed, fully indexed loop.
     */
    @Test
    void smallGridProducesClosedFullyIndexedFinalPatches() {
        ChimeraAlgorithmResult result =
                ChimeraAlgorithm.run(ChimeraGridPresets.smallDebugGrid());

        List<Patch> patches = result.getPatches();
        assertFalse(patches.isEmpty(), "phi splice must produce final patches");
        assertTrue(patches.size() >= result.getThetaPatchCount(),
                "phi splicing may preserve or subdivide a theta patch, never discard it");

        for (Patch patch : patches) {
            assertTrue(patch.isFullyIndexed(), patch::fullIndex);
            assertTrue(BasePatch.validateLoop(patch.curves), patch::fullIndex);
        }
    }

    /**
     * Verifies approximate area conservation across both splice stages. The
     * tolerance accommodates sampled curved-edge area estimates while still
     * detecting missing or multiply assembled loops.
     */
    @Test
    void spliceStagesApproximatelyPreserveArea() {
        ChimeraAlgorithmResult result =
                ChimeraAlgorithm.run(ChimeraGridPresets.smallDebugGrid());

        double prePatchArea = result.getPrePatches().stream()
                .mapToDouble(BasePatch::areaEstimate).sum();
        double thetaPatchArea = result.getThetaPatches().stream()
                .mapToDouble(BasePatch::areaEstimate).sum();
        double finalPatchArea = result.getPatches().stream()
                .mapToDouble(BasePatch::areaEstimate).sum();

        assertRelativeDifferenceBelow(prePatchArea, thetaPatchArea, 0.01,
                "theta splice area");
        assertRelativeDifferenceBelow(thetaPatchArea, finalPatchArea, 0.01,
                "phi splice area");
    }

    /**
     * Verifies the defining final-patch invariant: every non-grid-boundary sample
     * lies in the patch's assigned theta and phi bands. It also checks the
     * thesis's requirement that the five-tuple identify one final patch.
     */
    @Test
    void finalPatchesStayInsideTheirAssignedSphericalCell() {
        ChimeraGridSpec spec = ChimeraGridPresets.smallDebugGrid();
        assertFinalPatchCellInvariants(spec);
    }

    /**
     * Applies the same seam, pole, loop, and five-tuple checks to the denser grid
     * derived from the thesis test case.
     */
    @Test
    void paperGridFinalPatchesStayInsideTheirAssignedSphericalCell() {
        assertFinalPatchCellInvariants(ChimeraGridPresets.paperTestGrid());
    }

    /**
     * Checks the thesis's two global validation criteria on the paper grid:
     * sampled five-tuple classification and total spherical area. Samples that
     * land in intentionally unsupported kiss cells are excluded.
     */
    @Test
    void paperGridCoversSampledSphereAndPreservesArea() {
        ChimeraGridSpec spec = ChimeraGridPresets.paperTestGrid();
        ChimeraAlgorithmResult result = ChimeraAlgorithm.run(spec);

        double finalArea = result.getPatches().stream()
                .mapToDouble(BasePatch::areaEstimate).sum();
        assertRelativeDifferenceBelow(1.0, finalArea, 0.005,
                "paper-grid spherical coverage");

        Set<String> kissCells = new HashSet<>();
        for (Cell cell : result.getKissCells()) {
            kissCells.add(cartesianKey(cell.nx, cell.ny, cell.nz));
        }

        int[] indices = new int[5];
        int classified = 0;
        int sampleCount = 2_000;
        double goldenRatio = (1.0 + Math.sqrt(5.0)) / 2.0;
        double radius = spec.getSphericalGrid().getRadius();

        for (int i = 0; i < sampleCount; i++) {
            double zFraction = 1.0 - 2.0 * (i + 0.5) / sampleCount;
            double theta = Math.acos(zFraction);
            double phi = MathUtil.normalizeAngle(
                    2.0 * Math.PI * i / goldenRatio);

            spec.getPatchIndices(theta, phi, radius, indices);
            Patch patch = (Patch) BasePatch.fromSortedList(result.getPatches(),
                    indices[0], indices[1], indices[2], indices[3], indices[4]);

            if (patch == null) {
                assertTrue(kissCells.contains(
                        cartesianKey(indices[0], indices[1], indices[2])),
                        () -> "no final patch for sampled five-tuple "
                                + java.util.Arrays.toString(indices));
                continue;
            }

            double sinTheta = Math.sin(theta);
            assertTrue(patch.containsPoint(
                    radius * sinTheta * Math.cos(phi),
                    radius * sinTheta * Math.sin(phi),
                    radius * Math.cos(theta)),
                    () -> "five-tuple lookup returned a non-containing patch "
                            + patch.fullIndex());
            classified++;
        }

        assertTrue(classified >= 0.995 * sampleCount,
                "too many samples fell outside analytic final patches: "
                        + classified + "/" + sampleCount);
    }

    private static void assertFinalPatchCellInvariants(ChimeraGridSpec spec) {
        ChimeraAlgorithmResult result = ChimeraAlgorithm.run(spec);
        Grid1D thetaGrid = spec.getSphericalGrid().getThetaGrid();
        Grid1D phiGrid = spec.getSphericalGrid().getPhiGrid();
        Set<String> indices = new HashSet<>();

        for (Patch patch : result.getPatches()) {
            assertTrue(indices.add(patch.fullIndex()),
                    () -> "duplicate final-patch five-tuple " + patch.fullIndex());

            for (BaseCurve curve : patch.curves) {
                for (double t : new double[] { 0.25, 0.5, 0.75 }) {
                    double theta = curve.theta(t);
                    if (!onGridCut(thetaGrid, theta, false)) {
                        assertTrue(thetaGrid.cellIndex(theta) == patch.nTheta,
                                () -> patch.fullIndex() + " contains theta=" + theta);
                    }

                    double phi = MathUtil.normalizeAngle(curve.phi(t));
                    if (!onGridCut(phiGrid, phi, true)) {
                        assertTrue(phiGrid.cellIndex(phi) == patch.nPhi,
                                () -> patch.fullIndex() + " contains phi=" + phi
                                        + " on " + curve.getClass().getSimpleName()
                                        + " " + curve.shortString()
                                        + "\nloop=" + loopSummary(patch));
                    }
                }
            }
        }
    }

    private static boolean onGridCut(Grid1D grid, double value, boolean periodic) {
        for (int i = 0; i < grid.numPoints(); i++) {
            double difference = periodic
                    ? Math.abs(MathUtil.normalizeAngle(value - grid.valueAt(i)))
                    : Math.abs(value - grid.valueAt(i));
            if (difference <= 1.0e-8) {
                return true;
            }
        }
        return false;
    }

    private static String loopSummary(Patch patch) {
        StringBuilder summary = new StringBuilder();
        for (BaseCurve curve : patch.curves) {
            summary.append("\n  ")
                    .append(curve.getClass().getSimpleName())
                    .append(' ')
                    .append(curve.shortString());
        }
        return summary.toString();
    }

    private static String cartesianKey(int nx, int ny, int nz) {
        return nx + ":" + ny + ":" + nz;
    }

    private static void assertRelativeDifferenceBelow(
            double expected, double actual, double tolerance, String label) {
        double scale = Math.max(Math.abs(expected), 1.0e-12);
        double relativeDifference = Math.abs(actual - expected) / scale;
        assertTrue(relativeDifference <= tolerance,
                () -> String.format("%s changed by %.6f%% (%.12f -> %.12f)",
                        label, 100.0 * relativeDifference, expected, actual));
    }
}
