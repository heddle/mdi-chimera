package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.analysis.UnivariateFunction;

import edu.cnu.mdi.chimera.app.ChimeraApp;
import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.Crossing;
import edu.cnu.mdi.chimera.curve.GeneralCurve;
import edu.cnu.mdi.chimera.curve.ThetaCurve;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * A theta-patch is a {@link PrePatch} sliced by one band
 * {@code [theta_lo, theta_hi]} of the spherical theta grid.
 *
 * <p>The patch carries a 4-tuple {@code (nx, ny, nz, nTheta)} with
 * {@code nPhi = -1} until the phi splice.</p>
 */
public class ThetaPatch extends BasePatch {

    private ThetaPatch(List<BaseCurve> curves,
                       int nx, int ny, int nz, int nTheta) {
        super(curves, nx, ny, nz, nTheta, -1);
    }

    @Override
    public boolean containsPoint(double x, double y, double z) {
        CartesianGrid cartesianGrid = ChimeraApp.getInstance().getCartesianGrid();
        SphericalGrid sphericalGrid = ChimeraApp.getInstance().getSphericalGrid();

        int[] cart = cartesianGrid.getIndices(x, y, z, new int[3]);
        if (cart[0] != nx || cart[1] != ny || cart[2] != nz) {
            return false;
        }

        double theta = Math.acos(Math.max(-1.0, Math.min(1.0, z / radius)));
        int[] sph = sphericalGrid.getIndices(theta, Math.atan2(y, x), new int[2]);
        return sph[0] == nTheta;
    }

    // ---------------------------------------------------------------------
    // Constants
    // ---------------------------------------------------------------------

    private static final double THETA_TOL = 1.0e-9;
    private static final double ENDPOINT_TOL = 1.0e-8;
    private static final double ENDPOINT_ARTIFACT_TOL = 1.0e-4;

    /*
     * Loop-building tolerance must be larger than BasePatch.LOOP_TOL because
     * we sometimes project crossing endpoints onto exact theta grid lines.
     * The assembled curves themselves are constructed with matching projected
     * endpoints, so BasePatch validation should still pass with its tighter
     * tolerance.
     */
    private static final double LOOP_BUILD_TOL = 1.0e-6;

    /*
     * Used when deciding whether a curve fragment lies along a theta cut.
     */
    private static final double ON_CUT_TOL = 1.0e-8;

    /*
     * Used to snap a fragment endpoint to an exact theta grid line when the
     * endpoint is numerically very close to one.
     */
    private static final double PROJECT_THETA_TOL = 1.0e-6;

    // ---------------------------------------------------------------------
    // Public splice entry point
    // ---------------------------------------------------------------------

    /**
     * Slices a {@link PrePatch} along the spherical theta grid.
     *
     * @param pre the prepatch to slice
     * @return theta patches produced from the prepatch
     */
    public static List<ThetaPatch> splice(PrePatch pre) {
        SphericalGrid sphGrid = ChimeraApp.getInstance().getSphericalGrid();
        ArrayList<ThetaPatch> result = new ArrayList<>();

        List<Crossing> allCrossings = pre.getAllThetaCrossings();
        allCrossings = Crossing.removeDuplicates(allCrossings);
        allCrossings = removeEndpointCrossings(allCrossings, pre.curves);

        int numThetaCrossings = allCrossings.size();

        if (numThetaCrossings == 0) {
            int thetaIndex = bestThetaIndex(pre, sphGrid);
            result.add(new ThetaPatch(pre.curves, pre.nx, pre.ny, pre.nz, thetaIndex));
            return result;
        }

//        if (pre.polar()) {
//            System.err.println("[ThetaPatch] polar splice with "
//                    + numThetaCrossings + " theta crossings.");
//            return polarSplice(pre);
//        }

        if (numThetaCrossings % 2 != 0) {
            System.err.println("[ThetaPatch] Warning: odd crossing count after endpoint removal ("
                    + numThetaCrossings + ") for prepatch ("
                    + pre.nx + "," + pre.ny + "," + pre.nz + ").");

            System.err.println("Curves:");
            for (BaseCurve c : pre.curves) {
                System.err.println("  " + c.shortString());
            }

            System.err.println("\nCrossings:");
            for (Crossing c : allCrossings) {
                System.err.println("  " + c.summaryString());
            }

            throw new IllegalStateException("Odd number of theta crossings");
        }

        result.addAll(doSplice(pre, sphGrid, allCrossings));
        return result;
    }

    // ---------------------------------------------------------------------
    // Endpoint-crossing removal
    // ---------------------------------------------------------------------

    /**
     * Removes endpoint crossings that are only touches, not genuine theta-band
     * changes.
     *
     * <p>An endpoint crossing is genuine only if the first non-on-cut curve
     * before the junction and the first non-on-cut curve after the junction are
     * on opposite sides of the theta cut. Curves that lie along the theta cut
     * are skipped.</p>
     */
    private static List<Crossing> removeEndpointCrossings(
            List<Crossing> crossings, List<BaseCurve> curves) {

        final double probe = 0.02;
        List<Crossing> filtered = new ArrayList<>();

        for (Crossing c : crossings) {
            double t = c.t();
            boolean isEndpoint = (t < ENDPOINT_ARTIFACT_TOL
                    || t > 1.0 - ENDPOINT_ARTIFACT_TOL);

            if (!isEndpoint) {
                filtered.add(c);
                continue;
            }

            BaseCurve thisCurve = c.curve();
            int curveIdx = curveIndexOf(curves, thisCurve);

            if (curveIdx < 0) {
                filtered.add(c);
                continue;
            }

            double cut = c.value();
            boolean isAtEnd = (t > 1.0 - ENDPOINT_ARTIFACT_TOL);

            SideSample before;
            SideSample after;

            if (isAtEnd) {
                before = sampleBeforeJunction(curves, curveIdx, cut, probe);
                after = sampleAfterJunction(curves,
                        (curveIdx + 1) % curves.size(), cut, probe);
            } else {
                before = sampleBeforeJunction(curves,
                        (curveIdx - 1 + curves.size()) % curves.size(),
                        cut, probe);
                after = sampleAfterJunction(curves, curveIdx, cut, probe);
            }

            if (!before.valid || !after.valid) {
                filtered.add(c);
                continue;
            }

            if (before.above == after.above) {
                continue;
            }

            filtered.add(c);
        }

        return filtered;
    }

    private record SideSample(boolean valid, boolean above) {
    }

    private static SideSample sampleBeforeJunction(
            List<BaseCurve> curves, int startIdx, double cut, double probe) {

        int n = curves.size();

        for (int step = 0; step < n; step++) {
            int idx = (startIdx - step + n) % n;
            BaseCurve curve = curves.get(idx);

            if (curveLiesOnThetaCut(curve, cut)) {
                continue;
            }

            double theta = curve.theta(Math.max(0.0, 1.0 - probe));

            if (Math.abs(theta - cut) < THETA_TOL) {
                theta = curve.theta(0.5);
            }

            if (Math.abs(theta - cut) < THETA_TOL) {
                continue;
            }

            return new SideSample(true, theta > cut);
        }

        return new SideSample(false, false);
    }

    private static SideSample sampleAfterJunction(
            List<BaseCurve> curves, int startIdx, double cut, double probe) {

        int n = curves.size();

        for (int step = 0; step < n; step++) {
            int idx = (startIdx + step) % n;
            BaseCurve curve = curves.get(idx);

            if (curveLiesOnThetaCut(curve, cut)) {
                continue;
            }

            double theta = curve.theta(probe);

            if (Math.abs(theta - cut) < THETA_TOL) {
                theta = curve.theta(0.5);
            }

            if (Math.abs(theta - cut) < THETA_TOL) {
                continue;
            }

            return new SideSample(true, theta > cut);
        }

        return new SideSample(false, false);
    }

    private static boolean curveLiesOnThetaCut(BaseCurve curve, double cut) {
        return Math.abs(curve.theta(0.0) - cut) < THETA_TOL
            && Math.abs(curve.theta(0.5) - cut) < THETA_TOL
            && Math.abs(curve.theta(1.0) - cut) < THETA_TOL;
    }

    // ---------------------------------------------------------------------
    // Best theta index for no-splice prepatches
    // ---------------------------------------------------------------------

    private static int bestThetaIndex(PrePatch pre, SphericalGrid sphGrid) {
        Grid1D thetaGrid = sphGrid.getThetaGrid();
        Map<Integer, Integer> counts = new HashMap<>();

        for (BaseCurve curve : pre.curves) {
            double midTheta = curve.theta(0.5);
            int idx = thetaGrid.locateInterval(midTheta);
            if (idx >= 0) {
                counts.merge(idx, 1, Integer::sum);
            }
        }

        if (counts.isEmpty()) {
            return thetaGrid.locateInterval(pre.curves.get(0).sv0.theta);
        }

        if (counts.size() > 1) {
            System.err.printf(
                    "[ThetaPatch] bestThetaIndex: midpoints span %d theta cells "
                    + "for prepatch (%d,%d,%d); using majority%n",
                    counts.size(), pre.nx, pre.ny, pre.nz);
        }

        return counts.entrySet().stream()
                .max(Map.Entry.comparingByValue())
                .get()
                .getKey();
    }

    // ---------------------------------------------------------------------
    // Main theta splice implementation
    // ---------------------------------------------------------------------

    /**
     * Performs the theta splice using unordered edge assembly.
     */
    private static List<ThetaPatch> doSplice(
            PrePatch pre,
            SphericalGrid sphGrid,
            List<Crossing> crossings) {

        Grid1D thetaGrid = sphGrid.getThetaGrid();
        double radius = pre.radius;

        ThetaBandRange range = thetaBandRange(pre, thetaGrid);

        if (range.nThetaLo == range.nThetaHi) {
            int idx = bestThetaIndex(pre, sphGrid);
            return List.of(new ThetaPatch(pre.curves, pre.nx, pre.ny, pre.nz, idx));
        }

        int numBands = range.nThetaHi - range.nThetaLo + 1;
        int numThetaVertices = thetaGrid.numPoints();

        @SuppressWarnings("unchecked")
        List<BaseCurve>[] bandEdges = new List[numBands];

        for (int i = 0; i < numBands; i++) {
            bandEdges[i] = new ArrayList<>();
        }

        /*
         * Pending theta-boundary connectors are tracked by both band and cut.
         * This prevents accidentally connecting points on different theta lines.
         */
        Point3D.Double[][] pendingExit =
                new Point3D.Double[numBands][numThetaVertices];

        Point3D.Double[][] firstEntry =
                new Point3D.Double[numBands][numThetaVertices];

        @SuppressWarnings("unchecked")
        List<Crossing>[] crossingsByCurve = new List[pre.curves.size()];

        for (int i = 0; i < crossingsByCurve.length; i++) {
            crossingsByCurve[i] = new ArrayList<>();
        }

        for (Crossing crossing : crossings) {
            int curveIndex = curveIndexOf(pre.curves, crossing.curve());
            if (curveIndex >= 0) {
                crossingsByCurve[curveIndex].add(crossing);
            }
        }

        for (List<Crossing> list : crossingsByCurve) {
            list.sort((a, b) -> Double.compare(a.t(), b.t()));
        }

        for (int curveIndex = 0; curveIndex < pre.curves.size(); curveIndex++) {
            BaseCurve curve = pre.curves.get(curveIndex);
            List<Crossing> curveCrossings = crossingsByCurve[curveIndex];

            double tPrev = 0.0;
            Crossing prevCrossing = null;

            for (Crossing crossing : curveCrossings) {
                double tCut = clamp01(crossing.t());

                if (tCut - tPrev > ENDPOINT_TOL) {
                    addFragmentEdgesToBands(
                            bandEdges,
                            curve,
                            tPrev,
                            tCut,
                            prevCrossing,
                            crossing,
                            thetaGrid,
                            range.nThetaLo,
                            range.nThetaHi,
                            radius);
                }

                processThetaCrossing(
                        pre,
                        crossing,
                        curveIndex,
                        thetaGrid,
                        range.nThetaLo,
                        range.nThetaHi,
                        bandEdges,
                        pendingExit,
                        firstEntry,
                        radius);

                tPrev = tCut;
                prevCrossing = crossing;
            }

            if (1.0 - tPrev > ENDPOINT_TOL) {
                addFragmentEdgesToBands(
                        bandEdges,
                        curve,
                        tPrev,
                        1.0,
                        prevCrossing,
                        null,
                        thetaGrid,
                        range.nThetaLo,
                        range.nThetaHi,
                        radius);
            }
        }

        /*
         * Close chains that wrap around the artificial start of the original
         * prepatch boundary traversal.
         */
        for (int localBand = 0; localBand < numBands; localBand++) {
            for (int cutIndex = 0; cutIndex < numThetaVertices; cutIndex++) {
                Point3D.Double exit = pendingExit[localBand][cutIndex];
                Point3D.Double entry = firstEntry[localBand][cutIndex];

                if (exit != null && entry != null) {
                    double cutTheta = thetaGrid.valueAt(cutIndex);
                    bandEdges[localBand].add(
                            thetaConnector(exit, entry, cutTheta, radius));
                    pendingExit[localBand][cutIndex] = null;
                }
            }
        }

        List<ThetaPatch> result = new ArrayList<>();

        for (int localBand = 0; localBand < numBands; localBand++) {
            int globalBand = range.nThetaLo + localBand;

            List<List<BaseCurve>> loops = assembleClosedLoops(
                    bandEdges[localBand],
                    pre.nx, pre.ny, pre.nz, globalBand);

            for (List<BaseCurve> loop : loops) {
                if (loop.size() < 2) {
                    continue;
                }

                try {
                    result.add(new ThetaPatch(
                            loop, pre.nx, pre.ny, pre.nz, globalBand));
                } catch (IllegalArgumentException ex) {
                    System.err.printf(
                            "[ThetaSplice] Band %d of prepatch (%d,%d,%d): %s%n",
                            globalBand, pre.nx, pre.ny, pre.nz, ex.getMessage());
                }
            }
        }

        return result;
    }

    private record ThetaBandRange(int nThetaLo, int nThetaHi) {
    }

    private static ThetaBandRange thetaBandRange(PrePatch pre, Grid1D thetaGrid) {
        double thetaMin = Double.MAX_VALUE;
        double thetaMax = -Double.MAX_VALUE;

        for (BaseCurve curve : pre.curves) {
            for (int s = 0; s <= 20; s++) {
                double theta = curve.theta(s / 20.0);
                thetaMin = Math.min(thetaMin, theta);
                thetaMax = Math.max(thetaMax, theta);
            }
        }

        int lo = Math.max(0, thetaGrid.locateInterval(thetaMin));
        int hi = Math.min(thetaGrid.numCells() - 1,
                thetaGrid.locateInterval(thetaMax));

        return new ThetaBandRange(lo, hi);
    }

    // ---------------------------------------------------------------------
    // Fragment creation and band assignment
    // ---------------------------------------------------------------------

    private static void addFragmentEdgesToBands(
            List<BaseCurve>[] bandEdges,
            BaseCurve curve,
            double t0,
            double t1,
            Crossing startCrossing,
            Crossing endCrossing,
            Grid1D thetaGrid,
            int nThetaLo,
            int nThetaHi,
            double radius) {

        if (t1 - t0 <= ENDPOINT_TOL) {
            return;
        }

        Point3D.Double p0 = endpointPoint(
                curve, t0, startCrossing, thetaGrid, radius);

        Point3D.Double p1 = endpointPoint(
                curve, t1, endCrossing, thetaGrid, radius);

        int cutIndex = thetaCutIndexForInterval(curve, t0, t1, p0, p1, thetaGrid);

        BaseCurve fragment;

        if (cutIndex >= 0) {
            double cutTheta = thetaGrid.valueAt(cutIndex);
            fragment = thetaConnector(p0, p1, cutTheta, radius);

            addEdgeToGlobalBand(
                    bandEdges, fragment, cutIndex - 1, nThetaLo, nThetaHi);

            addEdgeToGlobalBand(
                    bandEdges, fragment, cutIndex, nThetaLo, nThetaHi);

            return;
        }

        fragment = new CurveFragment(curve, t0, t1, p0, p1);

        double thetaMid = curve.theta(0.5 * (t0 + t1));
        int globalBand = thetaGrid.locateInterval(thetaMid);

        addEdgeToGlobalBand(
                bandEdges, fragment, globalBand, nThetaLo, nThetaHi);
    }

    private static void addEdgeToGlobalBand(
            List<BaseCurve>[] bandEdges,
            BaseCurve edge,
            int globalBand,
            int nThetaLo,
            int nThetaHi) {

        if (globalBand < nThetaLo || globalBand > nThetaHi) {
            return;
        }

        bandEdges[globalBand - nThetaLo].add(edge);
    }

    private static Point3D.Double endpointPoint(
            BaseCurve curve,
            double t,
            Crossing crossing,
            Grid1D thetaGrid,
            double radius) {

        Point3D.Double raw = curve.getPoint(clamp01(t));

        if (crossing != null) {
            return projectToTheta(raw, crossing.value(), radius);
        }

        double theta = thetaOf(raw);
        int nearest = nearestThetaVertex(thetaGrid, theta, PROJECT_THETA_TOL);

        if (nearest >= 0) {
            return projectToTheta(raw, thetaGrid.valueAt(nearest), radius);
        }

        return raw;
    }

    private static int thetaCutIndexForInterval(
            BaseCurve curve,
            double t0,
            double t1,
            Point3D.Double p0,
            Point3D.Double p1,
            Grid1D thetaGrid) {

        double tm = 0.5 * (t0 + t1);
        double th0 = thetaOf(p0);
        double thm = curve.theta(tm);
        double th1 = thetaOf(p1);

        for (int i = 0; i < thetaGrid.numPoints(); i++) {
            double cut = thetaGrid.valueAt(i);

            if (Math.abs(th0 - cut) < ON_CUT_TOL
                    && Math.abs(thm - cut) < ON_CUT_TOL
                    && Math.abs(th1 - cut) < ON_CUT_TOL) {
                return i;
            }
        }

        return -1;
    }

    private static int nearestThetaVertex(
            Grid1D thetaGrid, double theta, double tolerance) {

        int bestIndex = -1;
        double bestDiff = Double.MAX_VALUE;

        for (int i = 0; i < thetaGrid.numPoints(); i++) {
            double diff = Math.abs(thetaGrid.valueAt(i) - theta);
            if (diff < bestDiff) {
                bestDiff = diff;
                bestIndex = i;
            }
        }

        return (bestDiff <= tolerance) ? bestIndex : -1;
    }

    // ---------------------------------------------------------------------
    // Crossing processing and theta connectors
    // ---------------------------------------------------------------------

    private static void processThetaCrossing(
            PrePatch pre,
            Crossing crossing,
            int curveIndex,
            Grid1D thetaGrid,
            int nThetaLo,
            int nThetaHi,
            List<BaseCurve>[] bandEdges,
            Point3D.Double[][] pendingExit,
            Point3D.Double[][] firstEntry,
            double radius) {

        int cutIndex = crossing.index();
        double cutTheta = crossing.value();

        SideSample before = sampleBeforeCrossing(
                pre.curves, curveIndex, crossing.t(), cutTheta);

        SideSample after = sampleAfterCrossing(
                pre.curves, curveIndex, crossing.t(), cutTheta);

        if (!before.valid || !after.valid) {
            return;
        }

        int beforeBand = sideToGlobalBand(before, cutIndex, thetaGrid);
        int afterBand = sideToGlobalBand(after, cutIndex, thetaGrid);

        if (beforeBand == afterBand) {
            return;
        }

        Point3D.Double point = projectToTheta(
                crossing.curve().getPoint(clamp01(crossing.t())),
                cutTheta,
                radius);

        recordExit(
                beforeBand,
                cutIndex,
                point,
                nThetaLo,
                nThetaHi,
                pendingExit);

        recordEntry(
                afterBand,
                cutIndex,
                point,
                nThetaLo,
                nThetaHi,
                bandEdges,
                pendingExit,
                firstEntry,
                thetaGrid,
                radius);
    }

    private static void recordExit(
            int globalBand,
            int cutIndex,
            Point3D.Double point,
            int nThetaLo,
            int nThetaHi,
            Point3D.Double[][] pendingExit) {

        if (globalBand < nThetaLo || globalBand > nThetaHi) {
            return;
        }

        pendingExit[globalBand - nThetaLo][cutIndex] = point;
    }

    private static void recordEntry(
            int globalBand,
            int cutIndex,
            Point3D.Double point,
            int nThetaLo,
            int nThetaHi,
            List<BaseCurve>[] bandEdges,
            Point3D.Double[][] pendingExit,
            Point3D.Double[][] firstEntry,
            Grid1D thetaGrid,
            double radius) {

        if (globalBand < nThetaLo || globalBand > nThetaHi) {
            return;
        }

        int localBand = globalBand - nThetaLo;

        if (pendingExit[localBand][cutIndex] != null) {
            Point3D.Double exit = pendingExit[localBand][cutIndex];
            double cutTheta = thetaGrid.valueAt(cutIndex);

            bandEdges[localBand].add(
                    thetaConnector(exit, point, cutTheta, radius));

            pendingExit[localBand][cutIndex] = null;
        } else if (firstEntry[localBand][cutIndex] == null) {
            firstEntry[localBand][cutIndex] = point;
        }
    }

    private static SideSample sampleBeforeCrossing(
            List<BaseCurve> curves,
            int curveIndex,
            double t,
            double cut) {

        final double probe = 0.02;
        BaseCurve curve = curves.get(curveIndex);

        if (t > ENDPOINT_ARTIFACT_TOL) {
            double dt = Math.min(probe, 0.5 * t);
            double theta = curve.theta(t - dt);

            if (Math.abs(theta - cut) >= THETA_TOL) {
                return new SideSample(true, theta > cut);
            }
        }

        return sampleBeforeJunction(curves, curveIndex, cut, probe);
    }

    private static SideSample sampleAfterCrossing(
            List<BaseCurve> curves,
            int curveIndex,
            double t,
            double cut) {

        final double probe = 0.02;
        BaseCurve curve = curves.get(curveIndex);

        if (t < 1.0 - ENDPOINT_ARTIFACT_TOL) {
            double dt = Math.min(probe, 0.5 * (1.0 - t));
            double theta = curve.theta(t + dt);

            if (Math.abs(theta - cut) >= THETA_TOL) {
                return new SideSample(true, theta > cut);
            }
        }

        return sampleAfterJunction(
                curves, (curveIndex + 1) % curves.size(), cut, probe);
    }

    private static int sideToGlobalBand(
            SideSample side, int cutIndex, Grid1D thetaGrid) {

        int band = side.above ? cutIndex : cutIndex - 1;
        return Math.max(0, Math.min(thetaGrid.numCells() - 1, band));
    }

    private static ThetaCurve thetaConnector(
            Point3D.Double p0,
            Point3D.Double p1,
            double theta,
            double radius) {

        Point3D.Double q0 = projectToTheta(p0, theta, radius);
        Point3D.Double q1 = projectToTheta(p1, theta, radius);
        return new ThetaCurve(q0, q1, radius);
    }

    private static Point3D.Double projectToTheta(
            Point3D.Double p, double theta, double radius) {

        double phi = Math.atan2(p.y, p.x);
        double sinTheta = Math.sin(theta);

        return new Point3D.Double(
                radius * sinTheta * Math.cos(phi),
                radius * sinTheta * Math.sin(phi),
                radius * Math.cos(theta));
    }

    private static double thetaOf(Point3D.Double p) {
        double r = Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        return Math.acos(Math.max(-1.0, Math.min(1.0, p.z / r)));
    }

    // ---------------------------------------------------------------------
    // Closed-loop assembly
    // ---------------------------------------------------------------------

    private static List<List<BaseCurve>> assembleClosedLoops(
            List<BaseCurve> edges,
            int nx,
            int ny,
            int nz,
            int nTheta) {

        List<List<BaseCurve>> loops = new ArrayList<>();
        List<BaseCurve> unused = new ArrayList<>(edges);

        while (!unused.isEmpty()) {
            List<BaseCurve> loop = new ArrayList<>();

            BaseCurve first = unused.remove(0);
            loop.add(first);

            Point3D.Double start = first.p0;
            Point3D.Double end = first.p1;

            boolean closed = false;

            while (!unused.isEmpty()) {
                if (pointsClose(end, start)) {
                    closed = true;
                    break;
                }

                int matchIndex = -1;
                boolean reverse = false;

                for (int i = 0; i < unused.size(); i++) {
                    BaseCurve candidate = unused.get(i);

                    if (pointsClose(end, candidate.p0)) {
                        matchIndex = i;
                        reverse = false;
                        break;
                    }

                    if (pointsClose(end, candidate.p1)) {
                        matchIndex = i;
                        reverse = true;
                        break;
                    }
                }

                if (matchIndex < 0) {
                    break;
                }

                BaseCurve next = unused.remove(matchIndex);
                if (reverse) {
                    next = next.reverse();
                }

                loop.add(next);
                end = next.p1;
            }

            if (!closed && pointsClose(end, start)) {
                closed = true;
            }

            if (closed) {
                loops.add(loop);
            } else {
                /*
                 * A single open edge can occur when a prepatch boundary rides along a theta
                 * cut and the edge was conservatively assigned to both adjacent bands. If
                 * that adjacent band has no actual area inside this prepatch, the edge is
                 * an orphan boundary fragment, not a failed theta patch.
                 */
                if (loop.size() > 1 || !unused.isEmpty()) {
                    System.err.printf(
                            "[ThetaSplice] Could not close loop for prepatch "
                            + "(%d,%d,%d) theta band %d; partial curves=%d "
                            + "unused=%d gap=%.3e%n",
                            nx, ny, nz, nTheta,
                            loop.size(), unused.size(),
                            Point3D.Double.distance(end, start));
                }
            }       }

        return loops;
    }

    private static boolean pointsClose(Point3D.Double a, Point3D.Double b) {
        return Point3D.Double.distance(a, b) < LOOP_BUILD_TOL;
    }

    private static int curveIndexOf(List<BaseCurve> curves, BaseCurve target) {
        for (int i = 0; i < curves.size(); i++) {
            if (curves.get(i) == target) {
                return i;
            }
        }
        return -1;
    }

    private static double clamp01(double t) {
        return Math.max(0.0, Math.min(1.0, t));
    }

    // ---------------------------------------------------------------------
    // Polar placeholder
    // ---------------------------------------------------------------------

    private static List<ThetaPatch> polarSplice(PrePatch pre) {
        /*
         * Polar theta splicing is a separate special case. Returning an empty
         * list preserves the current behavior better than attempting a general
         * splice through a pole.
         */
        return List.of();
    }

    // ---------------------------------------------------------------------
    // Internal curve fragment wrapper
    // ---------------------------------------------------------------------

    /**
     * A directed subrange of another curve, with explicitly supplied endpoints.
     *
     * <p>The explicit endpoints are important during theta splicing because
     * crossing points are projected onto exact theta grid lines. The interior
     * of the fragment follows the source curve, but its endpoints match the
     * projected connector points exactly, allowing robust loop assembly.</p>
     */
    private static final class CurveFragment extends BaseCurve {

        private final BaseCurve source;
        private final double t0;
        private final double t1;

        CurveFragment(BaseCurve source,
                      double t0,
                      double t1,
                      Point3D.Double p0,
                      Point3D.Double p1) {
            super(p0, p1, source.radius);
            this.source = source;
            this.t0 = t0;
            this.t1 = t1;
        }

        @Override
        public UnivariateFunction getThetaFunction() {
            return t -> thetaOf(getPoint(t));
        }

        @Override
        public UnivariateFunction getPhiFunction() {
            return t -> {
                Point3D.Double p = getPoint(t);
                return Math.atan2(p.y, p.x);
            };
        }

        @Override
        public Point3D.Double getPoint(double t) {
            if (t <= 0.0) {
                return p0;
            }
            if (t >= 1.0) {
                return p1;
            }

            double sourceT = t0 + t * (t1 - t0);
            return source.getPoint(sourceT);
        }

        @Override
        public BaseCurve reverse() {
            return new CurveFragment(source, t1, t0, p1, p0);
        }

        @Override
        public boolean isConstantTheta() {
            double th0 = theta(0.0);
            double thm = theta(0.5);
            double th1 = theta(1.0);
            double min = Math.min(th0, Math.min(thm, th1));
            double max = Math.max(th0, Math.max(thm, th1));
            return (max - min) < THETA_TOL;
        }

        @Override
        public String toString() {
            return String.format("CurveFragment[source=%d t0=%.6f t1=%.6f]",
                    source.curveId, t0, t1);
        }
    }
}