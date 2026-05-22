package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.List;

import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.GeneralCurve;
import edu.cnu.mdi.chimera.curve.ThetaCurve;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * A theta-patch is a {@link PrePatch} sliced by one band {@code [θ_lo, θ_hi]}
 * of the spherical theta grid (paper §6.4). It carries a 4-tuple
 * {@code (nx, ny, nz, nTheta)} with {@code nPhi = -1} until the phi splice.
 */
public class ThetaPatch extends BasePatch {

    private ThetaPatch(CartesianGrid cartGrid, SphericalGrid sphGrid,
                       List<BaseCurve> curves,
                       int nx, int ny, int nz, int nTheta) {
        super(cartGrid, sphGrid, curves, nx, ny, nz, nTheta, -1);
    }

    @Override
    public boolean containsPoint(double x, double y, double z) {
        int[] cart = cartesianGrid.getIndices(x, y, z, new int[3]);
        if (cart[0] != nx || cart[1] != ny || cart[2] != nz) return false;
        double theta = Math.acos(Math.max(-1.0, Math.min(1.0, z / radius)));
        int[] sph = sphericalGrid.getIndices(theta, Math.atan2(y, x), new int[2]);
        return sph[0] == nTheta;
    }

    private static final double THETA_TOL    = 1.0e-9;
    private static final double ENDPOINT_TOL = 1.0e-8;

    // -----------------------------------------------------------------------
    // Theta splice
    // -----------------------------------------------------------------------

    /**
     * Slices a {@link PrePatch} along all theta grid lines it spans, returning
     * one {@link ThetaPatch} per theta band.
     *
     * <h3>Key insight</h3>
     * <p>The boundary curves of each band must form a closed loop in traversal
     * order. ThetaCurve connectors between bands must therefore be inserted
     * at the exact position in the sequence where the crossing occurs —
     * not appended at the end. We achieve this with a single ordered walk:</p>
     * <ul>
     *   <li>We maintain one open curve-list per band, representing the
     *       partially-built boundary traversed so far.</li>
     *   <li>When we cross from band A to band B at point P, we record P as
     *       band A's pending exit point.</li>
     *   <li>When we later cross back into band A at point Q (on the same or
     *       a different cut line), we immediately prepend a
     *       {@code ThetaCurve(P→Q)} — wait, we APPEND it since we are
     *       building in traversal order and P was already emitted.</li>
     *   <li>After the full loop, each band's pending exit is connected back
     *       to its very first point by the final ThetaCurve.</li>
     * </ul>
     */
    public static List<ThetaPatch> splice(PrePatch pre) {
        CartesianGrid cartGrid  = pre.cartesianGrid;
        SphericalGrid sphGrid   = pre.sphericalGrid;
        Grid1D        thetaGrid = sphGrid.getThetaGrid();
        double        R         = sphGrid.getRadius();

        // Find the true theta range including curve interiors.
        double thetaMin = Double.MAX_VALUE, thetaMax = -Double.MAX_VALUE;
        for (BaseCurve c : pre.curves) {
            for (int s = 0; s <= 10; s++) {
                double th = (s == 0) ? c.sv0.theta
                          : (s == 10) ? c.sv1.theta
                          : c.theta(s / 10.0);
                thetaMin = Math.min(thetaMin, th);
                thetaMax = Math.max(thetaMax, th);
            }
        }

        int nThetaLo = Math.max(0, thetaGrid.locateInterval(thetaMin));
        int nThetaHi = Math.min(thetaGrid.numCells() - 1,
                                thetaGrid.locateInterval(thetaMax));

        // Single band — no cutting needed.
        if (nThetaLo == nThetaHi) {
            try {
                return List.of(new ThetaPatch(cartGrid, sphGrid,
                        new ArrayList<>(pre.curves),
                        pre.nx, pre.ny, pre.nz, nThetaLo));
            } catch (IllegalArgumentException ex) {
                System.err.printf("[ThetaSplice] Single-band (%d,%d,%d): %s%n",
                        pre.nx, pre.ny, pre.nz, ex.getMessage());
                return List.of();
            }
        }

        int numBands = nThetaHi - nThetaLo + 1;

        // cutTheta[k] is the gridline separating band k from band k+1
        // (0-indexed relative to nThetaLo).
        double[] cutTheta = new double[numBands - 1];
        for (int k = 0; k < numBands - 1; k++) {
            cutTheta[k] = thetaGrid.valueAt(nThetaLo + k + 1);
        }

        // Per-band state:
        //   bandCurves[i]   — curves accumulated so far (in traversal order)
        //   pendingExit[i]  — the last exit point from this band (awaiting
        //                     a ThetaCurve to be appended when re-entered)
        //   firstEntry[i]   — the very first entry point into this band
        //                     (needed to close the loop at the end)
        @SuppressWarnings("unchecked")
        List<BaseCurve>[] bandCurves = new List[numBands];
        Point3D.Double[]  pendingExit = new Point3D.Double[numBands];
        Point3D.Double[]  firstEntry  = new Point3D.Double[numBands];
        for (int i = 0; i < numBands; i++) bandCurves[i] = new ArrayList<>();

        // Determine the band of the very first curve's start point.
        int curBand = bandOf(pre.curves.get(0).sv0.theta, cutTheta, numBands);
        int initialBand = curBand; // remember — this band closes via the prepatch loop itself

        // firstEntry for the initial band is left null: its loop is already
        // closed by the prepatch boundary (last point connects to p0 naturally).
        // All other bands entered via a gridline crossing will have firstEntry
        // set to that crossing point (an exact gridline point) below.

        // ---------------------------------------------------------------
        // Single ordered walk over all curves.
        // ---------------------------------------------------------------
        for (BaseCurve curve : pre.curves) {

            // Handle band changes at curve boundaries (when sv0 is on a cut line).
            // findCrossings skips sv0-near-cut to avoid double-counting, so we
            // must handle the junction-point band transition here instead.
            //
            // When sv0.theta is exactly on a cut, bandOf(sv0) gives the band BELOW
            // the cut. But the curve body may be in the band ABOVE (if the curve
            // travels upward). Use an interior sample to decide the actual band.
            double sv0Theta = curve.sv0.theta;
            // Junction tolerance: sv0 must be this close to a cut to trigger a junction.
            // This must be larger than floating-point noise from getPoint() but smaller
            // than the minimum distance between adjacent theta gridlines.
            final double JUNCTION_TOL = 1.0e-4;
            boolean sv0OnCut = false;
            for (double ct : cutTheta) {
                if (Math.abs(sv0Theta - ct) < JUNCTION_TOL) { sv0OnCut = true; break; }
            }
            int curveBodyBand;
            if (sv0OnCut) {
                // sv0 is on a cut — use interior theta to determine which band the curve is in.
                double interiorTheta = curve.theta(0.1);
                curveBodyBand = bandOf(interiorTheta, cutTheta, numBands);
            } else {
                curveBodyBand = bandOf(sv0Theta, cutTheta, numBands);
            }
            int curveStartBand = curveBodyBand;
            if (sv0OnCut && curveStartBand != curBand) {
                // We crossed a gridline at the junction between the previous curve
                // and this one. The crossing point is curve.p0 (= previous curve's p1).
                Point3D.Double junctionPt = curve.p0;
                // Leave curBand.
                pendingExit[curBand] = junctionPt;
                // Enter curveStartBand.
                if (firstEntry[curveStartBand] == null) {
                    firstEntry[curveStartBand] = junctionPt;
                } else if (pendingExit[curveStartBand] != null) {
                    ThetaCurve tc = new ThetaCurve(pendingExit[curveStartBand], junctionPt, R);
                    bandCurves[curveStartBand].add(tc);
                    pendingExit[curveStartBand] = null;
                }
                curBand = curveStartBand;
            }

            // After a junction, findCrossings may return a spurious near-t=0 crossing
            // for the same cut the junction already handled (because sv0 is not exactly
            // on the cut gridline, just close to it). Filter it out by advancing tPrev
            // past any crossing at t < 0.1 that the junction already resolved.
            List<double[]> crossings = findCrossings(curve, cutTheta);
            double tPrev = 0.0;
            if (sv0OnCut && !crossings.isEmpty() && crossings.get(0)[0] < 0.15) {
                // This early crossing duplicates the junction — skip it by consuming
                // the fragment [0, tEarly] directly into curBand without a band change.
                double tEarly = crossings.get(0)[0];
                BaseCurve frag = subrange(curve, 0.0, tEarly);
                if (frag != null) bandCurves[curBand].add(frag);
                crossings.remove(0);
                tPrev = tEarly;
            }

            for (double[] cross : crossings) {
                double tCut = cross[0];

                // Emit fragment [tPrev, tCut] into curBand.
                if (tCut - tPrev > ENDPOINT_TOL) {
                    BaseCurve frag = subrange(curve, tPrev, tCut);
                    if (frag != null) bandCurves[curBand].add(frag);
                }

                Point3D.Double crossPt = curve.getPoint(tCut);

                // Leaving curBand: record the pending exit point.
                pendingExit[curBand] = crossPt;

                // Determine the next band by sampling theta just after the crossing.
                // If tCut is near t=1, sample just before instead and infer direction.
                int nextBand;
                if (tCut + 1e-4 < 1.0 - ENDPOINT_TOL) {
                    // Forward sample: gives the band we are entering.
                    nextBand = bandOf(curve.theta(tCut + 1e-4), cutTheta, numBands);
                } else {
                    // Near t=1: backward sample gives the band we came FROM (= curBand).
                    // The next band is simply the other side of the cut from curBand.
                    // Since crossings only move ±1 band, nextBand = 2*curBand_at_cut - curBand.
                    // Simpler: use sv1.theta to determine the final band.
                    nextBand = bandOf(curve.sv1.theta, cutTheta, numBands);
                }

                // Entering nextBand.
                if (pendingExit[nextBand] != null) {
                    // Re-entering a band we previously left — insert the ThetaCurve connector.
                    ThetaCurve tc = new ThetaCurve(pendingExit[nextBand], crossPt, R);
                    bandCurves[nextBand].add(tc);
                    pendingExit[nextBand] = null;
                } else if (firstEntry[nextBand] == null) {
                    // First time entering this band — record the entry point.
                    firstEntry[nextBand] = crossPt;
                }

                curBand = nextBand;
                tPrev   = tCut;
            }

            // Final fragment [tPrev, 1.0].
            if (1.0 - tPrev > ENDPOINT_TOL) {
                BaseCurve frag = subrange(curve, tPrev, 1.0);
                if (frag != null) bandCurves[curBand].add(frag);
            }
        }

        // ---------------------------------------------------------------
        // Close each band's loop with the remaining ThetaCurve connectors.
        //
        // After the full walk, each band that was exited but not re-entered
        // has a non-null pendingExit. The closing ThetaCurve runs from
        // pendingExit[i] back to firstEntry[i].
        //
        // The initial band is excluded: its closure is guaranteed by the
        // prepatch's own closed boundary (the last point of the last fragment
        // connects naturally back to p0, just as in the original prepatch).
        // ---------------------------------------------------------------
        for (int i = 0; i < numBands; i++) {
            if (i == initialBand) continue; // closed by prepatch loop itself
            if (pendingExit[i] != null && firstEntry[i] != null) {
                ThetaCurve closing = new ThetaCurve(pendingExit[i], firstEntry[i], R);
                bandCurves[i].add(closing);
            }
        }

        // ---------------------------------------------------------------
        // Wrap each band into a ThetaPatch.
        // ---------------------------------------------------------------
        List<ThetaPatch> result = new ArrayList<>();
        for (int i = 0; i < numBands; i++) {
            List<BaseCurve> bc = bandCurves[i];
            if (bc.size() < 2) continue;
            int bandIdx = nThetaLo + i;
            try {
                result.add(new ThetaPatch(cartGrid, sphGrid, bc,
                        pre.nx, pre.ny, pre.nz, bandIdx));
            } catch (IllegalArgumentException ex) {
                System.err.printf("[ThetaSplice] Band %d of prepatch (%d,%d,%d): %s%n",
                        bandIdx, pre.nx, pre.ny, pre.nz, ex.getMessage());
            }
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /**
     * Returns the 0-based band index (relative to nThetaLo) for the given
     * theta value.  Band k spans {@code [cutTheta[k-1], cutTheta[k]]}.
     */
    private static int bandOf(double theta, double[] cutTheta, int numBands) {
        int band = 0;
        for (double cut : cutTheta) {
            if (theta > cut + ENDPOINT_TOL) band++;
            else break;
        }
        return Math.max(0, Math.min(numBands - 1, band));
    }

    /**
     * Finds all t values where {@code curve.theta(t)} equals one of the
     * {@code cutTheta} values, sorted by ascending t, skipping near-endpoint
     * hits that are within {@link #ENDPOINT_TOL} of 0 or 1.
     */
    /**
     * Finds all t values in (0, 1) where {@code curve.theta(t)} equals one of
     * the {@code cutTheta} values, sorted by ascending t.
     *
     * <p>Endpoint hits (t ≈ 0 or t ≈ 1) are skipped intentionally. A crossing
     * at t=0 means the curve starts exactly on a gridline — the band change is
     * handled by the caller's {@code bandOf} logic at curve boundaries. A
     * near-endpoint hit at t≈1 must NOT be treated as a mid-curve crossing
     * because {@code curve.getPoint(t≈1)} would return the curve endpoint
     * (whose theta may differ slightly from the gridline) causing a ThetaCurve
     * to be constructed with mismatched theta values and a loop gap.</p>
     */
    /**
     * Finds ALL t values in (0,1) where {@code curve.theta(t)} equals one of
     * the cut thetas, sorted by ascending t. Searches repeatedly in the
     * remaining interval after each found root to catch multiple crossings
     * of the same gridline (e.g. a curve that bows above a cut and returns).
     */
    private static List<double[]> findCrossings(BaseCurve curve, double[] cutTheta) {
        List<double[]> crossings = new ArrayList<>();

        // True theta range including interior samples.
        double curveMin = curve.sv0.theta, curveMax = curve.sv0.theta;
        for (int s = 1; s <= 20; s++) {
            double th = (s == 20) ? curve.sv1.theta : curve.theta(s / 20.0);
            curveMin = Math.min(curveMin, th);
            curveMax = Math.max(curveMax, th);
        }

        for (int k = 0; k < cutTheta.length; k++) {
            double cut = cutTheta[k];
            if (cut < curveMin - ENDPOINT_TOL || cut > curveMax + ENDPOINT_TOL) continue;

            // Skip endpoints on the cut.
            if (Math.abs(curve.sv0.theta - cut) < ENDPOINT_TOL) continue;
            if (Math.abs(curve.sv1.theta - cut) < ENDPOINT_TOL) continue;

            // Search for ALL crossings of this cut by repeatedly scanning
            // the remaining interval after each found root.
            double tLo = ENDPOINT_TOL;
            double tHi = 1.0 - ENDPOINT_TOL;
            while (tLo < tHi - ENDPOINT_TOL) {
                double t = MathUtil.computeT(
                        curve.getThetaFunction(), cut, tLo, tHi, THETA_TOL, 200);
                if (Double.isNaN(t)) break;
                crossings.add(new double[]{t, k});
                // Advance past this root to find the next one.
                tLo = t + ENDPOINT_TOL * 100;
            }
        }

        crossings.sort((a, b) -> Double.compare(a[0], b[0]));
        return crossings;
    }

    /** Returns a subrange of {@code curve} over {@code [t0, t1]}. */
    private static BaseCurve subrange(BaseCurve curve, double t0, double t1) {
        if (t0 >= t1) return null;
        if (curve instanceof GeneralCurve gc) return gc.subrange(t0, t1);
        return new ThetaCurve(curve.getPoint(t0), curve.getPoint(t1), curve.radius);
    }
}