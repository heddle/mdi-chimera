package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    
    
    // doSplice: main splice logic for prepatches with crossings. Separated out from splice() for readability.
    private static List<ThetaPatch> doSplice(PrePatch pre,
            SphericalGrid sphGrid,
            List<Crossing> crossings) {
    	
    	boolean debug = (pre.nx == 29 && pre.ny == 15 && pre.nz == 20);
    	
    	if (debug) {
			System.out.println("[ThetaPatch] doSplice: prepatch (" + pre.nx + "," + pre.ny + "," + pre.nz + ") with " + crossings.size() + " crossings:");
			for (Crossing c : crossings) {
				System.out.println("  " + c.summaryString());
			}
		}
    	
    	List<ThetaPatch> result = new ArrayList<>();
    	List<BaseCurve> allCurves = new ArrayList<>();
    	
    	// step1: split curves at theta crossings
    	List<BaseCurve> curves = pre.curves;
    	for (BaseCurve curve : curves) {
    		GeneralCurve gcurve = (GeneralCurve) curve;
    		allCurves.addAll(gcurve.splitAtThetaCrossings());
    	}
    	
    	//create back and forth theta curves for each set of matching crossings
    	List<Crossing> remainingCrossings = new ArrayList<>(crossings);
    	List<BaseCurve> thetaCurves = new ArrayList<>();
    	
    	while (!remainingCrossings.isEmpty()) {
			Crossing c = remainingCrossings.remove(0);
			double theta = c.value();
			Crossing match = null;
			
			for (Crossing other : remainingCrossings) {
				if (Math.abs(other.value() - theta) < THETA_TOL) {
					match = other;
					break;
				}
			}
			
			if (match != null) {
				Point3D.Double p0 = c.curve().getPoint(clamp01(c.t()));
				Point3D.Double p1 = match.curve().getPoint(clamp01(match.t()));
				thetaCurves.add(new ThetaCurve(p0, p1, pre.radius));
				thetaCurves.add(new ThetaCurve(p1, p0, pre.radius));
				remainingCrossings.remove(match);
			}
		}
		
		//assemble all curves into loops and create theta patches from loops
		allCurves.addAll(thetaCurves);
		List<List<BaseCurve>> loops = assembleClosedLoops(allCurves, pre.nx, pre.ny, pre.nz, -1);
		
		for (List<BaseCurve> loop : loops) {
			if (loop.size() < 2) {
				continue;
			}
			
			try {
				int thetaIndex = bestThetaIndexForLoop(loop, sphGrid);
				result.add(new ThetaPatch(loop, pre.nx, pre.ny, pre.nz, thetaIndex));
			} catch (IllegalArgumentException ex) {
				System.err.printf(
						"[ThetaPatch] altDoSplice: %s%n",
						ex.getMessage());
			}
		}
    	return result;
    }
    
    /**
     * Chooses the theta-band index for a closed theta-splice loop.
     *
     * <p>
     * The loop should lie entirely within one theta band, except for any
     * {@link ThetaCurve} segments that run along theta-grid boundaries. This method
     * samples the midpoint of each non-boundary curve and votes for the theta cell
     * containing that midpoint. Boundary theta curves are ignored when possible
     * because their midpoint lies on a grid line and can be assigned ambiguously to
     * either adjacent band.
     * </p>
     *
     * <p>
     * If all curves are boundary-like, the method falls back to sampling all curve
     * midpoints and using {@link Grid1D#locateInterval(double)}. This should be rare,
     * but it keeps the method total.
     * </p>
     *
     * @param loop    closed loop produced by theta splicing
     * @param sphGrid spherical grid
     * @return theta cell index for the loop
     */
    private static int bestThetaIndexForLoop(List<BaseCurve> loop, SphericalGrid sphGrid) {
        Grid1D thetaGrid = sphGrid.getThetaGrid();
        Map<Integer, Integer> counts = new HashMap<>();

        /*
         * First pass: vote only with curves that are not lying on a theta cut.
         * These are the most reliable indicators of the interior theta band.
         */
        for (BaseCurve curve : loop) {
            if (curve == null) {
                continue;
            }

            double th0 = curve.theta(0.0);
            double thm = curve.theta(0.5);
            double th1 = curve.theta(1.0);

            boolean onSomeCut = false;
            for (int i = 0; i < thetaGrid.numPoints(); i++) {
                double cut = thetaGrid.valueAt(i);
                if (Math.abs(th0 - cut) < ON_CUT_TOL
                        && Math.abs(thm - cut) < ON_CUT_TOL
                        && Math.abs(th1 - cut) < ON_CUT_TOL) {
                    onSomeCut = true;
                    break;
                }
            }

            if (onSomeCut) {
                continue;
            }

            int idx = thetaGrid.locateInterval(thm);
            if (idx >= 0) {
                counts.merge(idx, 1, Integer::sum);
            }
        }

        if (!counts.isEmpty()) {
            return counts.entrySet().stream()
                    .max(Map.Entry.comparingByValue())
                    .get()
                    .getKey();
        }

        /*
         * Fallback: all curves looked like theta-boundary curves. Use midpoint
         * samples anyway, nudging exact grid-line values very slightly inward when
         * needed.
         */
        for (BaseCurve curve : loop) {
            if (curve == null) {
                continue;
            }

            double theta = curve.theta(0.5);
            int idx = thetaGrid.locateInterval(theta);

            if (idx < 0) {
                /*
                 * Handle exact theta-grid vertices or tiny roundoff excursions.
                 */
                int vertex = thetaGrid.valueIsAVertex(theta);
                if (vertex >= 0) {
                    if (vertex == 0) {
                        idx = 0;
                    } else if (vertex >= thetaGrid.numPoints() - 1) {
                        idx = thetaGrid.numCells() - 1;
                    } else {
                        /*
                         * Ambiguous interior cut. Pick the lower adjacent band as a
                         * deterministic fallback.
                         */
                        idx = vertex - 1;
                    }
                }
            }

            if (idx >= 0) {
                counts.merge(idx, 1, Integer::sum);
            }
        }

        if (!counts.isEmpty()) {
            return counts.entrySet().stream()
                    .max(Map.Entry.comparingByValue())
                    .get()
                    .getKey();
        }

        throw new IllegalStateException("Could not determine theta index for theta-splice loop.");
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


}