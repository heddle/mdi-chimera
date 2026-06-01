package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
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
    
    
    // doSplice: main splice logic for prepatches with crossings. 
    //Separated out from splice() for readability.
    private static List<ThetaPatch> doSplice(PrePatch pre,
            SphericalGrid sphGrid,
            List<Crossing> crossings) {
 
    	
   	boolean debug = (pre.nx == 16 && pre.ny == 16 && pre.nz == 30);
    	
    	if (debug) {
			System.out.println("[ThetaPatch] doSplice: prepatch (" + pre.nx + "," + pre.ny + "," + pre.nz + ") with " + crossings.size() + " crossings:");
			for (Crossing c : crossings) {
				System.out.println("  " + c.summaryString());
			}
			System.out.println();
		}
  
  	
    	List<ThetaPatch> result = new ArrayList<>();
    	List<BaseCurve> allCurves = new ArrayList<>();
    	
    	// step 1: split curves at theta crossings
    	List<BaseCurve> curves = pre.curves;
    	for (BaseCurve curve : curves) {
    		GeneralCurve gcurve = (GeneralCurve) curve;
    		allCurves.addAll(gcurve.splitAtThetaCrossings());
    	}
    	
    	if (debug) {
    	    System.out.println("[doSplice] After split, allCurves:");
    	    for (BaseCurve curve : allCurves) {
    	        System.out.printf("  %s  theta0=%.4f thetaMid=%.4f theta1=%.4f%n",
    	                curve.shortString(),
    	                Math.toDegrees(curve.theta(0.0)),
    	                Math.toDegrees(curve.theta(0.5)),
    	                Math.toDegrees(curve.theta(1.0)));
    	    }
    	    System.out.println("[doSplice] Crossings:");
    	    for (Crossing cx : crossings) {
    	        System.out.println("  " + cx.summaryString());
    	    }
    	}
    	
    	// step 2: create back and forth theta curves for each set of matching crossings.
    	// For each crossing, we find the closest unmatched crossing at the same theta value
    	// (within THETA_TOL), measured by great-circle distance. This ensures correct pairing
    	// when there are more than two crossings at the same theta (e.g. four crossings:
    	// we want nearest-neighbor pairs rather than arbitrary first-found pairs).
    	List<Crossing> remainingCrossings = new ArrayList<>(crossings);
    	List<BaseCurve> thetaCurves = new ArrayList<>();

    	while (!remainingCrossings.isEmpty()) {
    	    Crossing c = remainingCrossings.remove(0);
    	    double theta = c.value();

    	    // Find the closest unmatched crossing at the same theta value
    	    Crossing match = null;
    	    double bestDist = Double.MAX_VALUE;
    	    for (Crossing other : remainingCrossings) {
    	        if (Math.abs(other.value() - theta) < THETA_TOL) {
    	            double dist = c.distanceTo(other);
    	            if (dist < bestDist) {
    	                bestDist = dist;
    	                match = other;
    	            }
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

    	// step 3: assemble all curves into loops and create theta patches from loops
    	
    	
    	while (!allCurves.isEmpty()) {
    	    BaseCurve c = allCurves.remove(0);
    	    List<BaseCurve> patchCurves = new ArrayList<>();
    	    patchCurves.add(c);
    	    boolean madeLoop = false;

    	    while (!madeLoop) {

    	        //get the next connected curve, preferring theta curves
    	        boolean connected = false;
    	        
    	        boolean isTheta = (c instanceof ThetaCurve);

				// Check theta curves first (if this is not a theta curve)
				if (!isTheta) {
					Iterator<BaseCurve> thetaIt = thetaCurves.iterator();
					while (thetaIt.hasNext()) {
						BaseCurve theta = thetaIt.next();
						if (connects(c, theta)) {
							patchCurves.add(theta);
							thetaIt.remove();
							c = theta;
							connected = true;
							break;
						}
					}
				}

    	        // If not connected, check remaining curves
    	        if (!connected) {
    	            Iterator<BaseCurve> allIt = allCurves.iterator();
    	            while (allIt.hasNext()) {
    	                BaseCurve other = allIt.next();
    	                if (connects(c, other)) {
    	                    patchCurves.add(other);
    	                    allIt.remove();
    	                    c = other;
    	                    connected = true;
    	                    break;
    	                }
    	            }
    	        }

    	        if (!connected) {
    	            System.err.printf(
    	                    "[ThetaPatch] doSplice: could not connect curve %s in prepatch "
    	                            + "(%d,%d,%d) with %d crossings%n",
    	                    c.shortString(), pre.nx, pre.ny, pre.nz, crossings.size());
    	            break;
    	        }

    	        //have we completed a loop
    	        madeLoop = makesLoop(patchCurves);
    	        if (madeLoop) {
    	            try {
    	                int thetaIndex = bestThetaIndexForLoop(patchCurves, sphGrid);
    	                ThetaPatch patch = new ThetaPatch(patchCurves, pre.nx, pre.ny, pre.nz, thetaIndex);
                        ThetaPatch reversed = reverseCurves(patch);
    	                result.add(reversed);
    	            } catch (IllegalArgumentException ex) {
    	                System.err.printf(
    	                        "[ThetaPatch] doSplice: %s%n",
    	                        ex.getMessage());
    	            }
    	        }
    	    }
    	}

		return result;

	}
  
    // Checks whether the given curves form a closed loop by comparing the 
    //start of the first curve and the end of the last curve.
    private static boolean makesLoop(List<BaseCurve> curves) {
		if (curves.isEmpty()) {
			return false;
		}

		Point3D.Double start = curves.get(0).p0;
		Point3D.Double end = curves.get(curves.size() - 1).p1;

		return pointsClose(start, end);
    }
    
    // Checks whether the end of c1 and the start of c2 are close enough 
    // to be considered connected.
    private static boolean connects(BaseCurve c1, BaseCurve c2) {
 		Point3D.Double c1End = c1.p1;
		Point3D.Double c2Start = c2.p0;
		
		return  pointsClose(c1End, c2Start);
    	
    }
    
 
    // Reverses the direction of all curves in the patch and reverses their order to maintain connectivity.
    private static ThetaPatch reverseCurves(ThetaPatch patch) {
        if (patch == null) {
            throw new IllegalArgumentException("Cannot reverse a null ThetaPatch.");
        }

        List<BaseCurve> reversed = new ArrayList<>(patch.curves.size());

        /*
         * If the original loop is:
         *
         *   c0: p0 -> p1
         *   c1: p1 -> p2
         *   c2: p2 -> p3
         *   ...
         *   cn: pn -> p0
         *
         * then the reversed loop must be:
         *
         *   cn.reverse(): p0 -> pn
         *   ...
         *   c2.reverse(): p3 -> p2
         *   c1.reverse(): p2 -> p1
         *   c0.reverse(): p1 -> p0
         *
         * Reversing each curve without reversing the list order does not preserve
         * connectivity.
         */
        for (int i = patch.curves.size() - 1; i >= 0; i--) {
            reversed.add(patch.curves.get(i).reverse());
        }

        return new ThetaPatch(reversed, patch.nx, patch.ny, patch.nz, patch.nTheta);
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