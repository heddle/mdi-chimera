package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import edu.cnu.mdi.chimera.app.ChimeraApp;
import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.Crossing;
import edu.cnu.mdi.chimera.curve.PhiCurve;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;
import edu.cnu.mdi.chimera.util.SphericalVector;

/**
 * Represents a patch defined by the intersection of curves on a spherical grid.
 * A patch is defined by its Cartesian grid indices (nx, ny, nz) and its spherical
 * grid indices (nTheta, nPhi). The patch contains points that lie within the
 * corresponding grid cell in Cartesian coordinates and the corresponding cell in
 * spherical coordinates. This is the "final stage" of patch construction, where the
 * patch is fully defined and can be used for rendering or analysis.
 */
public class Patch extends BasePatch {


	private static final double PHI_TOL = 1.0e-9;
	private static final double ENDPOINT_ARTIFACT_TOL = 1.0e-4;
	private static final double ON_CUT_TOL = 1.0e-8;

	/**
	 * Constructs a Patch with the specified curves and grid indices.
	 *
	 * @param curves the list of curves that define the patch
	 * @param nx     the Cartesian grid index in the x direction
	 * @param ny     the Cartesian grid index in the y direction
	 * @param nz     the Cartesian grid index in the z direction
	 * @param nTheta the spherical grid index in the theta direction
	 * @param nPhi   the spherical grid index in the phi direction
	 */
	public Patch(List<BaseCurve> curves, int nx, int ny, int nz, int nTheta, int nPhi) {
		super(curves, nx, ny, nz, nTheta, nPhi);
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
		return sph[0] == nTheta && sph[1] == nPhi;
	}
	
	/**
	 * Slices a {@link ThetaPatch} along the spherical phi grid.
	 *
	 * @param theta the ThetaPatch to slice
	 * @return full patches produced from the ThetaPatch, sliced along the spherical
	 *         phi grid
	 */
	public static List<Patch> splice(ThetaPatch theta) {
	    SphericalGrid sphGrid = ChimeraApp.getInstance().getSphericalGrid();
	    ArrayList<Patch> result = new ArrayList<>();

	    List<Crossing> allCrossings = theta.getAllPhiCrossings();
	    allCrossings = Crossing.removeDuplicates(allCrossings);

	    int numPhiCrossings = allCrossings.size();

	    // If there are no phi crossings, we can create a single patch without splicing.
//	    if (numPhiCrossings == 0) {
//	        int phiIndex = bestPhiIndex(theta, sphGrid);
//	        result.add(new Patch(theta.curves, theta.nx, theta.ny, theta.nz,
//	                theta.nTheta, phiIndex));
//	        return result;
//	    }
	    
	    if (theta.polar()) {
	    	System.err.println("Polar patch: " + theta.nx + "," + theta.ny + "," + theta.nz + "," + theta.nTheta);
	    	for (Crossing c : allCrossings) {
	    		System.err.println("  " + c.summaryString());
	    	}
	    }

	    return result;
	}

	/**
	 * Slices a {@link ThetaPatch} along the spherical phi grid.
	 *
	 * @param theta the ThetaPatch to slice
	 * @return full patches produced from the ThetaPatch, sliced along the spherical
	 *         phi grid
	 */
	public static List<Patch> Xsplice(ThetaPatch theta) {
	    SphericalGrid sphGrid = ChimeraApp.getInstance().getSphericalGrid();
	    ArrayList<Patch> result = new ArrayList<>();

	    List<Crossing> allCrossings = theta.getAllPhiCrossings();
	    allCrossings = Crossing.removeDuplicates(allCrossings);
	    allCrossings = removeEndpointPhiCrossings(allCrossings, theta.curves);

	    int numPhiCrossings = allCrossings.size();

	    if (numPhiCrossings == 0) {
	        int phiIndex = bestPhiIndex(theta, sphGrid);
	        result.add(new Patch(theta.curves, theta.nx, theta.ny, theta.nz,
	                theta.nTheta, phiIndex));
	        return result;
	    }

	    if (!theta.polar() && (numPhiCrossings % 2 != 0)) {
	        System.err.println("[Patch] Warning: odd phi crossing count after endpoint removal ("
	                + numPhiCrossings + ") for non-polar theta patch ("
	                + theta.nx + "," + theta.ny + "," + theta.nz + ","
	                + theta.nTheta + ").");

	        System.err.println("Curves:");
	        for (BaseCurve c : theta.curves) {
	            System.err.println("  " + c.shortString());
	        }

	        System.err.println("\nCrossings:");
	        for (Crossing c : allCrossings) {
	            System.err.println("  " + c.summaryString());
	        }

	        throw new IllegalStateException("Odd number of phi crossings");
	    }
	    result.addAll(doSplice(theta, sphGrid, allCrossings));

	    return result;
	}

	private static List<Patch> doSplice(
	        ThetaPatch theta,
	        SphericalGrid sphGrid,
	        List<Crossing> crossings) {

	    List<BaseCurve> splitCurves = splitThetaPatchCurves(theta);

	    List<BaseCurve> phiCurves;

	    if (theta.polar()) {
	        phiCurves = buildPolarPhiConnectors(theta, crossings);
	    } else {
	        List<PhiConnectorPair> connectorPairs =
	                buildPhiConnectorPairs(theta, crossings);
	        phiCurves = selectPhiConnectors(connectorPairs);
	    }

	    /*
	     * Important:
	     * Some original theta-patch boundary curves already lie on phi grid cuts.
	     * In the polar case these are often the meridian curves that run into or
	     * out of the pole. They must behave as splice connectors, not ordinary
	     * boundary fragments, or a final patch can wrap across several phi wedges.
	     */
	    moveExistingPhiCutCurvesToConnectors(splitCurves, phiCurves, sphGrid);

	    return assemblePatches(theta, sphGrid, splitCurves, phiCurves, crossings);
	}
	
	/**
	 * Moves any existing boundary curve that lies on a phi grid line out of the
	 * ordinary boundary list and into the phi-connector list.
	 *
	 * <p>
	 * This is essential for polar theta patches. A theta patch may already contain
	 * a curve from the outer boundary to the pole along a meridian. That curve is a
	 * phi boundary of final patches, not an ordinary exterior edge to be walked
	 * through while assembling one final patch.
	 * </p>
	 *
	 * <p>
	 * Both directions are added so the loop assembler can use whichever direction
	 * closes the current final-patch loop.
	 * </p>
	 */
	private static void moveExistingPhiCutCurvesToConnectors(
	        List<BaseCurve> splitCurves,
	        List<BaseCurve> phiCurves,
	        SphericalGrid sphGrid) {

	    Grid1D phiGrid = sphGrid.getPhiGrid();

	    Iterator<BaseCurve> it = splitCurves.iterator();

	    while (it.hasNext()) {
	        BaseCurve curve = it.next();

	        int cutIndex = phiCutIndexIncludingPole(curve, phiGrid);
	        if (cutIndex < 0) {
	            continue;
	        }

	        double cutPhi = MathUtil.normalizeAngle(phiGrid.valueAt(cutIndex));

	        /*
	         * Rebuild as an explicit PhiCurve on the known meridian. Do not rely on
	         * the PhiCurve ordinary constructor when a pole is involved, because phi
	         * is singular at the pole.
	         */
	        Point3D.Double q0 = projectToPhiOrPole(curve.p0, cutPhi, curve.radius);
	        Point3D.Double q1 = projectToPhiOrPole(curve.p1, cutPhi, curve.radius);

	        if (!pointsClose(q0, q1)) {
	            PhiCurve forward = PhiCurve.onMeridian(q0, q1, curve.radius, cutPhi);
	            PhiCurve backward = PhiCurve.onMeridian(q1, q0, curve.radius, cutPhi);

	            phiCurves.add(forward);
	            phiCurves.add(backward);
	        }

	        it.remove();
	    }
	}
	
	/**
	 * Returns the phi-grid cut index if the curve lies on a phi grid line.
	 *
	 * <p>
	 * Pole samples are ignored because phi is singular there. This lets us detect
	 * meridian curves that run from an ordinary boundary point to a pole.
	 * </p>
	 */
	private static int phiCutIndexIncludingPole(BaseCurve curve, Grid1D phiGrid) {
	    for (int i = 0; i < phiGrid.numPoints(); i++) {
	        double cut = MathUtil.normalizeAngle(phiGrid.valueAt(i));

	        if (curveLiesOnPhiCutIncludingPole(curve, cut)) {
	            return i;
	        }
	    }

	    return -1;
	}

	/**
	 * Tests whether a curve lies on the given phi cut, allowing one or more samples
	 * to be at a pole.
	 */
	private static boolean curveLiesOnPhiCutIncludingPole(BaseCurve curve, double cut) {
	    double[] ts = { 0.0, 0.5, 1.0 };

	    boolean sawNonPoleSample = false;

	    for (double t : ts) {
	        Point3D.Double p = curve.getPoint(t);

	        if (atPole(p, curve.radius)) {
	            continue;
	        }

	        sawNonPoleSample = true;

	        double phi = MathUtil.normalizeAngle(Math.atan2(p.y, p.x));
	        double diff = Math.abs(MathUtil.normalizeAngle(phi - cut));

	        if (diff > ON_CUT_TOL) {
	            return false;
	        }
	    }

	    return sawNonPoleSample;
	}

	/**
	 * Projects an ordinary point to the given meridian, but leaves pole points
	 * exactly at the pole.
	 */
	private static Point3D.Double projectToPhiOrPole(
	        Point3D.Double p,
	        double phi,
	        double radius) {

	    if (atPole(p, radius)) {
	        return new Point3D.Double(
	                0.0,
	                0.0,
	                p.z >= 0.0 ? radius : -radius);
	    }

	    return projectToPhi(p, phi, radius);
	}
	
	/**
	 * Builds meridian connectors for a polar theta patch.
	 *
	 * <p>
	 * In a polar theta patch, a meridian generally intersects the ordinary boundary
	 * once and then ends at the pole. The pole is not reported as a phi crossing
	 * because phi is singular there. Therefore the ordinary non-polar strategy of
	 * pairing two crossings on the same meridian is wrong for polar theta patches.
	 * </p>
	 *
	 * <p>
	 * For each boundary crossing, add both directed meridian connectors:
	 * crossing -> pole and pole -> crossing. The loop assembler will use the
	 * direction needed for each final phi-cell wedge.
	 * </p>
	 */
	private static List<BaseCurve> buildPolarPhiConnectors(
	        ThetaPatch theta,
	        List<Crossing> crossings) {

	    List<BaseCurve> phiCurves = new ArrayList<>();

	    Point3D.Double pole = polarPoint(theta);

	    for (Crossing c : crossings) {
	        if (c == null) {
	            continue;
	        }

	        Point3D.Double boundary = projectToPhi(
	                c.curve().getPoint(clamp01(c.t())),
	                c.value(),
	                theta.radius);

	        /*
	         * Avoid a zero-length connector if numerical weirdness puts a crossing
	         * essentially at the pole.
	         */
	        if (pointsClose(boundary, pole)) {
	            continue;
	        }

	        double phi = MathUtil.normalizeAngle(c.value());

	        phiCurves.add(PhiCurve.onMeridian(boundary, pole, theta.radius, phi));
	        phiCurves.add(PhiCurve.onMeridian(pole, boundary, theta.radius, phi));	    }

	    return phiCurves;
	}
	
	/**
	 * Returns the pole enclosed by a polar theta patch.
	 */
	private static Point3D.Double polarPoint(ThetaPatch theta) {
	    if (theta.enclosesNorthPole()) {
	        return new Point3D.Double(0.0, 0.0, theta.radius);
	    }

	    if (theta.enclosesSouthPole()) {
	        return new Point3D.Double(0.0, 0.0, -theta.radius);
	    }

	    throw new IllegalArgumentException(
	            "polarPoint called for a non-polar theta patch.");
	}

	/**
	 * Splits every original theta-patch curve at its interior phi crossings.
	 *
	 * <p>
	 * The returned list preserves the original boundary order and direction.
	 * </p>
	 */
	private static List<BaseCurve> splitThetaPatchCurves(ThetaPatch theta) {
	    List<BaseCurve> splitCurves = new ArrayList<>();

	    for (BaseCurve curve : theta.curves) {
	        List<? extends BaseCurve> pieces = curve.splitAtPhiCrossings();

	        if (pieces == null || pieces.isEmpty()) {
	            splitCurves.add(curve);
	        } else {
	            splitCurves.addAll(pieces);
	        }
	    }

	    return splitCurves;
	}
	// ---------------------------------------------------------------------
	// Best phi index for no-splice theta patches
	// ---------------------------------------------------------------------

	private static int bestPhiIndex(ThetaPatch theta, SphericalGrid sphGrid) {
		Grid1D phiGrid = sphGrid.getPhiGrid();
		Map<Integer, Integer> counts = new HashMap<>();

		for (BaseCurve curve : theta.curves) {
			double midPhi = curve.phi(0.5);
			int idx = phiGrid.locateInterval(midPhi);
			if (idx >= 0) {
				counts.merge(idx, 1, Integer::sum);
			}
		}

		if (counts.isEmpty()) {
			return phiGrid.locateInterval(theta.curves.get(0).sv0.phi);
		}

		if (counts.size() > 1) {
			System.err.printf("[Patch] bestPhiIndex: midpoints span %d phi cells "
					+ "for ThetPatch (%d,%d,%d, %d); using majority%n",
					counts.size(), theta.nx, theta.ny, theta.nz, theta.nTheta);
		}

		return counts.entrySet().stream().max(Map.Entry.comparingByValue()).get().getKey();
	}

	private record PhiConnectorPair(
	        PhiCurve forward,
	        PhiCurve backward) {
	}

	/**
	 * Removes endpoint phi crossings that are only touches, not genuine phi-band
	 * changes.
	 *
	 * <p>
	 * An endpoint crossing is genuine only if the first non-on-cut curve before the
	 * junction and the first non-on-cut curve after the junction are on opposite
	 * sides of the phi cut. Curves that lie along the phi cut are skipped.
	 * </p>
	 */
	private static List<Crossing> removeEndpointPhiCrossings(
	        List<Crossing> crossings,
	        List<BaseCurve> curves) {

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

	        PhiSideSample before;
	        PhiSideSample after;

	        if (isAtEnd) {
	            before = sampleBeforePhiJunction(curves, curveIdx, cut, probe);
	            after = sampleAfterPhiJunction(curves,
	                    (curveIdx + 1) % curves.size(), cut, probe);
	        } else {
	            before = sampleBeforePhiJunction(curves,
	                    (curveIdx - 1 + curves.size()) % curves.size(),
	                    cut, probe);
	            after = sampleAfterPhiJunction(curves, curveIdx, cut, probe);
	        }

	        if (!before.valid || !after.valid) {
	            filtered.add(c);
	            continue;
	        }

	        if (before.positive == after.positive) {
	            continue;
	        }

	        filtered.add(c);
	    }

	    return filtered;
	}

	private record PhiSideSample(boolean valid, boolean positive) {
	}

	private static PhiSideSample sampleBeforePhiJunction(
	        List<BaseCurve> curves,
	        int startIdx,
	        double cut,
	        double probe) {

	    int n = curves.size();

	    for (int step = 0; step < n; step++) {
	        int idx = (startIdx - step + n) % n;
	        BaseCurve curve = curves.get(idx);

	        if (curveLiesOnPhiCut(curve, cut)) {
	            continue;
	        }

	        double phi = curve.phi(Math.max(0.0, 1.0 - probe));
	        double diff = MathUtil.normalizeAngle(phi - cut);

	        if (Math.abs(diff) < PHI_TOL) {
	            diff = MathUtil.normalizeAngle(curve.phi(0.5) - cut);
	        }

	        /*
	         * Avoid branch-cut ambiguity: this side test is only meaningful near the
	         * meridian being crossed.
	         */
	        if ((Math.abs(diff) < PHI_TOL) || (Math.abs(diff) > Math.PI / 2.0)) {
	            continue;
	        }

	        return new PhiSideSample(true, diff > 0.0);
	    }

	    return new PhiSideSample(false, false);
	}

	private static PhiSideSample sampleAfterPhiJunction(
	        List<BaseCurve> curves,
	        int startIdx,
	        double cut,
	        double probe) {

	    int n = curves.size();

	    for (int step = 0; step < n; step++) {
	        int idx = (startIdx + step) % n;
	        BaseCurve curve = curves.get(idx);

	        if (curveLiesOnPhiCut(curve, cut)) {
	            continue;
	        }

	        double phi = curve.phi(probe);
	        double diff = MathUtil.normalizeAngle(phi - cut);

	        if (Math.abs(diff) < PHI_TOL) {
	            diff = MathUtil.normalizeAngle(curve.phi(0.5) - cut);
	        }

	        if ((Math.abs(diff) < PHI_TOL) || (Math.abs(diff) > Math.PI / 2.0)) {
	            continue;
	        }

	        return new PhiSideSample(true, diff > 0.0);
	    }

	    return new PhiSideSample(false, false);
	}

	private static boolean curveLiesOnPhiCut(BaseCurve curve, double cut) {
	    /*
	     * Do not classify pole samples as lying on a phi cut. Phi is singular there.
	     */
	    if (Math.sin(curve.theta(0.0)) < ON_CUT_TOL
	            || Math.sin(curve.theta(0.5)) < ON_CUT_TOL
	            || Math.sin(curve.theta(1.0)) < ON_CUT_TOL) {
	        return false;
	    }

	    double d0 = Math.abs(MathUtil.normalizeAngle(curve.phi(0.0) - cut));
	    double dm = Math.abs(MathUtil.normalizeAngle(curve.phi(0.5) - cut));
	    double d1 = Math.abs(MathUtil.normalizeAngle(curve.phi(1.0) - cut));

	    return d0 < ON_CUT_TOL && dm < ON_CUT_TOL && d1 < ON_CUT_TOL;
	}

	/**
	 * Builds phi connector pairs from matching phi crossings.
	 *
	 * <p>
	 * Each crossing pair receives two directed meridian connectors: one in each
	 * direction. Unlike theta connectors, there is no short/long azimuthal choice.
	 * A phi connector is a meridian segment parameterized by theta.
	 * </p>
	 */
	private static List<PhiConnectorPair> buildPhiConnectorPairs(
	        ThetaPatch theta,
	        List<Crossing> crossings) {

	    List<Crossing> remainingCrossings = new ArrayList<>(crossings);
	    List<PhiConnectorPair> connectorPairs = new ArrayList<>();

	    while (!remainingCrossings.isEmpty()) {
	        Crossing c = remainingCrossings.remove(0);
	        double phi = c.value();

	        /*
	         * Find the closest unmatched crossing at the same phi value.
	         * This matters when four or more crossings occur on the same meridian.
	         */
	        Crossing match = null;
	        double bestDist = Double.MAX_VALUE;

	        for (Crossing other : remainingCrossings) {
	            if (Math.abs(MathUtil.normalizeAngle(other.value() - phi)) < PHI_TOL) {
	                double dist = c.distanceTo(other);
	                if (dist < bestDist) {
	                    bestDist = dist;
	                    match = other;
	                }
	            }
	        }

	        if (match == null) {
	            System.err.printf(
	                    "[Patch] buildPhiConnectorPairs: unmatched crossing "
	                            + "for theta patch (%d,%d,%d,%d): %s%n",
	                    theta.nx, theta.ny, theta.nz, theta.nTheta,
	                    c.summaryString());
	            continue;
	        }

	        Point3D.Double p0 = projectToPhi(
	                c.curve().getPoint(clamp01(c.t())),
	                phi,
	                theta.radius);

	        Point3D.Double p1 = projectToPhi(
	                match.curve().getPoint(clamp01(match.t())),
	                phi,
	                theta.radius);

	        double fixedPhi = MathUtil.normalizeAngle(phi);

	        PhiCurve forward = PhiCurve.onMeridian(p0, p1, theta.radius, fixedPhi);
	        PhiCurve backward = PhiCurve.onMeridian(p1, p0, theta.radius, fixedPhi);
	        connectorPairs.add(new PhiConnectorPair(forward, backward));

	        remainingCrossings.remove(match);
	    }

	    return connectorPairs;
	}

	/**
	 * Selects both directions of every phi connector pair.
	 */
	private static List<BaseCurve> selectPhiConnectors(
	        List<PhiConnectorPair> connectorPairs) {

	    List<BaseCurve> phiCurves = new ArrayList<>(2 * connectorPairs.size());

	    for (PhiConnectorPair pair : connectorPairs) {
	        phiCurves.add(pair.forward);
	        phiCurves.add(pair.backward);
	    }

	    return phiCurves;
	}

	/**
	 * Assembles split theta-patch curves plus meridian connectors into final closed
	 * patches.
	 */
	private static List<Patch> assemblePatches(
	        ThetaPatch theta,
	        SphericalGrid sphGrid,
	        List<BaseCurve> splitCurves,
	        List<BaseCurve> selectedPhiCurves,
	        List<Crossing> crossings) {

	    List<Patch> result = new ArrayList<>();

	    List<BaseCurve> allCurves = new ArrayList<>(splitCurves);
	    List<BaseCurve> phiCurves = new ArrayList<>(selectedPhiCurves);

	    while (!allCurves.isEmpty()) {
	        BaseCurve c = allCurves.remove(0);

	        List<BaseCurve> patchCurves = new ArrayList<>();
	        patchCurves.add(c);

	        boolean madeLoop = false;

	        while (!madeLoop) {
	            boolean connected = false;
	            boolean isPhi = (c instanceof PhiCurve);

	            /*
	             * Polar special case:
	             *
	             * A polar final patch often has the form:
	             *
	             *   boundary segment -> meridian to pole -> meridian from pole
	             *
	             * So after a PhiCurve ending at the pole, we must allow another
	             * PhiCurve starting at that same pole. Prefer the one that closes
	             * the current loop back to the start point.
	             */
	            if (theta.polar() && atPole(c.p1, theta.radius)) {
	                Point3D.Double loopStart = patchCurves.get(0).p0;

	                /*
	                 * First priority at a pole:
	                 * find an edge that connects from the pole back to the loop start.
	                 *
	                 * That edge may be a generated PhiCurve, but it may also be an original
	                 * theta-patch boundary curve. The previous version only searched phiCurves,
	                 * which allowed polar patches to walk through multiple wedges.
	                 */

	                Iterator<BaseCurve> phiIt = phiCurves.iterator();
	                while (phiIt.hasNext()) {
	                    BaseCurve phi = phiIt.next();

	                    if (connects(c, phi) && pointsClose(phi.p1, loopStart)) {
	                        patchCurves.add(phi);
	                        phiIt.remove();
	                        c = phi;
	                        connected = true;
	                        break;
	                    }
	                }

	                if (!connected) {
	                    Iterator<BaseCurve> allIt = allCurves.iterator();
	                    while (allIt.hasNext()) {
	                        BaseCurve other = allIt.next();

	                        if (connects(c, other) && pointsClose(other.p1, loopStart)) {
	                            patchCurves.add(other);
	                            allIt.remove();
	                            c = other;
	                            connected = true;
	                            break;
	                        }
	                    }
	                }

	                /*
	                 * Second priority:
	                 * if no immediate closure is possible, prefer an ordinary boundary curve
	                 * from the pole before using another generated phi connector.
	                 *
	                 * This matters for polar edge pieces whose boundary is partly an original
	                 * cell/sphere intersection curve rather than a spherical-grid meridian.
	                 */
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

	                /*
	                 * Last resort only: use another generated phi connector from the pole.
	                 * This should now be uncommon. If this happens often, it is worth logging.
	                 */
	                if (!connected) {
	                    phiIt = phiCurves.iterator();
	                    while (phiIt.hasNext()) {
	                        BaseCurve phi = phiIt.next();

	                        if (connects(c, phi)) {
	                            patchCurves.add(phi);
	                            phiIt.remove();
	                            c = phi;
	                            connected = true;
	                            break;
	                        }
	                    }
	                }
	            }

	            /*
	             * Normal case: prefer meridian connectors immediately after non-phi
	             * fragments. This is the analogue of preferring theta connectors
	             * during theta splice.
	             */
	            if (!connected && !isPhi) {
	                Iterator<BaseCurve> phiIt = phiCurves.iterator();
	                while (phiIt.hasNext()) {
	                    BaseCurve phi = phiIt.next();

	                    if (connects(c, phi)) {
	                        patchCurves.add(phi);
	                        phiIt.remove();
	                        c = phi;
	                        connected = true;
	                        break;
	                    }
	                }
	            }

	            /*
	             * Otherwise continue along ordinary split boundary curves.
	             */
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
	                        "[Patch] assemblePatches: could not connect curve %s "
	                                + "in theta patch (%d,%d,%d,%d) with %d crossings%n",
	                        c.shortString(),
	                        theta.nx, theta.ny, theta.nz, theta.nTheta,
	                        crossings.size());
	                break;
	            }

	            madeLoop = makesLoop(patchCurves);

	            if (madeLoop) {
	                try {
	                    int phiIndex = bestPhiIndexForLoop(patchCurves, sphGrid);
	                    Patch patch = new Patch(
	                            patchCurves,
	                            theta.nx, theta.ny, theta.nz,
	                            theta.nTheta, phiIndex);
	                    result.add(patch);
	                } catch (IllegalArgumentException ex) {
	                    System.err.printf("[Patch] assemblePatches: %s%n",
	                            ex.getMessage());
	                }
	            }
	        }
	    }

	    return result;
	}
	
	/**
	 * Returns true if a point is effectively at either pole.
	 */
	private static boolean atPole(Point3D.Double p, double radius) {
	    double rho = Math.hypot(p.x, p.y);
	    return rho < 1.0e-6 * radius
	            && Math.abs(Math.abs(p.z) - radius) < 1.0e-6 * radius;
	}

	private static boolean makesLoop(List<BaseCurve> curves) {
	    if (curves.isEmpty()) {
	        return false;
	    }

	    Point3D.Double start = curves.get(0).p0;
	    Point3D.Double end = curves.get(curves.size() - 1).p1;

	    return pointsClose(start, end);
	}

	private static boolean connects(BaseCurve c1, BaseCurve c2) {
	    return pointsClose(c1.p1, c2.p0);
	}

	/**
	 * Chooses the phi-band index for a closed final-patch loop.
	 *
	 * <p>
	 * The loop should lie inside one phi band, except for any {@link PhiCurve}
	 * segments that lie along meridian grid boundaries. Boundary phi curves are
	 * ignored when possible because their midpoint lies exactly on a grid line and
	 * can be assigned ambiguously to either adjacent phi band.
	 * </p>
	 */
	private static int bestPhiIndexForLoop(
	        List<BaseCurve> loop,
	        SphericalGrid sphGrid) {

	    Grid1D phiGrid = sphGrid.getPhiGrid();
	    Map<Integer, Integer> counts = new HashMap<>();

	    /*
	     * First pass: ignore curves that lie along a phi cut.
	     */
	    for (BaseCurve curve : loop) {
	        if (curve == null) {
	            continue;
	        }

	        double phm = MathUtil.normalizeAngle(curve.phi(0.5));

	        boolean onSomeCut = false;
	        for (int i = 0; i < phiGrid.numPoints(); i++) {
	            double cut = MathUtil.normalizeAngle(phiGrid.valueAt(i));
	            if (curveLiesOnPhiCut(curve, cut)) {
	                onSomeCut = true;
	                break;
	            }
	        }

	        if (onSomeCut) {
	            continue;
	        }

	        int idx = phiGrid.locateInterval(phm);
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
	     * Fallback: vote using all curve midpoints.
	     */
	    for (BaseCurve curve : loop) {
	        if (curve == null) {
	            continue;
	        }

	        double phi = MathUtil.normalizeAngle(curve.phi(0.5));
	        int idx = phiGrid.locateInterval(phi);

	        if (idx < 0) {
	            int vertex = phiGrid.valueIsAVertex(phi);
	            if (vertex >= 0) {
	                if (vertex == 0) {
	                    idx = 0;
	                } else if (vertex >= phiGrid.numPoints() - 1) {
	                    idx = phiGrid.numCells() - 1;
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

	    throw new IllegalStateException(
	            "Could not determine phi index for final patch loop.");
	}

	/**
	 * Projects a point to the exact meridian {@code phi}, preserving theta.
	 */
	private static Point3D.Double projectToPhi(
	        Point3D.Double p,
	        double phi,
	        double radius) {

	    SphericalVector sv = new SphericalVector(p);
	    double sinTheta = Math.sin(sv.theta);

	    return new Point3D.Double(
	            radius * sinTheta * Math.cos(phi),
	            radius * sinTheta * Math.sin(phi),
	            radius * Math.cos(sv.theta));
	}
}
