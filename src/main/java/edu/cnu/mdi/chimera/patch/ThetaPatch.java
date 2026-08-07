package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.Crossing;
import edu.cnu.mdi.chimera.curve.GeneralCurve;
import edu.cnu.mdi.chimera.curve.ThetaCurve;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.model.ChimeraGridContext;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * A theta-patch is a {@link PrePatch} sliced by one band
 * {@code [theta_lo, theta_hi]} of the spherical theta grid.
 *
 * <p>
 * The patch carries a 4-tuple {@code (nx, ny, nz, nTheta)} with
 * {@code nPhi = -1} until the phi splice.
 * </p>
 */
public class ThetaPatch extends BasePatch {

	/*
	 * Polar connector ambiguity is usually small. If there are many crossing pairs,
	 * exhaustive 2^N testing could grow too much, so fall back to trying all-short
	 * and all-long.
	 */
	private static final int MAX_POLAR_CONNECTOR_ENUMERATION = 12;

	/**
	 * Constructs a theta patch with the given curves and grid indices. The curves
	 * should already be split at theta crossings, and the theta index should be
	 * consistent with the curves' theta values.
	 * @param curves    boundary curves for the patch, including theta connectors
	 * @param nx        Cartesian grid x index
	 * @param ny        Cartesian grid y index
	 * @param nz        Cartesian grid z index
	 * @param nTheta    spherical grid theta index
	 *
	 */
	private ThetaPatch(List<BaseCurve> curves, int nx, int ny, int nz, int nTheta) {
		super(curves, nx, ny, nz, nTheta, -1);
		for (BaseCurve c : curves) {
			c.getPhiCrossings(); // precompute phi crossings for later use in phi splice
		}
	}

	@Override
	public boolean containsPoint(double x, double y, double z) {
		CartesianGrid cartesianGrid = ChimeraGridContext.cartesianGrid();
		SphericalGrid sphericalGrid = ChimeraGridContext.sphericalGrid();

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
		SphericalGrid sphGrid = ChimeraGridContext.sphericalGrid();
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
			System.err.println("[ThetaPatch] Warning: odd crossing count after endpoint removal (" + numThetaCrossings
					+ ") for prepatch (" + pre.nx + "," + pre.ny + "," + pre.nz + ").");

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
	 * <p>
	 * An endpoint crossing is genuine only if the first non-on-cut curve before the
	 * junction and the first non-on-cut curve after the junction are on opposite
	 * sides of the theta cut. Curves that lie along the theta cut are skipped.
	 * </p>
	 */
	private static List<Crossing> removeEndpointCrossings(List<Crossing> crossings, List<BaseCurve> curves) {

		final double probe = 0.02;
		List<Crossing> filtered = new ArrayList<>();

		for (Crossing c : crossings) {
			double t = c.t();
			boolean isEndpoint = (t < ENDPOINT_ARTIFACT_TOL || t > 1.0 - ENDPOINT_ARTIFACT_TOL);

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
				after = sampleAfterJunction(curves, (curveIdx + 1) % curves.size(), cut, probe);
			} else {
				before = sampleBeforeJunction(curves, (curveIdx - 1 + curves.size()) % curves.size(), cut, probe);
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

	private record ThetaConnectorPair(ThetaCurve shortForward, ThetaCurve shortBackward, ThetaCurve longForward,
			ThetaCurve longBackward) {
	}

	// Helper class to represent the result of sampling a curve near a junction.
	private record SideSample(boolean valid, boolean above) {
	}

	// Samples the curve list in the specified direction from the junction,
	// skipping curves that lie along the theta cut.
	private static SideSample sampleBeforeJunction(List<BaseCurve> curves, int startIdx, double cut, double probe) {

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

	// Similar to sampleBeforeJunction but samples in the forward direction.
	private static SideSample sampleAfterJunction(List<BaseCurve> curves, int startIdx, double cut, double probe) {

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
		return Math.abs(curve.theta(0.0) - cut) < THETA_TOL && Math.abs(curve.theta(0.5) - cut) < THETA_TOL
				&& Math.abs(curve.theta(1.0) - cut) < THETA_TOL;
	}

	// ---------------------------------------------------------------------
	// Best theta index for no-splice prepatches
	// ---------------------------------------------------------------------

	private static int bestThetaIndex(PrePatch pre, SphericalGrid sphGrid) {
		Grid1D thetaGrid = sphGrid.getThetaGrid();
		Map<Integer, Integer> counts = new HashMap<>();

		for (BaseCurve curve : pre.curves) {
			boolean onThetaCut = false;
			for (int i = 0; i < thetaGrid.numPoints(); i++) {
				if (curveLiesOnThetaCut(curve, thetaGrid.valueAt(i))) {
					onThetaCut = true;
					break;
				}
			}
			if (onThetaCut) {
				continue;
			}

			double midTheta = curve.theta(0.5);
			int idx = thetaGrid.locateInterval(midTheta);
			if (idx >= 0) {
				counts.merge(idx, 1, Integer::sum);
			}
		}

		if (counts.isEmpty()) {
			return thetaGrid.locateInterval(pre.curves.get(0).sv0.theta);
		}

		return counts.entrySet().stream().max(Map.Entry.comparingByValue()).get().getKey();
	}

	/*
	 * Core theta splice logic: builds candidate theta connectors, tries different
	 * combinations if the prepatch is polar, and assembles loops into theta
	 * patches.
	 */
	private static List<ThetaPatch> doSplice(PrePatch pre, SphericalGrid sphGrid, List<Crossing> crossings) {

		List<BaseCurve> splitCurves = splitPrePatchCurves(pre);

		List<ThetaConnectorPair> connectorPairs = buildThetaConnectorPairs(pre, crossings);

		/*
		 * Non-polar prepatches should use ordinary short theta connectors.
		 */
		if (!pre.polar() || connectorPairs.isEmpty()) {
			List<BaseCurve> thetaCurves = selectThetaConnectors(connectorPairs, 0L, false);
			return assembleThetaPatches(pre, sphGrid, splitCurves, thetaCurves, crossings);
		}

		/*
		 * Polar prepatches can be ambiguous because phi is singular at the pole. Try
		 * connector choices and keep the candidate whose total area best matches the
		 * original prepatch area.
		 */
		double targetArea = pre.areaEstimate();
		List<ThetaPatch> best = List.of();
		double bestError = Double.POSITIVE_INFINITY;

		int n = connectorPairs.size();

		if (n <= MAX_POLAR_CONNECTOR_ENUMERATION) {
			long numMasks = 1L << n;

			for (long mask = 0L; mask < numMasks; mask++) {
				List<BaseCurve> thetaCurves = selectThetaConnectors(connectorPairs, mask, false);

				List<ThetaPatch> candidate = assembleThetaPatches(pre, sphGrid, splitCurves, thetaCurves, crossings);

				if (candidate.isEmpty()) {
					continue;
				}

				double area = totalArea(candidate);
				double error = Math.abs(area - targetArea);

				if (error < bestError) {
					bestError = error;
					best = candidate;
				}
			}
		} else {
			/*
			 * Fallback for unexpectedly many connector pairs: try all-short and all-long.
			 */
			List<BaseCurve> shortThetaCurves = selectThetaConnectors(connectorPairs, 0L, false);

			List<ThetaPatch> shortCandidate = assembleThetaPatches(pre, sphGrid, splitCurves, shortThetaCurves,
					crossings);

			if (!shortCandidate.isEmpty()) {
				double error = Math.abs(totalArea(shortCandidate) - targetArea);
				if (error < bestError) {
					bestError = error;
					best = shortCandidate;
				}
			}

			List<BaseCurve> longThetaCurves = selectThetaConnectors(connectorPairs, 0L, true);

			List<ThetaPatch> longCandidate = assembleThetaPatches(pre, sphGrid, splitCurves, longThetaCurves,
					crossings);

			if (!longCandidate.isEmpty()) {
				double error = Math.abs(totalArea(longCandidate) - targetArea);
				if (error < bestError) {
					bestError = error;
					best = longCandidate;
				}
			}
		}

		return best;
	}

	/**
	 * Splits every original prepatch curve at its interior theta crossings.
	 *
	 * <p>
	 * The returned list preserves the original boundary order and direction.
	 * </p>
	 */
	private static List<BaseCurve> splitPrePatchCurves(PrePatch pre) {
		List<BaseCurve> splitCurves = new ArrayList<>();

		for (BaseCurve curve : pre.curves) {
			if (curve instanceof GeneralCurve) {
				GeneralCurve gcurve = (GeneralCurve) curve;
				splitCurves.addAll(gcurve.splitAtThetaCrossings());
			} else {
				splitCurves.add(curve);
			}
		}

		return splitCurves;
	}

	/**
	 * Builds candidate theta connector pairs from matching theta crossings.
	 *
	 * <p>
	 * Each crossing pair gets two possible geometries:
	 * </p>
	 *
	 * <ul>
	 * <li>short arc, both directions</li>
	 * <li>long arc, both directions</li>
	 * </ul>
	 *
	 * <p>
	 * The backward connector is included because loop assembly may need either
	 * direction.
	 * </p>
	 */
	private static List<ThetaConnectorPair> buildThetaConnectorPairs(PrePatch pre, List<Crossing> crossings) {

		List<Crossing> remainingCrossings = new ArrayList<>(crossings);
		List<ThetaConnectorPair> connectorPairs = new ArrayList<>();

		while (!remainingCrossings.isEmpty()) {
			Crossing c = remainingCrossings.remove(0);
			double theta = c.value();

			/*
			 * Find the closest unmatched crossing at the same theta value. This is
			 * important when there are four or more crossings on the same theta grid line.
			 */
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

			if (match == null) {
				throw new IllegalStateException(String.format(
						"Unmatched theta crossing for prepatch (%d,%d,%d): %s",
						pre.nx, pre.ny, pre.nz, c.summaryString()));
			}

			Point3D.Double p0 = c.curve().getPoint(clamp01(c.t()));
			Point3D.Double p1 = match.curve().getPoint(clamp01(match.t()));

			ThetaCurve shortForward = ThetaCurve.between(p0, p1, pre.radius, false);
			ThetaCurve shortBackward = ThetaCurve.between(p1, p0, pre.radius, false);

			ThetaCurve longForward = ThetaCurve.between(p0, p1, pre.radius, true);
			ThetaCurve longBackward = ThetaCurve.between(p1, p0, pre.radius, true);

			connectorPairs.add(new ThetaConnectorPair(shortForward, shortBackward, longForward, longBackward));

			remainingCrossings.remove(match);
		}

		return connectorPairs;
	}

	/**
	 * Selects one geometry for every theta connector pair.
	 *
	 * @param connectorPairs available connector-pair options
	 * @param longMask       bit mask: bit i set means pair i uses the long arc
	 * @param forceAllLong   if true, ignore {@code longMask} and make every pair
	 *                       long
	 * @return directed theta curves for loop assembly
	 */
	private static List<BaseCurve> selectThetaConnectors(List<ThetaConnectorPair> connectorPairs, long longMask,
			boolean forceAllLong) {

		List<BaseCurve> thetaCurves = new ArrayList<>(2 * connectorPairs.size());

		for (int i = 0; i < connectorPairs.size(); i++) {
			ThetaConnectorPair pair = connectorPairs.get(i);

			boolean useLong = forceAllLong || ((longMask & (1L << i)) != 0L);

			if (useLong) {
				thetaCurves.add(pair.longForward);
				thetaCurves.add(pair.longBackward);
			} else {
				thetaCurves.add(pair.shortForward);
				thetaCurves.add(pair.shortBackward);
			}
		}

		return thetaCurves;
	}

	/**
	 * Assembles split prepatch curves plus selected theta connectors into closed
	 * theta patches.
	 *
	 * <p>
	 * This is intentionally destructive on local copies only. It does not mutate
	 * the supplied {@code splitCurves} or {@code selectedThetaCurves}.
	 * </p>
	 */
	private static List<ThetaPatch> assembleThetaPatches(PrePatch pre, SphericalGrid sphGrid,
			List<BaseCurve> splitCurves, List<BaseCurve> selectedThetaCurves, List<Crossing> crossings) {

		List<ThetaPatch> result = new ArrayList<>();

		List<BaseCurve> allCurves = new ArrayList<>(splitCurves);
		List<BaseCurve> thetaCurves = new ArrayList<>(selectedThetaCurves);

		while (!allCurves.isEmpty()) {
			BaseCurve c = allCurves.remove(0);

			List<BaseCurve> patchCurves = new ArrayList<>();
			patchCurves.add(c);

			boolean madeLoop = false;

			while (!madeLoop) {
				boolean connected = false;
				boolean isTheta = (c instanceof ThetaCurve);

				/*
				 * Prefer theta connectors immediately after non-theta fragments. This tends to
				 * close a theta-sliced boundary before walking into the neighboring theta band.
				 */
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
					throw new IllegalStateException(String.format(
							"Could not connect curve %s in prepatch (%d,%d,%d) "
									+ "with %d crossings.",
							c.shortString(), pre.nx, pre.ny, pre.nz, crossings.size()));
				}

				madeLoop = makesLoop(patchCurves);

				if (madeLoop) {
					int thetaIndex = bestThetaIndexForLoop(patchCurves, sphGrid);
					ThetaPatch patch = new ThetaPatch(
							patchCurves, pre.nx, pre.ny, pre.nz, thetaIndex);
					result.add(patch);
				}
			}
		}

		return result;
	}

	/**
	 * Computes the total normalized area of a theta-patch list.
	 */
	private static double totalArea(List<ThetaPatch> patches) {
		double sum = 0.0;

		if (patches == null) {
			return sum;
		}

		for (ThetaPatch patch : patches) {
			if (patch != null) {
				sum += patch.areaEstimate();
			}
		}

		return sum;
	}

	// Checks whether the given curves form a closed loop by comparing the
	// start of the first curve and the end of the last curve.
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

		return pointsClose(c1End, c2Start);

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
	 * midpoints and using {@link Grid1D#locateInterval(double)}. This should be
	 * rare, but it keeps the method total.
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
		 * First pass: vote only with curves that are not lying on a theta cut. These
		 * are the most reliable indicators of the interior theta band.
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
				if (Math.abs(th0 - cut) < ON_CUT_TOL && Math.abs(thm - cut) < ON_CUT_TOL
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
			return counts.entrySet().stream().max(Map.Entry.comparingByValue()).get().getKey();
		}

		/*
		 * Fallback: all curves looked like theta-boundary curves. Use midpoint samples
		 * anyway, nudging exact grid-line values very slightly inward when needed.
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
						 * Ambiguous interior cut. Pick the lower adjacent band as a deterministic
						 * fallback.
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
			return counts.entrySet().stream().max(Map.Entry.comparingByValue()).get().getKey();
		}

		throw new IllegalStateException("Could not determine theta index for theta-splice loop.");
	}

	/**
	 * Retrieves a list of all phi crossings across all curves in this theta patch.
	 * Used to draw markers at phi crossings during the phi splice step.
	 * @return list of all theta crossings in this prepatch
	 */
	public List<Crossing> getAllPhiCrossings() {
		List<Crossing> crossings = new ArrayList<>();
		for (BaseCurve curve : curves) {
			crossings.addAll(curve.getPhiCrossings());
		}
		return Crossing.removeDuplicates(crossings);
	}

}
