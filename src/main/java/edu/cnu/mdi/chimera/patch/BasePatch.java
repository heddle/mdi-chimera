package edu.cnu.mdi.chimera.patch;

import java.util.ArrayList;
import java.util.List;

import edu.cnu.mdi.chimera.area.SphericalPolygonArea;
import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.BaseCurve.PoleStatus;
import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.model.ChimeraGridContext;
import edu.cnu.mdi.chimera.util.Point3D;
import edu.cnu.mdi.chimera.util.ThetaPhi;

/**
 * Abstract base class for all Mosaic patches.
 *
 * <p>
 * A patch is a region on the sphere surface enclosed by a closed loop of
 * {@link BaseCurve}s, uniquely identified by a 5-tuple {@code (nx, ny, nz)}
 * from the Cartesian grid and {@code (nTheta, nPhi)} from the spherical grid.
 * Prepatches carry only the Cartesian 3-tuple until the splice steps assign the
 * spherical indices.
 * </p>
 *
 * <h2>Area</h2>
 * <p>
 * Estimated via {@link SphericalPolygonArea} using fan triangulation with the
 * Oosterom-Strackee formula. Vertices are sampled directly as
 * {@link Point3D.Double} Cartesian points from {@link BaseCurve#getPoint},
 * bypassing {@link edu.cnu.mdi.chimera.util.ThetaPhi} entirely to avoid any
 * theta/latitude or degree/radian conversion errors.
 * </p>
 *
 * <h2>Pole enclosure</h2>
 * <p>
 * Computed once and cached using the winding-number algorithm from paper
 * Appendix A, including the inter-curve phi jumps at each junction.
 * </p>
 */
public abstract class BasePatch implements Comparable<BasePatch> {

    /*
     * Loop-building tolerance must be larger than BasePatch.LOOP_TOL because
     * we sometimes project crossing endpoints onto exact theta grid lines.
     * The assembled curves themselves are constructed with matching projected
     * endpoints, so BasePatch validation should still pass with its tighter
     * tolerance.
     */
    private static final double LOOP_BUILD_TOL = 1.0e-6;


	/** Samples per curve for area estimation. The paper uses n=5. */
	protected static final int DEFAULT_AREA_SAMPLES = 50;

	private static final double WINDING_TOL = 0.01; // accumulated rounding over 100 steps/curve

	// -----------------------------------------------------------------------
	// Fields
	// -----------------------------------------------------------------------

	public final List<BaseCurve> curves;
	public final double radius;
	public final int nx;
	public final int ny;
	public final int nz;
	public final int nTheta;
	public final int nPhi;

	private Boolean _enclosesNorthPole = null;
	private Boolean _enclosesSouthPole = null;


	// -----------------------------------------------------------------------
	// Construction
	// -----------------------------------------------------------------------

	public BasePatch(List<BaseCurve> curves, int nx, int ny, int nz, int nTheta, int nPhi) {

		if (curves == null || curves.size() < 2)
			throw new IllegalArgumentException("A patch requires at least two curves.");

		this.curves = curves;
		this.radius = ChimeraGridContext.radius();
		this.nx = nx;
		this.ny = ny;
		this.nz = nz;
		this.nTheta = nTheta;
		this.nPhi = nPhi;

		if (!validateLoop(curves)) {
			throw new IllegalArgumentException("Curves do not form a closed loop for patch (" + nx + "," + ny + "," + nz
					+ "," + nTheta + "," + nPhi + ").");
		}
	}

	/** Convenience constructor for prepatches (nTheta = nPhi = -1). */
	public BasePatch(List<BaseCurve> curves, int nx, int ny, int nz) {
		this(curves, nx, ny, nz, -1, -1);
	}
	
	// -----------------------------------------------------------------------
	// Abstract interface
	// -----------------------------------------------------------------------

	public abstract boolean containsPoint(double x, double y, double z);

	// -----------------------------------------------------------------------
	// Loop validation
	// -----------------------------------------------------------------------

	public static boolean validateLoop(List<BaseCurve> curves) {
		int n = curves.size();
		for (int i = 0; i < n; i++) {
			Point3D.Double end = curves.get(i).p1;
			Point3D.Double start = curves.get((i + 1) % n).p0;
			double dist = Point3D.Double.distance(end, start);
			if (dist > LOOP_BUILD_TOL) {
				System.err.printf("[BasePatch] Loop gap at curve %d→%d: distance=%.3e%n", i, (i + 1) % n, dist);
				return false;
			}
		}
		return true;
	}

	// -----------------------------------------------------------------------
	// Pole enclosure 
	// -----------------------------------------------------------------------
	/**
	 * Returns true if this patch encloses the north pole, false otherwise.
	 * The result is computed lazily and cached for future calls.
	 *
	 * @return true if this patch encloses the north pole, false otherwise
	 */
	public boolean enclosesNorthPole() {
		if (_enclosesNorthPole == null)
			computePoleEnclosure();
		return _enclosesNorthPole;
	}

	/**
	 * Returns true if this patch encloses the south pole, false otherwise.
	 * The result is computed lazily and cached for future calls.
	 *
	 * @return true if this patch encloses the south pole, false otherwise
	 */
	public boolean enclosesSouthPole() {
		if (_enclosesSouthPole == null)
			computePoleEnclosure();
		return _enclosesSouthPole;
	}

	/**
	 * Returns true if this patch encloses either pole, i.e. if it is a "polar patch" as defined in the paper.
	 * @return true if this patch encloses the north pole or the south pole, false otherwise
	 */
	public boolean polar() {
		return enclosesNorthPole() || enclosesSouthPole();
	}

	private void computePoleEnclosure() {
		double totalWinding = 0.0;
		double thetaSum = 0.0;
		int thetaCount = 0;

		int n = curves.size();
		for (int i = 0; i < n; i++) {
			BaseCurve curve = curves.get(i);

			PoleStatus onCurve = curve.poleOnCurve();
			if (onCurve == PoleStatus.NORTH_ON_CURVE) {
				_enclosesNorthPole = true;
				_enclosesSouthPole = false;
				return;
			}
			if (onCurve == PoleStatus.SOUTH_ON_CURVE) {
				_enclosesNorthPole = false;
				_enclosesSouthPole = true;
				return;
			}

			totalWinding += curve.windingContribution();

			// Inter-curve phi jump at each junction.
			BaseCurve next = curves.get((i + 1) % n);
			double phiEnd = MathUtil.normalizeAngle(curve.sv1.phi);
			double phiStart = MathUtil.normalizeAngle(next.sv0.phi);
			totalWinding += MathUtil.normalizeAngle(phiStart - phiEnd);

			thetaSum += curve.averageTheta();
			thetaCount++;
		}

		double avgTheta = (thetaCount > 0) ? thetaSum / thetaCount : Math.PI / 2;

		if (Math.abs(Math.abs(totalWinding) - 2.0 * Math.PI) < WINDING_TOL) {
			_enclosesNorthPole = avgTheta < Math.PI / 2;
			_enclosesSouthPole = !_enclosesNorthPole;
		} else {
			_enclosesNorthPole = false;
			_enclosesSouthPole = false;
		}
	}

	// -----------------------------------------------------------------------
	// Geometry
	// -----------------------------------------------------------------------

	/**
	 * Returns a sampled list of boundary vertices as {@link ThetaPhi} points,
	 * suitable for drawing.
	 *
	 * <p>
	 * Delegates to {@link #getBoundaryPoints(int)} and converts each Cartesian
	 * point to spherical coordinates.
	 * </p>
	 *
	 * @param n samples per curve (≥ 1)
	 * @return ordered list of {@link ThetaPhi} boundary points
	 */
	public List<edu.cnu.mdi.chimera.util.ThetaPhi> getSphericalVertices(int n) {
		List<Point3D.Double> pts = getBoundaryPoints(n);
		List<edu.cnu.mdi.chimera.util.ThetaPhi> result = new ArrayList<>(pts.size());
		for (Point3D.Double p : pts) {
			edu.cnu.mdi.chimera.util.SphericalVector sv = new edu.cnu.mdi.chimera.util.SphericalVector(p);
			result.add(new edu.cnu.mdi.chimera.util.ThetaPhi(radius, sv.theta, sv.phi));
		}
		return result;
	}

	/**
	 * Returns boundary vertices using the default sample count
	 * ({@value #DEFAULT_AREA_SAMPLES} samples per curve).
	 *
	 * @return ordered list of {@link edu.cnu.mdi.chimera.util.ThetaPhi} boundary
	 *         points
	 */
	public List<ThetaPhi> getSphericalVertices() {
		return getSphericalVertices(DEFAULT_AREA_SAMPLES);
	}

	/**
	 * Samples the boundary as Cartesian {@link Point3D.Double} points on the sphere
	 * surface — one point per step per curve, at {@code t = step, 2·step, ..., 1}
	 * (i.e. i = 1..n inclusive).
	 *
	 * <p>
	 * Starting at {@code t = step} rather than {@code t = 0} avoids duplicating the
	 * junction point shared with the previous curve's end. Each junction appears
	 * exactly once, as the {@code t = 1} sample of the curve that ends there.
	 * </p>
	 *
	 * @param n samples per curve (≥ 1)
	 * @return ordered list of boundary points
	 */
	public List<Point3D.Double> getBoundaryPoints(int n) {
		if (n < 1)
			throw new IllegalArgumentException("n must be at least 1.");
		List<Point3D.Double> pts = new ArrayList<>(curves.size() * n);
		double step = 1.0 / n;
		for (BaseCurve curve : curves) {
			for (int i = 1; i <= n; i++) {
				pts.add(curve.getPoint(i * step));
			}
		}
		return pts;
	}

	/**
	 * Estimates the normalized area of this patch as a fraction of {@code 4πR²}.
	 *
	 * @param n samples per curve
	 * @return normalized area in {@code [0, 1]}
	 */
	public double areaEstimate(int n) {
		List<Point3D.Double> pts = getBoundaryPoints(n);
		return SphericalPolygonArea.computeAreaFraction(pts);
	}

	/** 
	 * Area estimate using {@value #DEFAULT_AREA_SAMPLES} samples per curve.
	 * @return the normalized area estimate 
	 */
	public double areaEstimate() {
		return areaEstimate(DEFAULT_AREA_SAMPLES);
	}

	/**
	 * Computes the perimeter of this patch as the sum of the arc lengths of its
	 * boundary curves.
	 *
	 * @return the perimeter length in the same units as the curve arc lengths
	 */
	public double perimeter() {
		double total = 0.0;
		for (BaseCurve curve : curves)
			total += curve.arcLength();
		return total;
	}

	// -----------------------------------------------------------------------
	// Index accessors
	// -----------------------------------------------------------------------

	/**
	 * Returns true if this patch has valid spherical grid indices (nTheta and nPhi are non-negative).
	 *
	 * @return true if both nTheta and nPhi are non-negative, false otherwise
	 */
	public boolean isFullyIndexed() {
		return nTheta >= 0 && nPhi >= 0;
	}

	/**
	 * Returns a string representation of the Cartesian 3-tuple index of this patch.
	 *
	 * @return a string in the format "(nx, ny, nz)"
	 */
	public String cartesianIndex() {
		return String.format("(%d,%d,%d)", nx, ny, nz);
	}

	/**
	 * Returns a string representation of the full 5-tuple index of this patch, including
	 * both Cartesian and spherical indices. If the spherical indices are not assigned (i.e. nTheta or nPhi is negative),
	 * the method will still include the Cartesian indices and indicate that the spherical indices are unassigned.
	 * For example, if nTheta and nPhi are both -1, the method might return "(nx, ny, nz) (unassigned)". If nTheta and nPhi are assigned, it will return "(nx, ny, nz) (nTheta, nPhi)".
	 * 
	 * @return
	 */
	public String fullIndex() {
		return String.format("(%d,%d,%d) (%d,%d)", nx, ny, nz, nTheta, nPhi);
	}

	@Override
	public String toString() {
		return String.format("%s[%s curves=%d polar=N%s S%s]", getClass().getSimpleName(), fullIndex(), curves.size(),
				enclosesNorthPole() ? "✓" : "✗", enclosesSouthPole() ? "✓" : "✗");
	}

	@Override
	public int compareTo(BasePatch other) {
		int cmp = Integer.compare(this.nx, other.nx);
		if (cmp != 0)
			return cmp;
		cmp = Integer.compare(this.ny, other.ny);
		if (cmp != 0)
			return cmp;
		cmp = Integer.compare(this.nz, other.nz);
		if (cmp != 0)
			return cmp;
		cmp = Integer.compare(this.nTheta, other.nTheta);
		if (cmp != 0)
			return cmp;
		return Integer.compare(this.nPhi, other.nPhi);
	}

	/**
	 * Performs a binary search on the sorted list of patches to find the patch
	 * with the specified indices. The list must be sorted in ascending order by
	 * Cartesian indices (nx, ny, nz) and then by spherical indices (nTheta, nPhi).
	 *
	 * @param patches the sorted list of patches to search
	 * @param nx      the Cartesian x-index of the patch to find
	 * @param ny      the Cartesian y-index of the patch to find
	 * @param nz      the Cartesian z-index of the patch to find
	 * @param nTheta  the spherical theta-index of the patch to find
	 * @param nPhi    the spherical phi-index of the patch to find
	 * @return the patch with the specified indices, or null if not found
	 */
	public static BasePatch fromSortedList(List<? extends BasePatch> patches, int nx, int ny, int nz, int nTheta,
			int nPhi) {
		if (patches == null || patches.isEmpty()) {
			return null;
		}
		

		int lo = 0, hi = patches.size() - 1;
		while (lo <= hi) {
			int mid = (lo + hi) >>> 1;
			BasePatch p = patches.get(mid);

			int cmp = Integer.compare(p.nx, nx);
			if (cmp == 0)
				cmp = Integer.compare(p.ny, ny);
			if (cmp == 0)
				cmp = Integer.compare(p.nz, nz);
			if (cmp == 0)
				cmp = Integer.compare(p.nTheta, nTheta);
			if (cmp == 0)
				cmp = Integer.compare(p.nPhi, nPhi);

			if (cmp < 0)
				lo = mid + 1;
			else if (cmp > 0)
				hi = mid - 1;
			else
				return p;
		}
		return null;
	}
	

    /**
	 * Determines whether two points are close enough to be considered the same for loop-building purposes.
	 *
	 * <p>
	 * This method uses a tolerance defined by {@code LOOP_BUILD_TOL} to account for minor discrepancies
	 * in curve endpoints that may arise from numerical precision issues during curve construction and
	 * projection. If the distance between the two points is less than this tolerance, they are considered
	 * close enough to be treated as the same point when validating that curves form a closed loop.
	 * </p>
	 *
	 * @param a the first point to compare
	 * @param b the second point to compare
	 * @return true if the points are close enough to be considered the same, false otherwise
	 */
    static boolean pointsClose(Point3D.Double a, Point3D.Double b) {
        return Point3D.Double.distance(a, b) < LOOP_BUILD_TOL;
    }

    /**
	 * Returns the index of the target curve in the list, or -1 if not found.
	 *
	 * @param curves the list of curves to search
	 * @param target the curve to find
	 * @return the index of the target curve, or -1 if not found
	 */
    static int curveIndexOf(List<BaseCurve> curves, BaseCurve target) {
        for (int i = 0; i < curves.size(); i++) {
            if (curves.get(i) == target) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Clamps the input value to the range [0, 1]. If t is less than 0, returns 0. If t is greater than 1, returns 1. Otherwise, returns t.
     * @param t the input value to clamp
     * @return the clamped value in the range [0, 1]
     */
    static double clamp01(double t) {
        return Math.max(0.0, Math.min(1.0, t));
    }

}
