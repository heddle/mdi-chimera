package edu.cnu.mdi.chimera.curve;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.UnivariateFunction;

import edu.cnu.mdi.chimera.grid.Grid1D;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.model.ChimeraGridContext;
import edu.cnu.mdi.chimera.util.ChimeraPlane;
import edu.cnu.mdi.chimera.util.ChimeraRotation;
import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * A GENERAL curve: the intersection of the sphere surface with one face of a
 * rectangular grid cell (paper §5).
 *
 * <p>
 * A GENERAL curve lies simultaneously on the sphere and exactly one cell face.
 * Its endpoints are the two points where cell edges pierce the sphere surface.
 * The parameterization is constructed via the paper's rotation algorithm:
 * </p>
 * <ol>
 * <li>Build the face plane P from three face corners and obtain its normal n̂
 * and rotation matrix R that maps ẑ → n̂.</li>
 * <li>Rotate both endpoints into the primed frame. In that frame they share a
 * common z' = const, hence a common θ* = θ'.</li>
 * <li>Read off φ'₀ and φ'₁; normalize Δφ' to (−π, π] to take the shorter arc
 * (branch-cut fix).</li>
 * <li>Parameterize in the primed frame: θ'(t) = θ*, φ'(t) = φ'₀ + t·Δφ'.</li>
 * <li>Rotate each point back with R⁻¹ to recover (R, θ(t), φ(t)) in the
 * original frame.</li>
 * </ol>
 *
 * <p>
 * Steps 1–3 are performed once at construction time; the returned
 * {@link UnivariateFunction}s for θ and φ apply step 5 on every evaluation,
 * which is O(1) per sample (one matrix–vector multiply and two trig calls).
 * </p>
 */
public class GeneralCurve extends BaseCurve {
	
	/**
	 * Tolerance for suppressing theta crossings that occur at curve endpoints.
	 *
	 * <p>Endpoint hits are junction events, not curve-interior crossings. They are
	 * handled, if needed, by the patch-level splice logic that can inspect the
	 * neighboring curves.</p>
	 */
	private static final double ENDPOINT_CROSSING_TOL = 1.0e-8;

	/**
	 * Small offset used to test whether a candidate crossing actually changes
	 * theta band. This filters tangencies and endpoint grazes.
	 */
	private static final double SIDE_TEST_DT = 1.0e-6;
	
	/**
	 * The inverse rotation matrix (3×3) that maps the primed frame back to the
	 * original frame. Stored as a flat array for efficiency.
	 */
	private final double[][] invMatrix;

	/**
	 * The constant polar angle θ* in the primed (rotated) frame. Both endpoints
	 * share this value by construction.
	 */
	private final double thetaStar;

	/** φ'₀: azimuthal angle of the start point in the primed frame. */
	private final double primePhi0;

	/** Δφ': signed azimuthal sweep in the primed frame, normalised to (−π, π]. */
	private final double deltaPrimePhi;

	/** Cached list of crossings where the curve crosses theta grid lines. */
	private List<Crossing> thetaCrossings;
	
	/** Cached list of crossings where the curve crosses phi grid lines. */
	private List<Crossing> phiCrossings;

	// -----------------------------------------------------------------------
	// Construction
	// -----------------------------------------------------------------------

	/**
	 * Constructs a GENERAL curve using the paper's rotation algorithm.
	 *
	 * @param p0      start point on the sphere (edge–sphere intersection)
	 * @param p1      end point on the sphere (edge–sphere intersection)
	 * @param r       sphere radius
	 * @param corner0 any corner of the cell face containing this curve
	 * @param corner1 a second corner of the same face
	 * @param corner2 a third corner of the same face
	 * @throws IllegalArgumentException if the three corners are collinear
	 */
	public GeneralCurve(Point3D.Double p0, Point3D.Double p1, double r, Point3D.Double corner0, Point3D.Double corner1,
			Point3D.Double corner2) {
		super(p0, p1, r);

		// ------------------------------------------------------------------
		// Steps 3–5 (paper): build plane and rotation matrix R.
		// ChimeraPlane computes the face normal and the rotation that maps
		// ẑ onto that normal (equivalently, maps the normal to ẑ in the
		// inverse direction).
		// ------------------------------------------------------------------
		ChimeraPlane plane = new ChimeraPlane(corner0, corner1, corner2);

		// rmat.matrix maps original → primed (ẑ aligns with face normal).
		// rmat.invMatrix maps primed → original.
		double[][] rotMatrix = plane.rmat.matrix;
		invMatrix = plane.rmat.invMatrix;

		// ------------------------------------------------------------------
		// Steps 6–7 (paper): rotate both endpoints into the primed frame
		// and read off θ* and φ'₀, φ'₁.
		// ------------------------------------------------------------------
		double[] v0prime = ChimeraRotation.multiplyMatrixVector(rotMatrix, new double[] { p0.x, p0.y, p0.z });
		double[] v1prime = ChimeraRotation.multiplyMatrixVector(rotMatrix, new double[] { p1.x, p1.y, p1.z });

		// Both rotated points are on the sphere; their z' values should
		// agree. Use their average as θ* for numerical stability.
		double zPrime0 = v0prime[2];
		double zPrime1 = v1prime[2];

		if (Math.abs(zPrime0 - zPrime1) > TOL * r) {
			System.err.printf("[GeneralCurve] z' values differ after rotation: %.6f vs %.6f%n", zPrime0, zPrime1);
		}

		double zPrime = 0.5 * (zPrime0 + zPrime1);
		// Clamp to [-R, R] before acos to guard against tiny numerical overruns.
		thetaStar = Math.acos(Math.max(-1.0, Math.min(1.0, zPrime / r)));

		// ------------------------------------------------------------------
		// Step 8 (paper, with branch-cut fix): compute Δφ' and normalise.
		// ------------------------------------------------------------------
		primePhi0 = MathUtil.normalizeAngle(Math.atan2(v0prime[1], v0prime[0]));
		double phi1p = MathUtil.normalizeAngle(Math.atan2(v1prime[1], v1prime[0]));
		deltaPrimePhi = MathUtil.normalizeAngle(phi1p - primePhi0);

		thetaCrossings = getThetaCrossings(); // Precompute crossings for efficiency
	}

	/**
	 * Checks if this curve is effectively a constant-θ curve, which can happen when
	 * the endpoints are very close together or when the face plane is nearly
	 * tangent to the sphere. In such cases, the curve is essentially a small circle
	 * and can be treated as a ThetaCurve for efficiency.
	 *
	 * @return true if this curve is effectively constant-θ, false otherwise
	 */
	@Override
	public boolean isConstantTheta() {
		/*
		 * A GeneralCurve has constant original-frame theta iff z(t) is constant. Since
		 * theta = acos(z/R), this is equivalent to saying that the original-frame
		 * z-coordinate does not vary along the primed small-circle arc.
		 *
		 * Rather than relying only on deltaPrimePhi == 0, sample the actual theta
		 * function. This correctly catches horizontal face arcs, where theta is
		 * constant even though the primed azimuthal sweep is nonzero.
		 */
		double th0 = theta(0.0);
		double th1 = theta(1.0);
		double thm = theta(0.5);

		double min = Math.min(th0, Math.min(th1, thm));
		double max = Math.max(th0, Math.max(th1, thm));

		return (max - min) < TOL;
	}

	/**
	 * Private constructor used by {@link #reverse()}.
	 */
	private GeneralCurve(Point3D.Double p0, Point3D.Double p1, double r, double[][] invMatrix, double thetaStar,
			double primePhi0, double deltaPrimePhi) {
		super(p0, p1, r);
		this.invMatrix = invMatrix;
		this.thetaStar = thetaStar;
		this.primePhi0 = primePhi0;
		this.deltaPrimePhi = deltaPrimePhi;
	}

	/**
	 * Reconstructs a directed general curve from its canonical small-circle
	 * parameterization.
	 *
	 * <p>The vectors {@code basisU} and {@code basisV} form the first two columns
	 * of the primed-to-original orthonormal rotation. Their cross product supplies
	 * the circle-plane normal. The signed center displacement along that normal
	 * determines the primed colatitude.</p>
	 *
	 * @param p0 exported start point
	 * @param p1 exported end point
	 * @param sphereRadius radius of the containing sphere
	 * @param center center of the small circle
	 * @param basisU first in-plane unit basis vector
	 * @param basisV second in-plane unit basis vector
	 * @param alpha0 starting circle angle
	 * @param deltaAlpha signed angular sweep
	 * @return reconstructed directed general curve
	 */
	public static GeneralCurve fromCircleParameters(Point3D.Double p0,
			Point3D.Double p1, double sphereRadius, Point3D.Double center,
			Point3D.Double basisU, Point3D.Double basisV, double alpha0,
			double deltaAlpha) {
		double nx = basisU.y * basisV.z - basisU.z * basisV.y;
		double ny = basisU.z * basisV.x - basisU.x * basisV.z;
		double nz = basisU.x * basisV.y - basisU.y * basisV.x;
		double normalLength = Math.sqrt(nx * nx + ny * ny + nz * nz);
		if (!Double.isFinite(normalLength) || normalLength < TOL) {
			throw new IllegalArgumentException("General-curve basis vectors are degenerate.");
		}
		nx /= normalLength;
		ny /= normalLength;
		nz /= normalLength;

		double[][] inverse = {
			{ basisU.x, basisV.x, nx },
			{ basisU.y, basisV.y, ny },
			{ basisU.z, basisV.z, nz }
		};
		double zPrime = center.x * nx + center.y * ny + center.z * nz;
		double thetaPrime = Math.acos(Math.max(-1.0,
				Math.min(1.0, zPrime / sphereRadius)));
		return new GeneralCurve(p0, p1, sphereRadius, inverse, thetaPrime,
				alpha0, deltaAlpha);
	}

	// -----------------------------------------------------------------------
	// BaseCurve implementation
	// -----------------------------------------------------------------------

	/**
	 * Returns θ(t) by evaluating the curve in the primed frame and rotating back
	 * (paper step 9).
	 *
	 * <p>
	 * In the primed frame: θ'(t) = θ*, φ'(t) = φ'₀ + t·Δφ'. The corresponding
	 * Cartesian primed point is rotated back by R⁻¹ and the spherical coordinates
	 * of the result give θ(t).
	 * </p>
	 */
	@Override
	public UnivariateFunction getThetaFunction() {
		return t -> {
			double[] xyz = primeToOriginal(t);
			double r2 = Math.sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
			return Math.acos(Math.max(-1.0, Math.min(1.0, xyz[2] / r2)));
		};
	}

	/**
	 * Returns φ(t) by evaluating the curve in the primed frame and rotating back
	 * (paper step 9).
	 */
	@Override
	public UnivariateFunction getPhiFunction() {
		return t -> {
			double[] xyz = primeToOriginal(t);
			return MathUtil.normalizeAngle(Math.atan2(xyz[1], xyz[0]));
		};
	}
	
	@Override
	public List<GeneralCurve> splitAtPhiCrossings() {
			    List<Crossing> crossings = getPhiCrossings();

	    if (crossings == null || crossings.isEmpty()) {
	        return List.of(this);
	    }

	    List<Double> splitTs = new ArrayList<>();

	    for (Crossing crossing : crossings) {
	        if (crossing == null) {
	            continue;
	        }

	        double t = crossing.t();

	        if (!Double.isFinite(t)) {
	            continue;
	        }

	        t = Math.max(0.0, Math.min(1.0, t));

	        boolean duplicate = false;
	        for (double existing : splitTs) {
	            if (Math.abs(existing - t) < 1.0e-10) {
	                duplicate = true;
	                break;
	            }
	        }

	        if (!duplicate) {
	            splitTs.add(t);
	        }
	    }

	    if (splitTs.isEmpty()) {
	        return List.of(this);
	    }

	    splitTs.sort(Double::compare);

	    List<GeneralCurve> pieces = new ArrayList<>();

	    double t0 = 0.0;
	    for (double t1 : splitTs) {
	        if (t1 - t0 > TOL) {
	            pieces.add(subCurve(t0, t1));
	        }
	        t0 = t1;
	    }

	    if (1.0 - t0 > TOL) {
	        pieces.add(subCurve(t0, 1.0));
	    }

	    if (pieces.isEmpty()) {
	        return List.of(this);
	    }

	    return pieces;
	}
	
	/**
	 * Splits this GENERAL curve at its interior theta-grid crossings.
	 *
	 * <p>
	 * The returned curves follow the same geometric path, in the same order and
	 * direction, as this curve. If this curve has no interior theta crossings, or
	 * if the only theta crossings occur at the endpoints, the returned list contains
	 * only this curve.
	 * </p>
	 *
	 * <p>
	 * Endpoint crossings are deliberately ignored here. They are junction events
	 * between adjacent boundary curves, not reasons to split this curve internally.
	 * </p>
	 *
	 * @return an ordered list of curve pieces, preserving the direction of this curve
	 */
	public List<GeneralCurve> splitAtThetaCrossings() {
	    List<Crossing> crossings = getThetaCrossings();

	    if (crossings == null || crossings.isEmpty()) {
	        return List.of(this);
	    }

	    /*
	     * Collect only true interior split parameters. getThetaCrossings() may include
	     * endpoint candidate crossings so the patch-level splicer can reason about
	     * junctions. Those are not internal split points for this curve.
	     */
	    List<Double> splitTs = new ArrayList<>();

	    for (Crossing crossing : crossings) {
	        if (crossing == null) {
	            continue;
	        }

	        double t = crossing.t();

	        if (!Double.isFinite(t)) {
	            continue;
	        }

	        if (isEndpointT(t)) {
	            continue;
	        }

	        t = Math.max(0.0, Math.min(1.0, t));

	        boolean duplicate = false;
	        for (double existing : splitTs) {
	            if (Math.abs(existing - t) < 1.0e-10) {
	                duplicate = true;
	                break;
	            }
	        }

	        if (!duplicate) {
	            splitTs.add(t);
	        }
	    }

	    if (splitTs.isEmpty()) {
	        return List.of(this);
	    }

	    /*
	     * Sorting by t preserves the original direction of travel. This is true
	     * regardless of the sign of deltaPrimePhi because t is the curve's own
	     * parameter, increasing from p0 to p1.
	     */
	    splitTs.sort(Double::compare);

	    List<GeneralCurve> pieces = new ArrayList<>();

	    double t0 = 0.0;
	    for (double t1 : splitTs) {
	        if (t1 - t0 > TOL) {
	            pieces.add(subCurve(t0, t1));
	        }
	        t0 = t1;
	    }

	    if (1.0 - t0 > TOL) {
	        pieces.add(subCurve(t0, 1.0));
	    }

	    if (pieces.isEmpty()) {
	        return List.of(this);
	    }

	    return pieces;
	}
	
	/**
	 * Creates a sub-curve corresponding to the parameter interval [t0, t1] of this
	 * curve.
	 *
	 * <p>
	 * The new curve uses the same primed-frame small-circle representation as this
	 * curve, but with a restricted primed azimuthal sweep. This preserves the
	 * original orientation: increasing local parameter on the returned curve
	 * corresponds to increasing {@code t} on this curve.
	 * </p>
	 *
	 * @param t0 starting parameter on this curve
	 * @param t1 ending parameter on this curve
	 * @return a GENERAL curve representing this curve from {@code t0} to {@code t1}
	 */
	private GeneralCurve subCurve(double t0, double t1) {
	    if (t0 < -TOL || t0 > 1.0 + TOL || t1 < -TOL || t1 > 1.0 + TOL) {
	        throw new IllegalArgumentException(String.format(
	                "Sub-curve parameters out of bounds: t0=%.12f, t1=%.12f", t0, t1));
	    }

	    if (t1 <= t0 + TOL) {
	        throw new IllegalArgumentException(String.format(
	                "Sub-curve parameters are not increasing: t0=%.12f, t1=%.12f", t0, t1));
	    }

	    t0 = Math.max(0.0, Math.min(1.0, t0));
	    t1 = Math.max(0.0, Math.min(1.0, t1));

	    Point3D.Double q0 = (t0 <= TOL) ? p0 : getPoint(t0);
	    Point3D.Double q1 = (t1 >= 1.0 - TOL) ? p1 : getPoint(t1);

	    double newPrimePhi0 = MathUtil.normalizeAngle(primePhi0 + t0 * deltaPrimePhi);
	    double newDeltaPrimePhi = deltaPrimePhi * (t1 - t0);

	    return new GeneralCurve(q0, q1, radius, invMatrix, thetaStar,
	            newPrimePhi0, newDeltaPrimePhi);
	}
	/**
	 * Overrides {@link BaseCurve#getPoint} for efficiency: rather than computing
	 * (θ, φ) and re-converting to Cartesian, we rotate directly from the primed
	 * frame.
	 */
	@Override
	public Point3D.Double getPoint(double t) {
		double[] xyz = primeToOriginal(t);
		return new Point3D.Double(xyz[0], xyz[1], xyz[2]);
	}


	@Override
	public GeneralCurve reverse() {
		double newPhi0 = MathUtil.normalizeAngle(primePhi0 + deltaPrimePhi);
		double newDeltaPhi = -deltaPrimePhi;
		return new GeneralCurve(p1, p0, radius, invMatrix, thetaStar, newPhi0, newDeltaPhi);
	}

	// -----------------------------------------------------------------------
	// Exact geometric export accessors
	// -----------------------------------------------------------------------

	/**
	 * Returns the center of the small circle containing this arc.
	 *
	 * <p>The returned point is expressed in the original Cartesian coordinate
	 * system. Together with {@link #getCircleRadius()}, {@link #getCircleBasisU()},
	 * {@link #getCircleBasisV()}, {@link #getCircleAlpha0()}, and
	 * {@link #getCircleDeltaAlpha()}, it defines the exact directed arc by
	 * {@code c + rho * (u*cos(alpha) + v*sin(alpha))}.</p>
	 *
	 * @return a new point containing the circle center
	 */
	public Point3D.Double getCircleCenter() {
		double zPrime = radius * Math.cos(thetaStar);
		double[] center = ChimeraRotation.multiplyMatrixVector(invMatrix,
				new double[] { 0.0, 0.0, zPrime });
		return new Point3D.Double(center[0], center[1], center[2]);
	}

	/** @return the radius of the small circle containing this arc */
	public double getCircleRadius() {
		return radius * Math.sin(thetaStar);
	}

	/**
	 * Returns the first unit basis vector in the circle plane.
	 * @return a new unit vector represented as a point triple
	 */
	public Point3D.Double getCircleBasisU() {
		return new Point3D.Double(invMatrix[0][0], invMatrix[1][0], invMatrix[2][0]);
	}

	/**
	 * Returns the second unit basis vector in the circle plane.
	 * @return a new unit vector represented as a point triple
	 */
	public Point3D.Double getCircleBasisV() {
		return new Point3D.Double(invMatrix[0][1], invMatrix[1][1], invMatrix[2][1]);
	}

	/** @return the starting circle angle in radians */
	public double getCircleAlpha0() {
		return primePhi0;
	}

	/** @return the signed angular sweep of the directed arc in radians */
	public double getCircleDeltaAlpha() {
		return deltaPrimePhi;
	}

	/**
	 * Returns the exact length of this small-circle arc.
	 *
	 * <p>In the primed frame the curve has constant circle radius
	 * {@code R sin(thetaStar)} and signed angular sweep
	 * {@code deltaPrimePhi}. Rotation back to the original frame preserves
	 * length, so numerical integration is unnecessary.</p>
	 *
	 * @return arc length in the same units as the sphere radius
	 */
	@Override
	public double arcLength() {
		return getCircleRadius() * Math.abs(deltaPrimePhi);
	}

	// -----------------------------------------------------------------------
	// Private helpers
	// -----------------------------------------------------------------------

	/**
	 * Evaluates the curve at {@code t} in the primed frame and rotates back to the
	 * original frame, returning a Cartesian xyz triple.
	 *
	 * <p>
	 * In the primed frame the point is the small-circle arc (R, θ*, φ'(t)), which
	 * in Cartesian is:
	 * </p>
	 * 
	 * <pre>
	 *   x' = R·sin(θ*)·cos(φ'(t))
	 *   y' = R·sin(θ*)·sin(φ'(t))
	 *   z' = R·cos(θ*)
	 * </pre>
	 */
	private double[] primeToOriginal(double t) {
		double phiPrime = MathUtil.normalizeAngle(primePhi0 + t * deltaPrimePhi);
		double sinTheta = Math.sin(thetaStar);
		double xp = radius * sinTheta * Math.cos(phiPrime);
		double yp = radius * sinTheta * Math.sin(phiPrime);
		double zp = radius * Math.cos(thetaStar);
		return ChimeraRotation.multiplyMatrixVector(invMatrix, new double[] { xp, yp, zp });
	}

	// -----------------------------------------------------------------------
	// Accessors
	// -----------------------------------------------------------------------

	/** @return θ* — the constant polar angle in the primed (rotated) frame */
	public double getThetaStar() {
		return thetaStar;
	}

	/** @return φ'₀ — start azimuthal angle in the primed frame */
	public double getPrimePhi0() {
		return primePhi0;
	}

	/** @return Δφ' — signed azimuthal sweep in the primed frame */
	public double getDeltaPrimePhi() {
		return deltaPrimePhi;
	}

	/**
	 * Computes candidate crossings where this curve meets theta grid lines.
	 *
	 * <p>This method deliberately distinguishes three cases:</p>
	 *
	 * <ul>
	 *   <li>A curve that lies along a theta grid line is not a crossing. It is a
	 *       theta-boundary segment, so this method returns no crossings.</li>
	 *   <li>An interior pass-through crossing is a genuine crossing and is kept.</li>
	 *   <li>An endpoint hit is kept as a candidate crossing. Whether it is a real
	 *       splice crossing depends on the neighboring curve, so final endpoint
	 *       filtering belongs in {@code ThetaPatch}, not here.</li>
	 * </ul>
	 *
	 * @return list of candidate crossings with theta grid lines
	 */
	public List<Crossing> getThetaCrossings() {
		if (thetaCrossings != null) {
			return thetaCrossings;
		}
	    thetaCrossings = new ArrayList<>();

	    SphericalGrid grid = ChimeraGridContext.sphericalGrid();
	    Grid1D thetaGrid = grid.getThetaGrid();

	    /*
	     * Constant-theta case.
	     *
	     * If this curve lies along a theta grid line, it is a boundary segment,
	     * not a crossing. Reporting its two endpoints as crossings creates
	     * artificial odd counts such as 1 or 3 at the prepatch level.
	     */
	    if (isConstantTheta()) {
	        return thetaCrossings;
	    }

	    /*
	     * Add analytic interior crossings. Endpoint solutions are handled
	     * explicitly below so that endpoint policy is clear and stable.
	     */
	    for (int i = 0; i < thetaGrid.numPoints(); i++) {
	        double targetTheta = thetaGrid.valueAt(i);
	        addThetaCrossingsForTarget(thetaCrossings, targetTheta, i);
	    }

	    /*
	     * Add endpoint candidate crossings.
	     *
	     * These are not automatically genuine crossings; the neighboring curve
	     * determines that. But the patch-level splicer needs to see them so it can
	     * pair cases like:
	     *
	     *     curve 22 ends on θ_k, curve 23 starts on θ_k,
	     *     and the boundary crosses the θ_k line at that junction.
	     */
	    addEndpointThetaCrossing(thetaCrossings, thetaGrid, 0.0);
	    addEndpointThetaCrossing(thetaCrossings, thetaGrid, 1.0);

	    thetaCrossings = Crossing.removeDuplicates(thetaCrossings);
	    return thetaCrossings;
	}	
	
	/**
	 * Adds crossings with one target theta grid line.
	 *
	 * <p>
	 * The calculation uses the fact that a GeneralCurve is a small circle in the
	 * primed frame. In the original frame, the z-coordinate is
	 * </p>
	 *
	 * <pre>
	 *   z(phi') = A cos(phi') + B sin(phi') + C
	 * </pre>
	 *
	 * <p>
	 * where phi' = primePhi0 + t * deltaPrimePhi. Since theta = acos(z/R), solving
	 * theta = targetTheta is equivalent to solving z = R cos(targetTheta).
	 * </p>
	 *
	 * @param crossings   list to receive any crossings found
	 * @param targetTheta target theta value
	 * @param thetaIndex  index of the theta grid line
	 */
	private void addThetaCrossingsForTarget(List<Crossing> crossings, double targetTheta, int thetaIndex) {
		if (Math.abs(deltaPrimePhi) < TOL) {
			return;
		}

		/*
		 * z = row 2 of invMatrix dot [x', y', z'].
		 *
		 * x' = R sin(thetaStar) cos(phi') y' = R sin(thetaStar) sin(phi') z' = R
		 * cos(thetaStar)
		 */
		double sinThetaStar = Math.sin(thetaStar);
		double cosThetaStar = Math.cos(thetaStar);

		double a = radius * sinThetaStar * invMatrix[2][0];
		double b = radius * sinThetaStar * invMatrix[2][1];
		double c = radius * cosThetaStar * invMatrix[2][2];

		double zTarget = radius * Math.cos(targetTheta);

		/*
		 * Solve:
		 *
		 * a cos(phi') + b sin(phi') + c = zTarget
		 *
		 * or:
		 *
		 * amp cos(phi' - alpha) = zTarget - c
		 */
		double amp = Math.hypot(a, b);

		if (amp < TOL * radius) {
			/*
			 * z is effectively constant. The constant-theta branch should already have
			 * handled the useful case, so there is nothing to add here.
			 */
			return;
		}

		double q = (zTarget - c) / amp;

		/*
		 * Allow a little numerical forgiveness near tangency.
		 */
		double eps = 1.0e-12;
		if (q > 1.0 + eps || q < -1.0 - eps) {
			return;
		}

		q = Math.max(-1.0, Math.min(1.0, q));

		double alpha = Math.atan2(b, a);
		double gamma = Math.acos(q);

		addPhiPrimeSolution(crossings, alpha + gamma, targetTheta, thetaIndex);
		addPhiPrimeSolution(crossings, alpha - gamma, targetTheta, thetaIndex);
	}

	/**
	 * Adds a crossing corresponding to one periodic phi' solution, if that solution
	 * lies in the interior of this curve's primed azimuthal sweep and represents a
	 * genuine pass-through crossing.
	 *
	 * @param crossings   list to receive the crossing
	 * @param phiSolution one representative solution for phi'
	 * @param targetTheta target theta value
	 * @param thetaIndex  index of the theta grid line
	 */
	private void addPhiPrimeSolution(List<Crossing> crossings,
	                                 double phiSolution,
	                                 double targetTheta,
	                                 int thetaIndex) {
	    double start = primePhi0;
	    double end = primePhi0 + deltaPrimePhi;

	    double lo = Math.min(start, end);
	    double hi = Math.max(start, end);

	    double twoPi = 2.0 * Math.PI;

	    int kMin = (int) Math.floor((lo - phiSolution) / twoPi) - 1;
	    int kMax = (int) Math.ceil((hi - phiSolution) / twoPi) + 1;

	    double phiTol = Math.max(1.0e-12, 1.0e-10 * Math.abs(deltaPrimePhi));

	    for (int k = kMin; k <= kMax; k++) {
	        double phi = phiSolution + k * twoPi;

	        if (phi < lo - phiTol || phi > hi + phiTol) {
	            continue;
	        }

	        double t = (phi - start) / deltaPrimePhi;

	        if (t < -1.0e-10 || t > 1.0 + 1.0e-10) {
	            continue;
	        }

	        t = Math.max(0.0, Math.min(1.0, t));

	        /*
	         * Endpoint hits are added explicitly by addEndpointThetaCrossing().
	         * Keeping them here as well only makes the endpoint policy harder to
	         * reason about.
	         */
	        if (isEndpointT(t)) {
	            continue;
	        }

	        /*
	         * Interior tangency: theta touches the cut and turns around.
	         * That is not a splice crossing.
	         */
	        if (!changesThetaSideAt(t, targetTheta)) {
	            continue;
	        }

	        crossings.add(new Crossing(this, t, targetTheta, thetaIndex));
	    }
	}
	
	/**
	 * Adds a candidate endpoint theta crossing if the endpoint lies on a theta
	 * grid line and the curve immediately leaves that line.
	 *
	 * <p>This method does not decide whether the prepatch boundary changes theta
	 * band at the junction. That requires the neighboring curve and must be done
	 * at the patch level.</p>
	 *
	 * @param crossings list to receive the crossing
	 * @param thetaGrid theta grid
	 * @param tEndpoint endpoint parameter, expected to be 0.0 or 1.0
	 */
	private void addEndpointThetaCrossing(List<Crossing> crossings,
	                                      Grid1D thetaGrid,
	                                      double tEndpoint) {
	    double th = theta(tEndpoint);
	    int index = thetaGrid.valueIsAVertex(th);

	    if (index < 0) {
	        return;
	    }

	    /*
	     * Suppress the endpoint if this curve merely lies along the theta line.
	     * The global constant-theta case is already handled earlier, but this
	     * local test also protects against very short nearly-on-cut segments.
	     */
	    double probeT = (tEndpoint <= 0.5) ? SIDE_TEST_DT : 1.0 - SIDE_TEST_DT;
	    probeT = Math.max(0.0, Math.min(1.0, probeT));

	    double probeTheta = theta(probeT);
	    if (Math.abs(probeTheta - th) < TOL) {
	        return;
	    }

	    crossings.add(new Crossing(this, tEndpoint, th, index));
	}

	/**
	 * Tests whether {@code t} is effectively an endpoint of the curve.
	 *
	 * @param t curve parameter
	 * @return true if {@code t} is near 0 or 1
	 */
	private static boolean isEndpointT(double t) {
	    return t <= ENDPOINT_CROSSING_TOL || t >= 1.0 - ENDPOINT_CROSSING_TOL;
	}

	/**
	 * Tests whether theta changes sides across an interior candidate crossing.
	 *
	 * <p>This filters interior tangencies. A tangency can produce a mathematically
	 * valid solution of theta(t) = thetaGridLine, but it does not move the boundary
	 * from one theta band into the next.</p>
	 *
	 * @param t           candidate crossing parameter
	 * @param targetTheta theta grid-line value
	 * @return true if theta - targetTheta changes sign across {@code t}
	 */
	private boolean changesThetaSideAt(double t, double targetTheta) {
	    double dt = Math.min(SIDE_TEST_DT,
	            0.25 * Math.min(t, 1.0 - t));

	    if (dt <= 0.0) {
	        return false;
	    }

	    double before = theta(t - dt) - targetTheta;
	    double after  = theta(t + dt) - targetTheta;

	    if (Math.abs(before) < TOL || Math.abs(after) < TOL) {
	        dt = Math.min(1.0e-4, 0.25 * Math.min(t, 1.0 - t));
	        if (dt <= 0.0) {
	            return false;
	        }

	        before = theta(t - dt) - targetTheta;
	        after  = theta(t + dt) - targetTheta;
	    }

	    return before * after < 0.0;
	}
	
	/**
	 * Computes candidate crossings where this GENERAL curve meets phi grid lines.
	 *
	 * <p>
	 * A phi grid line is a meridian: the half-plane through the z-axis with azimuth
	 * {@code targetPhi}. For a point {@code (x,y,z)}, lying in that meridian plane is
	 * equivalent to
	 * </p>
	 *
	 * <pre>
	 *   y cos(targetPhi) - x sin(targetPhi) = 0.
	 * </pre>
	 *
	 * <p>
	 * In the primed parameter of a {@code GeneralCurve}, both {@code x} and
	 * {@code y} are linear combinations of {@code cos(phiPrime)} and
	 * {@code sin(phiPrime)}, so the crossing equation can be solved analytically.
	 * Candidate endpoint crossings are included; patch-level splice code can later
	 * decide whether they are genuine band changes.
	 * </p>
	 *
	 * @return list of candidate crossings with phi grid lines
	 */
	@Override
	public List<Crossing> getPhiCrossings() {
	    if (phiCrossings != null) {
	        return phiCrossings;
	    }

	    phiCrossings = new ArrayList<>();

	    /*
	     * If the curve lies along a meridian, it is a phi-boundary segment, not a
	     * crossing. Reporting its endpoints as crossings would create artificial
	     * odd counts at the patch level.
	     */
	    if (isConstantPhi()) {
	        return phiCrossings;
	    }

	    SphericalGrid grid = ChimeraGridContext.sphericalGrid();
	    Grid1D phiGrid = grid.getPhiGrid();

	    for (int i = 0; i < phiGrid.numPoints(); i++) {
	        double targetPhi = MathUtil.normalizeAngle(phiGrid.valueAt(i));
	        addPhiCrossingsForTarget(phiCrossings, targetPhi, i);
	    }

	    addEndpointPhiCrossing(phiCrossings, phiGrid, 0.0);
	    addEndpointPhiCrossing(phiCrossings, phiGrid, 1.0);

	    phiCrossings = Crossing.removeDuplicates(phiCrossings);
	    return phiCrossings;
	}
	
	/**
	 * Returns true if this curve lies effectively along a constant-phi meridian.
	 */
	private boolean isConstantPhi() {
	    /*
	     * Avoid treating a pole-touching curve as constant phi. At a pole, phi is
	     * singular and the sampled value is not meaningful enough for this test.
	     */
	    if (Math.sin(theta(0.0)) < TOL
	            || Math.sin(theta(0.5)) < TOL
	            || Math.sin(theta(1.0)) < TOL) {
	        return false;
	    }

	    double ph0 = phi(0.0);
	    double phm = phi(0.5);
	    double ph1 = phi(1.0);

	    double d0m = Math.abs(MathUtil.normalizeAngle(phm - ph0));
	    double d01 = Math.abs(MathUtil.normalizeAngle(ph1 - ph0));

	    return d0m < TOL && d01 < TOL;
	}

	/**
	 * Adds crossings with one target phi grid line.
	 */
	private void addPhiCrossingsForTarget(List<Crossing> crossings,
	                                      double targetPhi,
	                                      int phiIndex) {

	    if (Math.abs(deltaPrimePhi) < TOL) {
	        return;
	    }

	    /*
	     * In the primed frame:
	     *
	     *   x' = R sin(thetaStar) cos(u)
	     *   y' = R sin(thetaStar) sin(u)
	     *   z' = R cos(thetaStar)
	     *
	     * where u = phiPrime.
	     *
	     * In the original frame:
	     *
	     *   x = a0 cos(u) + b0 sin(u) + c0
	     *   y = a1 cos(u) + b1 sin(u) + c1
	     *
	     * The meridian condition is:
	     *
	     *   y cos(phi) - x sin(phi) = 0.
	     */
	    double sinThetaStar = Math.sin(thetaStar);
	    double cosThetaStar = Math.cos(thetaStar);

	    double sx = radius * sinThetaStar;
	    double sz = radius * cosThetaStar;

	    double cosPhi = Math.cos(targetPhi);
	    double sinPhi = Math.sin(targetPhi);

	    double a = sx * (invMatrix[1][0] * cosPhi - invMatrix[0][0] * sinPhi);
	    double b = sx * (invMatrix[1][1] * cosPhi - invMatrix[0][1] * sinPhi);
	    double c = sz * (invMatrix[1][2] * cosPhi - invMatrix[0][2] * sinPhi);

	    double amp = Math.hypot(a, b);

	    if (amp < TOL * radius) {
	        return;
	    }

	    /*
	     * Solve:
	     *
	     *   a cos(u) + b sin(u) + c = 0
	     *
	     * as:
	     *
	     *   amp cos(u - alpha) = -c
	     */
	    double q = -c / amp;

	    double eps = 1.0e-12;
	    if (q > 1.0 + eps || q < -1.0 - eps) {
	        return;
	    }

	    q = Math.max(-1.0, Math.min(1.0, q));

	    double alpha = Math.atan2(b, a);
	    double gamma = Math.acos(q);

	    addPhiPrimeSolutionForPhiCrossing(crossings,
	            alpha + gamma, targetPhi, phiIndex);

	    addPhiPrimeSolutionForPhiCrossing(crossings,
	            alpha - gamma, targetPhi, phiIndex);
	}

	/**
	 * Adds one periodic phi-prime solution if it lies in the curve sweep and is on
	 * the selected meridian ray rather than the opposite meridian ray.
	 */
	private void addPhiPrimeSolutionForPhiCrossing(List<Crossing> crossings,
	                                               double phiPrimeSolution,
	                                               double targetPhi,
	                                               int phiIndex) {

	    double start = primePhi0;
	    double end = primePhi0 + deltaPrimePhi;

	    double lo = Math.min(start, end);
	    double hi = Math.max(start, end);

	    double twoPi = 2.0 * Math.PI;

	    int kMin = (int) Math.floor((lo - phiPrimeSolution) / twoPi) - 1;
	    int kMax = (int) Math.ceil((hi - phiPrimeSolution) / twoPi) + 1;

	    double phiTol = Math.max(1.0e-12, 1.0e-10 * Math.abs(deltaPrimePhi));

	    for (int k = kMin; k <= kMax; k++) {
	        double phiPrime = phiPrimeSolution + k * twoPi;

	        if (phiPrime < lo - phiTol || phiPrime > hi + phiTol) {
	            continue;
	        }

	        double t = (phiPrime - start) / deltaPrimePhi;

	        if (t < -1.0e-10 || t > 1.0 + 1.0e-10) {
	            continue;
	        }

	        t = Math.max(0.0, Math.min(1.0, t));

	        /*
	         * Endpoint hits are added explicitly below.
	         */
	        if (isEndpointT(t)) {
	            continue;
	        }

	        /*
	         * The meridian plane equation is satisfied both by phi and phi+pi.
	         * Keep only the point on the target meridian ray:
	         *
	         *   x cos(phi) + y sin(phi) >= 0.
	         */
	        if (!onTargetPhiRay(t, targetPhi)) {
	            continue;
	        }

	        if (!changesPhiSideAt(t, targetPhi)) {
	            continue;
	        }

	        crossings.add(new Crossing(this, t, targetPhi, phiIndex));
	    }
	}

	/**
	 * Adds a candidate endpoint phi crossing if the endpoint lies on a phi grid line
	 * and the curve immediately leaves that line.
	 */
	private void addEndpointPhiCrossing(List<Crossing> crossings,
	                                    Grid1D phiGrid,
	                                    double tEndpoint) {

	    /*
	     * Do not create endpoint phi crossings at a pole. Phi is singular there.
	     */
	    if (Math.sin(theta(tEndpoint)) < TOL) {
	        return;
	    }

	    double ph = MathUtil.normalizeAngle(phi(tEndpoint));

	    for (int i = 0; i < phiGrid.numPoints(); i++) {
	        double targetPhi = MathUtil.normalizeAngle(phiGrid.valueAt(i));

	        if (Math.abs(MathUtil.normalizeAngle(ph - targetPhi)) > TOL) {
	            continue;
	        }

	        double probeT = (tEndpoint <= 0.5) ? SIDE_TEST_DT : 1.0 - SIDE_TEST_DT;
	        probeT = Math.max(0.0, Math.min(1.0, probeT));

	        double probePhi = MathUtil.normalizeAngle(phi(probeT));

	        if (Math.abs(MathUtil.normalizeAngle(probePhi - targetPhi)) < TOL) {
	            continue;
	        }

	        crossings.add(new Crossing(this, tEndpoint, targetPhi, i));
	    }
	}

	/**
	 * Tests whether the point at {@code t} is on the target meridian ray rather than
	 * the opposite meridian ray.
	 */
	private boolean onTargetPhiRay(double t, double targetPhi) {
	    Point3D.Double p = getPoint(t);

	    double radialProjection = p.x * Math.cos(targetPhi)
	            + p.y * Math.sin(targetPhi);

	    return radialProjection >= -1.0e-10 * radius;
	}

	/**
	 * Tests whether phi changes sides across an interior candidate crossing.
	 *
	 * <p>
	 * This filters tangencies and numerical grazes. The side is measured by the
	 * signed normalized angular difference {@code phi(t) - targetPhi}.
	 * </p>
	 */
	private boolean changesPhiSideAt(double t, double targetPhi) {
	    double dt = Math.min(SIDE_TEST_DT,
	            0.25 * Math.min(t, 1.0 - t));

	    if (dt <= 0.0) {
	        return false;
	    }

	    double before = MathUtil.normalizeAngle(phi(t - dt) - targetPhi);
	    double after  = MathUtil.normalizeAngle(phi(t + dt) - targetPhi);

	    if (Math.abs(before) < TOL || Math.abs(after) < TOL) {
	        dt = Math.min(1.0e-4, 0.25 * Math.min(t, 1.0 - t));
	        if (dt <= 0.0) {
	            return false;
	        }

	        before = MathUtil.normalizeAngle(phi(t - dt) - targetPhi);
	        after  = MathUtil.normalizeAngle(phi(t + dt) - targetPhi);
	    }

	    /*
	     * Near the target meridian, a genuine crossing changes the sign of the
	     * wrapped difference. Guard against branch-cut artifacts by requiring both
	     * samples to be reasonably close to the target meridian.
	     */
	    if (Math.abs(before) > Math.PI / 2.0 || Math.abs(after) > Math.PI / 2.0) {
	        return false;
	    }

	    return before * after < 0.0;
	}

	@Override
	public String toString() {
		return String.format("GeneralCurve[thetaStar=%.4f primePhi0=%.4f deltaPhi=%.4f]", thetaStar, primePhi0,
				deltaPrimePhi);
	}
}
