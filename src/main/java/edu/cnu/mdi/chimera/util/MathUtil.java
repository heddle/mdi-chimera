package edu.cnu.mdi.chimera.util;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.NoBracketingException;

public class MathUtil {

	// Private constructor to prevent instantiation
	private MathUtil() {
	}

	/**
	 * Normalizes an azimuthal angle to be within the range [-π, π].
	 *
	 * @param angle the angle in radians to normalize
	 * @return the normalized angle in the range [-π, π]
	 */
	public static double normalizeAngle(double angle) {

		while (angle > Math.PI) {
			angle -= 2 * Math.PI;
		}
		while (angle <= -Math.PI) {
			angle += 2 * Math.PI;
		}

		return angle;
	}

	/**
	 * Normalizes an angle difference to be within the range [0, 2π].
	 *
	 * @param angle1 the first angle in radians
	 * @param angle2 the second angle in radians
	 * @return the normalized angle difference in the range [0, 2π]
	 */
	public static double normalizedAngleDifference(double angle1, double angle2) {
		angle1 = normalizeAngle(angle1);
		angle2 = normalizeAngle(angle2);
		double diff = normalizeAngle(angle1 - angle2);

		if (diff < 0) {
			diff += 2 * Math.PI;
		}
		return diff;
	}

	/**
	 * Finds the parameter value {@code t} in {@code [tmin, tmax]} such that
	 * {@code func(t) == value}, using Brent's root-finding algorithm.
	 *
	 * <p>This is a numerical function inversion: given a target output {@code value},
	 * the method returns the input {@code t} that produces it. It is particularly
	 * useful for arc-length parameterisation, where the cumulative arc-length
	 * function is monotonic but has no closed-form inverse.</p>
	 *
	 * <h2>Algorithm</h2>
	 * <p>The search range is divided into {@code numIntervals} equal subintervals.
	 * Each subinterval is tested for a sign change in {@code func(t) - value}; the
	 * first subinterval where a sign change is detected is handed to
	 * {@link BrentSolver}, which refines the root to within {@code tolerance}.
	 * The scan proceeds left-to-right, so the <em>smallest</em> qualifying {@code t}
	 * is returned when multiple solutions exist.</p>
	 *
	 * <h2>Choosing {@code numIntervals}</h2>
	 * <ul>
	 *   <li>For smooth, monotonic functions a small value (even {@code 1}) is
	 *       sufficient and fastest.</li>
	 *   <li>For functions that cross {@code value} multiple times, or that have
	 *       rapid local variation, increase {@code numIntervals} so that no crossing
	 *       falls entirely inside a single subinterval and is missed.</li>
	 *   <li>Each interval requires two function evaluations during the scan, so
	 *       runtime scales linearly with {@code numIntervals}.</li>
	 * </ul>
	 *
	 * @param func         the function to invert; must be defined and finite on
	 *                     {@code [tmin, tmax]}
	 * @param value        the target function value to match
	 * @param tmin         the lower bound of the search range (inclusive)
	 * @param tmax         the upper bound of the search range (inclusive);
	 *                     must be strictly greater than {@code tmin}
	 * @param tolerance    the convergence tolerance passed to {@link BrentSolver};
	 *                     the returned {@code t} satisfies
	 *                     {@code |func(t) - value| < tolerance}
	 * @param numIntervals the number of equal subintervals used during the
	 *                     bracketing scan; must be &ge; 1
	 * @return the first {@code t} in {@code [tmin, tmax]} such that
	 *         {@code func(t) ≈ value}, or {@link Double#NaN} if no bracketed root
	 *         was found in any subinterval
	 * @throws IllegalArgumentException if {@code tmax <= tmin} or
	 *                                  {@code numIntervals < 1}
	 */
	public static double computeT(UnivariateFunction func, double value,
	        double tmin, double tmax, double tolerance, int numIntervals) {

	    if (tmax <= tmin) {
	        throw new IllegalArgumentException(
	                "tmax must be greater than tmin, got [" + tmin + ", " + tmax + "]");
	    }
	    if (numIntervals < 1) {
	        throw new IllegalArgumentException(
	                "numIntervals must be >= 1, got " + numIntervals);
	    }

	    final int MAX_SOLVER_ITERATIONS = 1000;

	    BrentSolver solver = new BrentSolver(tolerance);
	    UnivariateFunction rootFunc = t -> func.value(t) - value;

	    double step = (tmax - tmin) / numIntervals;

	    for (int i = 0; i < numIntervals; i++) {
	        double t1 = tmin + i * step;
	        double t2 = tmin + (i + 1) * step;
	        double f1 = rootFunc.value(t1);
	        double f2 = rootFunc.value(t2);

	        // Short-circuit: grid point is exactly the root.
	        if (f1 == 0.0) return t1;
	        if (f2 == 0.0) return t2;

	        // Safe sign-change test: avoids overflow/underflow from f1*f2.
	        if (Math.signum(f1) != Math.signum(f2)) {
	            try {
	                return solver.solve(MAX_SOLVER_ITERATIONS, rootFunc, t1, t2);
	            } catch (NoBracketingException e) {
	                // Numerically, the sign change vanished by the time the solver
	                // evaluated the endpoints (e.g. a tangential root). Continue
	                // scanning; a better bracket may exist in a later subinterval.
	                continue;
	            }
	        }
	    }

	    return Double.NaN;
	}
	/**
	 * Finds the point where a function transitions from satisfying a test to not
	 * @param phi the function
	 * @param test the test function
	 * @param tol the tolerance
	 * @return the point where the function transitions
	 */
	public static double findTogglePoint(UnivariateFunction phi, java.util.function.Predicate<Double> test,
			double tol) {
		double tLow = 0.0;
		double tHigh = 1.0;
		boolean initialTest = test.test(phi.value(tLow));

// Binary search for the transition point
		while (tHigh - tLow > tol) {
			double tMid = (tLow + tHigh) / 2.0;
			boolean midTest = test.test(phi.value(tMid));

			if (midTest == initialTest) {
				tLow = tMid; // Move right
			} else {
				tHigh = tMid; // Move left
			}
		}

		return (tLow + tHigh) / 2.0; // Best estimate of the toggle point
	}

	// Test the normalizeAngle method
	public static void main(String[] args) {
		double[] testAngles = { 0, Math.PI, -Math.PI, 2 * Math.PI, -2 * Math.PI, 3 * Math.PI, -3 * Math.PI, 10, -10 };

		for (double angle : testAngles) {
			System.out.printf("Original: %f, Normalized: %f%n", angle, normalizeAngle(angle));
		}

		for (int i = 0; i < testAngles.length - 1; i++) {
			System.out.printf("Angle difference between %f and %f: %f%n", testAngles[i], testAngles[i + 1],
					normalizedAngleDifference(testAngles[i], testAngles[i + 1]));
		}
	}
}
