package edu.cnu.mdi.chimera.curve;

import org.apache.commons.math3.analysis.UnivariateFunction;

import edu.cnu.mdi.chimera.util.MathUtil;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * Represents a parameter value t such that func(t) is approximately equal to a target value. This is used for
 * arc-length parameterisation, where we want to find the t value that corresponds to a given arc length along the curve.
 */
public class TValue {
	
	private static int MAX_INTERVALS = 100; // Maximum number of intervals for root finding

	// The curve that this TValue is associated with. This is needed to compute the point on the curve at the computed t value.
	private final BaseCurve curve;
	
	// The parameter value t such that func(t) is approximately equal to the target value. This is computed in the constructor and stored as a final field.
	public final double t;
	
	// The target value that we want func(t) to be approximately equal to. This is stored as a final field for reference.
	public final double value;
	
	// The function that we are trying to invert. This is stored as a final field for reference and for computing the interpolated value at t.
	private final UnivariateFunction _func;


	public TValue(BaseCurve curve, UnivariateFunction func, double targetValue, double tmin, double tmax,
			double tolerance) {
		this.curve = curve;
		this.value = targetValue;
		this._func = func;
		t = MathUtil.computeT(func, value, tmin, tmax, tolerance, MAX_INTERVALS); // Assign computed value to `this.t`
	}

	/**
	 * Get the point on the curve at the computed t value.
	 *
	 * @return the point on the curve at the computed t value.
	 */
	public Point3D.Double point() {
		return curve.getPoint(t);
	}


	@Override
	public String toString() {
		String out = String.format("t = %.3f, targVal = %.3f  interpVal = %.3f",
				t, value, _func.value(t));
		return out;
	}

}
