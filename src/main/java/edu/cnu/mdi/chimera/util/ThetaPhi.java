package edu.cnu.mdi.chimera.util;

import java.awt.geom.Point2D;
import java.util.Random;

public class ThetaPhi extends Point2D.Double {

	/** The small theta character */
	public static final String SMALL_THETA = "\u03B8";

	/** The small phi character */
	public static final String SMALL_PHI = "\u03C6";

	/** The degree character */
	public static final String DEGREE = "\u00B0";

	private double radius;

	/**
	 * Constructor for ThetaPhi.
	 *
	 * @param theta The polar angle in the range [0, π].
	 * @param phi   The azimuthal angle in the range [-π, π].
	 */
	public ThetaPhi(double radius, double theta, double phi) {
		this.radius = radius;
		setTheta(theta);
		setPhi(phi);
	}

	/**
	 * Get the radius
	 *
	 * @return the radius
	 */
	public double getRadius() {
		return radius;
	}

	/**
	 * Gets the polar angle theta in radians.
	 *
	 * @return The polar angle theta.
	 */
	public double getTheta() {
		return x; // Theta is stored in x
	}

	/**
	 * Gets the polar angle theta in degrees.
	 *
	 * @return The polar angle theta.
	 */
	public double getThetaDegrees() {
		return Math.toDegrees(getTheta());
	}

	/**
	 * Sets the polar angle theta, ensuring it is in the range [0, π].
	 *
	 * @param theta The polar angle to set.
	 * @throws IllegalArgumentException if theta is out of range [0, π].
	 */
	public void setTheta(double theta) {
		if (theta < 0 || theta > Math.PI) {
			throw new IllegalArgumentException("Theta must be in the range [0, π]. value = " + theta + " radians.");
		}
		x = theta; // Use Point2D's x to store theta
	}

	/**
	 * Gets the azimuthal angle phi in radians.
	 *
	 * @return The azimuthal angle phi.
	 */
	public double getPhi() {
		return y; // Phi is stored in y
	}

	/**
	 * Gets the azimuthal angle phi in degrees.
	 *
	 * @return The azimuthal angle phi.
	 */
	public double getPhiDegrees() {
		return Math.toDegrees(getPhi());
	}

	/**
	 * Sets the azimuthal angle phi, normalizing it to the range [-π, π].
	 *
	 * @param phi The azimuthal angle to set.
	 */
	public void setPhi(double phi) {
		y = MathUtil.normalizeAngle(phi); // Use MathUtil to normalize phi
	}

	public void azimuthShift(double deltaPhi) {
		setPhi(getPhi() + deltaPhi);
	}

	/**
	 * Converts the spherical coordinates to Cartesian coordinates.
	 *
	 * @param cartesian The Point3D.Double to store the Cartesian
	 */
	public void toCartesian(Point3D.Double cartesian) {
		double sinTheta = Math.sin(getTheta());
		cartesian.x = radius * sinTheta * Math.cos(getPhi());
		cartesian.y = radius * sinTheta * Math.sin(getPhi());
		cartesian.z = radius * Math.cos(getTheta());
	}

	/**
	 * Converts the spherical coordinates to Cartesian coordinates.
	 *
	 * @return The Cartesian coordinates as a Point
	 */
	public Point3D.Double toCartesian() {
		Point3D.Double cartesian = new Point3D.Double();
		toCartesian(cartesian);
		return cartesian;
	}

	/**
	 * Converts polar angle theta (0 to 180) to latitude (-90 to 90)
	 *
	 * @return the latitude in radians
	 */
	public double getLatitude() {
		return Math.PI / 2 - getTheta();
	}

	/**
	 * Static method to generate random values for theta and phi.
	 *
	 * @param random   The Random object to use for generating random numbers.
	 * @param thetaPhi The ThetaPhi object to set random values.
	 */
	public static void setRandomThetaPhi(Random random, ThetaPhi thetaPhi) {
		// Generate random theta using theta = arccos(2v - 1)
		double v = random.nextDouble(); // Uniform random [0, 1)
		double theta = Math.acos(2 * v - 1);

		// Generate random phi in the range [-π, π]
		double phi = MathUtil.normalizeAngle(random.nextDouble() * 2 * Math.PI - Math.PI);

		thetaPhi.setTheta(theta);
		thetaPhi.setPhi(phi);
	}

	@Override
	public String toString() {
		return String.format("[%s = %.4f%s, %s = %.4f%s] ", SMALL_THETA, getThetaDegrees(), DEGREE, SMALL_PHI,
				getPhiDegrees(), DEGREE);
	}

}
