package edu.cnu.mdi.chimera.curve;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.cnu.mdi.chimera.util.SphericalVector;

/**
 * Represents a crossing of a curve with a grid line, storing the parameter value
 * at which the crossing occurs and the index of the grid line.
 *
 * @param curve the curve that crosses the grid line
 * @param t     the parameter value along the curve where the crossing occurs
 * @param value the value of the curve's coordinate at the crossing,
 *              e.g., theta or phi in radians
 * @param index the index of the grid line that is crossed
 */
public record Crossing(BaseCurve curve, double t, double value, int index) {

    private static final double TOL = 1.0e-9;

    /**
     * Returns a string representation of the crossing, including the parameter value,
     * the coordinate value in degrees, and the grid line index.
     *
     * @return a formatted string summarizing the crossing
     */
    public String summaryString() {
        return String.format("%s, t: %.3f, value: %.3f, index: %d",
                curve.shortString(), t, Math.toDegrees(value), index);
    }

    /**
     * Computes the great-circle distance between this crossing and another crossing.
     * The distance is calculated based on the spherical coordinates theta and phi
     * of the two crossings.
     *
     * @param other the other crossing to compare against
     * @return the great-circle distance in radians between the two crossings
     * @throws IllegalArgumentException if {@code other} is {@code null}
     */
    public double distanceTo(Crossing other) {
        if (other == null) {
            throw new IllegalArgumentException("Cannot compute distance to a null crossing.");
        }

        BaseCurve c1 = this.curve;
        BaseCurve c2 = other.curve;

        double theta1 = c1.getThetaFunction().value(this.t);
        double phi1 = c1.getPhiFunction().value(this.t);
        double theta2 = c2.getThetaFunction().value(other.t);
        double phi2 = c2.getPhiFunction().value(other.t);

        return SphericalVector.gcDistance(theta1, phi1, theta2, phi2);
    }

    /**
     * Determines whether this crossing is effectively the same as another crossing,
     * based on whether the great-circle distance between them is less than a
     * specified tolerance.
     *
     * @param other the other crossing to compare against
     * @return {@code true} if the crossings are considered the same,
     *         {@code false} otherwise
     */
    public boolean duplicate(Crossing other) {
        return (other != null) && (distanceTo(other) < TOL);
    }

    /**
     * Determines whether two crossings are effectively duplicates.
     *
     * <p>This static convenience method is useful when duplicate testing is needed
     * from collection utilities or from code that may not have a preferred
     * receiver object.</p>
     *
     * @param c1 the first crossing
     * @param c2 the second crossing
     * @return {@code true} if both crossings are non-null and are considered
     *         duplicates
     */
    public static boolean areDuplicates(Crossing c1, Crossing c2) {
        return (c1 != null) && c1.duplicate(c2);
    }

    /**
     * Returns a new list with duplicate crossings removed.
     *
     * <p>The input list is not modified. The first occurrence of each distinct
     * crossing is retained, and later crossings that duplicate an earlier retained
     * crossing are omitted.</p>
     *
     * <p>This method is safe to call on a list that may be concurrently modified
     * only if the caller synchronizes on the same lock used by all writers, or if
     * the caller passes a stable snapshot. The method itself makes an immediate
     * snapshot of the supplied list before processing, so later modifications to
     * the original list do not affect the returned result.</p>
     *
     * @param crossings the crossings to deduplicate; may be {@code null}
     * @return an unmodifiable list containing the first occurrence of each
     *         distinct crossing; never {@code null}
     */
    public static List<Crossing> removeDuplicates(List<Crossing> crossings) {
        if (crossings == null || crossings.isEmpty()) {
            return Collections.emptyList();
        }

        /*
         * Snapshot first. This protects the deduplication pass from later changes
         * to the caller's list. Note, however, that this copy constructor still
         * requires the caller not to mutate a plain ArrayList concurrently during
         * the copy itself unless external synchronization is used.
         */
        List<Crossing> snapshot = new ArrayList<>(crossings);
        List<Crossing> unique = new ArrayList<>();

        for (Crossing crossing : snapshot) {
            if (crossing == null) {
                continue;
            }

            boolean duplicate = false;

            for (Crossing existing : unique) {
                if (Crossing.areDuplicates(crossing, existing)) {
                    duplicate = true;
                    break;
                }
            }

            if (!duplicate) {
                unique.add(crossing);
            }
        }

        return Collections.unmodifiableList(unique);
    }

    /**
     * Returns a new list with duplicate crossings removed, using caller-provided
     * synchronization while taking the snapshot.
     *
     * <p>Use this overload when {@code crossings} is a mutable list that may be
     * accessed by multiple threads. Every other read or write of that list should
     * use the same lock object.</p>
     *
     * @param crossings the crossings to deduplicate; may be {@code null}
     * @param lock      the lock protecting access to {@code crossings}; must not
     *                  be {@code null} when {@code crossings} is non-null
     * @return an unmodifiable list containing the first occurrence of each
     *         distinct crossing; never {@code null}
     */
    public static List<Crossing> removeDuplicates(List<Crossing> crossings, Object lock) {
        if (crossings == null || crossings.isEmpty()) {
            return Collections.emptyList();
        }
        if (lock == null) {
            throw new IllegalArgumentException("Lock must not be null.");
        }

        final List<Crossing> snapshot;
        synchronized (lock) {
            snapshot = new ArrayList<>(crossings);
        }

        return removeDuplicates(snapshot);
    }
}