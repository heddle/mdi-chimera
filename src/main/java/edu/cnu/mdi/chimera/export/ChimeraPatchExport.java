package edu.cnu.mdi.chimera.export;

import java.util.List;

import com.fasterxml.jackson.annotation.JsonSubTypes;
import com.fasterxml.jackson.annotation.JsonTypeInfo;

/**
 * Language-neutral, versioned representation of one completed Chimera run.
 *
 * <p>The records in this class are persistence DTOs, deliberately separate from
 * Chimera's computational object graph. Array-valued coordinates are always
 * ordered {@code [x, y, z]}; all angles are radians. Patch boundaries are
 * ordered directed loops.</p>
 */
public final class ChimeraPatchExport {

    private ChimeraPatchExport() {
    }

    /** Root object of the {@code mdi-chimera-patches} version 1 format. */
    public record Document(String format, int version, Conventions conventions,
            Grid grid, Summary summary, List<PatchData> patches) { }

    /** Coordinate and parameter conventions needed by independent consumers. */
    public record Conventions(String lengthUnit, String coordinateSystem,
            String angleUnit, String theta, String phi, double[] phiRange,
            double[] curveParameterRange, String boundaryOrientation) { }

    /** Complete input grid definition. */
    public record Grid(String name, Sphere sphere, CartesianGrid cartesian,
            SphericalGrid spherical) { }

    /** Sphere center and radius. */
    public record Sphere(double[] center, double radius) { }

    /** Cartesian vertex arrays and the translation applied to them. */
    public record CartesianGrid(double[] xVertices, double[] yVertices,
            double[] zVertices, double[] offset) { }

    /** Spherical-grid vertex arrays. */
    public record SphericalGrid(double[] thetaVertices, double[] phiVertices) { }

    /** Export-wide counts and numerical checks. */
    public record Summary(int intersectingCellCount, int kissCellCount,
            int patchCount, double normalizedAreaFraction,
            double physicalAreaEstimate, double totalPerimeter) { }

    /** One final patch and its ordered, directed boundary. */
    public record PatchData(int id, Indices indices, Diagnostics diagnostics,
            List<CurveData> boundary) { }

    /** Cartesian and spherical cell indices identifying a final patch. */
    public record Indices(int nx, int ny, int nz, int nTheta, int nPhi) { }

    /** Values useful for validating an independent implementation. */
    public record Diagnostics(double normalizedAreaFraction,
            double physicalAreaEstimate, double perimeter,
            boolean enclosesNorthPole, boolean enclosesSouthPole,
            double maximumClosureError) { }

    /**
     * A directed analytic boundary curve. The {@code type} discriminator is
     * written by Jackson and is part of the public interchange format.
     */
    @JsonTypeInfo(use = JsonTypeInfo.Id.NAME, property = "type")
    @JsonSubTypes({
        @JsonSubTypes.Type(value = ThetaCurveData.class, name = "THETA"),
        @JsonSubTypes.Type(value = PhiCurveData.class, name = "PHI"),
        @JsonSubTypes.Type(value = GeneralCurveData.class, name = "GENERAL")
    })
    public sealed interface CurveData permits ThetaCurveData, PhiCurveData,
            GeneralCurveData { }

    /** Directed constant-colatitude arc: phi(t) = phi0 + t*deltaPhi. */
    public record ThetaCurveData(double[] start, double[] end, double theta,
            double phi0, double deltaPhi) implements CurveData { }

    /** Directed constant-azimuth arc: theta(t) = theta0 + t*deltaTheta. */
    public record PhiCurveData(double[] start, double[] end, double phi,
            double theta0, double deltaTheta) implements CurveData { }

    /**
     * Directed small-circle arc in Cartesian coordinates.
     *
     * <p>For {@code alpha = alpha0 + t*deltaAlpha}, reconstruct the curve as
     * {@code center + radius*(basisU*cos(alpha) + basisV*sin(alpha))}.</p>
     */
    public record GeneralCurveData(double[] start, double[] end,
            double[] center, double radius, double[] basisU, double[] basisV,
            double alpha0, double deltaAlpha) implements CurveData { }
}
