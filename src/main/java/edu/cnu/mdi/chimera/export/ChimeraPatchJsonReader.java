package edu.cnu.mdi.chimera.export;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.ObjectMapper;

import edu.cnu.mdi.chimera.alg.ChimeraAlgorithmResult;
import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.GeneralCurve;
import edu.cnu.mdi.chimera.curve.PhiCurve;
import edu.cnu.mdi.chimera.curve.ThetaCurve;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.CurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Document;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.GeneralCurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.PatchData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.PhiCurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.ThetaCurveData;
import edu.cnu.mdi.chimera.grid.CartesianGrid;
import edu.cnu.mdi.chimera.grid.SphericalGrid;
import edu.cnu.mdi.chimera.model.ChimeraGridContext;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.model.CoordinateSystem;
import edu.cnu.mdi.chimera.model.LengthUnit;
import edu.cnu.mdi.chimera.patch.BasePatch;
import edu.cnu.mdi.chimera.patch.Patch;
import edu.cnu.mdi.chimera.util.Point3D;

/** Reads, validates, and reconstructs version 1 final-patch JSON files. */
public final class ChimeraPatchJsonReader {

    private static final double GEOMETRY_TOL = 1.0e-7;
    private static final double UNIT_TOL = 1.0e-8;

    private static final ObjectMapper MAPPER = new ObjectMapper()
            .enable(DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES);

    private ChimeraPatchJsonReader() {
    }

    /** A validated grid specification paired with reconstructed final patches. */
    public record ImportedRun(ChimeraGridSpec gridSpec,
            ChimeraAlgorithmResult algorithmResult) { }

    /**
     * Reads and structurally validates a JSON file without modifying application
     * state or constructing live patch objects.
     *
     * @param input JSON file to read
     * @return validated interchange document
     * @throws IOException if parsing or validation fails
     */
    public static Document validate(Path input) throws IOException {
        if (input == null) {
            throw new IllegalArgumentException("Input path is required.");
        }
        Document document = MAPPER.readValue(input.toFile(), Document.class);
        validateDocument(document);
        return document;
    }

    /**
     * Reads, validates, and reconstructs the exported grid and final patches.
     * No application model is changed by this method.
     *
     * @param input JSON file to import
     * @return reconstructed grid and algorithm result
     * @throws IOException if parsing, validation, or reconstruction fails
     */
    public static ImportedRun read(Path input) throws IOException {
        Document document = validate(input);
        try {
            ChimeraGridSpec spec = createGridSpec(document);
            ChimeraGridContext.activate(spec);
            List<Patch> patches = new ArrayList<>(document.patches().size());
            for (PatchData data : document.patches()) {
                List<BaseCurve> boundary = data.boundary().stream()
                        .map(curve -> reconstructCurve(curve, spec.getSphericalGrid().getRadius()))
                        .toList();
                var indices = data.indices();
                patches.add(new Patch(boundary, indices.nx(), indices.ny(),
                        indices.nz(), indices.nTheta(), indices.nPhi()));
            }
            ChimeraAlgorithmResult result = ChimeraAlgorithmResult.empty();
            result.setPatches(patches);
            return new ImportedRun(spec, result);
        } catch (RuntimeException exception) {
            throw invalid("Could not reconstruct patch geometry: "
                    + exception.getMessage(), exception);
        }
    }

    private static void validateDocument(Document document) throws IOException {
        require(document != null, "Document is null.");
        require(ChimeraPatchJsonExporter.FORMAT.equals(document.format()),
                "Unsupported format: " + document.format());
        require(document.version() == ChimeraPatchJsonExporter.VERSION,
                "Unsupported format version: " + document.version());
        require(document.conventions() != null, "Missing conventions.");
        require("RADIANS".equals(document.conventions().angleUnit()),
                "Only radian angles are supported.");
        require(document.grid() != null, "Missing grid.");
        require(document.summary() != null, "Missing summary.");
        require(document.patches() != null, "Missing patches.");
        require(document.summary().patchCount() == document.patches().size(),
                "Summary patchCount does not match patches array.");

        ChimeraGridSpec spec;
        try {
            spec = createGridSpec(document);
        } catch (RuntimeException exception) {
            throw invalid("Invalid grid: " + exception.getMessage(), exception);
        }

        var cartesian = spec.getCartesianGrid();
        var spherical = spec.getSphericalGrid();
        Set<Integer> ids = new HashSet<>();
        Set<String> fiveTuples = new HashSet<>();
        for (PatchData patch : document.patches()) {
            require(patch != null, "Patch entry is null.");
            require(ids.add(patch.id()), "Duplicate patch id " + patch.id());
            require(patch.indices() != null, "Patch " + patch.id() + " has no indices.");
            var index = patch.indices();
            require(inRange(index.nx(), cartesian.getNumXCells())
                    && inRange(index.ny(), cartesian.getNumYCells())
                    && inRange(index.nz(), cartesian.getNumZCells())
                    && inRange(index.nTheta(), spherical.getNumThetaCells())
                    && inRange(index.nPhi(), spherical.getNumPhiCells()),
                    "Patch " + patch.id() + " has an out-of-range index.");
            String key = index.nx() + ":" + index.ny() + ":" + index.nz()
                    + ":" + index.nTheta() + ":" + index.nPhi();
            require(fiveTuples.add(key), "Duplicate patch five-tuple " + key);
            require(patch.boundary() != null && patch.boundary().size() >= 2,
                    "Patch " + patch.id() + " has an incomplete boundary.");

            for (int i = 0; i < patch.boundary().size(); i++) {
                CurveData curve = patch.boundary().get(i);
                validateCurve(curve, spherical.getRadius(), patch.id(), i);
                double[] end = end(curve);
                double[] nextStart = start(patch.boundary()
                        .get((i + 1) % patch.boundary().size()));
                require(distance(end, nextStart) <= GEOMETRY_TOL,
                        "Patch " + patch.id() + " boundary is not closed at curve " + i + ".");
            }
        }
    }

    private static ChimeraGridSpec createGridSpec(Document document) {
        var grid = document.grid();
        var c = grid.cartesian();
        var s = grid.spherical();
        if (c == null || s == null || grid.sphere() == null) {
            throw new IllegalArgumentException("Incomplete grid definition.");
        }
        double[] offset = requireVector(c.offset(), "Cartesian offset");
        double[] center = requireVector(grid.sphere().center(), "Sphere center");
        if (norm(center) > GEOMETRY_TOL) {
            throw new IllegalArgumentException("Only origin-centered spheres are supported.");
        }
        CartesianGrid cartesian = new CartesianGrid(c.xVertices(), c.yVertices(),
                c.zVertices(), offset[0], offset[1], offset[2]);
        SphericalGrid spherical = new SphericalGrid(s.thetaVertices(),
                s.phiVertices(), grid.sphere().radius());
        return new ChimeraGridSpec(grid.name(),
                CoordinateSystem.valueOf(document.conventions().coordinateSystem()),
                LengthUnit.valueOf(document.conventions().lengthUnit()),
                cartesian, spherical);
    }

    private static void validateCurve(CurveData curve, double sphereRadius,
            int patchId, int curveIndex) throws IOException {
        require(curve != null, "Patch " + patchId + " curve " + curveIndex + " is null.");
        double[] start = start(curve);
        double[] end = end(curve);
        requireVector(start, "curve start");
        requireVector(end, "curve end");
        require(Math.abs(norm(start) - sphereRadius) <= GEOMETRY_TOL,
                "Curve start is not on the sphere.");
        require(Math.abs(norm(end) - sphereRadius) <= GEOMETRY_TOL,
                "Curve end is not on the sphere.");

        double[] evaluatedStart = evaluate(curve, 0.0, sphereRadius);
        double[] evaluatedEnd = evaluate(curve, 1.0, sphereRadius);
        require(distance(start, evaluatedStart) <= GEOMETRY_TOL,
                "Curve parameters do not reproduce its start point.");
        require(distance(end, evaluatedEnd) <= GEOMETRY_TOL,
                "Curve parameters do not reproduce its end point.");

        if (curve instanceof GeneralCurveData general) {
            require(Math.abs(norm(requireVector(general.basisU(), "basisU")) - 1.0) <= UNIT_TOL,
                    "General-curve basisU is not a unit vector.");
            require(Math.abs(norm(requireVector(general.basisV(), "basisV")) - 1.0) <= UNIT_TOL,
                    "General-curve basisV is not a unit vector.");
            require(Math.abs(dot(general.basisU(), general.basisV())) <= UNIT_TOL,
                    "General-curve basis vectors are not orthogonal.");
            requireFinite(general.radius(), "general-curve radius");
            require(general.radius() >= 0.0, "General-curve radius is negative.");
        }
    }

    private static BaseCurve reconstructCurve(CurveData curve, double radius) {
        Point3D.Double p0 = point(start(curve));
        Point3D.Double p1 = point(end(curve));
        if (curve instanceof ThetaCurveData theta) {
            return ThetaCurve.parameterized(p0, p1, radius, theta.theta(),
                    theta.phi0(), theta.deltaPhi());
        }
        if (curve instanceof PhiCurveData phi) {
            return PhiCurve.onMeridian(p0, p1, radius, phi.phi());
        }
        GeneralCurveData general = (GeneralCurveData) curve;
        return GeneralCurve.fromCircleParameters(p0, p1, radius,
                point(general.center()), point(general.basisU()),
                point(general.basisV()), general.alpha0(), general.deltaAlpha());
    }

    private static double[] evaluate(CurveData curve, double t, double radius)
            throws IOException {
        if (curve instanceof ThetaCurveData theta) {
            requireFinite(theta.theta(), "theta");
            requireFinite(theta.phi0(), "phi0");
            requireFinite(theta.deltaPhi(), "deltaPhi");
            return spherical(radius, theta.theta(), theta.phi0() + t * theta.deltaPhi());
        }
        if (curve instanceof PhiCurveData phi) {
            requireFinite(phi.phi(), "phi");
            requireFinite(phi.theta0(), "theta0");
            requireFinite(phi.deltaTheta(), "deltaTheta");
            return spherical(radius, phi.theta0() + t * phi.deltaTheta(), phi.phi());
        }
        GeneralCurveData general = (GeneralCurveData) curve;
        double[] center = requireVector(general.center(), "center");
        double[] u = requireVector(general.basisU(), "basisU");
        double[] v = requireVector(general.basisV(), "basisV");
        requireFinite(general.alpha0(), "alpha0");
        requireFinite(general.deltaAlpha(), "deltaAlpha");
        double alpha = general.alpha0() + t * general.deltaAlpha();
        return new double[] {
            center[0] + general.radius() * (u[0] * Math.cos(alpha) + v[0] * Math.sin(alpha)),
            center[1] + general.radius() * (u[1] * Math.cos(alpha) + v[1] * Math.sin(alpha)),
            center[2] + general.radius() * (u[2] * Math.cos(alpha) + v[2] * Math.sin(alpha))
        };
    }

    private static double[] start(CurveData curve) {
        if (curve instanceof ThetaCurveData c) return c.start();
        if (curve instanceof PhiCurveData c) return c.start();
        if (curve instanceof GeneralCurveData c) return c.start();
        throw new IllegalArgumentException("Unsupported curve type.");
    }

    private static double[] end(CurveData curve) {
        if (curve instanceof ThetaCurveData c) return c.end();
        if (curve instanceof PhiCurveData c) return c.end();
        if (curve instanceof GeneralCurveData c) return c.end();
        throw new IllegalArgumentException("Unsupported curve type.");
    }

    private static double[] spherical(double r, double theta, double phi) {
        double sinTheta = Math.sin(theta);
        return new double[] { r * sinTheta * Math.cos(phi),
                r * sinTheta * Math.sin(phi), r * Math.cos(theta) };
    }

    private static Point3D.Double point(double[] v) {
        return new Point3D.Double(v[0], v[1], v[2]);
    }

    private static double[] requireVector(double[] value, String label) {
        if (value == null || value.length != 3) {
            throw new IllegalArgumentException(label + " must contain exactly three numbers.");
        }
        for (double coordinate : value) {
            if (!Double.isFinite(coordinate)) {
                throw new IllegalArgumentException(label + " contains a non-finite number.");
            }
        }
        return value;
    }

    private static boolean inRange(int value, int upperExclusive) {
        return value >= 0 && value < upperExclusive;
    }

    private static double norm(double[] v) {
        return Math.sqrt(dot(v, v));
    }

    private static double dot(double[] a, double[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    private static double distance(double[] a, double[] b) {
        double dx = a[0] - b[0];
        double dy = a[1] - b[1];
        double dz = a[2] - b[2];
        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }

    private static void requireFinite(double value, String label) throws IOException {
        require(Double.isFinite(value), label + " is not finite.");
    }

    private static void require(boolean condition, String message) throws IOException {
        if (!condition) throw invalid(message, null);
    }

    private static IOException invalid(String message, Throwable cause) {
        return new IOException("Invalid Chimera patch JSON: " + message, cause);
    }
}
