package edu.cnu.mdi.chimera.export;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;

import edu.cnu.mdi.chimera.alg.ChimeraAlgorithmResult;
import edu.cnu.mdi.chimera.curve.BaseCurve;
import edu.cnu.mdi.chimera.curve.GeneralCurve;
import edu.cnu.mdi.chimera.curve.PhiCurve;
import edu.cnu.mdi.chimera.curve.ThetaCurve;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.CartesianGrid;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Conventions;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.CurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Diagnostics;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Document;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.GeneralCurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Grid;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Indices;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.PatchData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.PhiCurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Sphere;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.SphericalGrid;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.Summary;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.ThetaCurveData;
import edu.cnu.mdi.chimera.model.ChimeraGridSpec;
import edu.cnu.mdi.chimera.patch.Patch;
import edu.cnu.mdi.chimera.util.Point3D;

/**
 * Converts final patches to the versioned Chimera JSON interchange format.
 *
 * <p>The exporter writes analytic curve parameters rather than sampled
 * polylines. Consequently, another implementation can reconstruct every
 * directed boundary and independently compute area or perimeter without
 * inheriting Chimera's sampling choices.</p>
 */
public final class ChimeraPatchJsonExporter {

    /** Public format identifier written into every file. */
    public static final String FORMAT = "mdi-chimera-patches";

    /** Current interchange-format version. */
    public static final int VERSION = 1;

    private static final ObjectMapper MAPPER = new ObjectMapper()
            .enable(SerializationFeature.INDENT_OUTPUT);

    private ChimeraPatchJsonExporter() {
    }

    /**
     * Builds an immutable export document from an input grid and completed run.
     *
     * @param spec exact grid specification supplied to the algorithm
     * @param result completed algorithm result containing final patches
     * @return export document ready for JSON serialization
     * @throws IllegalArgumentException if either argument is null or the result
     *         does not contain final patches
     */
    public static Document createDocument(ChimeraGridSpec spec,
            ChimeraAlgorithmResult result) {
        if (spec == null || result == null) {
            throw new IllegalArgumentException("Grid specification and result are required.");
        }
        if (result.getPatches() == null) {
            throw new IllegalArgumentException("The algorithm result has no final patches.");
        }

        var cartesian = spec.getCartesianGrid();
        var spherical = spec.getSphericalGrid();
        Grid grid = new Grid(spec.getName(),
                new Sphere(new double[] { 0.0, 0.0, 0.0 }, spherical.getRadius()),
                new CartesianGrid(
                        cartesian.getXGrid().getPoints(),
                        cartesian.getYGrid().getPoints(),
                        cartesian.getZGrid().getPoints(),
                        new double[] { cartesian.getXOffset(), cartesian.getYOffset(),
                                cartesian.getZOffset() }),
                new SphericalGrid(spherical.getThetaGrid().getPoints(),
                        spherical.getPhiGrid().getPoints()));

        List<PatchData> patches = new ArrayList<>(result.getPatchCount());
        double areaFraction = 0.0;
        double perimeter = 0.0;
        for (int id = 0; id < result.getPatches().size(); id++) {
            Patch patch = result.getPatches().get(id);
            PatchData data = exportPatch(id, patch);
            patches.add(data);
            areaFraction += data.diagnostics().normalizedAreaFraction();
            perimeter += data.diagnostics().perimeter();
        }

        double sphereArea = 4.0 * Math.PI * spherical.getRadius()
                * spherical.getRadius();
        Summary summary = new Summary(result.getIntersectingCellCount(),
                result.getKissCellCount(), patches.size(), areaFraction,
                areaFraction * sphereArea, perimeter);

        Conventions conventions = new Conventions(
                spec.getLengthUnit().name(), spec.getCoordinateSystem().name(),
                "RADIANS", "COLATITUDE_FROM_POSITIVE_Z",
                "AZIMUTH_FROM_POSITIVE_X_TOWARD_POSITIVE_Y",
                new double[] { -Math.PI, Math.PI }, new double[] { 0.0, 1.0 },
                "ORDERED_DIRECTED_CLOSED_LOOP");

        return new Document(FORMAT, VERSION, conventions, grid, summary,
                List.copyOf(patches));
    }

    /**
     * Writes a pretty-printed UTF-8 JSON interchange file.
     *
     * @param output destination file; its parent directory must already exist
     * @param spec exact input grid specification
     * @param result completed algorithm result
     * @throws IOException if the output cannot be written
     */
    public static void write(Path output, ChimeraGridSpec spec,
            ChimeraAlgorithmResult result) throws IOException {
        if (output == null) {
            throw new IllegalArgumentException("Output path is required.");
        }

		/*
		 * Build and validate the complete document before opening the destination.
		 * If geometry conversion fails, an existing file remains untouched and a new
		 * destination is not left behind as a misleading zero-byte file.
		 */
		Document document = createDocument(spec, result);
        try (var writer = Files.newBufferedWriter(output)) {
            MAPPER.writeValue(writer, document);
        }
    }

    private static PatchData exportPatch(int id, Patch patch) {
        double areaFraction = patch.areaEstimate();
        double sphereArea = 4.0 * Math.PI * patch.radius * patch.radius;
        double perimeter = patch.perimeter();
        Diagnostics diagnostics = new Diagnostics(areaFraction,
                areaFraction * sphereArea, perimeter, patch.enclosesNorthPole(),
                patch.enclosesSouthPole(), maximumClosureError(patch));

        List<CurveData> boundary = patch.curves.stream()
                .map(ChimeraPatchJsonExporter::exportCurve)
                .toList();
        return new PatchData(id,
                new Indices(patch.nx, patch.ny, patch.nz, patch.nTheta, patch.nPhi),
                diagnostics, boundary);
    }

    private static CurveData exportCurve(BaseCurve curve) {
        double[] start = xyz(curve.p0);
        double[] end = xyz(curve.p1);
        if (curve instanceof ThetaCurve theta) {
            return new ThetaCurveData(start, end, theta.getThetaStar(),
                    theta.getPhi0(), theta.getDeltaPhi());
        }
        if (curve instanceof PhiCurve phi) {
            return new PhiCurveData(start, end, phi.getPhiStar(), phi.sv0.theta,
                    phi.getDeltaTheta());
        }
        if (curve instanceof GeneralCurve general) {
            return new GeneralCurveData(start, end,
                    xyz(general.getCircleCenter()), general.getCircleRadius(),
                    xyz(general.getCircleBasisU()), xyz(general.getCircleBasisV()),
                    general.getCircleAlpha0(), general.getCircleDeltaAlpha());
        }
        throw new IllegalArgumentException(
                "Unsupported curve type: " + curve.getClass().getName());
    }

    private static double maximumClosureError(Patch patch) {
        double maximum = 0.0;
        for (int i = 0; i < patch.curves.size(); i++) {
            Point3D.Double end = patch.curves.get(i).p1;
            Point3D.Double nextStart = patch.curves
                    .get((i + 1) % patch.curves.size()).p0;
            maximum = Math.max(maximum, Point3D.Double.distance(end, nextStart));
        }
        return maximum;
    }

    private static double[] xyz(Point3D.Double point) {
        return new double[] { point.x, point.y, point.z };
    }
}
