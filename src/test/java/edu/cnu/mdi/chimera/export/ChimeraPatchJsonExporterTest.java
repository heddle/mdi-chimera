package edu.cnu.mdi.chimera.export;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTimeout;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.util.ArrayList;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import edu.cnu.mdi.chimera.alg.ChimeraAlgorithm;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.CurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.GeneralCurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.PhiCurveData;
import edu.cnu.mdi.chimera.export.ChimeraPatchExport.ThetaCurveData;
import edu.cnu.mdi.chimera.curve.GeneralCurve;
import edu.cnu.mdi.chimera.model.ChimeraGridPresets;

/** End-to-end checks for the public patch interchange format. */
class ChimeraPatchJsonExporterTest {

    private static final double GEOMETRY_TOLERANCE = 1.0e-12;

    @TempDir
    Path temporaryDirectory;

    /** Verifies metadata, discriminator fields, and physical-area units. */
    @Test
    void writesSelfDescribingJson() throws Exception {
        var spec = ChimeraGridPresets.smallDebugGrid();
        var result = ChimeraAlgorithm.run(spec);
        Path output = temporaryDirectory.resolve("patches.json");

        ChimeraPatchJsonExporter.write(output, spec, result);

        JsonNode root = new ObjectMapper().readTree(output.toFile());
        assertEquals("mdi-chimera-patches", root.path("format").asText());
        assertEquals(1, root.path("version").asInt());
        assertEquals(result.getPatchCount(), root.path("summary").path("patchCount").asInt());
        assertEquals(result.getPatchCount(), root.path("patches").size());
        assertTrue(Files.size(output) > 0);

        JsonNode firstCurve = root.path("patches").get(0).path("boundary").get(0);
        assertTrue(firstCurve.hasNonNull("type"));
        assertEquals(3, firstCurve.path("start").size());
        assertEquals(3, firstCurve.path("end").size());
    }

    /**
     * Reconstructs every exported analytic curve at t=0 and t=1 and compares it
     * with the independently exported endpoints.
     */
    @Test
    void analyticCurveParametersReconstructEndpoints() {
        var spec = ChimeraGridPresets.smallDebugGrid();
        var document = ChimeraPatchJsonExporter.createDocument(spec,
                ChimeraAlgorithm.run(spec));

        for (var patch : document.patches()) {
            for (CurveData curve : patch.boundary()) {
                assertVectorClose(start(curve), exportedStart(curve));
                assertVectorClose(end(curve), exportedEnd(curve));
            }
            assertTrue(patch.diagnostics().maximumClosureError() <= 1.0e-6);
        }
    }

	/** Ensures general-curve perimeter uses its exact small-circle formula. */
	@Test
    void generalCurveArcLengthMatchesExportedCircleParameters() {
		var spec = ChimeraGridPresets.smallDebugGrid();
		var result = ChimeraAlgorithm.run(spec);

		for (var patch : result.getPatches()) {
			for (var curve : patch.curves) {
				if (curve instanceof GeneralCurve general) {
					double expected = general.getCircleRadius()
							* Math.abs(general.getCircleDeltaAlpha());
					assertEquals(expected, general.arcLength(), GEOMETRY_TOLERANCE);
				}
			}
		}
	}

	/** Guards against accidentally reintroducing per-curve numerical integration. */
	@Test
	void paperGridExportCompletesPromptly() {
		assertTimeout(Duration.ofSeconds(30), () -> {
			var spec = ChimeraGridPresets.paperTestGrid();
			Path output = temporaryDirectory.resolve("paper-grid-patches.json");
			var original = ChimeraAlgorithm.run(spec);
			ChimeraPatchJsonExporter.write(output, spec, original);
			assertTrue(Files.size(output) > 1_000_000,
					"paper-grid export should contain complete patch geometry");

			var imported = ChimeraPatchJsonReader.read(output);
			assertEquals(original.getPatchCount(),
					imported.algorithmResult().getPatchCount());
			assertEquals(spec.getSphericalGrid().getRadius(),
					imported.gridSpec().getSphericalGrid().getRadius());
			for (int patchIndex = 0; patchIndex < original.getPatchCount(); patchIndex++) {
				var expectedPatch = original.getPatches().get(patchIndex);
				var actualPatch = imported.algorithmResult().getPatches().get(patchIndex);
				assertEquals(expectedPatch.fullIndex(), actualPatch.fullIndex());
				assertEquals(expectedPatch.curves.size(), actualPatch.curves.size());
				for (int curveIndex = 0; curveIndex < expectedPatch.curves.size(); curveIndex++) {
					var expectedMidpoint = expectedPatch.curves.get(curveIndex).getPoint(0.5);
					var actualMidpoint = actualPatch.curves.get(curveIndex).getPoint(0.5);
					assertEquals(expectedMidpoint.x, actualMidpoint.x, GEOMETRY_TOLERANCE);
					assertEquals(expectedMidpoint.y, actualMidpoint.y, GEOMETRY_TOLERANCE);
					assertEquals(expectedMidpoint.z, actualMidpoint.z, GEOMETRY_TOLERANCE);
				}
			}
		});
	}

	/** Ensures incompatible files are rejected before reconstruction. */
	@Test
	void rejectsUnsupportedFormatVersion() throws Exception {
		var spec = ChimeraGridPresets.smallDebugGrid();
		Path output = temporaryDirectory.resolve("unsupported.json");
		ChimeraPatchJsonExporter.write(output, spec, ChimeraAlgorithm.run(spec));
		String json = Files.readString(output).replaceFirst("\"version\" : 1",
				"\"version\" : 999");
		Files.writeString(output, json);

		assertThrows(java.io.IOException.class,
				() -> ChimeraPatchJsonReader.read(output));
	}

	/** Imported results report final-patch area without intermediate cell data. */
	@Test
	void importedResultFeedbackIncludesFinalPatchArea() throws Exception {
		var spec = ChimeraGridPresets.smallDebugGrid();
		Path output = temporaryDirectory.resolve("feedback-round-trip.json");
		ChimeraPatchJsonExporter.write(output, spec, ChimeraAlgorithm.run(spec));
		var imported = ChimeraPatchJsonReader.read(output).algorithmResult();

		var feedback = new ArrayList<String>();
		imported.feedbackSummary("", feedback);

		assertTrue(feedback.stream().anyMatch(line ->
				line.matches("final patches: 672 normalized area: [0-9.]+")),
				() -> "missing final-patch area feedback: " + feedback);
	}

    private static double[] start(CurveData curve) {
        return evaluate(curve, 0.0);
    }

    private static double[] end(CurveData curve) {
        return evaluate(curve, 1.0);
    }

    private static double[] evaluate(CurveData curve, double t) {
        if (curve instanceof ThetaCurveData c) {
            double phi = c.phi0() + t * c.deltaPhi();
            double radius = norm(c.start());
            return spherical(radius, c.theta(), phi);
        }
        if (curve instanceof PhiCurveData c) {
            double theta = c.theta0() + t * c.deltaTheta();
            double radius = norm(c.start());
            return spherical(radius, theta, c.phi());
        }
        GeneralCurveData c = (GeneralCurveData) curve;
        double alpha = c.alpha0() + t * c.deltaAlpha();
        double[] point = new double[3];
        for (int i = 0; i < 3; i++) {
            point[i] = c.center()[i] + c.radius()
                    * (c.basisU()[i] * Math.cos(alpha)
                    + c.basisV()[i] * Math.sin(alpha));
        }
        return point;
    }

    private static double[] exportedStart(CurveData curve) {
        if (curve instanceof ThetaCurveData c) return c.start();
        if (curve instanceof PhiCurveData c) return c.start();
        return ((GeneralCurveData) curve).start();
    }

    private static double[] exportedEnd(CurveData curve) {
        if (curve instanceof ThetaCurveData c) return c.end();
        if (curve instanceof PhiCurveData c) return c.end();
        return ((GeneralCurveData) curve).end();
    }

    private static double[] spherical(double radius, double theta, double phi) {
        double sinTheta = Math.sin(theta);
        return new double[] {
            radius * sinTheta * Math.cos(phi),
            radius * sinTheta * Math.sin(phi),
            radius * Math.cos(theta)
        };
    }

    private static double norm(double[] point) {
        return Math.sqrt(point[0] * point[0] + point[1] * point[1]
                + point[2] * point[2]);
    }

    private static void assertVectorClose(double[] expected, double[] actual) {
        for (int i = 0; i < 3; i++) {
            assertEquals(expected[i], actual[i], GEOMETRY_TOLERANCE);
        }
    }
}
