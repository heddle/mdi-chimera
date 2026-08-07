package edu.cnu.mdi.chimera.model;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTimeout;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.time.Duration;
import java.util.Arrays;
import java.util.stream.Stream;

import org.junit.jupiter.api.DynamicTest;
import org.junit.jupiter.api.TestFactory;

import edu.cnu.mdi.chimera.alg.ChimeraAlgorithm;
import edu.cnu.mdi.chimera.model.ChimeraGridPresets.GalleryPreset;
import edu.cnu.mdi.chimera.patch.BasePatch;

/** End-to-end checks for every grid exposed by the test-grid gallery menu. */
class ChimeraGridGalleryTest {

    /**
     * Runs every gallery grid through the complete algorithm and verifies that it
     * produces closed, fully indexed final patches within an interactive timeout.
     *
     * @return one independently reported dynamic test per gallery entry
     */
    @TestFactory
    Stream<DynamicTest> everyGalleryGridCompletes() {
        return Arrays.stream(GalleryPreset.values()).map(preset ->
                DynamicTest.dynamicTest(preset.toString(), () ->
                        assertTimeout(Duration.ofSeconds(30), () -> {
                            var result = ChimeraAlgorithm.run(preset.createGridSpec());
                            assertFalse(result.getPatches().isEmpty());
                            result.getPatches().forEach(patch -> {
                                assertTrue(patch.isFullyIndexed(), patch::fullIndex);
                                assertTrue(BasePatch.validateLoop(patch.curves),
                                        patch::fullIndex);
                            });
                        })));
    }
}
