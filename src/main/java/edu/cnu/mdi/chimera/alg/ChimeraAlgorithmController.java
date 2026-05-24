package edu.cnu.mdi.chimera.alg;

import edu.cnu.mdi.chimera.model.ChimeraModel;
import edu.cnu.mdi.log.Log;

/**
 * Controller that runs the algorithm and stores the result in the shared
 * model.
 * <p>
 * This class is the bridge between UI actions, such as a menu item, and the
 * stateless algorithm implementation.
 * </p>
 */
public final class ChimeraAlgorithmController {

    /** Shared model. */
    private final ChimeraModel model;

    /**
     * Creates an algorithm controller.
     *
     * @param model shared model
     */
    public ChimeraAlgorithmController(ChimeraModel model) {
        if (model == null) {
            throw new IllegalArgumentException("model must not be null.");
        }

        this.model = model;
    }

    /**
     * Runs the currently implemented algorithm synchronously.
     * <p>
     * This is fine for the current step-1 implementation. When later stages
     * become heavier, this method can be replaced or complemented by a
     * SwingWorker-based asynchronous run method without changing the algorithm
     * result/model structure.
     * </p>
     *
     * @return {@code true} if the run succeeded
     */
    public boolean runAlgorithm() {
        try {
            model.algorithmStarted();

            ChimeraAlgorithmResult result =
                    ChimeraAlgorithm.run(model.getGridSpec());
            model.setAlgorithmResult(result);

            model.algorithmCompleted("Algorithm completed.");

            return true;
        } catch (RuntimeException ex) {
            Log.getInstance().error("Algorithm failed: " + ex.getMessage());
            ex.printStackTrace();
            model.algorithmFailed(ex);
            return false;
        }
    }

    /**
     * Clears the current algorithm result from the model.
     */
    public void clearAlgorithmResult() {
        model.clearAlgorithmResult();
    }
}