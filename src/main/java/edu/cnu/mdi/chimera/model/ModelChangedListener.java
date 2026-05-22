package edu.cnu.mdi.chimera.model;

/**
 * Listener notified when the shared model changes.
 * <p>
 * This is intentionally small and domain-specific. MDI remains the application
 * framework; this listener mechanism only lets views and controllers
 * react to changes in the data model.
 * </p>
 */
@FunctionalInterface
public interface ModelChangedListener {

    /**
     * Called when the model changes.
     *
     * @param event the model change event
     */
    void modelChanged(ModelChangedEvent event);
}