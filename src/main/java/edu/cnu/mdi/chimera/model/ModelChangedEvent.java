package edu.cnu.mdi.chimera.model;

import java.time.Instant;
import java.util.Objects;

/**
 * Describes a change to the shared model.
 * <p>
 * A model changed event has a {@link Type}, the source {@link ChimeraModel},
 * optional old and new values, and a timestamp. Most views will only need the
 * event type and the model reference.
 * </p>
 */
public final class ModelChangedEvent {

    /**
     * Type of model change.
     */
    public enum Type {
    	
    	/**
    	 * Algorithm options changed.
    	 */
    	ALGORITHM_OPTIONS_CHANGED,
    	
    	/**
    	 * The complete algorithm result changed.
    	 */
    	ALGORITHM_RESULT_CHANGED,

        /**
         * The active grid specification changed.
         */
        GRID_SPEC_CHANGED,

        /**
         * The model was reset.
         */
        MODEL_RESET,

        /**
         * A algorithm run started.
         */
        ALGORITHM_STARTED,

        /**
         * A algorithm run completed.
         */
        ALGORITHM_COMPLETED,

        /**
         * A algorithm run failed.
         */
        ALGORITHM_FAILED,

        /**
         * The list of Cartesian cells intersecting the sphere changed.
         */
        INTERSECTING_CELLS_CHANGED,

        /**
         * The list of prepatches changed.
         */
        PREPATCHES_CHANGED,

        /**
         * The list of theta-spliced patches changed.
         */
        THETA_PATCHES_CHANGED,

        /**
         * The list of final patches changed.
         */
        PATCHES_CHANGED,

        /**
         * The current selected object changed.
         */
        SELECTION_CHANGED,

        /**
         * Display options changed, without changing the underlying computed data.
         */
        DISPLAY_OPTIONS_CHANGED,

        /**
         * A status message changed.
         */
        STATUS_CHANGED,
        
        /**
         * Monte Carlo point data changed.
         */
        MONTE_CARLO_CHANGED,

        /**
         * Monte Carlo point data was cleared.
         */
        MONTE_CARLO_CLEARED
    }

    /** The model that changed. */
    private final ChimeraModel model;

    /** The type of change. */
    private final Type type;

    /** Optional old value. */
    private final Object oldValue;

    /** Optional new value. */
    private final Object newValue;

    /** Optional human-readable message. */
    private final String message;

    /** Time at which the event was created. */
    private final Instant timestamp;

    /**
     * Creates a model changed event.
     *
     * @param model the model that changed
     * @param type the event type
     * @param oldValue optional old value
     * @param newValue optional new value
     * @param message optional status or diagnostic message
     */
    public ModelChangedEvent(ChimeraModel model, Type type,
            Object oldValue, Object newValue, String message) {

        this.model = Objects.requireNonNull(model, "model must not be null");
        this.type = Objects.requireNonNull(type, "type must not be null");
        this.oldValue = oldValue;
        this.newValue = newValue;
        this.message = message;
        this.timestamp = Instant.now();
    }

    /**
     * Gets the model that changed.
     *
     * @return the model
     */
    public ChimeraModel getModel() {
        return model;
    }

    /**
     * Gets the event type.
     *
     * @return the event type
     */
    public Type getType() {
        return type;
    }

    /**
     * Gets the optional old value.
     *
     * @return the old value, or {@code null}
     */
    public Object getOldValue() {
        return oldValue;
    }

    /**
     * Gets the optional new value.
     *
     * @return the new value, or {@code null}
     */
    public Object getNewValue() {
        return newValue;
    }

    /**
     * Gets the optional message.
     *
     * @return the message, or {@code null}
     */
    public String getMessage() {
        return message;
    }

    /**
     * Gets the event timestamp.
     *
     * @return the timestamp
     */
    public Instant getTimestamp() {
        return timestamp;
    }

    /**
     * Gets the old value cast to the requested type.
     *
     * @param <T> the requested type
     * @param clazz the requested class
     * @return the old value cast to {@code clazz}, or {@code null}
     */
    public <T> T getOldValue(Class<T> clazz) {
        if (clazz == null || oldValue == null) {
            return null;
        }
        return clazz.isInstance(oldValue) ? clazz.cast(oldValue) : null;
    }

    /**
     * Gets the new value cast to the requested type.
     *
     * @param <T> the requested type
     * @param clazz the requested class
     * @return the new value cast to {@code clazz}, or {@code null}
     */
    public <T> T getNewValue(Class<T> clazz) {
        if (clazz == null || newValue == null) {
            return null;
        }
        return clazz.isInstance(newValue) ? clazz.cast(newValue) : null;
    }

    @Override
    public String toString() {
        return "ModelChangedEvent[type=" + type
                + ", message=" + message
                + ", timestamp=" + timestamp + "]";
    }
}