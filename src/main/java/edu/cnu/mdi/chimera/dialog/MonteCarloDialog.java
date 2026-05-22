package edu.cnu.mdi.chimera.dialog;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.Window;
import java.util.List;
import java.util.Random;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;

import edu.cnu.mdi.chimera.mc.MonteCarloGenerator;
import edu.cnu.mdi.chimera.mc.MonteCarloPoint;
import edu.cnu.mdi.chimera.model.ChimeraModel;


/**
 * Dialog for generating Monte Carlo sample points on the current
 * spherical shell.
 */
@SuppressWarnings("serial")
public class MonteCarloDialog extends JDialog {

    /** Default number of points to add. */
    private static final int DEFAULT_COUNT = 3_000_000;

    /** Shared model. */
    private final ChimeraModel model;

    /** Point count spinner. */
    private JSpinner countSpinner;

    /** Optional seed field. */
    private JTextField seedField;

    /** Append checkbox. */
    private JCheckBox appendCheckBox;

    /** Progress bar. */
    private JProgressBar progressBar;

    /** Generate button. */
    private JButton generateButton;

    /** Clear button. */
    private JButton clearButton;

    /** Close/cancel button. */
    private JButton closeButton;

    /** Current worker, if generation is running. */
    private SwingWorker<List<MonteCarloPoint>, Void> worker;

    /**
     * Creates the Monte Carlo dialog.
     *
     * @param owner owner window
     * @param model shared model
     */
    public MonteCarloDialog(Window owner, ChimeraModel model) {
        super(owner, "Generate Monte Carlo Points", ModalityType.MODELESS);

        if (model == null) {
            throw new IllegalArgumentException("model must not be null.");
        }

        this.model = model;

        setLayout(new BorderLayout(8, 8));
        add(createInputPanel(), BorderLayout.CENTER);
        add(createButtonPanel(), BorderLayout.SOUTH);

        pack();
        setMinimumSize(getPreferredSize());
        setLocationRelativeTo(owner);
    }

    /**
     * Shows the Monte Carlo dialog.
     *
     * @param parent parent component
     * @param model shared model
     */
    public static void showDialog(Component parent, ChimeraModel model) {
        Window owner = (parent == null) ? null : SwingUtilities.getWindowAncestor(parent);
        MonteCarloDialog dialog = new MonteCarloDialog(owner, model);
        dialog.setVisible(true);
    }

    /**
     * Creates the input panel.
     *
     * @return input panel
     */
    private Component createInputPanel() {
        JPanel main = new JPanel(new BorderLayout(8, 8));
        main.setBorder(BorderFactory.createEmptyBorder(10, 10, 4, 10));

        JPanel fields = new JPanel(new GridLayout(0, 2, 8, 6));

        countSpinner = new JSpinner(new SpinnerNumberModel(DEFAULT_COUNT, 1, 50_000_000, 10_000));
        seedField = new JTextField();
        appendCheckBox = new JCheckBox("Append to existing points", true);

        fields.add(new JLabel("Points to generate"));
        fields.add(countSpinner);

        fields.add(new JLabel("Random seed, blank = random"));
        fields.add(seedField);

        fields.add(new JLabel("Current point count"));
        fields.add(new JLabel(String.format("%,d", model.getMonteCarloPointCount())));

        fields.add(new JLabel("Mode"));
        fields.add(appendCheckBox);

        progressBar = new JProgressBar(0, 100);
        progressBar.setStringPainted(true);
        progressBar.setValue(0);

        main.add(fields, BorderLayout.CENTER);
        main.add(progressBar, BorderLayout.SOUTH);

        return main;
    }

    /**
     * Creates the button panel.
     *
     * @return button panel
     */
    private Component createButtonPanel() {
        JPanel panel = new JPanel();

        generateButton = new JButton("Generate");
        clearButton = new JButton("Clear Existing");
        closeButton = new JButton("Close");

        generateButton.addActionListener(e -> startGeneration());
        clearButton.addActionListener(e -> clearExisting());
        closeButton.addActionListener(e -> closeOrCancel());

        panel.add(generateButton);
        panel.add(clearButton);
        panel.add(closeButton);

        return panel;
    }

    /**
     * Starts Monte Carlo generation on a worker thread.
     */
    private void startGeneration() {
        if (worker != null && !worker.isDone()) {
            return;
        }

        int count = ((Number) countSpinner.getValue()).intValue();
        boolean append = appendCheckBox.isSelected();
        Random random = createRandom();

        setRunning(true);
        progressBar.setValue(0);

        worker = new SwingWorker<>() {
            @Override
            protected List<MonteCarloPoint> doInBackground() {
                return MonteCarloGenerator.generate(
                        model.getGridSpec(),
                        count,
                        random,
                        this::setProgress,
                        this::isCancelled);
            }

            @Override
            protected void done() {
                try {
                    if (!isCancelled()) {
                        List<MonteCarloPoint> points = get();

                        if (append) {
                            model.addMonteCarloPoints(points);
                        } else {
                            model.setMonteCarloPoints(points);
                        }
                    }
                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(
                            MonteCarloDialog.this,
                            "Monte Carlo generation failed:\n" + ex.getMessage(),
                            "Monte Carlo Error",
                            JOptionPane.ERROR_MESSAGE);
                } finally {
                    setRunning(false);
                    progressBar.setValue(100);
                }
            }
        };

        worker.addPropertyChangeListener(evt -> {
            if ("progress".equals(evt.getPropertyName())) {
                progressBar.setValue((Integer) evt.getNewValue());
            }
        });

        worker.execute();
    }

    /**
     * Creates the random number generator from the seed field.
     *
     * @return random number generator
     */
    private Random createRandom() {
        String seedText = seedField.getText();

        if (seedText == null || seedText.isBlank()) {
            return new Random();
        }

        try {
            long seed = Long.parseLong(seedText.trim());
            return new Random(seed);
        } catch (NumberFormatException e) {
            /*
             * Allow a text seed by hashing it. This makes names like "debug1"
             * reproducible while keeping the UI simple.
             */
            return new Random(seedText.trim().hashCode());
        }
    }

    /**
     * Clears existing Monte Carlo points after confirmation.
     */
    private void clearExisting() {
        int count = model.getMonteCarloPointCount();
        if (count == 0) {
            return;
        }

        int choice = JOptionPane.showConfirmDialog(
                this,
                "Clear " + String.format("%,d", count) + " Monte Carlo points?",
                "Clear Monte Carlo Points",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.WARNING_MESSAGE);

        if (choice == JOptionPane.OK_OPTION) {
            model.clearMonteCarloPoints();
        }
    }

    /**
     * Closes the dialog or cancels the running worker.
     */
    private void closeOrCancel() {
        if (worker != null && !worker.isDone()) {
            worker.cancel(true);
            return;
        }

        dispose();
    }

    /**
     * Updates UI state while generation is running.
     *
     * @param running true while running
     */
    private void setRunning(boolean running) {
        generateButton.setEnabled(!running);
        clearButton.setEnabled(!running);
        countSpinner.setEnabled(!running);
        seedField.setEnabled(!running);
        appendCheckBox.setEnabled(!running);
        closeButton.setText(running ? "Cancel" : "Close");
    }
}