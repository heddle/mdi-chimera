package edu.cnu.mdi.chimera.map;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import edu.cnu.mdi.component.CommonBorder;
import edu.cnu.mdi.mapping.MapView2D;
import edu.cnu.mdi.mapping.projection.EProjection;
import edu.cnu.mdi.mapping.theme.MapTheme;
import edu.cnu.mdi.ui.fonts.Fonts;

/**
 * app-specific map control panel.
 *
 * <p>
 * This panel demonstrates how an MDI application can replace the standard
 * {@code MapControlPanel} and include an application-supplied projection
 * without adding that projection to MDI's built-in {@code EProjection} enum.
 * </p>
 */
@SuppressWarnings("serial")
public class ChimeraMapControlPanel extends JPanel {

    /** Controlled view. */
    private final MapView2D mapView;

    /** Control font. */
    private final Font font = Fonts.plainFontDelta(-2);

    /** Current theme. */
    private MapTheme currentTheme = MapTheme.light();

    /**
     * Creates the map control panel.
     *
     * @param mapView the controlled map view
     */
    public ChimeraMapControlPanel(MapView2D mapView) {
        if (mapView == null) {
            throw new IllegalArgumentException("mapView must not be null.");
        }

        this.mapView = mapView;

        setAlignmentX(Component.CENTER_ALIGNMENT);
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        setBorder(BorderFactory.createEmptyBorder(4, 4, 4, 4));

        createProjectionSelector();
    }
    
    @Override
    public Dimension getMaximumSize() {
        return new Dimension(Integer.MAX_VALUE, getPreferredSize().height);
    }

    /**
     * Gets the current theme.
     *
     * @return the current theme
     */
    public MapTheme getCurrentTheme() {
        return currentTheme;
    }

    /**
     * Creates the projection selector.
     */
    private void createProjectionSelector() {
        JPanel panel = new JPanel();
        panel.setLayout(new GridLayout(3, 2, 4, 2));
        panel.setAlignmentX(Component.LEFT_ALIGNMENT);
        panel.setBorder(new CommonBorder("Projection"));

        ButtonGroup group = new ButtonGroup();

        JRadioButton archimedes = radio("Archimedes/Lambert", true);
        JRadioButton mercator = radio("Mercator", false);
        JRadioButton mollweide = radio("Mollweide", false);
        JRadioButton orthographic = radio("Orthographic", false);
        JRadioButton lambertEqualArea = radio("Lambert Equal-Area", false);

        archimedes.addActionListener(e ->
                mapView.setProjection(new ArchimedesLambertCylindricalProjection(currentTheme)));

        mercator.addActionListener(e -> mapView.setProjection(EProjection.MERCATOR));
        mollweide.addActionListener(e -> mapView.setProjection(EProjection.MOLLWEIDE));
        orthographic.addActionListener(e -> mapView.setProjection(EProjection.ORTHOGRAPHIC));
        lambertEqualArea.addActionListener(e -> mapView.setProjection(EProjection.LAMBERT_EQUAL_AREA));

        group.add(archimedes);
        group.add(mercator);
        group.add(mollweide);
        group.add(orthographic);
        group.add(lambertEqualArea);

        panel.add(archimedes);
        panel.add(mercator);
        panel.add(mollweide);
        panel.add(orthographic);
        panel.add(lambertEqualArea);

        add(panel);
    }


    /**
     * Creates a radio button with standard map styling.
     *
     * @param label label text
     * @param selected initial selected state
     * @return the radio button
     */
    private JRadioButton radio(String label, boolean selected) {
        JRadioButton button = new JRadioButton(label, selected);
        button.setFont(font);
        button.setAlignmentX(Component.LEFT_ALIGNMENT);
        return button;
    }
}