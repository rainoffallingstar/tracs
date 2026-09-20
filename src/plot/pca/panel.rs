//! PCA panel rendering: the sample scatter plot and the variance scree plot.
//!
//! Split out from [`crate::plot::pca`] so the numerical core stays free of
//! rendering dependencies. The oracle test for `prcomp()` compatibility includes
//! `pca.rs` directly and therefore cannot resolve `crate::plot::*`; keeping the
//! drawing here lets that test stay a faithful, dependency-free comparison.
//!
//! Layout follows `pca_plot()`: a scatter panel of samples on two components,
//! optionally beside a scree panel showing each component's share of the
//! variance with the cumulative curve over it.

use crate::plot::pca::Pca;
use crate::plot::svg::PanelWriter;

/// One plotted sample: where it falls on the two chosen components, and how it
/// should look.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatterPoint {
    pub label: String,
    pub x: f64,
    pub y: f64,
    pub color: String,
}

/// Builds the scatter points for one pair of components.
///
/// Mirrors `pca_plot()`'s `pca_dat[, xpc]` / `[, ypc]` extraction and its merge
/// against the colour (and optional shape) lookup. Samples whose score is not
/// finite on either axis are dropped, matching R's `complete.cases()` behaviour
/// after the merge.
pub fn scatter_points(
    pca: &Pca,
    x_component: usize,
    y_component: usize,
    sample_names: &[String],
    colors: &[String],
) -> Vec<ScatterPoint> {
    let mut points = Vec::with_capacity(pca.n_samples);
    for sample in 0..pca.n_samples {
        let x = pca
            .components
            .get(x_component)
            .and_then(|component| component.scores.get(sample))
            .copied()
            .unwrap_or(f64::NAN);
        let y = pca
            .components
            .get(y_component)
            .and_then(|component| component.scores.get(sample))
            .copied()
            .unwrap_or(f64::NAN);
        if !x.is_finite() || !y.is_finite() {
            continue;
        }
        points.push(ScatterPoint {
            label: sample_names
                .get(sample)
                .cloned()
                .unwrap_or_else(|| format!("sample{}", sample + 1)),
            x,
            y,
            color: colors
                .get(sample)
                .cloned()
                .unwrap_or_else(|| DEFAULT_POINT_COLOR.to_string()),
        });
    }
    points
}

/// Default point colour when no grouping is supplied, matching `pca_plot()`'s
/// `group_df$color = "black"` branch.
pub const DEFAULT_POINT_COLOR: &str = "black";

/// Formats an axis title the way `pca_plot()` does:
/// `paste0(xpc, " [", round(var_explained, 2), "]")`.
///
/// R's `round()` is "round half to even", and the share is already a fraction,
/// so the value is printed with two decimals rather than as a percentage.
pub fn axis_title(component_name: &str, variance_explained: f64) -> String {
    format!(
        "{} [{}]",
        component_name,
        format_two_decimals(variance_explained)
    )
}

/// Formats a fraction with two decimals, matching R's `round(x, 2)` printing.
///
/// Rust's `{:.2}` rounds half away from zero while R rounds half to even, so the
/// tie case is handled explicitly. This only affects exact halves, which for a
/// variance share means values like 0.125.
fn format_two_decimals(value: f64) -> String {
    if !value.is_finite() {
        return "NA".to_string();
    }
    let scaled = value * 100.0;
    let floor = scaled.floor();
    let fraction = scaled - floor;
    let rounded = if (fraction - 0.5).abs() < 1e-9 {
        // Half-to-even, as R does.
        if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        scaled.round()
    };
    // `-0.00` would print for a tiny negative share; R normalises that to `0`.
    let rounded = if rounded == 0.0 { 0.0 } else { rounded };
    format!("{:.2}", rounded / 100.0)
}

/// Draws the PCA scatter panel: gridlines, points and their sample labels,
/// axis titles carrying the variance explained, and the legend for `color_by`.
///
/// Ported from `pca_plot()`. The axis ranges come from `pretty()` on the scores,
/// so the gridlines land on the same round numbers R would choose.
#[allow(clippy::too_many_arguments)]
pub fn draw_scatter_panel(
    panel: &mut PanelWriter,
    points: &[ScatterPoint],
    x_title: &str,
    y_title: &str,
    legend_entries: &[(String, String)],
    show_axis: bool,
    label_size: f64,
    point_size: f64,
    font_size: f64,
) {
    if points.is_empty() {
        return;
    }

    // R's margins: `par(mar = c(3, 4, 2, 1))`. The bottom and left margins also
    // have to hold the axis titles, which sit two lines out.
    let line_height = font_size * 1.2;
    let plot_left = 5.0 * line_height;
    let plot_right = panel.width() - 1.0 * line_height;
    let plot_top = 2.0 * line_height;
    let plot_bottom = panel.height() - 3.5 * line_height;
    if plot_right <= plot_left || plot_bottom <= plot_top {
        return;
    }

    // `xlim = range(pretty(scores))` keeps zero inside the view so the origin
    // cross is always visible, which is what makes the sign of a component
    // readable.
    let x_pretty = pretty_of(points.iter().map(|point| point.x));
    let y_pretty = pretty_of(points.iter().map(|point| point.y));
    if !x_pretty.low().is_finite() || !y_pretty.low().is_finite() {
        return;
    }

    let map_x = |value: f64| -> f64 {
        let span = x_pretty.high() - x_pretty.low();
        if span <= 0.0 {
            return (plot_left + plot_right) / 2.0;
        }
        plot_left + (value - x_pretty.low()) / span * (plot_right - plot_left)
    };
    let map_y = |value: f64| -> f64 {
        let span = y_pretty.high() - y_pretty.low();
        if span <= 0.0 {
            return (plot_bottom + plot_top) / 2.0;
        }
        plot_bottom - (value - y_pretty.low()) / span * (plot_bottom - plot_top)
    };

    // Gridlines at the pretty ticks, then a slightly stronger line through the
    // origin: `abline(h = pretty(x), v = pretty(y), lwd = 0.1)` followed by
    // `abline(h = 0, v = 0, lwd = 0.8)`.
    for value in &x_pretty.values {
        let x = map_x(*value);
        panel.line(x, plot_top, x, plot_bottom, "gray90", 1.0);
    }
    for value in &y_pretty.values {
        let y = map_y(*value);
        panel.line(plot_left, y, plot_right, y, "gray90", 1.0);
    }
    if x_pretty.low() <= 0.0 && x_pretty.high() >= 0.0 {
        let x = map_x(0.0);
        panel.line(x, plot_top, x, plot_bottom, "gray70", 1.6);
    }
    if y_pretty.low() <= 0.0 && y_pretty.high() >= 0.0 {
        let y = map_y(0.0);
        panel.line(plot_left, y, plot_right, y, "gray70", 1.6);
    }

    // Scatter points with their labels just above (`pos = 3`).
    let radius = (2.5 * point_size).max(1.5);
    for point in points {
        let x = map_x(point.x);
        let y = map_y(point.y);
        panel.circle(x, y, radius, &point.color);
        if label_size > 0.0 {
            panel.text(
                x,
                y - radius - 2.0,
                &point.label,
                font_size * label_size,
                "middle",
                &point.color,
            );
        }
    }

    if show_axis {
        for value in &x_pretty.values {
            let x = map_x(*value);
            panel.line(x, plot_bottom, x, plot_bottom + 4.0, "black", 1.0);
            panel.text(
                x,
                plot_bottom + 6.0 + font_size,
                &crate::plot::svg::format_tick(*value),
                font_size * 0.8,
                "middle",
                "black",
            );
        }
        for value in &y_pretty.values {
            let y = map_y(*value);
            panel.line(plot_left - 4.0, y, plot_left, y, "black", 1.0);
            panel.text(
                plot_left - 6.0,
                y + font_size * 0.3,
                &crate::plot::svg::format_tick(*value),
                font_size * 0.8,
                "end",
                "black",
            );
        }
    }

    // Axis titles carrying the variance explained, matching `mtext(line = 2)`.
    panel.text(
        (plot_left + plot_right) / 2.0,
        panel.height() - font_size * 0.6,
        x_title,
        font_size * 0.8,
        "middle",
        "black",
    );
    panel.text_rotated(
        font_size * 0.9,
        (plot_top + plot_bottom) / 2.0,
        y_title,
        font_size * 0.8,
        "middle",
        "black",
        -90.0,
    );

    // `color_by` legend at the requested corner, drawn inside the panel because
    // R's `xpd = TRUE` still clips to the figure.
    for (index, (label, color)) in legend_entries.iter().enumerate() {
        let y = plot_top + index as f64 * (font_size + 1.0) - font_size * 0.3;
        let text_left = plot_right - 6.0;
        panel.circle(text_left - 10.0, y, 2.5, color);
        panel.text(text_left, y, label, font_size, "end", color);
    }
}

/// Draws the variance-explained scree panel.
///
/// Ported from `pca_plot()`'s `show_cree` block: bars of decreasing height with
/// a cumulative line over them. The right axis carries the cumulative curve
/// because both quantities live on 0..1, which is what lets R use a single
/// y-range for the pair.
pub fn draw_scree_panel(
    panel: &mut PanelWriter,
    component_names: &[String],
    variance_explained: &[f64],
    show_axis: bool,
    font_size: f64,
) {
    if component_names.is_empty() {
        return;
    }

    // R's `par(mar = c(3, 4, 2, 4))`: both sides have room for their own axis.
    let line_height = font_size * 1.2;
    let plot_left = 4.5 * line_height;
    let plot_right = panel.width() - 4.5 * line_height;
    let plot_top = 2.0 * line_height;
    let plot_bottom = panel.height() - 4.0 * line_height;
    if plot_right <= plot_left || plot_bottom <= plot_top {
        return;
    }

    // `barplot(ylim = c(0, 1))` pins the range, so the bars stay comparable
    // between figures even when PC1 explains little of the variance.
    let map_y = |value: f64| -> f64 { plot_bottom - value.clamp(0.0, 1.0) * (plot_bottom - plot_top) };

    let count = component_names.len();
    let slot = (plot_right - plot_left) / count as f64;
    // R's `barplot` leaves a gap between bars; 0.8 is the familiar default.
    let bar_width = slot * 0.8;

    for (index, value) in variance_explained.iter().enumerate() {
        let left = plot_left + index as f64 * slot + (slot - bar_width) / 2.0;
        let top = map_y(*value);
        panel.rect(
            left,
            top,
            bar_width,
            plot_bottom - top,
            SCREE_BAR_COLOR,
        );
    }

    // Cumulative curve, matching `points(..., type = "l", lty = 2)` plus the
    // abline-style markers at each component.
    let mut cumulative = 0.0;
    let mut path = String::new();
    for (index, value) in variance_explained.iter().enumerate() {
        cumulative += value;
        let x = plot_left + index as f64 * slot + slot / 2.0;
        let y = map_y(cumulative);
        if index == 0 {
            let _ = std::fmt::Write::write_fmt(&mut path, format_args!("M {x:.3} {y:.3}"));
        } else {
            let _ = std::fmt::Write::write_fmt(&mut path, format_args!(" L {x:.3} {y:.3}"));
        }
    }
    panel.path(&path, SCREE_CUMULATIVE_COLOR, 1.2);

    let mut cumulative = 0.0;
    for (index, value) in variance_explained.iter().enumerate() {
        cumulative += value;
        let x = plot_left + index as f64 * slot + slot / 2.0;
        panel.circle(x, map_y(cumulative), 2.0, SCREE_CUMULATIVE_COLOR);
    }

    if show_axis {
        for step in 0..=10 {
            let value = step as f64 / 10.0;
            let y = map_y(value);
            panel.line(plot_left - 4.0, y, plot_left, y, "black", 1.0);
            panel.text(
                plot_left - 6.0,
                y + font_size * 0.3,
                &format!("{value:.1}"),
                font_size * 0.8,
                "end",
                "black",
            );
            // The cumulative axis on the right, as `axis(side = 4)`.
            panel.line(plot_right, y, plot_right + 4.0, y, "black", 1.0);
            panel.text(
                plot_right + 6.0,
                y + font_size * 0.3,
                &format!("{value:.1}"),
                font_size * 0.8,
                "start",
                "black",
            );
        }

        // Component labels along the bottom, rotated like `las = 2`.
        for (index, name) in component_names.iter().enumerate() {
            let x = plot_left + index as f64 * slot + slot / 2.0;
            panel.text_rotated(
                x,
                plot_bottom + font_size * 0.4,
                name,
                font_size * 0.8,
                "end",
                "black",
                -90.0,
            );
        }
    }

    // Axis titles, matching the two `mtext()` calls.
    panel.text_rotated(
        font_size * 0.8,
        (plot_top + plot_bottom) / 2.0,
        "var. explained",
        font_size * 0.8,
        "middle",
        "black",
        -90.0,
    );
    panel.text_rotated(
        panel.width() - font_size * 0.6,
        (plot_top + plot_bottom) / 2.0,
        "cumulative var. explained",
        font_size * 0.8,
        "middle",
        "black",
        -90.0,
    );
}

/// Bar colour in the scree panel (`pca_plot()`'s `col = "#2c3e50"`).
pub const SCREE_BAR_COLOR: &str = "#2c3e50";
/// Cumulative-curve colour (`pca_plot()`'s `col = "#c0392b"`).
pub const SCREE_CUMULATIVE_COLOR: &str = "#c0392b";

/// Collects the finite values and runs them through R's `pretty()`.
fn pretty_of(values: impl Iterator<Item = f64>) -> crate::plot::pretty::Pretty {
    let collected: Vec<f64> = values.filter(|value| value.is_finite()).collect();
    if collected.is_empty() {
        return crate::plot::pretty::pretty_range(0.0, 1.0, 5);
    }
    let min = collected.iter().copied().fold(f64::INFINITY, f64::min);
    let max = collected.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    crate::plot::pretty::pretty_range(min, max, 5)
}
