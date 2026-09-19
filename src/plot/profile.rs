//! Profile plots: signal aggregated around a focal point (usually a TSS).
//!
//! Ported from `trackplot.R`'s `profile_summarize()` / `profile_plot()`.
//!
//! Unlike the track stack, a profile plot is a single XY panel: every sample
//! becomes one line, x runs from `up` bases upstream to `down` bases
//! downstream, and each column is the mean (or median) across all regions.
//!
//! Strand handling lives in the *region* set rather than here: R builds a BED
//! whose start column is the TSS for plus-strand transcripts and whose end
//! column is the TSS for minus-strand ones, then always anchors with
//! `matrix -starts`.

use std::fmt::Write as _;

use crate::plot::svg::{format_tick, PanelWriter, YAxis};

/// Default line colours, taken from `profile_plot()`.
pub const PROFILE_COLORS: [&str; 12] = [
    "#2f4f4f", "#8b4513", "#228b22", "#00008b", "#ff0000", "#ffd700", "#7fff00", "#00ffff",
    "#ff00ff", "#6495ed", "#ffe4b5", "#ff69b4",
];

/// How replicate values are collapsed into one profile line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SummaryStat {
    Mean,
    Median,
}

impl SummaryStat {
    /// Parses the CLI spelling, matching `profile_summarize()`'s `stat`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "mean" => Some(SummaryStat::Mean),
            "median" => Some(SummaryStat::Median),
            _ => None,
        }
    }
}

/// Collapses a numeric column, skipping non-finite values like R's `na.rm`.
fn collapse(values: &[f64], stat: SummaryStat) -> f64 {
    let mut finite: Vec<f64> = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect();
    if finite.is_empty() {
        return f64::NAN;
    }
    match stat {
        SummaryStat::Mean => finite.iter().sum::<f64>() / finite.len() as f64,
        SummaryStat::Median => {
            finite.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mid = finite.len() / 2;
            if finite.len().is_multiple_of(2) {
                (finite[mid - 1] + finite[mid]) / 2.0
            } else {
                finite[mid]
            }
        }
    }
}

/// Summarizes matrix rows column-wise, as `profile_summarize()` does.
///
/// Each input row is one region and each column one bin; the output has one
/// value per bin. Rows may differ in length, in which case the longest row
/// defines the output width and shorter rows contribute only where they have
/// data.
pub fn summarize_columns(rows: &[Vec<f64>], stat: SummaryStat) -> Vec<f64> {
    let width = rows.iter().map(Vec::len).max().unwrap_or(0);
    (0..width)
        .map(|column| {
            let values: Vec<f64> = rows
                .iter()
                .filter_map(|row| row.get(column).copied())
                .collect();
            collapse(&values, stat)
        })
        .collect()
}

/// Summarizes several samples column-wise, preserving order.
pub fn summarize_samples(
    samples: &[(String, Vec<Vec<f64>>)],
    stat: SummaryStat,
) -> Vec<(String, Vec<f64>)> {
    samples
        .iter()
        .map(|(name, rows)| (name.clone(), summarize_columns(rows, stat)))
        .collect()
}

/// Summarizes replicates that share a condition.
///
/// Matches `profile_summarize(condition = ...)`: all rows from every sample in a
/// group are pooled before collapsing, so the result has one line per distinct
/// condition in first-seen order.
pub fn summarize_by_condition(
    samples: &[(String, Vec<Vec<f64>>)],
    conditions: &[String],
    stat: SummaryStat,
) -> Vec<(String, Vec<f64>)> {
    let mut order: Vec<String> = Vec::new();
    for condition in conditions {
        if !order.contains(condition) {
            order.push(condition.clone());
        }
    }

    order
        .into_iter()
        .map(|condition| {
            let pooled: Vec<Vec<f64>> = samples
                .iter()
                .zip(conditions.iter())
                .filter(|(_, sample_condition)| **sample_condition == condition)
                .flat_map(|((_, rows), _)| rows.clone())
                .collect();
            (condition, summarize_columns(&pooled, stat))
        })
        .collect()
}

/// One line in a profile plot.
#[derive(Clone, Debug, PartialEq)]
pub struct ProfileSeries {
    pub name: String,
    pub values: Vec<f64>,
}

/// Resolves the line colours, cycling the built-in palette.
pub fn resolve_profile_colors(requested: &[String], count: usize) -> Vec<String> {
    if count == 0 {
        return Vec::new();
    }
    if requested.is_empty() {
        return (0..count)
            .map(|index| PROFILE_COLORS[index % PROFILE_COLORS.len()].to_string())
            .collect();
    }
    (0..count)
        .map(|index| requested[index % requested.len()].clone())
        .collect()
}

/// Axis ticks for the profile x axis.
///
/// R draws three: the upstream end, the focal point, and the downstream end,
/// labelled with the `up`/`0`/`down` distances. The focal point sits at
/// `len * up / (up + down)`.
pub fn profile_xticks(nbins: usize, up: u32, down: u32) -> (Vec<f64>, Vec<String>) {
    let total = up as f64 + down as f64;
    let focal = if total > 0.0 {
        nbins as f64 * up as f64 / total
    } else {
        0.0
    };
    (
        vec![0.0, focal, nbins as f64],
        vec![up.to_string(), "0".to_string(), down.to_string()],
    )
}

/// Draws the profile panel: gridlines, one polyline per sample, axes and legend.
///
/// Ported from `profile_plot()`. `font_size` replaces R's `cex` scaling; R's
/// `par(mar = c(4, 4, 2, 1))` becomes the inset below.
#[allow(clippy::too_many_arguments)]
pub fn draw_profile_panel(
    panel: &mut PanelWriter,
    series: &[ProfileSeries],
    colors: &[String],
    up: u32,
    down: u32,
    show_axis: bool,
    font_size: f64,
    xlab: &str,
    ylab: &str,
) {
    if series.is_empty() {
        return;
    }

    // R's margins: mar = c(4, 4, 2, 1) lines, at ~1.2x the font size per line.
    // The left margin must additionally house the y tick labels and the rotated
    // y axis label, so it is widened by one line when a label is present.
    let line_height = font_size * 1.2;
    let plot_left = if ylab.is_empty() {
        4.0 * line_height
    } else {
        5.5 * line_height
    };
    let plot_right = if show_axis {
        // The last x tick sits at the plot edge and is centred on it, so it needs
        // half its width of clearance to stay inside the figure.
        panel.width() - 2.5 * line_height
    } else {
        panel.width() - 1.0 * line_height
    };
    let plot_top = 2.0 * line_height;
    let plot_bottom = panel.height() - 4.0 * line_height;
    if plot_right <= plot_left || plot_bottom <= plot_top {
        return;
    }

    let nbins = series.iter().map(|s| s.values.len()).max().unwrap_or(0);
    if nbins == 0 {
        return;
    }

    // R sets ylim from pretty(c(y_min, y_max), n = 5) so the axis ends on a tick.
    let y_min = series
        .iter()
        .flat_map(|s| s.values.iter())
        .copied()
        .filter(|value| value.is_finite())
        .fold(f64::INFINITY, f64::min);
    let y_max = series
        .iter()
        .flat_map(|s| s.values.iter())
        .copied()
        .filter(|value| value.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    if !y_min.is_finite() || !y_max.is_finite() {
        return;
    }
    let y_ticks = crate::plot::pretty::pretty_range(y_min, y_max, 5);
    let axis = YAxis::new(plot_top, plot_bottom, y_ticks.low(), y_ticks.high());

    let map_x = |position: f64| -> f64 {
        plot_left + position * (plot_right - plot_left) / ((nbins - 1).max(1) as f64)
    };

    // Gridlines: horizontal at the y ticks, vertical at the x ticks, matching
    // R's `abline(h = ylabs, v = pretty(xticks), col = "gray90", lty = 2)`.
    let (x_tick_positions, x_tick_labels) = profile_xticks(nbins, up, down);
    for value in &y_ticks.values {
        let y = axis.map(*value);
        panel.dashed_line(plot_left, y, plot_right, "gray90");
    }
    for position in &x_tick_positions {
        let x = map_x(*position);
        panel.line(x, plot_top, x, plot_bottom, "gray90", 1.0);
    }

    // Frame.
    panel.line(plot_left, plot_bottom, plot_right, plot_bottom, "black", 1.0);
    panel.line(plot_left, plot_top, plot_left, plot_bottom, "black", 1.0);

    // One polyline per sample, skipping gaps rather than breaking the line.
    for (index, sample) in series.iter().enumerate() {
        let color = colors
            .get(index)
            .map(String::as_str)
            .unwrap_or(PROFILE_COLORS[0]);
        let mut path = String::new();
        let mut started = false;
        for (bin, value) in sample.values.iter().enumerate() {
            if !value.is_finite() {
                continue;
            }
            let x = map_x(bin as f64);
            let y = axis.map(*value);
            if started {
                let _ = write!(path, " L {x:.3} {y:.3}");
            } else {
                let _ = write!(path, "M {x:.3} {y:.3}");
                started = true;
            }
        }
        if started {
            panel.path(&path, color, 1.5);
        }
    }

    if show_axis {
        // x axis: R's `axis(side = 1, at = xticks, labels = xlabs)`.
        for (position, label) in x_tick_positions.iter().zip(x_tick_labels.iter()) {
            let x = map_x(*position);
            panel.line(x, plot_bottom, x, plot_bottom + 4.0, "black", 1.0);
            panel.text(
                x,
                plot_bottom + 6.0 + font_size,
                label,
                font_size,
                "middle",
                "black",
            );
        }
        // y axis: every pretty() tick, matching `axis(side = 2, at = ylabs)`.
        for value in &y_ticks.values {
            let y = axis.map(*value);
            panel.line(plot_left - 4.0, y, plot_left, y, "black", 1.0);
            panel.text(
                plot_left - 6.0,
                y + font_size * 0.35,
                &format_tick(*value),
                font_size,
                "end",
                "black",
            );
        }
    }

    // Axis labels, matching R's `mtext(..., line = 2.5)`.
    if !xlab.is_empty() {
        panel.text(
            (plot_left + plot_right) / 2.0,
            panel.height() - font_size * 0.5,
            xlab,
            font_size,
            "middle",
            "black",
        );
    }
    if !ylab.is_empty() {
        // R draws this with `mtext(side = 2)`, i.e. rotated 90 degrees, centred
        // in the left margin. A horizontal label would be clipped there.
        panel.text_rotated(
            font_size * 0.8,
            (plot_top + plot_bottom) / 2.0,
            ylab,
            font_size,
            "middle",
            "black",
            -90.0,
        );
    }

    // Legend at top right, matching R's `legend("topright", bty = "n")`. R
    // allows this to overflow the plot region (`xpd = TRUE`) but still draw
    // inside the figure, so the text is right-aligned at the panel edge to
    // avoid running off the canvas.
    for (index, sample) in series.iter().enumerate() {
        let color = colors
            .get(index)
            .map(String::as_str)
            .unwrap_or(PROFILE_COLORS[0]);
        let y = plot_top + index as f64 * (font_size + 1.0) - font_size * 0.3;
        let text_left = plot_right - 6.0;
        panel.line(text_left - 14.0, y, text_left - 2.0, y, color, 2.0);
        panel.text(text_left, y, &sample.name, font_size, "end", color);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_and_median_skip_non_finite_values() {
        let rows = vec![
            vec![1.0, 2.0, f64::NAN],
            vec![3.0, f64::NAN, 30.0],
            vec![5.0, 4.0, 60.0],
        ];
        let mean = summarize_columns(&rows, SummaryStat::Mean);
        assert!((mean[0] - 3.0).abs() < 1e-9, "mean col 0: {}", mean[0]);
        assert!((mean[1] - 3.0).abs() < 1e-9, "mean col 1: {}", mean[1]);
        assert!((mean[2] - 45.0).abs() < 1e-9, "mean col 2: {}", mean[2]);

        let median = summarize_columns(&rows, SummaryStat::Median);
        assert!((median[0] - 3.0).abs() < 1e-9, "median col 0: {}", median[0]);
        assert!((median[1] - 3.0).abs() < 1e-9, "median col 1: {}", median[1]);
        assert!((median[2] - 45.0).abs() < 1e-9, "median col 2: {}", median[2]);
    }

    #[test]
    fn median_averages_the_middle_pair_for_even_counts() {
        let rows = vec![vec![1.0], vec![2.0], vec![10.0], vec![20.0]];
        let median = summarize_columns(&rows, SummaryStat::Median);
        // R's median(c(1, 2, 10, 20)) = (2 + 10) / 2
        assert!((median[0] - 6.0).abs() < 1e-9, "got {}", median[0]);
    }

    #[test]
    fn all_missing_column_yields_nan() {
        let rows = vec![vec![f64::NAN], vec![f64::NAN]];
        assert!(summarize_columns(&rows, SummaryStat::Mean)[0].is_nan());
        assert!(summarize_columns(&[], SummaryStat::Mean).is_empty());
    }

    #[test]
    fn conditions_pool_replicates_in_first_seen_order() {
        let samples = vec![
            ("s1".to_string(), vec![vec![1.0, 1.0]]),
            ("s2".to_string(), vec![vec![3.0, 3.0]]),
            ("s3".to_string(), vec![vec![10.0, 20.0]]),
        ];
        let conditions = vec!["A".to_string(), "A".to_string(), "B".to_string()];
        let summarized = summarize_by_condition(&samples, &conditions, SummaryStat::Mean);

        assert_eq!(summarized.len(), 2);
        assert_eq!(summarized[0].0, "A");
        // A pools s1 and s2: (1+3)/2 = 2 per column.
        assert!((summarized[0].1[0] - 2.0).abs() < 1e-9);
        assert!((summarized[0].1[1] - 2.0).abs() < 1e-9);
        // B is s3 alone.
        assert_eq!(summarized[1].0, "B");
        assert!((summarized[1].1[0] - 10.0).abs() < 1e-9);
    }

    #[test]
    fn xticks_put_the_focal_point_at_the_up_fraction() {
        // Symmetric flanks put the focal point in the middle.
        let (positions, labels) = profile_xticks(100, 2500, 2500);
        assert_eq!(labels, vec!["2500", "0", "2500"]);
        assert_eq!(positions, vec![0.0, 50.0, 100.0]);

        // Asymmetric: up = 1000, down = 3000 => focal at 25 of 100 bins.
        let (positions, labels) = profile_xticks(100, 1000, 3000);
        assert_eq!(labels, vec!["1000", "0", "3000"]);
        assert!((positions[1] - 25.0).abs() < 1e-9, "focal: {}", positions[1]);
    }

    #[test]
    fn colors_cycle_the_builtin_palette() {
        let defaults = resolve_profile_colors(&[], 3);
        assert_eq!(
            defaults,
            vec![PROFILE_COLORS[0], PROFILE_COLORS[1], PROFILE_COLORS[2]]
        );

        // More series than palette entries wraps around.
        let many = resolve_profile_colors(&[], PROFILE_COLORS.len() + 1);
        assert_eq!(many[PROFILE_COLORS.len()], PROFILE_COLORS[0]);

        let requested = resolve_profile_colors(&["#111".to_string()], 2);
        assert_eq!(requested, vec!["#111", "#111"]);
        assert!(resolve_profile_colors(&[], 0).is_empty());
    }

    #[test]
    fn summary_stat_parses_cli_spellings() {
        assert_eq!(SummaryStat::from_name("mean"), Some(SummaryStat::Mean));
        assert_eq!(SummaryStat::from_name(" MEDIAN "), Some(SummaryStat::Median));
        assert_eq!(SummaryStat::from_name("sum"), None);
    }
}
