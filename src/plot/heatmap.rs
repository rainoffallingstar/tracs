//! Profile heatmaps: signal around a focal point, one panel per sample.
//!
//! Ported from `trackplot.R`'s `profile_heatmap()`. Each panel is one bigWig's
//! regions-by-bins matrix, sorted by row mean/median so the strongest regions
//! cluster together, then coloured by a sequential palette.
//!
//! Two behaviours from R are worth calling out because they are easy to get
//! wrong:
//!
//! - R sorts rows **per sample** (`.order_matrix()` is called on each matrix
//!   independently inside the panel loop), so panels are not row-aligned with
//!   each other. That is reproduced here.
//! - The colour ramp runs light-to-dark as the value rises, because R reverses
//!   `hcl.colors()` before handing it to `image()`. [`colors::resolve_ramp`]
//!   keeps that orientation.

pub mod colors;

use crate::plot::svg::PanelWriter;

/// How rows are ordered inside each panel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SortBy {
    Mean,
    Median,
}

impl SortBy {
    /// Parses the CLI spelling, matching `profile_heatmap()`'s `sortBy`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "mean" => Some(SortBy::Mean),
            "median" => Some(SortBy::Median),
            _ => None,
        }
    }
}

/// Row-wise statistic used for sorting.
fn row_stat(row: &[f64], sort_by: SortBy) -> f64 {
    let mut finite: Vec<f64> = row
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect();
    if finite.is_empty() {
        return f64::NAN;
    }
    match sort_by {
        SortBy::Mean => finite.iter().sum::<f64>() / finite.len() as f64,
        SortBy::Median => {
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

/// Sorts rows by descending row statistic, as `.order_matrix()` does.
///
/// R's `order(x, decreasing = TRUE)` places `NA` last, which `sort_by` on a
/// `NaN`-aware key reproduces by treating missing rows as smallest.
pub fn order_rows(matrix: &[Vec<f64>], sort_by: SortBy) -> Vec<Vec<f64>> {
    let mut rows: Vec<(f64, Vec<f64>)> = matrix
        .iter()
        .map(|row| (row_stat(row, sort_by), row.clone()))
        .collect();
    rows.sort_by(|a, b| {
        let a_key = if a.0.is_nan() { f64::NEG_INFINITY } else { a.0 };
        let b_key = if b.0.is_nan() { f64::NEG_INFINITY } else { b.0 };
        b_key
            .partial_cmp(&a_key)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    rows.into_iter().map(|(_, row)| row).collect()
}

/// One heatmap panel: a sample name plus its rows-by-bins matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct HeatmapPanel {
    pub name: String,
    /// Sorted matrix, rows are regions, columns are bins.
    pub matrix: Vec<Vec<f64>>,
    /// Value mapped to the light end of the ramp.
    pub z_min: f64,
    /// Value mapped to the dark end; values above are clamped, as R does.
    pub z_max: f64,
}

/// Resolves the colour limits for one panel.
///
/// `profile_heatmap()` defaults `zmin` to the matrix minimum and `zmax` to the
/// **maximum row mean** rather than the matrix maximum, which keeps a few
/// extreme regions from washing out the rest of the panel. Explicit overrides
/// win.
pub fn resolve_limits(
    matrix: &[Vec<f64>],
    z_min: Option<f64>,
    z_max: Option<f64>,
    sort_by: SortBy,
) -> (f64, f64) {
    let z_min = z_min.unwrap_or_else(|| {
        crate::plot::svg::round_two_of(
            matrix
                .iter()
                .flat_map(|row| row.iter())
                .copied()
                .filter(|value| value.is_finite())
                .fold(f64::INFINITY, f64::min),
        )
    });

    let z_max = z_max.unwrap_or_else(|| {
        // R: max(rowMeans(hm.dat, na.rm = TRUE)) after the transpose/rev dance.
        let row_means: Vec<f64> = matrix.iter().map(|row| row_stat(row, sort_by)).collect();
        row_means
            .iter()
            .copied()
            .filter(|value| value.is_finite())
            .fold(f64::NEG_INFINITY, f64::max)
    });

    // Guard a degenerate range so the ramp lookup cannot divide by zero. Using
    // `partial_cmp` keeps NaN limits (an all-missing matrix) on the guard path
    // rather than silently producing a zero-width scale.
    let ordered = z_max
        .partial_cmp(&z_min)
        .map(|order| order.is_gt())
        .unwrap_or(false);
    if !ordered {
        return (z_min, z_min + 1.0);
    }
    (z_min, z_max)
}

/// Builds sorted panels for every sample, resolving limits per panel.
pub fn build_panels(
    samples: &[(String, Vec<Vec<f64>>)],
    sort_by: SortBy,
    z_mins: Option<&[f64]>,
    z_maxs: Option<&[f64]>,
) -> Vec<HeatmapPanel> {
    samples
        .iter()
        .enumerate()
        .map(|(index, (name, matrix))| {
            let sorted = order_rows(matrix, sort_by);
            let (z_min, z_max) = resolve_limits(
                &sorted,
                z_mins.map(|values| values[index % values.len()]),
                z_maxs.map(|values| values[index % values.len()]),
                sort_by,
            );
            HeatmapPanel {
                name: name.clone(),
                matrix: sorted,
                z_min,
                z_max,
            }
        })
        .collect()
}

/// Draws one heatmap panel with its colour bar and axis labels.
///
/// Mirrors R's `par(mar = c(3, 2, 1.5, 1))` layout: the heatmap occupies the
/// data area, a slim colour bar sits to its left margin, the scale's numeric
/// labels are to the right of that bar, the sample name is the title, and the
/// x axis is labelled with the `-up / 0 / down` distances.
#[allow(clippy::too_many_arguments)]
pub fn draw_heatmap_panel(
    panel: &mut PanelWriter,
    heatmap: &HeatmapPanel,
    ramp: &[String],
    up: u32,
    down: u32,
    show_axis: bool,
    font_size: f64,
) {
    let rows = heatmap.matrix.len();
    let columns = heatmap.matrix.iter().map(Vec::len).max().unwrap_or(0);
    if rows == 0 || columns == 0 {
        return;
    }

    let line_height = font_size * 1.2;
    // The left margin holds the colour bar *and* its numeric labels outside it
    // (R draws the labels to the left of the bar via `mtext(side = 2)`), so it
    // is wider than a track panel's.
    let plot_left = 6.0 * line_height;
    let plot_right = panel.width() - 1.0 * line_height;
    let plot_top = 1.5 * line_height;
    let plot_bottom = panel.height() - 3.0 * line_height;
    if plot_right <= plot_left || plot_bottom <= plot_top {
        return;
    }

    let cell_width = (plot_right - plot_left) / columns as f64;
    let cell_height = (plot_bottom - plot_top) / rows as f64;
    let span = heatmap.z_max - heatmap.z_min;

    for (row_index, row) in heatmap.matrix.iter().enumerate() {
        for (column_index, value) in row.iter().enumerate() {
            if !value.is_finite() {
                continue;
            }
            // R clamps above zmax: `hm.dat[hm.dat >= zmax] = zmax`.
            let fraction = ((value - heatmap.z_min) / span).clamp(0.0, 1.0);
            let color = ramp_color(ramp, fraction);
            panel.rect(
                plot_left + column_index as f64 * cell_width,
                plot_top + row_index as f64 * cell_height,
                cell_width + 0.5,
                cell_height + 0.5,
                color,
            );
        }
    }

    // Frame, matching R's `rect(xleft = 0, ybottom = 0, xright = 1, ytop = 1)`.
    panel.rect_stroked(
        plot_left,
        plot_top,
        plot_right - plot_left,
        plot_bottom - plot_top,
        "none",
        "black",
        1.0,
    );

    // Title: the sample name.
    panel.text(
        (plot_left + plot_right) / 2.0,
        plot_top - font_size * 0.4,
        &heatmap.name,
        font_size,
        "middle",
        "black",
    );

    if show_axis {
        // Colour bar in the left margin, drawn dark-at-top like R's `image()`
        // over `seq(0, 1, length.out = length(hmcols) - 1)`.
        let bar_right = plot_left - 3.0;
        let bar_left = bar_right - 0.5 * line_height;
        let steps = ramp.len().max(2);
        for step in 0..steps {
            let fraction = step as f64 / (steps - 1) as f64;
            let color = ramp_color(ramp, fraction);
            let y_top = plot_bottom - fraction * (plot_bottom - plot_top);
            let y_bottom = plot_bottom - (step as f64 + 1.0) / (steps - 1) as f64
                * (plot_bottom - plot_top);
            panel.rect(
                bar_left,
                y_bottom.max(plot_top),
                bar_right - bar_left,
                (y_top - y_bottom).abs().max(1.0),
                color,
            );
        }

        // Five numeric labels spanning zmin..zmax, as R's `mtext(round(seq(...)))`.
        for step in 0..5 {
            let fraction = step as f64 / 4.0;
            let value = heatmap.z_min + fraction * span;
            let y = plot_bottom - fraction * (plot_bottom - plot_top);
            panel.text(
                bar_left - 3.0,
                y + font_size * 0.35,
                &crate::plot::svg::format_tick(crate::plot::svg::round_two_of(value)),
                font_size,
                "end",
                "black",
            );
        }

        // x axis: three ticks at the up / 0 / down distances, matching R's
        // `mtext(c(paste0("-", xlabs[1]), xlabs[2], xlabs[3]))`.
        let ticks = crate::plot::profile::profile_xticks(columns, up, down);
        let labels = [
            format!("-{}", ticks.1[0]),
            ticks.1[1].clone(),
            ticks.1[2].clone(),
        ];
        for (position, label) in ticks.0.iter().zip(labels.iter()) {
            let x = plot_left + position / (columns as f64) * (plot_right - plot_left);
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
    }
}

/// Samples the ramp at `fraction` in `0..=1`, matching R's `image()` behaviour
/// of mapping the bottom of the range to the first colour.
fn ramp_color(ramp: &[String], fraction: f64) -> &str {
    if ramp.is_empty() {
        return "black";
    }
    let index = (fraction.clamp(0.0, 1.0) * (ramp.len() - 1) as f64).round() as usize;
    ramp.get(index).map(String::as_str).unwrap_or("black")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_sort_by_descending_statistic() {
        let matrix = vec![
            vec![1.0, 1.0],   // mean 1
            vec![10.0, 10.0], // mean 10
            vec![5.0, 5.0],   // mean 5
        ];
        let sorted = order_rows(&matrix, SortBy::Mean);
        assert_eq!(sorted[0], vec![10.0, 10.0]);
        assert_eq!(sorted[1], vec![5.0, 5.0]);
        assert_eq!(sorted[2], vec![1.0, 1.0]);
    }

    #[test]
    fn median_sorting_differs_from_mean_when_a_row_has_an_outlier() {
        // Row A: mean 34, median 2. Row B: mean 10, median 10.
        let matrix = vec![vec![2.0, 2.0, 98.0], vec![10.0, 10.0, 10.0]];
        let by_mean = order_rows(&matrix, SortBy::Mean);
        assert_eq!(by_mean[0], vec![2.0, 2.0, 98.0], "mean puts the outlier row first");

        let by_median = order_rows(&matrix, SortBy::Median);
        assert_eq!(
            by_median[0],
            vec![10.0, 10.0, 10.0],
            "median ranks the consistently-high row first"
        );
    }

    #[test]
    fn missing_rows_sort_last() {
        let matrix = vec![vec![f64::NAN, f64::NAN], vec![3.0, 3.0]];
        let sorted = order_rows(&matrix, SortBy::Mean);
        assert_eq!(sorted[0], vec![3.0, 3.0]);
        assert!(sorted[1][0].is_nan(), "all-missing row should sink to the bottom");
    }

    #[test]
    fn limits_default_to_minimum_and_max_row_mean() {
        // R uses the matrix min for zmin but the max row *mean* for zmax, so an
        // extreme single row does not wash out the panel.
        let matrix = vec![vec![0.0, 0.0], vec![5.0, 5.0], vec![100.0, 0.0]];
        let (z_min, z_max) = resolve_limits(&matrix, None, None, SortBy::Mean);
        assert!((z_min - 0.0).abs() < 1e-9, "zmin: {z_min}");
        // Row means are 0, 5 and 50 -> zmax is 50, not the matrix max of 100.
        assert!((z_max - 50.0).abs() < 1e-9, "zmax should be the max row mean: {z_max}");
    }

    #[test]
    fn explicit_limits_win_and_degenerate_ranges_are_widened() {
        let matrix = vec![vec![1.0, 2.0]];
        let (z_min, z_max) = resolve_limits(&matrix, Some(0.0), Some(9.0), SortBy::Mean);
        assert_eq!((z_min, z_max), (0.0, 9.0));

        // A constant matrix would otherwise divide by zero.
        let flat = vec![vec![3.0, 3.0]];
        let (z_min, z_max) = resolve_limits(&flat, None, None, SortBy::Mean);
        assert!(z_max > z_min, "degenerate range not widened: {z_min}..{z_max}");
    }

    #[test]
    fn ramp_runs_light_to_dark_as_values_rise() {
        let ramp = vec![
            "#F4FAFE".to_string(), // light
            "#7FABD3".to_string(),
            "#273871".to_string(), // dark
        ];
        assert_eq!(ramp_color(&ramp, 0.0), "#F4FAFE", "low values are light");
        assert_eq!(ramp_color(&ramp, 1.0), "#273871", "high values are dark");
        assert_eq!(ramp_color(&ramp, 0.5), "#7FABD3");
        // Out-of-range fractions clamp rather than panicking.
        assert_eq!(ramp_color(&ramp, -1.0), "#F4FAFE");
        assert_eq!(ramp_color(&ramp, 2.0), "#273871");
        assert_eq!(ramp_color(&[], 0.5), "black");
    }

    #[test]
    fn build_panels_sorts_each_sample_independently() {
        // Different rows are strongest per sample, so the two panels must not
        // end up in the same order (R sorts inside the panel loop).
        let samples = vec![
            ("s1".to_string(), vec![vec![1.0], vec![9.0]]),
            ("s2".to_string(), vec![vec![9.0], vec![1.0]]),
        ];
        let panels = build_panels(&samples, SortBy::Mean, None, None);
        assert_eq!(panels.len(), 2);
        assert_eq!(panels[0].matrix[0], vec![9.0], "s1's strongest row first");
        assert_eq!(panels[1].matrix[0], vec![9.0], "s2's strongest row first");
        // ...which means row 0 corresponds to different regions across panels.
        assert_eq!(panels[0].matrix[1], vec![1.0]);
        assert_eq!(panels[1].matrix[1], vec![1.0]);
    }

    #[test]
    fn build_panels_applies_per_sample_overrides() {
        let samples = vec![
            ("a".to_string(), vec![vec![1.0, 2.0]]),
            ("b".to_string(), vec![vec![3.0, 4.0]]),
        ];
        // A single override is cycled across samples, as R warns and does.
        let panels = build_panels(&samples, SortBy::Mean, Some(&[0.0]), Some(&[10.0]));
        assert_eq!(panels[0].z_min, 0.0);
        assert_eq!(panels[0].z_max, 10.0);
        assert_eq!(panels[1].z_min, 0.0);
        assert_eq!(panels[1].z_max, 10.0);
    }

    #[test]
    fn sort_by_parses_cli_spellings() {
        assert_eq!(SortBy::from_name("mean"), Some(SortBy::Mean));
        assert_eq!(SortBy::from_name(" MEDIAN "), Some(SortBy::Median));
        assert_eq!(SortBy::from_name("hclust"), None);
    }
}
