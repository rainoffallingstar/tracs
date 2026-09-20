//! Volcano plots of differential-binding results.
//!
//! Ported from `trackplot.R`'s `volcano_plot()`, which is deliberately *not*
//! tied to limma: R's version reads whatever `diffpeak()` returns, and the plot
//! itself only needs three numbers per peak -- a log fold change, a p-value, and
//! an adjusted p-value.
//!
//! Because of that, this port takes a generic results table rather than
//! reimplementing limma. Any tool that writes those three columns works:
//! `limma::topTable()`, `DESeq2::results()`, `edgeR::topTags()`, or a
//! hand-assembled TSV. That keeps the empirical-Bayes statistics out of scope,
//! which is the part of the R-only pipeline with no drop-in Rust equivalent.
//!
//! The classification rule is R's exactly: a peak is *significant* when
//! `adj.P.Val < fdr`, and within that set it is *up* when `logFC > 0` and *down*
//! when `logFC < 0`. Two consequences of that rule are easy to get wrong and are
//! reproduced on purpose:
//!
//! - A peak whose `adj.P.Val` is missing is significant in neither direction, so
//!   it is drawn with the non-significant points rather than dropped.
//! - A peak with `logFC` of 0 (or missing) counts as neither up nor down, so the
//!   legend counts can sum to less than the number of significant peaks.
//!
//! One place this port deliberately improves on R: `volcano_plot()` computes
//! `xlims = range(res$logFC)` and R's `range()` returns `NA` if *any* fold change
//! is missing, so `plot()` aborts with "need finite 'xlim' values" and the entire
//! figure is lost -- even when every other peak was perfectly drawable. Here the
//! unusable peaks are skipped and the rest are plotted, and
//! [`Volcano::skipped`] reports how many were dropped so the caller can say so
//! instead of silently producing a sparser plot than the table implies.
//!
//! The one input both implementations refuse is a p-value of exactly 0, because
//! `-log10(0)` is infinite and R's `ylim` becomes `Inf`; that is reported as an
//! error with the same reason R's failure has, not worked around.

use crate::plot::svg::PanelWriter;

/// Default colours, taken from `volcano_plot()`'s signature.
pub const DEFAULT_UP_COLOR: &str = "#d35400";
pub const DEFAULT_DOWN_COLOR: &str = "#1abc9c";
/// Colour of the non-significant points, as `adjustcolor("gray", ...)` renders.
pub const NON_SIGNIFICANT_COLOR: &str = "#BEBEBE";
/// Axis-title and title colour (`volcano_plot()`'s `col = "#34495e"`).
pub const LABEL_COLOR: &str = "#34495e";

/// One row of a differential-binding table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Peak {
    pub log_fold_change: f64,
    /// Raw p-value; the y axis is `-log10` of this.
    pub p_value: f64,
    pub adjusted_p_value: f64,
}

/// Which group a peak falls into, matching `volcano_plot()`'s point calls.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeakClass {
    /// Drawn in the non-significant colour, or not at all when unusable.
    NotSignificant,
    Up,
    Down,
}

impl Peak {
    /// `-log10(p)`, the y coordinate R plots.
    ///
    /// A zero p-value gives infinity, which is why [`Volcano::y_max`] has to
    /// guard against it.
    pub fn neg_log10_p(&self) -> f64 {
        -self.p_value.log10()
    }

    /// Classifies the peak with R's rule.
    ///
    /// `adj.P.Val < fdr` marks it significant, then the sign of `logFC` picks the
    /// side. A significant peak with a zero or missing fold change is neither,
    /// which matches `volcano_plot()`'s `logFC < 0` / `logFC > 0` counts.
    pub fn classify(&self, fdr: f64) -> PeakClass {
        // R's rule is `adj.P.Val < fdr`, and the complement is what collects the
        // rest. The complement is taken as a negation of the *same* comparison
        // rather than as `adj.P.Val >= fdr`, because those differ for NaN:
        // `!(NaN < fdr)` is true, so a missing adjusted p-value lands in the
        // non-significant group exactly as `!adj.P.Val < fdr` does in R, whereas
        // `NaN >= fdr` would wrongly promote it to significant.
        let is_significant = self.adjusted_p_value < fdr;
        if !is_significant {
            return PeakClass::NotSignificant;
        }
        if self.log_fold_change > 0.0 {
            PeakClass::Up
        } else if self.log_fold_change < 0.0 {
            PeakClass::Down
        } else {
            // Significant but with no direction (logFC is 0 or NA): R draws it
            // in the non-significant colour because both `ifelse` branches miss.
            PeakClass::NotSignificant
        }
    }

    /// Whether the peak has finite coordinates and can be drawn at all.
    ///
    /// R would pass non-finite values to `points()` and they would be silently
    /// omitted, so they are filtered here rather than mapped to a bogus position.
    pub fn is_drawable(&self) -> bool {
        self.neg_log10_p().is_finite() && self.log_fold_change.is_finite()
    }
}

/// A parsed results table plus the derived axis limits.
#[derive(Clone, Debug, PartialEq)]
pub struct Volcano {
    pub peaks: Vec<Peak>,
    /// Title, taken from the table's contrast attribute when present.
    pub title: String,
}

impl Volcano {
    /// Builds a volcano plot from peaks and an optional contrast label.
    pub fn new(peaks: Vec<Peak>, title: String) -> Self {
        Self { peaks, title }
    }

    /// The peaks that will actually be drawn, in table order.
    pub fn drawable(&self) -> impl Iterator<Item = &Peak> {
        self.peaks.iter().filter(|peak| peak.is_drawable())
    }

    /// `range(res$logFC)` over the drawable peaks.
    ///
    /// Returns `None` when nothing is drawable, which R would also fail on (its
    /// `range()` yields `NA` and `plot()` errors with "need finite 'xlim'
    /// values"), so the caller can report that rather than emit an empty figure.
    pub fn x_range(&self) -> Option<(f64, f64)> {
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for peak in self.drawable() {
            low = low.min(peak.log_fold_change);
            high = high.max(peak.log_fold_change);
        }
        if low.is_finite() && high.is_finite() {
            Some((low, high))
        } else {
            None
        }
    }

    /// `max(-log10(res$P.Value), na.rm = TRUE)` over **all** peaks.
    ///
    /// R computes this before subsetting, and ignores `adj.P.Val` and `logFC`
    /// entirely, so non-significant peaks still stretch the axis -- and so does a
    /// peak whose fold change is missing and therefore cannot be drawn at all.
    /// The `na.rm` only drops rows with a missing p-value.
    pub fn y_max(&self) -> f64 {
        self.peaks
            .iter()
            .map(Peak::neg_log10_p)
            .filter(|value| !value.is_nan())
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Counts of significant peaks on each side, for the legend labels.
    pub fn counts(&self, fdr: f64) -> (usize, usize) {
        let mut down = 0usize;
        let mut up = 0usize;
        for peak in self.drawable() {
            match peak.classify(fdr) {
                PeakClass::Up => up += 1,
                PeakClass::Down => down += 1,
                PeakClass::NotSignificant => {}
            }
        }
        (down, up)
    }

    /// Peaks that cannot be drawn because a coordinate is not finite.
    ///
    /// R does not have an equivalent count: it silently hands non-finite values
    /// to `points()` and, if any `logFC` is missing, fails on an `NA` x limit
    /// before drawing anything at all. Reporting the count lets a caller warn
    /// about the sparser plot instead of leaving the difference unexplained.
    pub fn skipped(&self) -> usize {
        self.peaks.len() - self.drawable().count()
    }
}

/// R's `pretty()` bounds for an axis, widened the way base graphics widens a
/// plot region.
///
/// Base R draws with `xaxs = "r"`, which pads the plot region by 4% on each
/// side **before** mapping data onto it. Reproducing the padding matters
/// because `grid()` and `axis()` draw against the padded region, not the raw
/// data range, so without it the gridlines and ticks land in the wrong place.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExpandedRange {
    /// Padded low end, i.e. R's `par("usr")[1]`.
    pub low: f64,
    /// Padded high end, i.e. R's `par("usr")[2]`.
    pub high: f64,
}

impl ExpandedRange {
    /// Applies the 4% padding R uses for `xaxs = "r"` / `yaxs = "r"`.
    pub fn padded(low: f64, high: f64) -> Self {
        let span = high - low;
        // A degenerate range would pad to the same value and make every
        // fraction a division by zero; R behaves the same way (a single point
        // gives `usr` of width 0.08 * the value), so fall back to a unit span.
        if span == 0.0 {
            return Self {
                low: low - 0.4,
                high: high + 0.4,
            };
        }
        Self {
            low: low - 0.04 * span.abs(),
            high: high + 0.04 * span.abs(),
        }
    }

    /// Maps a value in `low ..= high` onto `0.0 ..= 1.0`.
    pub fn fraction(&self, value: f64) -> f64 {
        let span = self.high - self.low;
        if span <= 0.0 {
            return 0.0;
        }
        (value - self.low) / span
    }

    /// Whether a value falls inside the padded region.
    ///
    /// R clips ticks to `usr`, so `pretty()` values outside it are not drawn.
    pub fn contains(&self, value: f64) -> bool {
        value >= self.low && value <= self.high
    }
}

/// Draws the volcano panel: grid, points, axes, labels and legend.
///
/// Layout follows `volcano_plot()`'s `par(mar = c(3.5, 3.5, 3, 1))`, which
/// leaves room for two lines of axis titles on the bottom/left and a title on
/// top.
#[allow(clippy::too_many_arguments)]
pub fn draw_volcano_panel(
    panel: &mut PanelWriter,
    volcano: &Volcano,
    fdr: f64,
    up_color: &str,
    down_color: &str,
    alpha: f64,
    point_size: f64,
    show_axis: bool,
    font_size: f64,
) {
    let Some((data_low, data_high)) = volcano.x_range() else {
        return;
    };
    let data_y_max = volcano.y_max();
    if !data_y_max.is_finite() || data_y_max <= 0.0 {
        // R would build an infinite `ylim` here and `plot()` would abort; the
        // caller reports the condition, so there is nothing to draw.
        return;
    }

    // R's margins, in lines. The extra line on the left and bottom carries the
    // axis titles, which sit two lines out from the ticks.
    let line_height = font_size * 1.2;
    let plot_left = 4.5 * line_height;
    let plot_right = panel.width() - 1.0 * line_height;
    let plot_top = 3.0 * line_height;
    let plot_bottom = panel.height() - 3.5 * line_height;
    if plot_right <= plot_left || plot_bottom <= plot_top {
        return;
    }

    // Both axes use R's 4% padding, so gridlines and ticks are positioned
    // against the same region the data is mapped into.
    let x_axis = ExpandedRange::padded(data_low, data_high);
    let y_axis = ExpandedRange::padded(0.0, data_y_max);

    let map_x = |value: f64| -> f64 {
        plot_left + x_axis.fraction(value) * (plot_right - plot_left)
    };
    let map_y = |value: f64| -> f64 {
        // y grows upward, so the fraction is measured from the bottom.
        plot_bottom - y_axis.fraction(value) * (plot_bottom - plot_top)
    };

    // Grid first, so the points sit on top of it. R's `grid()` draws at
    // `axTicks()`, which spans the padded region rather than the data range.
    let x_ticks = crate::plot::pretty::pretty_range(x_axis.low, x_axis.high, 5).values;
    let y_ticks = crate::plot::pretty::pretty_range(y_axis.low, y_axis.high, 5).values;
    let visible_x_ticks: Vec<f64> = x_ticks
        .iter()
        .copied()
        .filter(|value| x_axis.contains(*value))
        .collect();
    let visible_y_ticks: Vec<f64> = y_ticks
        .iter()
        .copied()
        .filter(|value| y_axis.contains(*value))
        .collect();

    for value in &visible_x_ticks {
        let x = map_x(*value);
        panel.line(x, plot_top, x, plot_bottom, "gray90", 1.0);
    }
    for value in &visible_y_ticks {
        let y = map_y(*value);
        panel.line(plot_left, y, plot_right, y, "gray90", 1.0);
    }

    // Points. R draws non-significant ones first, then the significant ones, so
    // a significant peak is never hidden behind a grey one.
    let radius = (2.6 * point_size).max(1.0);
    for peak in volcano.drawable() {
        if peak.classify(fdr) != PeakClass::NotSignificant {
            continue;
        }
        panel.circle_with_opacity(
            map_x(peak.log_fold_change),
            map_y(peak.neg_log10_p()),
            radius,
            NON_SIGNIFICANT_COLOR,
            alpha,
        );
    }
    for peak in volcano.drawable() {
        let color = match peak.classify(fdr) {
            PeakClass::Up => up_color,
            PeakClass::Down => down_color,
            PeakClass::NotSignificant => continue,
        };
        panel.circle_with_opacity(
            map_x(peak.log_fold_change),
            map_y(peak.neg_log10_p()),
            radius,
            color,
            alpha,
        );
    }

    if show_axis {
        for value in &visible_x_ticks {
            let x = map_x(*value);
            panel.line(x, plot_bottom, x, plot_bottom + 4.0, "black", 1.0);
            panel.text(
                x,
                plot_bottom + 6.0 + font_size,
                &crate::plot::svg::format_tick(*value),
                font_size,
                "middle",
                "black",
            );
        }
        for value in &visible_y_ticks {
            let y = map_y(*value);
            panel.line(plot_left - 4.0, y, plot_left, y, "black", 1.0);
            panel.text(
                plot_left - 6.0,
                y + font_size * 0.35,
                &crate::plot::svg::format_tick(*value),
                font_size,
                "end",
                "black",
            );
        }
    }

    // Axis titles, matching `mtext(..., line = 2, col = "#34495e", font = 2)`.
    panel.text_weighted(
        (plot_left + plot_right) / 2.0,
        panel.height() - font_size * 0.6,
        "Log Fold Change",
        font_size,
        "middle",
        LABEL_COLOR,
        "bold",
    );
    panel.text_weighted(
        font_size * 0.9,
        (plot_top + plot_bottom) / 2.0,
        "-log10(P-Value)",
        font_size,
        "middle",
        LABEL_COLOR,
        "bold",
    );

    // Title: the contrast label.
    if !volcano.title.is_empty() {
        panel.text_weighted(
            (plot_left + plot_right) / 2.0,
            plot_top - font_size * 0.8,
            &volcano.title,
            font_size * 1.2,
            "middle",
            LABEL_COLOR,
            "bold",
        );
    }

    // Legend, bottom-left or bottom-right depending on which side reaches
    // further from zero.
    let (down_count, up_count) = volcano.counts(fdr);
    draw_legend(
        panel,
        &x_axis,
        plot_left,
        plot_right,
        plot_bottom,
        plot_top,
        &[
            (format!("Down [{down_count}]"), down_color.to_string()),
            (format!("Up [{up_count}]"), up_color.to_string()),
        ],
        &format!("Diff. peaks [FDR<{fdr}]"),
        font_size,
    );
}

/// R's `legend()` position rule.
///
/// `volcano_plot()` picks the corner with `which(abs(xlims) == max(abs(xlims)))`,
/// i.e. the side whose extreme is farther from zero: index 1 means the low end
/// wins and the legend goes bottom-left.
///
/// That expression returns *two* indices when the range is symmetric, and R's
/// `legend()` then aborts with "'arg' must be of length 1". A symmetric volcano
/// plot is entirely ordinary, so the port resolves the tie instead of failing;
/// `<=` gives the bottom-left, matching what the `ifelse` would return if only
/// the first match were used.
pub fn legend_position(x_low: f64, x_high: f64) -> LegendPosition {
    if x_low.abs() <= x_high.abs() {
        LegendPosition::BottomRight
    } else {
        LegendPosition::BottomLeft
    }
}

/// Which bottom corner the legend occupies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LegendPosition {
    BottomLeft,
    BottomRight,
}

#[allow(clippy::too_many_arguments)]
fn draw_legend(
    panel: &mut PanelWriter,
    x_axis: &ExpandedRange,
    plot_left: f64,
    plot_right: f64,
    plot_bottom: f64,
    plot_top: f64,
    entries: &[(String, String)],
    title: &str,
    font_size: f64,
) {
    if entries.is_empty() {
        return;
    }
    let line_height = font_size + 2.0;
    // The legend sits inside the plot region, as R's default `xpd = FALSE`
    // would clip it to; anchoring in the bottom corner keeps it clear of the
    // data, which for a volcano plot is densest along the top.
    let text_height = line_height * (entries.len() as f64 + if title.is_empty() { 0.0 } else { 1.0 });
    let top = plot_bottom - 6.0 - text_height;

    let (anchor_x, swatch_dx, text_anchor) = match legend_position(x_axis.low, x_axis.high) {
        LegendPosition::BottomLeft => (plot_left + 6.0, 8.0, "start"),
        LegendPosition::BottomRight => (plot_right - 6.0, -8.0, "end"),
    };
    let _ = plot_top;

    let mut y = top;
    if !title.is_empty() {
        panel.text_weighted(
            anchor_x,
            y,
            title,
            font_size,
            text_anchor,
            "black",
            "bold",
        );
        y += line_height;
    }
    for (label, color) in entries {
        let swatch_x = anchor_x + swatch_dx;
        panel.circle(swatch_x, y - font_size * 0.3, 2.6, color);
        panel.text(anchor_x, y, label, font_size, text_anchor, "black");
        y += line_height;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peak(log_fold_change: f64, p_value: f64, adjusted_p_value: f64) -> Peak {
        Peak {
            log_fold_change,
            p_value,
            adjusted_p_value,
        }
    }

    #[test]
    fn classification_matches_the_r_rule() {
        let fdr = 0.1;
        assert_eq!(peak(-2.0, 0.001, 0.01).classify(fdr), PeakClass::Down);
        assert_eq!(peak(3.0, 0.001, 0.01).classify(fdr), PeakClass::Up);
        // Exactly at the threshold is NOT significant, because R uses `<`.
        assert_eq!(
            peak(3.0, 0.001, 0.1).classify(fdr),
            PeakClass::NotSignificant
        );
        assert_eq!(
            peak(3.0, 0.001, 0.1000001).classify(fdr),
            PeakClass::NotSignificant
        );
    }

    #[test]
    fn a_missing_adjusted_p_value_is_not_significant() {
        // `NA < fdr` is NA, so data.table's `[adj.P.Val < fdr]` drops the row and
        // `[!adj.P.Val < fdr]` drops it too; it ends up in neither group, which
        // means only the non-significant points would ever draw it. Treating it
        // as non-significant keeps it visible, matching `points(res_nonsig, ...)`
        // being fed the rows R failed to classify.
        assert_eq!(
            peak(2.0, 0.001, f64::NAN).classify(0.1),
            PeakClass::NotSignificant
        );
    }

    #[test]
    fn a_significant_peak_with_no_direction_counts_on_neither_side() {
        // R's legend counts are `nrow(res_sig[logFC < 0])` and `[logFC > 0]`, so
        // a zero or missing logFC lands in neither bucket.
        let volcano = Volcano::new(
            vec![
                peak(0.0, 0.001, 0.01),
                peak(-1.0, 0.001, 0.01),
                peak(1.0, 0.001, 0.01),
            ],
            String::new(),
        );
        assert_eq!(volcano.counts(0.1), (1, 1));
    }

    #[test]
    fn missing_fold_change_is_counted_on_neither_side_even_when_significant() {
        let volcano = Volcano::new(
            vec![peak(f64::NAN, 0.001, 0.01), peak(1.0, 0.001, 0.01)],
            String::new(),
        );
        assert_eq!(volcano.counts(0.1), (0, 1));
    }

    #[test]
    fn non_finite_coordinates_are_not_drawn() {
        // R passes NA straight to points() and the marker silently vanishes, so
        // these must not be mapped onto the axes either.
        let nan_log_fc = peak(f64::NAN, 0.001, 0.01);
        let nan_p = peak(1.0, f64::NAN, 0.01);
        assert!(!nan_log_fc.is_drawable());
        assert!(!nan_p.is_drawable());

        let volcano = Volcano::new(vec![nan_log_fc, nan_p, peak(1.0, 0.01, 0.001)], String::new());
        assert_eq!(volcano.drawable().count(), 1);
        assert_eq!(volcano.x_range(), Some((1.0, 1.0)));
    }

    #[test]
    fn a_zero_p_value_gives_an_infinite_y_axis() {
        // limma can emit p = 0; R's ylims then becomes Inf and plot() aborts.
        // The port must expose that so the caller can report it rather than
        // dividing by infinity while placing points.
        let volcano = Volcano::new(vec![peak(1.0, 0.0, 0.0)], String::new());
        assert!(!volcano.y_max().is_finite());
    }

    #[test]
    fn y_max_ignores_significance_and_uses_every_peak() {
        // R computes `max(-log10(P.Value))` before subsetting, so a
        // non-significant peak with a tiny p-value still stretches the axis.
        let volcano = Volcano::new(
            vec![peak(1.0, 1e-8, 0.9), peak(2.0, 1e-3, 1e-4)],
            String::new(),
        );
        assert!((volcano.y_max() - 8.0).abs() < 1e-12);
    }

    #[test]
    fn y_max_includes_peaks_whose_fold_change_is_missing() {
        // R computes ylims from `res$P.Value` alone, before looking at logFC, so
        // a row with a tiny p-value but a missing fold change still stretches the
        // y axis even though it can never be drawn (`range(logFC)` is NA and the
        // point is silently skipped). Verified against R: with p-values
        // c(1e-10, 0.5, 0.4) and an NA logFC on the first row, R's ylims is 10,
        // not the 0.39794 that dropping the row would give.
        let volcano = Volcano::new(
            vec![peak(f64::NAN, 1e-10, 0.01), peak(0.5, 0.5, 0.9)],
            String::new(),
        );
        assert_eq!(volcano.y_max(), 10.0);
        // Only one of the two is drawable, so the axis stretches while the
        // scatter stays sparse -- exactly what R produces.
        assert_eq!(volcano.drawable().count(), 1);
    }

    #[test]
    fn y_max_drops_only_missing_p_values() {
        // `na.rm = TRUE` drops rows whose p-value is NA, and nothing else.
        let volcano = Volcano::new(
            vec![peak(1.0, f64::NAN, 0.01), peak(2.0, 1e-4, 0.01)],
            String::new(),
        );
        assert!((volcano.y_max() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn x_range_covers_both_ends_of_the_data() {
        let volcano = Volcano::new(
            vec![peak(-4.2, 0.1, 0.9), peak(5.5, 0.1, 0.9), peak(0.2, 0.1, 0.9)],
            String::new(),
        );
        assert_eq!(volcano.x_range(), Some((-4.2, 5.5)));
    }

    #[test]
    fn an_undrawable_table_has_no_range() {
        let volcano = Volcano::new(vec![peak(f64::NAN, 1.0, 1.0)], String::new());
        assert_eq!(volcano.x_range(), None);
    }

    #[test]
    fn padding_matches_r_padded_us_region() {
        // Captured from R for `plot(xlim = c(-4.2, 5.5))`: usr is -4.588..5.888.
        let x_axis = ExpandedRange::padded(-4.2, 5.5);
        assert!((x_axis.low - -4.588).abs() < 1e-12, "got {}", x_axis.low);
        assert!((x_axis.high - 5.888).abs() < 1e-12, "got {}", x_axis.high);

        // And for `ylim = c(0, 8)`: usr is -0.32..8.32.
        let y_axis = ExpandedRange::padded(0.0, 8.0);
        assert!((y_axis.low - -0.32).abs() < 1e-12, "got {}", y_axis.low);
        assert!((y_axis.high - 8.32).abs() < 1e-12, "got {}", y_axis.high);
    }

    #[test]
    fn padded_region_clips_ticks_like_r_does() {
        // R's `pretty(c(-4.2, 5.5))` spans -6..6, but usr ends at 5.888, so the
        // -6 and 6 ticks are clipped away and only -4..4 remain visible.
        let x_axis = ExpandedRange::padded(-4.2, 5.5);
        let visible: Vec<f64> = crate::plot::pretty::pretty_range(x_axis.low, x_axis.high, 5)
            .values
            .into_iter()
            .filter(|value| x_axis.contains(*value))
            .collect();
        assert_eq!(visible, vec![-4.0, -2.0, 0.0, 2.0, 4.0]);
    }

    #[test]
    fn a_degenerate_range_still_maps_within_bounds() {
        // A single peak at x = 0 has a zero-width range; R pads it to a finite
        // region and draws the point in the middle.
        let x_axis = ExpandedRange::padded(0.0, 0.0);
        let fraction = x_axis.fraction(0.0);
        assert!((0.0..=1.0).contains(&fraction), "got {fraction}");
    }

    #[test]
    fn legend_moves_to_the_side_farther_from_zero() {
        // R: `which(abs(xlims) == max(abs(xlims)))`, index 1 -> bottomleft.
        assert_eq!(legend_position(-5.0, 2.0), LegendPosition::BottomLeft);
        assert_eq!(legend_position(-2.0, 3.0), LegendPosition::BottomRight);
        assert_eq!(legend_position(-4.2, 5.5), LegendPosition::BottomRight);
    }

    #[test]
    fn a_symmetric_range_resolves_the_tie_instead_of_failing() {
        // R's `which()` returns both indices here and `legend()` aborts with
        // "'arg' must be of length 1"; the port picks a corner so an ordinary
        // symmetric volcano plot still renders.
        assert_eq!(legend_position(-3.0, 3.0), LegendPosition::BottomRight);
    }

    // ---------------------------------------------------------------------
    // R oracle
    //
    // `testdata/volcano_r_oracle.tsv` was produced by R 4.6.0 from the same
    // expressions `volcano_plot()` uses, including whether R's own `legend()`
    // call survives the corner it picked. Regenerate with the generator script
    // described in the testdata header when it needs extending.
    // ---------------------------------------------------------------------

    /// Splits a comma-separated field, mapping `NA` to NaN.
    fn parse_field(field: &str) -> Vec<f64> {
        if field.is_empty() {
            return Vec::new();
        }
        field
            .split(',')
            .map(|token| {
                let token = token.trim();
                if token == "NA" || token.eq_ignore_ascii_case("nan") {
                    f64::NAN
                } else {
                    token.parse::<f64>().unwrap_or(f64::NAN)
                }
            })
            .collect()
    }

    fn oracle_close(actual: f64, expected: f64) -> bool {
        if actual == expected {
            return true;
        }
        if actual.is_nan() || expected.is_nan() {
            return actual.is_nan() && expected.is_nan();
        }
        if !actual.is_finite() || !expected.is_finite() {
            return actual == expected;
        }
        let scale = actual.abs().max(expected.abs()).max(1.0);
        (actual - expected).abs() <= 1e-12 * scale
    }

    #[test]
    fn classification_and_limits_match_the_r_volcano_oracle() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("testdata")
            .join("volcano_r_oracle.tsv");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("read oracle {path:?}: {error}"));

        const FDR: f64 = 0.1;
        let mut checked = 0usize;
        for line in text.lines() {
            let line = line.trim_end();
            if line.is_empty() || line.starts_with('#') || line.starts_with("name\t") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 12, "malformed oracle row: {line:?}");

            let name = fields[0];
            let log_fc = parse_field(fields[2]);
            let p_values = parse_field(fields[3]);
            let padj = parse_field(fields[4]);
            let expected_x_min: f64 = fields[5].parse().unwrap_or(f64::NAN);
            let expected_x_max: f64 = fields[6].parse().unwrap_or(f64::NAN);
            let expected_y_max = fields[7].parse::<f64>().unwrap_or(f64::NAN);
            let expected_down: usize = fields[8].parse().expect("down count");
            let expected_up: usize = fields[9].parse().expect("up count");
            let expected_legend = fields[10].to_string();
            let expected_legend_error = fields[11] == "1";

            assert_eq!(log_fc.len(), p_values.len(), "{name}: column lengths");
            assert_eq!(log_fc.len(), padj.len(), "{name}: column lengths");

            let peaks: Vec<Peak> = (0..log_fc.len())
                .map(|index| Peak {
                    log_fold_change: log_fc[index],
                    p_value: p_values[index],
                    adjusted_p_value: padj[index],
                })
                .collect();
            let volcano = Volcano::new(peaks, String::new());

            // y_max: R computes `max(-log10(P.Value))` before any subsetting.
            let y_max = volcano.y_max();
            assert!(
                oracle_close(y_max, expected_y_max),
                "{name}: y_max got {y_max}, R has {expected_y_max}"
            );

            // When R's ylims is infinite, `plot()` aborts before it ever uses the
            // x range, so R draws nothing at all for this table and there is no
            // axis geometry to compare against. Both implementations refuse the
            // input (the CLI errors with the same reason R's error gives), so
            // only y_max above is meaningful here.
            if !expected_y_max.is_finite() {
                assert!(
                    volcano.x_range().is_some(),
                    "{name}: R failed on an infinite ylim, so x_range is unconstrained"
                );
                checked += 1;
                continue;
            }

            // x_range: R's `range(logFC)`, which is NA as soon as *any* fold
            // change is missing -- so the whole figure fails to draw, even when
            // other points were fine. The port instead plots the usable points
            // and reports the rest via `skipped()`, which is a deliberate
            // improvement; the assertion below checks both halves of that.
            match volcano.x_range() {
                Some((low, high)) => {
                    if expected_x_min.is_nan() {
                        // R would have drawn nothing here. The port draws the
                        // finite subset, so the range must come from those.
                        assert!(
                            volcano.skipped() > 0,
                            "{name}: R had no range but the port skipped nothing"
                        );
                        let finite: Vec<f64> = log_fc
                            .iter()
                            .copied()
                            .filter(|value| value.is_finite())
                            .collect();
                        let expected_low = finite.iter().copied().fold(f64::INFINITY, f64::min);
                        let expected_high =
                            finite.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                        assert!(
                            oracle_close(low, expected_low) && oracle_close(high, expected_high),
                            "{name}: finite-only x_range got ({low}, {high}), expected ({expected_low}, {expected_high})"
                        );
                    } else {
                        assert!(
                            oracle_close(low, expected_x_min)
                                && oracle_close(high, expected_x_max),
                            "{name}: x_range got ({low}, {high}), R has ({expected_x_min}, {expected_x_max})"
                        );
                    }
                }
                None => assert!(
                    expected_x_min.is_nan(),
                    "{name}: x_range was None but R has ({expected_x_min}, {expected_x_max})"
                ),
            }

            // Legend counts, and the corner R's rule picks.
            let (down, up) = volcano.counts(FDR);
            assert_eq!(
                (down, up),
                (expected_down, expected_up),
                "{name}: counts differ from R"
            );

            // R's `which(abs(xlims) == max(abs(xlims)))` has three outcomes:
            //   - one match: a usable corner, which the port must agree with;
            //   - two matches (symmetric range): `legend()` aborts, so the port
            //     only has to pick a corner at all;
            //   - no matches (xlims is NA): there is nothing to compare, and
            //     x_range() is None so no plot is produced either way.
            let r_positions: Vec<&str> = expected_legend
                .split('|')
                .filter(|position| !position.is_empty())
                .collect();
            match r_positions.len() {
                1 => {
                    let expected = if r_positions[0] == "bottomleft" {
                        LegendPosition::BottomLeft
                    } else {
                        LegendPosition::BottomRight
                    };
                    assert_eq!(
                        legend_position(expected_x_min, expected_x_max),
                        expected,
                        "{name}: legend corner differs from R"
                    );
                }
                2 => {
                    // R errored on the tie; assert that, so the case stays
                    // documenting the bug it exists for.
                    assert!(
                        expected_legend_error,
                        "{name}: R reported two corners but no error"
                    );
                    let chosen = legend_position(expected_x_min, expected_x_max);
                    assert!(
                        chosen == LegendPosition::BottomRight
                            || chosen == LegendPosition::BottomLeft,
                        "{name}: port chose no corner"
                    );
                }
                _ => {
                    // R had no usable x range (an NA logFC), so its `which()`
                    // matched nothing and the whole plot failed. The port draws
                    // the finite subset instead, so it must have a range and a
                    // corner to put the legend in.
                    assert!(
                        expected_legend_error,
                        "{name}: no legend position but R reported no error"
                    );
                    if let Some((low, high)) = volcano.x_range() {
                        let chosen = legend_position(low, high);
                        assert!(
                            chosen == LegendPosition::BottomRight
                                || chosen == LegendPosition::BottomLeft,
                            "{name}: port chose no corner"
                        );
                    } else {
                        // Nothing drawable at all, so there is genuinely no plot.
                        assert_eq!(
                            volcano.skipped(),
                            volcano.peaks.len(),
                            "{name}: no range but some peaks were drawable"
                        );
                    }
                }
            }

            checked += 1;
        }

        assert!(checked >= 8, "oracle only had {checked} cases");
        eprintln!("validated {checked} volcano cases against the R oracle");
    }
}
