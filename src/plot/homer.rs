//! Summaries of HOMER `annotatePeaks.pl` output.
//!
//! Ported from `trackplot.R`'s `summarize_homer_annots()`, which reads one
//! annotation file per sample, counts peaks per genomic feature, and draws one
//! horizontal stacked bar per sample.
//!
//! The interesting part is not the counting but the *filtering*. R computes each
//! category's fraction over **all** annotations:
//!
//! ```r
//! h[, .N, Annotation][, fract := N / sum(N)]
//! ```
//!
//! but then draws only the categories that appear in a fixed nine-entry palette:
//!
//! ```r
//! homer.anno.stats[names(pie.cols)[names(pie.cols) %in% rownames(homer.anno.stats)], ]
//! ```
//!
//! So a category outside that palette is dropped from the figure while still
//! counting toward the denominator, and the drawn segments of a bar sum to less
//! than 1. That is reproduced here ([`Summary::column_sum`] exposes the total so
//! a caller can see it), because a drop-in replacement has to produce the same
//! picture.
//!
//! One case deserves naming, because it looks like a bug and behaves like one.
//! The palette has an entry `'NA' = 'gray70'`, but `fread()` reads the literal
//! text `NA` that HOMER writes for unannotated peaks as a *missing value*, and
//! `%in%` never matches `NA`. The palette entry is therefore unreachable and
//! those peaks are always dropped.
//!
//! The port keeps those peaks **counted** (so they show up under
//! [`DroppedCategory`] and the loss is visible) and lets
//! [`build_summary`]'s `draw_unannotated` argument make the palette entry
//! reachable, which is what `--keep-unannotated` turns on.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use anyhow::{anyhow, Context, Result};

use crate::plot::svg::PanelWriter;

/// The fixed category palette, in draw order, from `summarize_homer_annots()`.
///
/// The order is load-bearing twice over: it decides the order of the stacked
/// segments *and* the order of the legend, both of which come from indexing
/// this vector rather than from the data.
pub const ANNOTATION_PALETTE: [(&str, &str); 9] = [
    ("3pUTR", "#E7298A"),
    ("5pUTR", "#D95F02"),
    ("Intergenic", "#BEBADA"),
    ("TTS", "#FB8072"),
    ("exon", "#80B1D3"),
    ("intron", "#FDB462"),
    ("non-coding", "#FFFFB3"),
    ("NA", "gray70"),
    ("promoter-TSS", "#1B9E77"),
];

/// The category HOMER writes for a peak it could not annotate.
///
/// R's palette names a colour for this, but the lookup can never reach it (see
/// the module docs), so by default the category is counted and then dropped.
pub const UNMAPPED_NA_CATEGORY: &str = "NA";

/// Whether a category's palette entry is reachable, given how the `NA` category
/// is treated.
///
/// R's palette lookup is `names(pie.cols)[names(pie.cols) %in% rownames(x)]`, and
/// `%in%` never matches `NA`, so the `'NA'` entry is unreachable there.
/// `draw_unannotated` models the caller asking for that entry to work.
pub fn palette_reachable(name: &str, draw_unannotated: bool) -> bool {
    if name == UNMAPPED_NA_CATEGORY && !draw_unannotated {
        return false;
    }
    palette_color(name).is_some()
}

/// Looks up a palette colour by category name.
///
/// This is the raw table lookup; use [`palette_reachable`] to ask whether the
/// entry can actually be matched for a given input.
pub fn palette_color(name: &str) -> Option<&'static str> {
    ANNOTATION_PALETTE
        .iter()
        .find(|(category, _)| *category == name)
        .map(|(_, color)| *color)
}

/// One sample's parsed annotation file.
#[derive(Clone, Debug, PartialEq)]
pub struct Sample {
    /// Sample name, from `--sample` or the file's first dot-separated field.
    pub name: String,
    /// Total peaks read, i.e. R's `npeaks` and the `sum(N)` denominator.
    pub n_peaks: usize,
    /// Peak count per normalized category, including
    /// [`UNMAPPED_NA_CATEGORY`] when the file has unannotated peaks.
    pub counts: BTreeMap<String, usize>,
}

impl Sample {
    /// `fract` for one category, i.e. `N / sum(N)`.
    ///
    /// Returns 0 for a category the sample has none of, matching the `fill = 0`
    /// that `dcast()` applies when a category is missing from one sample.
    pub fn fraction(&self, category: &str) -> f64 {
        if self.n_peaks == 0 {
            return 0.0;
        }
        self.counts.get(category).copied().unwrap_or(0) as f64 / self.n_peaks as f64
    }

    /// Peaks HOMER could not annotate.
    pub fn unannotated(&self) -> usize {
        self.counts
            .get(UNMAPPED_NA_CATEGORY)
            .copied()
            .unwrap_or(0)
    }
}

/// A category that exists in the data but is not drawn.
#[derive(Clone, Debug, PartialEq)]
pub struct DroppedCategory {
    pub name: String,
    /// Peak count per sample, in `Summary::samples` order.
    pub counts: Vec<usize>,
}

impl DroppedCategory {
    /// Peak count summed across samples.
    pub fn total(&self) -> usize {
        self.counts.iter().sum()
    }

    /// Whether this is the `NA` category, whose drop follows from R's `fread()`
    /// coercion rather than from a decision about that category, so a caller can
    /// explain it specifically.
    pub fn is_unmapped_na(&self) -> bool {
        self.name == UNMAPPED_NA_CATEGORY
    }
}

/// Normalizes one annotation cell the way the R code does.
///
/// HOMER writes the feature followed by the nearest feature in parentheses, e.g.
/// `promoter-TSS (NM_000001)`, and R keeps everything before the first `" ("`
/// via `strsplit(..., split = ' (', fixed = TRUE)[[1]]`. Categories such as
/// `Intergenic` have no suffix and pass through untouched.
///
/// HOMER's literal `NA` text is returned as [`UNMAPPED_NA_CATEGORY`] rather than
/// treated as missing, because R's grouping keeps it as a category before the
/// palette filter removes it; keeping it here is what lets the port report it as
/// dropped. An empty cell is genuinely missing and yields `None`.
pub fn normalize_annotation(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }

    let feature = match trimmed.find(" (") {
        Some(index) => &trimmed[..index],
        None => trimmed,
    };
    Some(rename_utr(feature))
}

/// Applies R's two `gsub()` renames.
///
/// The palette keys are `3pUTR`/`5pUTR` while HOMER writes `3' UTR`/`5' UTR`, so
/// these two substitutions are what make the UTR categories plottable at all.
/// `gsub()` replaces every occurrence, but these labels contain at most one, so a
/// plain replacement is equivalent.
fn rename_utr(feature: &str) -> String {
    feature.replace("3' UTR", "3pUTR").replace("5' UTR", "5pUTR")
}

/// Finds the index of the `Annotation` column in a header row.
///
/// R selects columns by name, so the header is required; there is no positional
/// fallback that would silently read the wrong column.
fn annotation_index(header: &[&str]) -> Result<usize> {
    header
        .iter()
        .position(|field| field.trim().eq_ignore_ascii_case("Annotation"))
        .ok_or_else(|| {
            anyhow!(
                "no `Annotation` column in the header; found: {}. \
                 This does not look like `annotatePeaks.pl` output.",
                header.join(", ")
            )
        })
}

/// Parses one `annotatePeaks.pl` output file into a [`Sample`].
///
/// `name` is the sample label; when empty it is derived from the file name the
/// way R derives it, by taking everything before the first dot.
pub fn parse_sample(path: &Path, name: &str) -> Result<Sample> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("read HOMER annotation file: {path:?}"))?;

    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let header_line = lines
        .next()
        .ok_or_else(|| anyhow!("HOMER annotation file is empty: {path:?}"))?;
    let header: Vec<&str> = header_line.split('\t').collect();
    let annotation_column = annotation_index(&header)
        .with_context(|| format!("in {path:?}"))?;

    let sample_name = if name.is_empty() {
        derive_sample_name(path)
    } else {
        name.to_string()
    };

    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut n_peaks = 0usize;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        n_peaks += 1;
        // A row shorter than the header cannot carry the annotation; R would
        // have failed reading the file at all, so treat it as unannotated rather
        // than guessing.
        let raw = fields.get(annotation_column).copied().unwrap_or("");
        if let Some(category) = normalize_annotation(raw) {
            *counts.entry(category).or_insert(0) += 1;
        }
    }

    if n_peaks == 0 {
        return Err(anyhow!(
            "HOMER annotation file has a header but no peaks: {path:?}"
        ));
    }

    Ok(Sample {
        name: sample_name,
        n_peaks,
        counts,
    })
}

/// Derives a sample name from a file path the way R's
/// `tstrsplit(basename(anno), "\\.", keep = 1)` does.
fn derive_sample_name(path: &Path) -> String {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name| name.split('.').next().unwrap_or(name).to_string())
        .unwrap_or_else(|| "sample".to_string())
}

/// A matrix of fractions with one row per drawn category and one column per
/// sample, ready to render.
#[derive(Clone, Debug, PartialEq)]
pub struct Summary {
    /// Category names in palette order; only categories present in the palette.
    pub categories: Vec<String>,
    /// Sample names in the order they were supplied.
    pub samples: Vec<String>,
    /// Total peaks per sample, for the right-hand labels.
    pub n_peaks: Vec<usize>,
    /// `fractions[category][sample]`, indexed like `categories` and `samples`.
    pub fractions: Vec<Vec<f64>>,
    /// Categories dropped because the palette has no entry, with their counts
    /// per sample. Empty when nothing was dropped.
    pub dropped: Vec<DroppedCategory>,
}

impl Summary {
    /// Fraction of a sample's bar that is actually drawn.
    ///
    /// R's bars sum to less than 1 whenever a category was dropped; this reports
    /// the drawn total so a caller can warn about it.
    pub fn column_sum(&self, sample_index: usize) -> f64 {
        self.fractions
            .iter()
            .map(|row| row.get(sample_index).copied().unwrap_or(0.0))
            .sum()
    }

    /// The colour for a drawn category, which is always available because
    /// categories come from the palette.
    pub fn category_color(&self, category_index: usize) -> &'static str {
        self.categories
            .get(category_index)
            .and_then(|name| palette_color(name))
            .unwrap_or("gray70")
    }

    /// The `Annotation [N]` legend text R computes into its `leg` column.
    ///
    /// `summarize_homer_annots()` only uses these in its commented-out pie-chart
    /// branch; the drawn bar chart labels the legend with bare category names.
    /// They are exposed here because they are the per-sample counts a caller
    /// usually wants, and recomputing them from `fractions` would introduce
    /// rounding.
    pub fn legend_labels(&self, counts: &[Sample]) -> Vec<Vec<String>> {
        self.categories
            .iter()
            .map(|category| {
                counts
                    .iter()
                    .map(|sample| {
                        let count = sample.counts.get(category).copied().unwrap_or(0);
                        format!("{category} [{count}]")
                    })
                    .collect()
            })
            .collect()
    }
}

/// Builds the plot matrix from parsed samples.
///
/// This mirrors R's chain of `rbindlist()` → `dcast(Annotation ~ Sample)` →
/// row-filter on the palette, including three details that are easy to miss:
///
/// - Category order comes from the **palette**, not from the data, so the rows
///   are always in the same order regardless of which sample is largest.
/// - A category missing from a sample contributes 0, matching `dcast(fill = 0)`.
/// - The `NA` entry is dropped unless `draw_unannotated` is set, because R's
///   `%in%` lookup can never match it (see [`palette_reachable`]).
///
/// `draw_unannotated` is spelled out rather than defaulted so that no call site
/// gets R's surprising behaviour by accident.
pub fn build_summary(samples: &[Sample], draw_unannotated: bool) -> Summary {
    let categories: Vec<String> = ANNOTATION_PALETTE
        .iter()
        .map(|(name, _)| (*name).to_string())
        .filter(|name| {
            // R keeps a palette entry only when some sample's row names contain
            // it, i.e. when at least one sample actually has that category...
            samples.iter().any(|sample| sample.counts.contains_key(name))
                // ...and the lookup can actually match it.
                && palette_reachable(name, draw_unannotated)
        })
        .collect();

    let fractions: Vec<Vec<f64>> = categories
        .iter()
        .map(|category| {
            samples
                .iter()
                .map(|sample| sample.fraction(category))
                .collect()
        })
        .collect();

    // Categories present in the data but not drawn. Collected in first-seen
    // order, which matches the row order `dcast()` would produce, and then
    // sorted by name so the report is stable across input orderings.
    let mut dropped_names: Vec<String> = Vec::new();
    for sample in samples {
        for category in sample.counts.keys() {
            if !categories.contains(category) && !dropped_names.contains(category) {
                dropped_names.push(category.clone());
            }
        }
    }
    dropped_names.sort();

    let dropped = dropped_names
        .into_iter()
        .map(|name| DroppedCategory {
            counts: samples
                .iter()
                .map(|sample| sample.counts.get(&name).copied().unwrap_or(0))
                .collect(),
            name,
        })
        .collect();

    Summary {
        categories,
        samples: samples.iter().map(|s| s.name.clone()).collect(),
        n_peaks: samples.iter().map(|s| s.n_peaks).collect(),
        fractions,
        dropped,
    }
}

/// Draws the horizontal stacked bars: one bar per sample, one segment per
/// category, with the tick axis below and the legend above.
///
/// Ported from `summarize_homer_annots()`'s `barplot()` call. R's margins are
/// `par(mar = c(2, 4, 5, 3))`: two lines below for the fraction axis, four to the
/// left for the sample names, five above for the three-column legend, and three
/// to the right for the peak counts.
#[allow(clippy::too_many_arguments)]
pub fn draw_homer_panel(
    panel: &mut PanelWriter,
    summary: &Summary,
    legend_font_size: f64,
    show_axis: bool,
    font_size: f64,
) {
    if summary.samples.is_empty() || summary.categories.is_empty() {
        return;
    }

    let line_height = font_size * 1.2;
    // R's `mar` sizes the left margin in lines and measures the labels itself.
    // The SVG writer has no text metrics, so the width is estimated from the
    // character count and the font size. The estimate is deliberately generous:
    // the labels are right-anchored, so an underestimated margin would push the
    // leading characters off the canvas entirely, whereas an overestimate only
    // wastes a little space. "H3K27ac" and "H3K4me3" are the same length but not
    // the same width in a proportional font, which is why a per-character budget
    // alone is not enough.
    let longest_name = summary
        .samples
        .iter()
        .map(|name| name.chars().count())
        .max()
        .unwrap_or(1);
    // ~0.62 em per character is a comfortable upper bound for the sans-serif
    // stack the writer emits, plus a gap so the text does not touch the bars.
    let label_width = longest_name as f64 * font_size * 0.62 + font_size;
    let plot_left = (4.0 * line_height).max(label_width);
    let plot_right = panel.width() - 3.0 * line_height;
    let plot_top = 5.0 * line_height;
    let plot_bottom = panel.height() - 2.0 * line_height;
    if plot_right <= plot_left || plot_bottom <= plot_top {
        return;
    }

    // `barplot()` uses the default `xaxs = "r"`, so the fraction axis is padded
    // by 4% just like a scatter plot's. The bars are drawn from 0, and the
    // padding is why the axis extends slightly past 1.
    let x_axis = crate::plot::svg::ExpandedRange::padded(0.0, 1.0);
    let map_x = |value: f64| -> f64 {
        plot_left + x_axis.fraction(value) * (plot_right - plot_left)
    };

    // R's `barplot` places `n` bars of width 1 with a 0.2 gap, then centres the
    // group on the device. Reproducing the spacing keeps the bars from touching
    // and matches the returned positions used for the labels.
    let n_samples = summary.samples.len();
    let plot_height = plot_bottom - plot_top;
    // `barplot` reserves one slot per bar plus one gap, all scaled to fit.
    let slot = plot_height / (n_samples as f64 * 1.2);
    let bar_height = slot;

    let bar_top = |index: usize| -> f64 {
        // Bars run top-to-bottom in sample order, each centred in its slot.
        plot_top + slot * (index as f64) + (slot - bar_height) / 2.0 + bar_height * 0.1
    };

    for (sample_index, _) in summary.samples.iter().enumerate() {
        let top = bar_top(sample_index);
        let height = bar_height * 0.8;
        let mut cursor = 0.0f64;
        for (category_index, _) in summary.categories.iter().enumerate() {
            let fraction = summary.fractions[category_index][sample_index];
            // A zero-width rectangle is invisible, so zeros are skipped rather
            // than emitted as degenerate elements. Testing `is_finite` first also
            // drops a NaN fraction, which a bare `> 0.0` comparison would let
            // through.
            if !fraction.is_finite() || fraction <= 0.0 {
                continue;
            }
            let x_start = map_x(cursor);
            let x_end = map_x(cursor + fraction);
            panel.rect_stroked(
                x_start,
                top,
                (x_end - x_start).max(0.0),
                height,
                summary.category_color(category_index),
                "#34495e",
                1.0,
            );
            cursor += fraction;
        }
    }

    if show_axis {
        // R: `axis(side = 1, at = seq(0, 1, 0.25), font = 2, lwd = 2)`.
        for step in 0..=4 {
            let value = step as f64 * 0.25;
            let x = map_x(value);
            panel.line(x, plot_bottom, x, plot_bottom - 4.0, "black", 2.0);
            panel.text_weighted(
                x,
                plot_bottom + 6.0 + font_size,
                &crate::plot::svg::format_tick(value),
                font_size,
                "middle",
                "black",
                "bold",
            );
        }
    }

    // Sample names along the left, and the peak count on the right, both against
    // the bar's own vertical centre as R's `mtext(at = b, ...)` does.
    for (sample_index, name) in summary.samples.iter().enumerate() {
        let top = bar_top(sample_index);
        let centre = top + bar_height * 0.4;
        panel.text_weighted(
            plot_left - 6.0,
            centre + font_size * 0.35,
            name,
            font_size,
            "end",
            "black",
            "bold",
        );
        // `mtext(side = 4, at = b, font = 4)` -- italic in R, approximated here
        // with the normal weight since the SVG writer has no italic primitive.
        panel.text(
            plot_right + 6.0,
            centre + font_size * 0.35,
            &summary.n_peaks[sample_index].to_string(),
            font_size,
            "start",
            "black",
        );
    }

    draw_category_legend(
        panel,
        summary,
        plot_left,
        plot_top,
        legend_font_size * font_size,
    );
}

/// Draws the category legend above the bars, in `ncol = 3` order.
///
/// R fills a multi-column legend by column, so the first `ceil(n/3)` entries form
/// the left column. Reproducing that matters for matching the layout.
fn draw_category_legend(
    panel: &mut PanelWriter,
    summary: &Summary,
    plot_left: f64,
    plot_top: f64,
    font_size: f64,
) {
    let entries = summary.categories.len();
    if entries == 0 {
        return;
    }
    let columns = 3usize;
    let rows = entries.div_ceil(columns);
    let line_height = font_size * 1.5;
    // R anchors the legend at "topleft" with `xjust = 0, yjust = 0`, i.e. its
    // top-left corner sits at the plot's top-left corner, and the rows grow
    // downward from there.
    let swatch = font_size * 1.2;
    let column_width = (panel.width() - plot_left) / columns as f64;

    for index in 0..entries {
        // Column-major fill, matching R's default `ncol` behaviour.
        let column = index / rows;
        let row = index % rows;
        let x = plot_left + column as f64 * column_width;
        let y = plot_top - 2.0 - (rows - 1 - row) as f64 * line_height;
        panel.rect(x, y - swatch * 0.7, swatch, swatch, summary.category_color(index));
        panel.text_weighted(
            x + swatch + 4.0,
            y,
            &summary.categories[index],
            font_size,
            "start",
            "black",
            "bold",
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotation_suffix_is_stripped_at_the_first_paren() {
        // HOMER's usual form: feature followed by the nearest feature.
        assert_eq!(
            normalize_annotation("promoter-TSS (NM_000001)"),
            Some("promoter-TSS".to_string())
        );
        assert_eq!(
            normalize_annotation("exon (NM_000002)"),
            Some("exon".to_string())
        );
        // No suffix at all, which HOMER writes for intergenic peaks.
        assert_eq!(
            normalize_annotation("Intergenic"),
            Some("Intergenic".to_string())
        );
    }

    #[test]
    fn utr_labels_are_renamed_to_palette_spellings() {
        // R applies these two `gsub()` calls because the palette keys are
        // "3pUTR"/"5pUTR" while HOMER writes "3' UTR"/"5' UTR".
        assert_eq!(
            normalize_annotation("3' UTR (NM_000003)"),
            Some("3pUTR".to_string())
        );
        assert_eq!(
            normalize_annotation("5' UTR (NM_000004)"),
            Some("5pUTR".to_string())
        );
    }

    #[test]
    fn a_literal_na_annotation_is_kept_as_a_category() {
        // HOMER's literal NA text is a real annotation value, so it is kept as a
        // category. R then drops it at the palette filter (see the next test),
        // which is what makes the drop reportable.
        assert_eq!(normalize_annotation("NA"), Some("NA".to_string()));
    }

    #[test]
    fn the_na_palette_entry_is_unreachable_unless_asked_for() {
        // The palette names a colour, and the raw lookup finds it...
        assert_eq!(palette_color(UNMAPPED_NA_CATEGORY), Some("gray70"));
        // ...but R's `%in%` lookup can never match it, so by default the port
        // reports it as unreachable too.
        assert!(!palette_reachable(UNMAPPED_NA_CATEGORY, false));
        // Opting in makes it reachable, which is what lets the bar reach 1.
        assert!(palette_reachable(UNMAPPED_NA_CATEGORY, true));
        // Other categories are unaffected either way.
        assert!(palette_reachable("exon", false));
        assert!(palette_reachable("exon", true));
    }

    #[test]
    fn an_empty_annotation_cell_is_missing() {
        assert_eq!(normalize_annotation(""), None);
        assert_eq!(normalize_annotation("   "), None);
    }

    #[test]
    fn fractions_divide_by_every_peak_not_just_the_kept_ones() {
        // The R code computes `fract := N / sum(N)` before filtering, so a
        // category outside the palette still shrinks the others. With 8 kept and
        // 2 dropped peaks, a 4-peak category is 0.4, not 0.5.
        let sample = Sample {
            name: "s".to_string(),
            n_peaks: 10,
            counts: BTreeMap::from([
                ("exon".to_string(), 4usize),
                ("intron".to_string(), 4),
                ("ncRNA".to_string(), 2),
            ]),
        };
        assert!((sample.fraction("exon") - 0.4).abs() < 1e-15);
        // The dropped category is reported through the summary rather than on
        // the sample, so a test-only accessor is not needed.
        let summary = build_summary(&[sample], false);
        assert_eq!(summary.dropped.len(), 1);
        assert_eq!(summary.dropped[0].name, "ncRNA");
        assert_eq!(summary.dropped[0].total(), 2);
        assert!(!summary.dropped[0].is_unmapped_na());
    }

    #[test]
    fn a_category_absent_from_a_sample_contributes_zero() {
        // `dcast(fill = 0)` gives a missing category a zero rather than dropping
        // the whole row, so a category present in only one sample still gets a
        // (zero-width) segment in the other.
        let first = Sample {
            name: "a".to_string(),
            n_peaks: 10,
            counts: BTreeMap::from([("exon".to_string(), 5usize), ("intron".to_string(), 5)]),
        };
        let second = Sample {
            name: "b".to_string(),
            n_peaks: 4,
            counts: BTreeMap::from([("exon".to_string(), 4usize)]),
        };
        let summary = build_summary(&[first, second], false);

        assert_eq!(summary.categories, vec!["exon", "intron"]);
        assert_eq!(summary.fractions[0], vec![0.5, 1.0]);
        assert_eq!(summary.fractions[1], vec![0.5, 0.0]);
    }

    #[test]
    fn row_order_follows_the_palette_not_the_data() {
        // Sample "b" is dominated by exon, sample "a" by promoter-TSS, and the
        // rows still come out in palette order.
        let first = Sample {
            name: "a".to_string(),
            n_peaks: 10,
            counts: BTreeMap::from([("promoter-TSS".to_string(), 9usize), ("exon".to_string(), 1)]),
        };
        let second = Sample {
            name: "b".to_string(),
            n_peaks: 10,
            counts: BTreeMap::from([("exon".to_string(), 9usize), ("promoter-TSS".to_string(), 1)]),
        };
        let summary = build_summary(&[first, second], false);

        // Palette order puts exon before promoter-TSS regardless of which sample
        // favours which.
        assert_eq!(summary.categories, vec!["exon", "promoter-TSS"]);
    }

    #[test]
    fn only_palette_categories_reachable_from_the_data_are_kept() {
        // A category R knows about but the data lacks must not appear, or every
        // bar would gain an empty segment and the legend a spurious entry.
        let sample = Sample {
            name: "s".to_string(),
            n_peaks: 3,
            counts: BTreeMap::from([("exon".to_string(), 3usize)]),
        };
        let summary = build_summary(&[sample], false);
        assert_eq!(summary.categories, vec!["exon"]);
    }

    #[test]
    fn dropped_categories_are_reported_with_their_counts() {
        let sample = Sample {
            name: "s".to_string(),
            n_peaks: 4,
            counts: BTreeMap::from([
                ("exon".to_string(), 1usize),
                ("ncRNA".to_string(), 2),
                ("UTR".to_string(), 1),
            ]),
        };
        let summary = build_summary(&[sample], false);

        assert_eq!(summary.categories, vec!["exon"]);
        // Dropped categories are listed in byte order, which is deterministic
        // and locale-independent (unlike R's locale-sensitive `sort()`), so the
        // report does not change with the environment. Uppercase therefore sorts
        // before lowercase.
        assert_eq!(
            summary
                .dropped
                .iter()
                .map(|d| d.name.as_str())
                .collect::<Vec<_>>(),
            vec!["UTR", "ncRNA"],
            "dropped categories should be listed in byte order"
        );
        // The drawn portion is only a quarter of the bar, which is exactly the
        // signal a caller needs to explain why the bar stops short.
        assert!((summary.column_sum(0) - 0.25).abs() < 1e-15);
    }

    #[test]
    fn legend_labels_report_the_raw_count_per_sample() {
        let first = Sample {
            name: "a".to_string(),
            n_peaks: 10,
            counts: BTreeMap::from([("exon".to_string(), 4usize)]),
        };
        let second = Sample {
            name: "b".to_string(),
            n_peaks: 4,
            counts: BTreeMap::from([("exon".to_string(), 4usize)]),
        };
        let samples = vec![first, second];
        let summary = build_summary(&samples, false);
        let labels = summary.legend_labels(&samples);

        // R's `leg` column is `paste0(Annotation, " [", N, "]")` per sample, so
        // the same category can carry different counts per sample.
        assert_eq!(labels, vec![vec!["exon [4]".to_string(), "exon [4]".to_string()]]);
    }

    #[test]
    fn a_sample_with_no_peaks_reports_zero_rather_than_nan() {
        let sample = Sample {
            name: "s".to_string(),
            n_peaks: 0,
            counts: BTreeMap::new(),
        };
        assert_eq!(sample.fraction("exon"), 0.0);
    }

    #[test]
    fn sample_names_come_from_the_first_dot_separated_field() {
        // R: `tstrsplit(basename(anno), "\\.", keep = 1)`, so "H3K27ac.homer.txt"
        // becomes "H3K27ac".
        assert_eq!(
            derive_sample_name(Path::new("/data/H3K27ac.homer.txt")),
            "H3K27ac"
        );
        assert_eq!(derive_sample_name(Path::new("plain.txt")), "plain");
        // No dot at all: keep the whole name.
        assert_eq!(derive_sample_name(Path::new("nodot")), "nodot");
    }

    #[test]
    fn palette_colors_resolve_for_every_entry() {
        for (name, color) in ANNOTATION_PALETTE {
            assert_eq!(palette_color(name), Some(color));
        }
        assert_eq!(palette_color("not-a-category"), None);
    }

    // ---------------------------------------------------------------------
    // R oracle
    //
    // `testdata/homer_r_oracle.tsv` was produced by R 4.6.0 by running the real
    // `summarize_homer_annots()` from trackplot.R against the committed fixtures
    // `testdata/H3K27ac.homer.txt` and `testdata/H3K4me3.homer.txt`, then
    // recomputing its internal stats expressions:
    //
    //   homer <- summarize_homer_annots(anno = paths, sample_names = NULL)
    //   s <- h[, .N, Annotation][, fract := N/sum(N)][order(N, decreasing = TRUE)]
    //   wide <- dcast(rbindlist(s, idcol = "Sample"), Annotation ~ Sample,
    //                 value.var = "fract", fill = 0)
    //   kept <- wide[names(pie.cols)[names(pie.cols) %in% rownames(wide)], , drop = FALSE]
    //
    // Regenerate with those expressions when the oracle needs extending.
    // ---------------------------------------------------------------------

    /// Relative tolerance for values R prints with 17 significant digits.
    fn oracle_close(actual: f64, expected: f64) -> bool {
        if actual == expected {
            return true;
        }
        let scale = actual.abs().max(expected.abs()).max(1.0);
        (actual - expected).abs() <= 1e-12 * scale
    }

    fn fixture(name: &str) -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("testdata")
            .join(name)
    }

    #[test]
    fn aggregates_match_the_r_homer_oracle() {
        let path = fixture("homer_r_oracle.tsv");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("read oracle {path:?}: {error}"));

        // R derives names from the file base up to the first dot.
        let first = parse_sample(&fixture("H3K27ac.homer.txt"), "").expect("parse H3K27ac");
        let second = parse_sample(&fixture("H3K4me3.homer.txt"), "").expect("parse H3K4me3");
        assert_eq!(first.name, "H3K27ac");
        assert_eq!(second.name, "H3K4me3");

        let samples = vec![first.clone(), second.clone()];
        let summary = build_summary(&samples, false);

        // Peak totals are R's `npeaks`, and the fraction denominator.
        assert_eq!(summary.n_peaks, vec![24, 14]);
        assert_eq!(summary.samples, vec!["H3K27ac", "H3K4me3"]);

        let mut checked = 0usize;
        for line in text.lines() {
            let line = line.trim_end();
            if line.is_empty() || line.starts_with('#') || line.starts_with("annotation\t") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 5, "malformed oracle row: {line:?}");
            let category = fields[0];
            let sample = fields[1];
            let expected_fraction: f64 = fields[2].parse().expect("fraction");
            let expected_legend = fields[4];

            let sample_index = summary
                .samples
                .iter()
                .position(|name| name == sample)
                .unwrap_or_else(|| panic!("oracle names unknown sample {sample}"));
            let category_index = summary
                .categories
                .iter()
                .position(|name| name == category)
                .unwrap_or_else(|| panic!("{category} missing from the built summary"));

            let actual = summary.fractions[category_index][sample_index];
            assert!(
                oracle_close(actual, expected_fraction),
                "{category}/{sample}: fraction got {actual}, R has {expected_fraction}"
            );

            // The `Annotation [N]` legend text, for samples that have the
            // category. R writes "NA" there for a category absent from a sample,
            // because the `dcast` fill is 0 and the join finds no row.
            let labels = summary.legend_labels(&samples);
            if expected_legend != "NA" {
                assert_eq!(
                    labels[category_index][sample_index], expected_legend,
                    "{category}/{sample}: legend text differs from R"
                );
            }

            checked += 1;
        }
        assert!(checked >= 16, "oracle only had {checked} rows");
        eprintln!("validated {checked} HOMER annotation rows against the R oracle");
    }

    #[test]
    fn the_committed_fixture_reproduces_rs_dropped_categories() {
        // The oracle records these two as dropped: `ncRNA` because the palette
        // has no entry, and `NA` because R's palette lookup can never match it.
        let first = parse_sample(&fixture("H3K27ac.homer.txt"), "").expect("parse H3K27ac");
        let second = parse_sample(&fixture("H3K4me3.homer.txt"), "").expect("parse H3K4me3");
        let summary = build_summary(&[first.clone(), second.clone()], false);

        let dropped: Vec<&str> = summary.dropped.iter().map(|d| d.name.as_str()).collect();
        assert_eq!(dropped, vec!["NA", "ncRNA"]);

        // R's kept fractions sum to 0.875 and 0.928571...; the remainder is
        // exactly the dropped share (1/24 + 2/24, and 1/14).
        assert!(oracle_close(summary.column_sum(0), 0.875));
        assert!(oracle_close(summary.column_sum(1), 13.0 / 14.0));

        // And drawing the `NA` segment recovers only that 1-peak share while
        // leaving every other fraction untouched -- the denominator does not
        // change, so the bars simply extend by the previously invisible sliver.
        let kept = build_summary(&[first, second], true);

        assert!(
            !kept.dropped.iter().any(|d| d.name == UNMAPPED_NA_CATEGORY),
            "the NA category should be present once it is parsed as a string"
        );
        // H3K4me3 has no unannotated peaks, so only H3K27ac gains the segment.
        assert!(oracle_close(kept.column_sum(0), 0.875 + 1.0 / 24.0));
        assert!(oracle_close(kept.column_sum(1), 13.0 / 14.0));

        // A shared category keeps the same fraction in both builds.
        let exon = |s: &Summary| {
            let index = s.categories.iter().position(|c| c == "exon").unwrap();
            s.fractions[index][0]
        };
        assert!(oracle_close(exon(&summary), exon(&kept)));
    }
}
